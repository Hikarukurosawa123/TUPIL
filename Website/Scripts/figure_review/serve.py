"""Serve the Results_Main figure review workspace using Python's standard library."""

import argparse
import errno
import hmac
import ipaddress
import json
import os
from pathlib import Path
from datetime import datetime, timezone
from queue import Empty, Queue
import re
import secrets
import shutil
import socket
import sqlite3
import subprocess
import threading
import time
from http import HTTPStatus
from http.cookies import SimpleCookie
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlsplit


DRIVE_FOLDER = Path(r"G:\Shared drives\TUPIL\Hikaru Kurosawa - MBP PhD Student - 2026\Results Storage\Results_Main")
ASSETS = Path(__file__).with_name("assets")
TYPES = {
    ".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg",
    ".webp": "image/webp", ".gif": "image/gif", ".avif": "image/avif",
    ".svg": "image/svg+xml", ".bmp": "image/bmp", ".pdf": "application/pdf",
}
STATIC = {
    "/": ("notebook.html", "text/html; charset=utf-8"),
    "/annotations": ("index.html", "text/html; charset=utf-8"),
    "/app.js": ("app.js", "text/javascript; charset=utf-8"),
    "/style.css": ("style.css", "text/css; charset=utf-8"),
    "/notebook.js": ("notebook.js", "text/javascript; charset=utf-8"),
    "/notebook.css": ("notebook.css", "text/css; charset=utf-8"),
    "/markdown.js": ("markdown.js", "text/javascript; charset=utf-8"),
    "/anchors.js": ("anchors.js", "text/javascript; charset=utf-8"),
}
PAGE_TITLES = {"results-1": "Results 1", "results-2": "Results 2", "project-1": "Project 1", "project-2": "Project 2"}
PAGE_SECTIONS = ("results", "notebook")
INITIAL_NAVIGATION = (("results-1", "results", 0), ("results-2", "results", 1),
                      ("project-1", "notebook", 0), ("project-2", "notebook", 1))
BLOCK_ID = re.compile(r"(?:block-[a-f0-9]{32}|legacy-\d{1,3})", re.ASCII)
CLOUDFLARE_URL = re.compile(r"https://[a-z0-9-]+\.trycloudflare\.com\b")
PINGGY_URL = re.compile(r"https://[a-z0-9-]+\.run\.pinggy-free\.link\b")
PINGGY_ID = re.compile(r"Tunnel started \(ID: ([a-f0-9-]+)\)")
HOSTNAME = re.compile(r"(?:[a-z0-9](?:[a-z0-9-]{0,61}[a-z0-9])?\.)+[a-z][a-z0-9-]{1,62}", re.ASCII)
NGROK_ERROR = re.compile(r"\bERR_NGROK_\d+\b")
CONTENT_POLICY = (
    "default-src 'none'; script-src 'self'; style-src 'self'; img-src 'self'; "
    "frame-src 'self'; connect-src 'self'; base-uri 'none'; frame-ancestors 'none'; form-action 'self'"
)


class ReviewServer(ThreadingHTTPServer):
    daemon_threads = True
    allow_reuse_address = os.name != "nt"

    def __init__(self, args):
        self.root = args.results.resolve()
        self.notes = args.notes.resolve()
        self.remote = args.remote
        self.remote_host = ""
        self.access_code = secrets.token_urlsafe(32 if args.remote else 12) if args.remote or args.share else None
        self.session_token = secrets.token_urlsafe(32)
        self.notes.parent.mkdir(parents=True, exist_ok=True)
        with self.database() as db:
            db.execute("""CREATE TABLE IF NOT EXISTS notes (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                figure TEXT NOT NULL, text TEXT NOT NULL,
                x REAL NOT NULL, y REAL NOT NULL
            )""")
            db.execute("CREATE INDEX IF NOT EXISTS notes_figure ON notes (figure)")
            db.execute("""CREATE TABLE IF NOT EXISTS notebook_pages (
                page TEXT PRIMARY KEY, revision INTEGER NOT NULL, blocks TEXT NOT NULL
            )""")
            db.execute("CREATE TABLE IF NOT EXISTS notebook_titles (page TEXT PRIMARY KEY, title TEXT NOT NULL)")
            db.execute("""CREATE TABLE IF NOT EXISTS notebook_navigation (
                page TEXT PRIMARY KEY, section TEXT NOT NULL, position INTEGER NOT NULL
            )""")
            db.executemany("INSERT OR IGNORE INTO notebook_navigation VALUES (?, ?, ?)", INITIAL_NAVIGATION)
            db.execute("CREATE TABLE IF NOT EXISTS notebook_trash (page TEXT PRIMARY KEY)")
            db.execute("""CREATE TABLE IF NOT EXISTS notebook_comments (
                id INTEGER PRIMARY KEY, page TEXT NOT NULL, author TEXT NOT NULL,
                text TEXT NOT NULL, created TEXT NOT NULL
            )""")
            db.execute("CREATE INDEX IF NOT EXISTS notebook_comments_page ON notebook_comments (page, id)")
            db.execute("""CREATE TABLE IF NOT EXISTS notebook_comment_anchors (
                id INTEGER PRIMARY KEY, block TEXT NOT NULL, mode TEXT NOT NULL,
                offset INTEGER NOT NULL, quote TEXT NOT NULL
            )""")
        super().__init__(("0.0.0.0" if args.share else "127.0.0.1", args.port), Handler)

    def server_bind(self):
        if os.name == "nt":
            self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_EXCLUSIVEADDRUSE, 1)
        super().server_bind()

    def database(self):
        return sqlite3.connect(self.notes, factory=ClosingConnection)

    def restore_invite(self, hostname):
        with self.database() as db:
            db.execute("CREATE TABLE IF NOT EXISTS remote_invites (hostname TEXT PRIMARY KEY, code TEXT NOT NULL)")
            db.execute("INSERT OR IGNORE INTO remote_invites VALUES (?, ?)", (hostname, self.access_code))
            self.access_code = db.execute("SELECT code FROM remote_invites WHERE hostname = ?", (hostname,)).fetchone()[0]

    def figure(self, key):
        relative = Path(key)
        if not key or relative.anchor or any(part == '..' for part in relative.parts):
            raise ValueError("Invalid figure path.")
        path = (self.root / relative).resolve()
        if not path.is_relative_to(self.root) or path.suffix.lower() not in TYPES or not path.is_file():
            raise FileNotFoundError
        return path


class ClosingConnection(sqlite3.Connection):
    def __exit__(self, *args):
        try:
            return super().__exit__(*args)
        finally:
            self.close()


class PageUnavailable(Exception):
    pass


class Handler(BaseHTTPRequestHandler):
    server_version = "FigureReview"
    sys_version = ""

    def send_body(self, status, body, mime, *, cookie=None, policy=CONTENT_POLICY):
        self.send_headers(status, len(body), mime, cookie=cookie, policy=policy)
        self.wfile.write(body)

    def send_headers(self, status, size, mime, *, cookie=None, policy=CONTENT_POLICY):
        self.send_response(status)
        self.send_header("Content-Type", mime)
        self.send_header("Content-Length", str(size))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        if self.server.remote:
            self.send_header("X-Robots-Tag", "noindex, nofollow")
        if policy:
            self.send_header("Content-Security-Policy", policy)
        if cookie:
            self.send_header("Set-Cookie", cookie)
        self.end_headers()

    def json(self, data, status=HTTPStatus.OK, **kwargs):
        self.send_body(status, json.dumps(data).encode(), "application/json; charset=utf-8", **kwargs)

    def valid_host(self):
        host = self.headers.get("Host", "")
        if self.server.remote:
            forwarded = self.headers.get_all("X-Forwarded-Proto", [])
            protocol = forwarded[-1].rsplit(",", 1)[-1].strip() if forwarded else ""
            return bool(host) and host == self.server.remote_host and protocol == "https"
        name, separator, port = host.rpartition(":")
        if not separator or port != str(self.server.server_port):
            return False
        if name.lower() == "localhost" or name == "127.0.0.1":
            return True
        if self.server.access_code:
            try:
                return ipaddress.ip_address(name).version == 4
            except ValueError:
                pass
        return False

    def authorized(self):
        if self.server.access_code is None:
            return True
        cookies = SimpleCookie(self.headers.get("Cookie", ""))
        token = cookies.get(self.session_cookie_name())
        return token is not None and hmac.compare_digest(token.value.encode(), self.server.session_token.encode())

    def session_cookie_name(self):
        return "__Host-figure_session" if self.server.remote else "figure_session"

    def do_GET(self):
        self.dispatch()

    def do_POST(self):
        self.dispatch()

    def dispatch(self):
        try:
            if not self.valid_host():
                return self.json({"error": "Use the address printed by the launcher."}, HTTPStatus.BAD_REQUEST)
            url = urlsplit(self.path)
            query = parse_qs(url.query)
            if self.command == "POST":
                protocol = "https" if self.server.remote else "http"
                if (self.headers.get("Origin") != f"{protocol}://{self.headers['Host']}"
                        or self.headers.get_content_type() != "application/json"):
                    return self.json({"error": "Open the website to make changes."}, HTTPStatus.FORBIDDEN)
                length = int(self.headers.get("Content-Length", 0))
                limit = 6_000_000 if url.path == "/api/notebook" else 65536
                if not 0 < length <= limit:
                    raise ValueError("The request exceeds the allowed size.")
                body = json.loads(self.rfile.read(length))
                if not isinstance(body, dict):
                    raise ValueError("Expected a JSON object.")
                if url.path == "/api/login":
                    code = body.get("code")
                    if (self.server.access_code and isinstance(code, str)
                            and hmac.compare_digest(code.strip().encode(), self.server.access_code.encode())):
                        cookie = f"{self.session_cookie_name()}={self.server.session_token}; Path=/; HttpOnly; SameSite=Strict"
                        if self.server.remote:
                            cookie += "; Secure"
                        return self.json({"ok": True}, cookie=cookie)
                    return self.json({"error": "That access code did not match."}, HTTPStatus.UNAUTHORIZED)
            elif url.path in STATIC:
                filename, mime = STATIC[url.path]
                return self.send_body(HTTPStatus.OK, (ASSETS / filename).read_bytes(), mime)
            if not self.authorized():
                return self.json({"error": "Enter the access code to open the figures."}, HTTPStatus.UNAUTHORIZED)
            if self.command == "POST" and url.path == "/api/notes":
                return self.save_note(body)
            if url.path == "/api/notebook":
                if self.command == "POST":
                    return self.save_page(body)
                return self.get_page(query.get("page", [""])[0])
            if url.path == "/api/notebook/title" and self.command == "POST":
                return self.rename_page(body)
            if url.path == "/api/notebook/pages":
                if self.command == "POST":
                    return self.create_page(body)
                with self.server.database() as db:
                    return self.json(self.page_index(db))
            if url.path == "/api/notebook/trash" and self.command == "POST":
                return self.trash_page(body)
            if url.path == "/api/notebook/move" and self.command == "POST":
                return self.move_page(body)
            if url.path == "/api/notebook/comments":
                if self.command == "POST":
                    return self.post_comment(body)
                return self.get_comments(query.get("page", [""])[0])
            if self.command == "GET":
                key = query.get("figure", [""])[0]
                if url.path == "/api/figures":
                    return self.list_figures()
                if url.path == "/api/notes":
                    self.server.figure(key)
                    with self.server.database() as db:
                        db.row_factory = sqlite3.Row
                        rows = db.execute("SELECT id, text, x, y FROM notes WHERE figure = ? ORDER BY id", (key,))
                        return self.json({"notes": [dict(row) for row in rows]})
                if url.path == "/figure":
                    path = self.server.figure(key)
                    mime = TYPES[path.suffix.lower()]
                    with path.open("rb") as file:
                        policy = "sandbox; default-src 'none'; style-src 'unsafe-inline'; img-src data:" if mime == "image/svg+xml" else None
                        self.send_headers(HTTPStatus.OK, os.fstat(file.fileno()).st_size, mime, policy=policy)
                        shutil.copyfileobj(file, self.wfile)
                    return
            self.json({"error": "Page not found."}, HTTPStatus.NOT_FOUND)
        except (BrokenPipeError, ConnectionResetError, ConnectionAbortedError):
            return
        except PageUnavailable as exc:
            self.json({"error": str(exc)}, HTTPStatus.NOT_FOUND)
        except FileNotFoundError:
            self.json({"error": "The figure or folder is unavailable."}, HTTPStatus.NOT_FOUND)
        except ValueError as exc:
            self.json({"error": str(exc)}, HTTPStatus.BAD_REQUEST)
        except (OSError, sqlite3.Error) as exc:
            self.log_error("%s", exc)
            self.json({"error": "The server could not read the figures or save the notes."}, HTTPStatus.SERVICE_UNAVAILABLE)

    def list_figures(self):
        if not self.server.root.is_dir():
            return self.json({"error": "Results_Main is unavailable. Check your Google Drive connection and refresh."}, HTTPStatus.SERVICE_UNAVAILABLE)

        def on_error(error):
            raise error

        figures = []
        for folder, dirs, files in os.walk(self.server.root, onerror=on_error):
            dirs[:] = [name for name in dirs if not name.startswith('.')
                       and (Path(folder) / name).resolve().is_relative_to(self.server.root)]
            for name in files:
                path = Path(folder) / name
                if not name.startswith('.') and path.suffix.lower() in TYPES and path.resolve().is_relative_to(self.server.root):
                    figures.append(path.relative_to(self.server.root).as_posix())
        figures.sort(key=str.casefold)
        self.json({"figures": figures})

    def get_page(self, page):
        with self.server.database() as db:
            self.require_page(db, page)
            row = db.execute("SELECT revision, blocks FROM notebook_pages WHERE page = ?", (page,)).fetchone()
            index = self.page_index(db)
        revision, blocks = (row[0], self.blocks_with_ids(row[1])) if row else (0, [])
        self.json({"revision": revision, "blocks": blocks, **index,
                   "titles": {item["id"]: item["title"] for item in index["pages"]}})

    @staticmethod
    def require_page(db, page):
        if not isinstance(page, str) or not db.execute("""SELECT 1 FROM notebook_navigation AS n
                WHERE n.page = ? AND NOT EXISTS (SELECT 1 FROM notebook_trash AS t WHERE t.page = n.page)""", (page,)).fetchone():
            raise PageUnavailable("The page is unavailable or in Trash.")

    @staticmethod
    def page_index(db):
        rows = db.execute("""SELECT n.page, n.section, t.title, trash.page IS NOT NULL
                             FROM notebook_navigation AS n
                             LEFT JOIN notebook_titles AS t ON t.page = n.page
                             LEFT JOIN notebook_trash AS trash ON trash.page = n.page
                             ORDER BY n.section, n.position, n.page""")
        result = {"pages": [], "trash": []}
        for page, section, title, trashed in rows:
            result["trash" if trashed else "pages"].append({"id": page, "section": section,
                "title": title or PAGE_TITLES.get(page, "Untitled page")})
        return result

    @staticmethod
    def blocks_with_ids(serialized):
        return [dict(block, id=block.get("id", f"legacy-{index}")) for index, block in enumerate(json.loads(serialized))]

    @staticmethod
    def clean_title(title):
        if not isinstance(title, str) or not 1 <= len(title.strip()) <= 120 or len(title.splitlines()) != 1:
            raise ValueError("Enter a page title between 1 and 120 characters on one line.")
        return title.strip()

    def create_page(self, body):
        section = body.get("section")
        if section not in PAGE_SECTIONS:
            raise ValueError("Choose a section for the new page.")
        title = self.clean_title(body.get("title"))
        page = f"page-{secrets.token_hex(12)}"
        with self.server.database() as db:
            db.execute("BEGIN IMMEDIATE")
            position = db.execute("SELECT COALESCE(MAX(position) + 1, 0) FROM notebook_navigation WHERE section = ?",
                                  (section,)).fetchone()[0]
            db.execute("INSERT INTO notebook_navigation VALUES (?, ?, ?)", (page, section, position))
            db.execute("INSERT INTO notebook_titles VALUES (?, ?)", (page, title))
            index = self.page_index(db)
        self.json({"id": page, **index}, HTTPStatus.CREATED)

    def move_page(self, body):
        page, section, after = body.get("page"), body.get("section"), body.get("after")
        if section not in PAGE_SECTIONS or (after is not None and not isinstance(after, str)):
            raise ValueError("Choose a section and position for the page.")
        with self.server.database() as db:
            db.execute("BEGIN IMMEDIATE")
            self.require_page(db, page)
            order = [row[0] for row in db.execute("""SELECT n.page FROM notebook_navigation AS n
                      WHERE n.section = ? AND n.page <> ? AND NOT EXISTS (SELECT 1 FROM notebook_trash AS t WHERE t.page = n.page)
                      ORDER BY n.position, n.page""", (section, page))]
            if after is not None and after not in order:
                raise ValueError("That position changed. Reopen the move form and try again.")
            position = 0 if after is None else order.index(after) + 1
            order.insert(position, page)
            db.executemany("UPDATE notebook_navigation SET section = ?, position = ? WHERE page = ?",
                           ((section, index, item) for index, item in enumerate(order)))
            index = self.page_index(db)
        self.json(index)

    def trash_page(self, body):
        page, action = body.get("page"), body.get("action")
        if not isinstance(page, str) or action not in ("delete", "restore"):
            raise ValueError("Choose a page and a valid Trash action.")
        with self.server.database() as db:
            db.execute("BEGIN IMMEDIATE")
            if action == "delete":
                self.require_page(db, page)
                db.execute("INSERT INTO notebook_trash VALUES (?)", (page,))
            else:
                restored = db.execute("DELETE FROM notebook_trash WHERE page = ?", (page,))
                if not restored.rowcount:
                    raise PageUnavailable("The page is no longer in Trash.")
                db.execute("""UPDATE notebook_navigation SET position = (
                              SELECT COALESCE(MAX(peer.position) + 1, 0) FROM notebook_navigation AS peer
                              WHERE peer.section = notebook_navigation.section
                            ) WHERE page = ?""", (page,))
            index = self.page_index(db)
        self.json(index)

    def rename_page(self, body):
        page, title = body.get("page"), self.clean_title(body.get("title"))
        with self.server.database() as db:
            db.execute("BEGIN IMMEDIATE")
            self.require_page(db, page)
            db.execute("""INSERT INTO notebook_titles (page, title) VALUES (?, ?)
                          ON CONFLICT(page) DO UPDATE SET title = excluded.title""", (page, title))
        self.json({"title": title})

    def get_comments(self, page):
        with self.server.database() as db:
            self.require_page(db, page)
            db.row_factory = sqlite3.Row
            rows = db.execute("""SELECT c.id, c.author, c.text, c.created, a.block, a.mode, a.offset, a.quote
                FROM notebook_comments AS c LEFT JOIN notebook_comment_anchors AS a ON a.id = c.id
                WHERE c.page = ? ORDER BY c.id DESC""", (page,))
            comments = []
            for row in rows:
                comment = {key: row[key] for key in ("id", "author", "text", "created")}
                if row["block"] is not None:
                    comment["anchor"] = {key: row[key] for key in ("block", "mode", "offset", "quote")}
                comments.append(comment)
        self.json({"comments": comments})

    def validate_anchor(self, db, page, anchor):
        if not isinstance(anchor, dict):
            raise ValueError("Select some document text before commenting.")
        block, mode, offset, quote = (anchor.get(key) for key in ("block", "mode", "offset", "quote"))
        if (not isinstance(block, str) or not BLOCK_ID.fullmatch(block) or mode not in ("source", "preview")
                or type(offset) is not int or not 0 <= offset <= 40_000
                or not isinstance(quote, str) or not quote.strip() or len(quote) > 2_000):
            raise ValueError("Select up to 2,000 characters in one text block.")
        row = db.execute("SELECT blocks FROM notebook_pages WHERE page = ?", (page,)).fetchone()
        target = next((item for item in self.blocks_with_ids(row[0]) if item["id"] == block), None) if row else None
        if target is None or (mode == "source" and quote not in target["markdown"].replace("\r\n", "\n").replace("\r", "\n")):
            raise ValueError("The selected text has changed. Select it again before posting.")
        return {"block": block, "mode": mode, "offset": offset, "quote": quote}

    def post_comment(self, body):
        page, author, comment = body.get("page"), body.get("author", ""), body.get("text")
        if not isinstance(author, str) or len(author.strip()) > 80 or len(author.splitlines()) > 1:
            raise ValueError("Enter a name of up to 80 characters on one line.")
        if not isinstance(comment, str) or not 1 <= len(comment.strip()) <= 4_000:
            raise ValueError("Enter a comment between 1 and 4,000 characters.")
        author, comment = author.strip() or "Reviewer", comment.strip()
        created = datetime.now(timezone.utc).isoformat(timespec="seconds")
        with self.server.database() as db:
            db.execute("BEGIN IMMEDIATE")
            self.require_page(db, page)
            anchor = self.validate_anchor(db, page, body["anchor"]) if body.get("anchor") is not None else None
            row = db.execute("INSERT INTO notebook_comments (page, author, text, created) VALUES (?, ?, ?, ?)",
                             (page, author, comment, created))
            comment = {"id": row.lastrowid, "author": author, "text": comment, "created": created}
            if anchor is not None:
                db.execute("INSERT INTO notebook_comment_anchors VALUES (?, ?, ?, ?, ?)",
                           (row.lastrowid, anchor["block"], anchor["mode"], anchor["offset"], anchor["quote"]))
                comment["anchor"] = anchor
        self.json({"comment": comment}, HTTPStatus.CREATED)

    def save_page(self, body):
        page, revision, blocks = body.get("page"), body.get("revision"), body.get("blocks")
        if (type(revision) is not int
                or not 0 <= revision < 2**63 - 1 or not isinstance(blocks, list) or len(blocks) > 200):
            raise ValueError("Invalid notebook page. A page can hold up to 200 blocks.")
        cleaned = []
        seen = set()
        for index, block in enumerate(blocks):
            if not isinstance(block, dict) or block.get("kind") not in ("text", "figure"):
                raise ValueError("Invalid notebook block.")
            kind = block["kind"]
            block_id = block.get("id", f"legacy-{index}")
            if not isinstance(block_id, str) or not BLOCK_ID.fullmatch(block_id) or block_id in seen:
                raise ValueError("Invalid or duplicate notebook block ID.")
            seen.add(block_id)
            markdown = block.get("markdown")
            if not isinstance(markdown, str) or len(markdown) > 20_000:
                raise ValueError("A text block can hold up to 20,000 characters.")
            item = {"id": block_id, "kind": kind, "markdown": markdown}
            if kind == "figure":
                figure, layout = block.get("figure"), block.get("layout")
                if not isinstance(figure, str) or not 1 <= len(figure) <= 4096 or layout not in ("full", "left", "right"):
                    raise ValueError("Invalid notebook figure.")
                relative = Path(figure)
                if relative.anchor or ".." in relative.parts or relative.suffix.lower() not in TYPES or "\x00" in figure:
                    raise ValueError("Invalid notebook figure path.")
                item.update(figure=figure, layout=layout)
            cleaned.append(item)
        serialized = json.dumps(cleaned, ensure_ascii=True)
        with self.server.database() as db:
            db.execute("BEGIN IMMEDIATE")
            self.require_page(db, page)
            if revision == 0:
                result = db.execute("INSERT OR IGNORE INTO notebook_pages VALUES (?, 1, ?)", (page, serialized))
            else:
                result = db.execute("UPDATE notebook_pages SET revision = revision + 1, blocks = ? WHERE page = ? AND revision = ?",
                                    (serialized, page, revision))
            if result.rowcount == 0:
                return self.json({"error": "Someone else saved this page. Copy your changes before reloading it."}, HTTPStatus.CONFLICT)
        self.json({"revision": revision + 1})

    def save_note(self, body):
        key = body.get("figure")
        if not isinstance(key, str):
            raise ValueError("Choose a figure first.")
        path = self.server.figure(key)
        key = path.relative_to(self.server.root).as_posix()
        action, note_id = body.get("action"), body.get("id")
        if action != "add" and (type(note_id) is not int or not 0 < note_id <= 2**63 - 1):
            raise ValueError("Invalid note ID.")
        if action in ("add", "edit"):
            note = body.get("text")
            if not isinstance(note, str) or not 1 <= len(note.strip()) <= 2000:
                raise ValueError("Enter a note between 1 and 2,000 characters.")
            note = note.strip()
        if action in ("add", "move"):
            x, y = body.get("x"), body.get("y")
            if any(type(value) not in (int, float) or not 0 <= value <= 1 for value in (x, y)):
                raise ValueError("Choose a position in the figure area.")
        with self.server.database() as db:
            if action == "add":
                result = db.execute("INSERT INTO notes (figure, text, x, y) VALUES (?, ?, ?, ?)", (key, note, x, y))
                note_id = result.lastrowid
            elif action == "edit":
                result = db.execute("UPDATE notes SET text = ? WHERE id = ? AND figure = ?", (note, note_id, key))
            elif action == "move":
                result = db.execute("UPDATE notes SET x = ?, y = ? WHERE id = ? AND figure = ?", (x, y, note_id, key))
            elif action == "delete":
                result = db.execute("DELETE FROM notes WHERE id = ? AND figure = ?", (note_id, key))
            else:
                raise ValueError("Unknown note action.")
            if result.rowcount == 0:
                return self.json({"error": "This note was already removed. Refresh the notes."}, HTTPStatus.NOT_FOUND)
        self.json({"id": note_id})


def serve_invite(server, url, persistent=False, on_ready=None):
    server.remote_host = urlsplit(url).hostname
    if persistent:
        server.restore_invite(server.remote_host)
    if on_ready is None:
        print(f"Private invite link: {url}/#invite={server.access_code}", flush=True)
    else:
        on_ready(url)
    server.serve_forever()


def serve_pinggy(server, executable):
    executable = Path(executable).resolve()

    def pinggy(args):
        return subprocess.run([str(executable), *args], cwd=executable.parent, capture_output=True,
                              encoding="utf-8", errors="replace", timeout=30)

    result = pinggy(["-l", f"http://127.0.0.1:{server.server_port}", "--b", "x:https"])
    tunnel_id = PINGGY_ID.search(result.stdout)
    url = PINGGY_URL.search(result.stdout)
    expiry = threading.Timer(3600, server.shutdown)
    expiry.daemon = True
    try:
        if result.returncode or not tunnel_id or not url:
            raise RuntimeError("Pinggy could not start an HTTPS tunnel.\n" + (result.stderr or result.stdout).strip())
        expiry.start()
        serve_invite(server, url.group())
    finally:
        expiry.cancel()
        if tunnel_id:
            pinggy(["stop", tunnel_id.group(1)])


def serve_cloudflare(server, executable):
    origin = f"http://127.0.0.1:{server.server_port}"
    command = [executable, "tunnel", "--no-autoupdate", "--url", origin]

    url = None

    def endpoint(line):
        nonlocal url
        match = CLOUDFLARE_URL.search(line)
        if match:
            url = match.group()
        return url if "Registered tunnel connection" in line else None

    run_tunnel_process(server, command, endpoint, "Cloudflare could not connect the Quick Tunnel. Check outbound port 7844.", True)


def serve_ngrok(server, executable, domain, on_ready=None):
    origin = f"http://127.0.0.1:{server.server_port}"
    command = [executable, "http", origin, "--inspect=false", "--log", "stdout", "--log-format", "json"]
    if domain:
        command.extend(["--url", f"https://{domain}"])
    error_code = None

    def endpoint(line):
        nonlocal error_code
        match = NGROK_ERROR.search(line)
        if match:
            error_code = match.group()
        try:
            record = json.loads(line)
            candidate = record.get("url") if isinstance(record, dict) else None
            parsed = urlsplit(candidate) if isinstance(candidate, str) else None
        except ValueError:
            return None
        if (parsed and parsed.scheme == "https" and HOSTNAME.fullmatch(parsed.netloc)
                and parsed.path in ("", "/") and not parsed.query and not parsed.fragment
                and (not domain or parsed.netloc == domain)):
            return f"https://{parsed.netloc}"
        return None

    def failure(_):
        if error_code == "ERR_NGROK_4018":
            return "ngrok needs an authenticated account. Run ngrok config add-authtoken with the token from your dashboard."
        message = "ngrok could not start. Check your ngrok authentication and domain settings."
        return f"{message} ({error_code})" if error_code else message

    run_tunnel_process(server, command, endpoint, failure, persistent=True, on_ready=on_ready)


def run_tunnel_process(server, command, endpoint, failure_message, show_logs=False, persistent=False, on_ready=None):
    events = Queue()
    ready = threading.Event()
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, encoding="utf-8", errors="replace") as tunnel:
        def read_output():
            for line in tunnel.stdout:
                if not ready.is_set():
                    events.put(line.rstrip())
            events.put(None)

        threading.Thread(target=read_output, daemon=True).start()
        try:
            deadline = time.monotonic() + 45
            recent = []
            while time.monotonic() < deadline:
                try:
                    line = events.get(timeout=max(0, deadline - time.monotonic()))
                except Empty:
                    break
                if line is None:
                    break
                recent = [*recent[-5:], line]
                url = endpoint(line)
                if url:
                    ready.set()
                    break
            if not ready.is_set():
                detail = "\n" + "\n".join(recent) if show_logs else ""
                message = failure_message(recent) if callable(failure_message) else failure_message
                raise RuntimeError(message + detail)

            def stop_when_tunnel_exits():
                tunnel.wait()
                server.shutdown()

            threading.Thread(target=stop_when_tunnel_exits, daemon=True).start()
            serve_invite(server, url, persistent, on_ready)
            raise RuntimeError("The tunnel agent stopped. Restart the launcher to bring the invite back online.")
        finally:
            if tunnel.poll() is None:
                tunnel.terminate()
                try:
                    tunnel.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    tunnel.kill()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, help="Local port. Defaults to 8765 locally and an available port in remote mode.")
    parser.add_argument("--results", type=Path, default=DRIVE_FOLDER)
    state = Path(os.environ.get("LOCALAPPDATA", Path.home())) / "TUPIL" / "FigureReview"
    parser.add_argument("--notes", type=Path, default=state / "notes.sqlite3")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--share", action="store_true", help="Enable access from the local network with a code.")
    mode.add_argument("--remote", action="store_true", help="Create a private HTTPS invite link.")
    parser.add_argument("--provider", choices=("pinggy", "cloudflare", "ngrok"), default="pinggy", help="HTTPS tunnel provider.")
    parser.add_argument("--domain", help="Override ngrok's automatically assigned hostname.")
    parser.add_argument("--tunnel-executable", type=Path, help="Path to the tunnel provider's executable.")
    args = parser.parse_args()
    if args.port is None:
        args.port = 0 if args.remote else 8765
    if not 0 <= args.port <= 65535:
        parser.error("--port must be between 0 and 65535.")
    if args.remote and args.provider == "ngrok" and args.domain is not None:
        args.domain = args.domain.lower().strip()
        if len(args.domain) > 253 or not HOSTNAME.fullmatch(args.domain):
            parser.error("Pass your assigned ngrok hostname to --domain, without https:// or a path.")
    if args.remote:
        cli = {"pinggy": "pinggy", "cloudflare": "cloudflared", "ngrok": "ngrok"}[args.provider]
        cached = state / "bin" / (cli + ".exe")
        executable = str(args.tunnel_executable.resolve()) if args.tunnel_executable else shutil.which(cli)
        if not executable and cached.is_file():
            executable = str(cached)
        if not executable:
            parser.error(f"Install {cli} or provide its path with --tunnel-executable.")
    try:
        server = ReviewServer(args)
    except OSError as exc:
        if exc.errno != errno.EADDRINUSE and getattr(exc, "winerror", None) != 10048:
            raise
        if not args.remote or args.port == 0:
            parser.exit(1, f"Port {args.port} is already in use. Run with --port 0 to choose an available port.\n")
        print(f"Port {args.port} is already in use. Using an available port for this invite.", flush=True)
        args.port = 0
        server = ReviewServer(args)
    with server:
        port = server.server_port
        if args.remote:
            try:
                if args.provider == "ngrok":
                    serve_ngrok(server, executable, args.domain)
                else:
                    start_tunnel = serve_pinggy if args.provider == "pinggy" else serve_cloudflare
                    start_tunnel(server, executable)
            except KeyboardInterrupt:
                pass
            except (OSError, RuntimeError, subprocess.TimeoutExpired) as exc:
                parser.exit(1, f"{exc}\n")
            return
        print(f"Local URL: http://127.0.0.1:{port}", flush=True)
        print(f"Results folder: {server.root}\nNotes database: {server.notes}", flush=True)
        if server.access_code:
            addresses = socket.getaddrinfo(socket.gethostname(), None, socket.AF_INET)
            for address in sorted({result[4][0] for result in addresses} - {"127.0.0.1"}):
                print(f"Network URL: http://{address}:{port}", flush=True)
            print(f"Access code: {server.access_code}", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()
