"""Publish a LaTeX results write-up as pages on the TUPIL figure review website.

Each \\section{...} in the .tex file becomes one page in the website's
"Results" section. Text is converted to the website's Markdown, and every
\\includegraphics PNG is copied from Results_Main into the Drive folder that
the website serves, then inserted as a figure block with its \\caption.

Rerunning the script updates pages with the same title instead of creating
duplicates. Pages in Trash are ignored, so a trashed page is recreated.

    py .\\Website\\Scripts\\figure_format\\publish_results.py            # publish
    py .\\Website\\Scripts\\figure_format\\publish_results.py --dry-run  # preview only
"""

import argparse
import hashlib
import json
import os
import re
import secrets
import shutil
import sqlite3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "figure_review"))
import serve  # noqa: E402  (website server: provides DRIVE_FOLDER and DB helpers)

DEFAULT_TEX = HERE / "results_pages.tex"
DEFAULT_SOURCE = HERE.parents[2] / "Results_Main"
DEFAULT_NOTES = Path(os.environ.get("LOCALAPPDATA", Path.home())) / "TUPIL" / "FigureReview" / "notes.sqlite3"
SECTION = "results"

MATH_SYMBOLS = {
    r"\kappa": "κ", r"\rho": "ρ", r"\alpha": "α", r"\beta": "β", r"\sigma": "σ", r"\mu": "μ",
    r"\Delta": "Δ", r"\geq": "≥", r"\ge": "≥", r"\leq": "≤", r"\le": "≤", r"\times": "×",
    r"\pm": "±", r"\approx": "≈", r"\to": "→", r"\cdot": "·", r"\neq": "≠",
}
SUPERSCRIPT = str.maketrans("0123456789-+()n", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺⁽⁾ⁿ")
NL = "\x01"  # newline that survives paragraph re-flowing


# ---------------------------------------------------------------- LaTeX parsing
def braced(text, start):
    """Return (content, end) for the {...} group that starts at text[start]."""
    if start >= len(text) or text[start] != "{":
        raise ValueError(f"Expected '{{' near: {text[start:start + 40]!r}")
    depth = 0
    for index in range(start, len(text)):
        if text[index] == "{" and text[index - 1] != "\\":
            depth += 1
        elif text[index] == "}" and text[index - 1] != "\\":
            depth -= 1
            if depth == 0:
                return text[start + 1:index], index + 1
    raise ValueError(f"Unbalanced braces near: {text[start:start + 40]!r}")


def replace_command(text, name, render):
    """Replace every \\name{arg} (optionally \\name[opt]{arg}) with render(arg)."""
    pattern = re.compile(r"\\" + name + r"(?![A-Za-z])\s*(\[[^\]]*\])?\s*")
    out, position = [], 0
    while match := pattern.search(text, position):
        arg, end = braced(text, match.end())
        out += [text[position:match.start()], render(arg)]
        position = end
    return "".join(out + [text[position:]])


def math_to_text(math):
    def superscript(match):
        value = match.group(1) or match.group(2)
        return value.translate(SUPERSCRIPT) if all(c in "0123456789-+()n" for c in value) else "^" + value

    for macro in sorted(MATH_SYMBOLS, key=len, reverse=True):
        math = re.sub(re.escape(macro) + r"(?![A-Za-z])", MATH_SYMBOLS[macro], math)
    math = re.sub(r"\^\{([^{}]*)\}|\^(\w)", superscript, math)
    math = math.replace("-", "−").replace("{", "").replace("}", "")
    return re.sub(r"\s+", " ", math)


def inline(text):
    text = re.sub(r"(?<!\\)\$(.+?)(?<!\\)\$", lambda m: math_to_text(m.group(1)), text, flags=re.S)
    text = replace_command(text, "textbf", lambda s: f"**{inline(s).strip()}**")
    for name in ("emph", "textit"):
        text = replace_command(text, name, lambda s: f"*{inline(s).strip()}*")
    text = replace_command(text, "label", lambda s: "")
    text = re.sub(r"\\centering\b", "", text)
    text = text.replace("---", "—").replace("--", "–").replace("``", "“").replace("''", "”")
    text = re.sub(r"\\([%&_#$])", r"\1", text)
    return text.replace("~", " ").replace(r"\\", NL)


def convert_lists(text):
    def render(match):
        bullet = "1. " if match.group(1) == "enumerate" else "- "
        items = [re.sub(r"\s+", " ", item).strip() for item in re.split(r"\\item\b", match.group(2))[1:]]
        return "\n\n" + NL.join(bullet + item for item in items if item) + "\n\n"

    return re.sub(r"\\begin\{(itemize|enumerate)\}(.*?)\\end\{\1\}", render, text, flags=re.S)


def to_markdown(text):
    text = convert_lists(inline(text))
    text = replace_command(text, "subsubsection", lambda s: f"\n\n### {s.strip()}\n\n")
    text = replace_command(text, "subsection", lambda s: f"\n\n## {s.strip()}\n\n")
    text = replace_command(text, "paragraph", lambda s: f"\n\n**{s.strip()}**{NL}")
    leftover = re.findall(r"\\[A-Za-z]+", text)
    if leftover:
        raise ValueError(f"Unsupported LaTeX command(s) for the website: {sorted(set(leftover))}")
    paragraphs = [re.sub(r"\s*\n\s*", " ", p).strip() for p in re.split(r"\n\s*\n", text)]
    paragraphs = [p.replace(NL + " ", NL).replace(NL, "\n").strip() for p in paragraphs]
    return "\n\n".join(p for p in paragraphs if p)


def parse_pages(tex):
    tex = re.sub(r"(?<!\\)%.*", "", tex)
    body = re.search(r"\\begin\{document\}(.*?)\\end\{document\}", tex, re.S)
    body = (body.group(1) if body else tex).replace(r"\maketitle", "")
    starts = [m for m in re.finditer(r"\\section\*?\s*(?=\{)", body)]
    pages = []
    for number, match in enumerate(starts):
        title, content_start = braced(body, match.end())
        content = body[content_start:starts[number + 1].start() if number + 1 < len(starts) else len(body)]
        blocks, position = [], 0
        for figure in re.finditer(r"\\begin\{figure\*?\}(.*?)\\end\{figure\*?\}", content, re.S):
            blocks.append({"kind": "text", "markdown": to_markdown(content[position:figure.start()])})
            graphic = re.search(r"\\includegraphics\s*(\[[^\]]*\])?\s*(?=\{)", figure.group(1))
            if not graphic:
                raise ValueError(f"A figure in '{title}' has no \\includegraphics.")
            path, _ = braced(figure.group(1), graphic.end())
            caption = re.search(r"\\caption\s*(?=\{)", figure.group(1))
            caption = to_markdown(braced(figure.group(1), caption.end())[0]) if caption else ""
            blocks.append({"kind": "figure", "figure": path.strip().replace("\\", "/"), "layout": "full",
                           "markdown": caption})
            position = figure.end()
        blocks.append({"kind": "text", "markdown": to_markdown(content[position:])})
        blocks = [b for b in blocks if b["kind"] == "figure" or b["markdown"]]
        pages.append({"title": re.sub(r"\s+", " ", inline(title)).strip(), "blocks": blocks})
    if not pages:
        raise ValueError("No \\section{...} found; each section becomes one website page.")
    return pages


def validate(pages, source):
    for page in pages:
        serve.Handler.clean_title(page["title"])
        if len(page["blocks"]) > 200:
            raise ValueError(f"'{page['title']}' has more than 200 blocks.")
        for index, block in enumerate(page["blocks"]):
            if len(block["markdown"]) > 20_000:
                raise ValueError(f"A text block in '{page['title']}' exceeds 20,000 characters.")
            if block["kind"] == "figure":
                relative = Path(block["figure"])
                if relative.suffix.lower() != ".png":
                    raise ValueError(f"Only .png figures are uploaded: {block['figure']}")
                if relative.anchor or ".." in relative.parts:
                    raise ValueError(f"Figure paths must be relative to Results_Main: {block['figure']}")
                if not (source / relative).is_file():
                    raise FileNotFoundError(f"Figure not found in {source}: {block['figure']}")
            digest = hashlib.md5(f"{page['title']}\0{index}\0{block['kind']}".encode()).hexdigest()
            block["id"] = f"block-{digest}"


# ---------------------------------------------------------------- publishing
def upload_figures(pages, source, dest):
    copied = skipped = 0
    for figure in sorted({b["figure"] for p in pages for b in p["blocks"] if b["kind"] == "figure"}):
        src, dst = source / figure, dest / figure
        if dst.is_file() and dst.stat().st_size == src.stat().st_size and dst.stat().st_mtime >= src.stat().st_mtime:
            skipped += 1
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        copied += 1
    return copied, skipped


def publish_pages(pages, notes):
    if not notes.is_file():
        raise FileNotFoundError(f"Website database not found: {notes}. Start the website once first.")
    results = []
    with sqlite3.connect(notes, factory=serve.ClosingConnection) as db:
        db.execute("BEGIN IMMEDIATE")
        for page in pages:
            serialized = json.dumps(page["blocks"], ensure_ascii=True)
            existing = db.execute("""SELECT n.page FROM notebook_navigation AS n
                JOIN notebook_titles AS t ON t.page = n.page
                WHERE n.section = ? AND t.title = ?
                  AND NOT EXISTS (SELECT 1 FROM notebook_trash AS x WHERE x.page = n.page)
                ORDER BY n.position LIMIT 1""", (SECTION, page["title"])).fetchone()
            if existing:
                page_id = existing[0]
                updated = db.execute("UPDATE notebook_pages SET revision = revision + 1, blocks = ? WHERE page = ?",
                                     (serialized, page_id))
                if not updated.rowcount:
                    db.execute("INSERT INTO notebook_pages VALUES (?, 1, ?)", (page_id, serialized))
                results.append(("updated", page["title"]))
            else:
                page_id = f"page-{secrets.token_hex(12)}"
                position = db.execute("SELECT COALESCE(MAX(position) + 1, 0) FROM notebook_navigation WHERE section = ?",
                                      (SECTION,)).fetchone()[0]
                db.execute("INSERT INTO notebook_navigation VALUES (?, ?, ?)", (page_id, SECTION, position))
                db.execute("INSERT INTO notebook_titles VALUES (?, ?)", (page_id, page["title"]))
                db.execute("INSERT INTO notebook_pages VALUES (?, 1, ?)", (page_id, serialized))
                results.append(("created", page["title"]))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tex", type=Path, default=DEFAULT_TEX, help="LaTeX write-up to publish.")
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE, help="Results_Main folder the figures come from.")
    parser.add_argument("--dest", type=Path, default=serve.DRIVE_FOLDER, help="Results_Main folder the website serves.")
    parser.add_argument("--notes", type=Path, default=DEFAULT_NOTES, help="Website SQLite database.")
    parser.add_argument("--dry-run", action="store_true", help="Print the converted pages without changing anything.")
    args = parser.parse_args()

    pages = parse_pages(args.tex.read_text(encoding="utf-8"))
    validate(pages, args.source)
    if args.dry_run:
        for page in pages:
            print(f"\n=== PAGE: {page['title']}")
            for block in page["blocks"]:
                if block["kind"] == "figure":
                    print(f"[FIGURE {block['figure']}]\n  caption: {block['markdown']}")
                else:
                    print(block["markdown"])
        return
    if not args.dest.is_dir():
        parser.exit(1, f"Website figure folder is unavailable: {args.dest}\n")
    copied, skipped = upload_figures(pages, args.source, args.dest)
    print(f"Figures: {copied} copied, {skipped} already up to date -> {args.dest}")
    for action, title in publish_pages(pages, args.notes):
        print(f"Page {action}: {title}")
    print("Reload the website to see the pages.")


if __name__ == "__main__":
    main()
