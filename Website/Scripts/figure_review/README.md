# TUPIL Figure Review

The workspace reads `Results_Main` and its subfolders directly from the Google Drive for desktop mount configured in this repository. It displays PNG, JPEG, WebP, GIF, AVIF, SVG, BMP, and PDF files. Click **Refresh** after adding or moving figures in Drive.

## Open on this computer

In PowerShell, from the TUPIL repository root, run the script with Python 3.10 or newer. It uses Python's standard library.

```powershell
py .\Website\Scripts\figure_review\serve.py
```

Open [http://127.0.0.1:8765](http://127.0.0.1:8765). Keep the PowerShell window running while you use the viewer. Press `Ctrl+C` in that window to stop it. Use `--port 8766` to choose a different port or `--results 'D:\path\to\Results_Main'` if Google Drive is mounted elsewhere.

## Write a results or project page

- **Results** contains **Results 1** and **Results 2**. **Research Notebook** contains **Project 1** and **Project 2**. Click a page in the navigation to open it.
- Drag a figure from **Drive figures** onto a dashed line in the document, or click **Insert**. Drag the figure's **Drag to move** handle to reposition it. **Move up** and **Move down** work with a keyboard or touch screen.
- Use **Add text** to write above, between, or below figures. Choose **Figure left** or **Figure right** to write beside a figure. **Text below** stacks a figure and its accompanying text.
- Write Markdown in any text field and use **Preview Markdown** to see the formatted page. The preview supports headings, bold, italic, lists, quotes, code, links, and horizontal rules.
- Pages save automatically to the host's SQLite database. A figure block links to the file in Drive, so moving or renaming that file can break the link. If two reviewers edit the same page at once, the second save displays a conflict so they can copy their changes before reloading.

## Place and edit figure annotations

Click **Figure annotations** to open the figure pin viewer.

1. Select a figure and click **Add note**.
2. Click on the figure or in the surrounding margin, type your note, and click **Save note**.
3. Click a pin to highlight its note. Drag it to move it, or use **Move** in the note list. Use **Edit** to change its text.

The host stores notes in `%LOCALAPPDATA%\TUPIL\FigureReview\notes.sqlite3` by default. Notes survive a restart and are shared by everyone connected to the same running host. Notes are associated with a figure's path inside `Results_Main`. Renaming or moving a figure changes the path the viewer uses to find its notes. PDF pins refer to the preview area. Use **Open original** to read a PDF separately.

## Invite someone on the same network

Stop the local session first, or specify another port when starting the shared session.

```powershell
py .\Website\Scripts\figure_review\serve.py --share
```

Give each reviewer the printed **Network URL** for the network you share and the **Access code**. The code changes each time the server starts. Reviewers can open figures, edit notebook pages, and add, edit, or delete notes. The host reads the Google Drive files on their behalf, so only share the session with people who should be able to see every figure in `Results_Main`.

Shared mode uses HTTP on your local network. Use it on a trusted network and keep the host computer on. If the network URL cannot be reached, check that reviewers are on the same network and that your firewall or network allows the selected port. You can also send the `figure_review` folder to anyone who wants to run their own copy. That copy needs access to the Drive files and has a separate notes database.

## Invite someone outside your network

Install the [Pinggy CLI](https://pinggy.io/cli/) on the host. The launcher looks on `PATH` and for `%LOCALAPPDATA%\TUPIL\FigureReview\bin\pinggy.exe`. You can also pass `--tunnel-executable 'C:\path\to\pinggy.exe'`.

```powershell
py .\Website\Scripts\figure_review\serve.py --remote
```

Remote mode picks an available local port automatically. Starting a new remote session leaves any existing invites active until they expire or are stopped.

Copy the full **Private invite link**, including the `#invite=…` part. On a reviewer's first visit, Pinggy may show a [screening page](https://pinggy.io/docs/http_tunnels/screening/). They can click **Enter** and reopen the full link if the workspace asks for a code. Anyone with the complete link can open the figures and edit workspace notes and notebook pages. Reviewers need a browser and an internet connection. The host still needs the workspace files and its Google Drive connection. The launcher uses the same local database as the local viewer.

Keep the host awake and the launcher running. A [free Pinggy tunnel](https://pinggy.io/) lasts up to one hour. Press `Ctrl+C` to close it early. Restarting the launcher generates a new URL and invite key.

`--provider cloudflare` uses a [Cloudflare Quick Tunnel](https://developers.cloudflare.com/tunnel/get-started/quick-tunnels/) instead. Install `cloudflared` on the host and ensure its network allows outbound connections on [port 7844](https://developers.cloudflare.com/cloudflare-one/networks/connectors/cloudflare-tunnel/configure-tunnels/tunnel-with-firewall/).

## Keep the same private invite link

Install the [ngrok CLI for Windows](https://ngrok.com/download/windows), sign in to ngrok, and follow the [dashboard instructions](https://dashboard.ngrok.com/get-started/your-authtoken) to add your agent authtoken on the host.

```powershell
py .\Website\Scripts\figure_review\serve.py --remote --provider ngrok
```

Share the full **Private invite link** printed by the script. The hostname comes from ngrok, and the workspace saves the invite key locally in the same database as your notebook, so the link stays the same when you run the command again with the same hostname and database. Reviewers may see a one-time **Visit Site** screen on ngrok's free plan. [ngrok dev domains](https://ngrok.com/docs/gateway/domains)

## Start or quit the background website

On the Windows host, double-click [`start.cmd`](./start.cmd) to start the private website and enable automatic startup at Windows sign-in. The window prints the private invite link and can be closed once the background site starts. Run [`quit.cmd`](./quit.cmd) to stop the site and disable automatic startup. Run **Start** again to re-enable it. The same hostname and database give you the same invite link each time.

The background session asks Windows to keep the computer awake while it runs. It requires the host to be powered on, online, and signed in. The task uses the signed-in user's configured ngrok account and reads figures from that user's Google Drive mount. Startup errors are recorded in `%LOCALAPPDATA%\TUPIL\FigureReview\background-site.log`.
