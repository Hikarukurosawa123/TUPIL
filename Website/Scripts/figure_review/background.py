"""Start, stop, and run TUPIL's per-user background website on Windows."""

import argparse
import ctypes
import json
import os
from argparse import Namespace
from http.client import HTTPConnection, HTTPException
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import threading
import time
from urllib.parse import urlsplit

import serve


TASK_NAME = "TUPIL Private Website"
PORT = 8766
ENTRY = Path(__file__).resolve()
STATE = Path(os.environ.get("LOCALAPPDATA", Path.home())) / "TUPIL" / "FigureReview"
STATUS = STATE / "background-site.json"
STOP = STATE / "background-site.stop"
LOG = STATE / "background-site.log"
NOTES = STATE / "notes.sqlite3"

TASK_HEADER = r"""
$ErrorActionPreference = 'Stop'
$name = $env:TUPIL_TASK_NAME
$task = Get-ScheduledTask -TaskName $name -TaskPath '\' -ErrorAction SilentlyContinue
if ($task -and ($task.Actions.Count -ne 1 -or
    $task.Actions[0].Execute -ne $env:TUPIL_PYTHONW -or
    $task.Actions[0].Arguments -ne $env:TUPIL_TASK_ARGUMENTS)) {
    throw 'A different task already uses the TUPIL Private Website name.'
}
"""


def task_command(script, timeout=55):
    pythonw = Path(sys.executable).with_name("pythonw.exe")
    if not pythonw.is_file():
        raise RuntimeError("Python for Windows with pythonw.exe is required.")
    env = {
        **os.environ,
        "TUPIL_TASK_NAME": TASK_NAME,
        "TUPIL_PYTHONW": str(pythonw),
        "TUPIL_TASK_ARGUMENTS": subprocess.list2cmdline(["-B", str(ENTRY), "run"]),
    }
    result = subprocess.run(["powershell.exe", "-NoProfile", "-NonInteractive", "-Command", TASK_HEADER + script],
                            env=env, capture_output=True, text=True, errors="replace", timeout=timeout)
    if result.returncode:
        raise RuntimeError((result.stderr or result.stdout).strip())
    return result.stdout.strip()


def active_invite():
    try:
        info = json.loads(STATUS.read_text(encoding="utf-8"))
        url = info["url"]
        host = urlsplit(url).hostname
        client = HTTPConnection("127.0.0.1", PORT, timeout=1)
        try:
            client.request("GET", "/", headers={"Host": host, "X-Forwarded-Proto": "https"})
            response = client.getresponse()
            if response.status != 200 or b"<title>TUPIL Research Workspace</title>" not in response.read():
                return None
        finally:
            client.close()
        with sqlite3.connect(NOTES, factory=serve.ClosingConnection) as db:
            key = db.execute("SELECT code FROM remote_invites WHERE hostname = ?", (host,)).fetchone()
        if key:
            return f"{url}/#invite={key[0]}"
    except (OSError, ValueError, KeyError, HTTPException, sqlite3.Error):
        return None
    return None


def start():
    if not shutil.which("ngrok"):
        raise RuntimeError("Install and authenticate ngrok on this computer before starting the background website.")
    invite = active_invite()
    if invite and task_command("if ($task -and $task.State -eq 'Running') { 'running' }") == "running":
        print(f"TUPIL is already running in the background.\nPrivate invite link: {invite}")
        return
    STOP.unlink(missing_ok=True)
    task_command(r"""
if (-not $task) {
    $account = [System.Security.Principal.WindowsIdentity]::GetCurrent().Name
    $action = New-ScheduledTaskAction -Execute $env:TUPIL_PYTHONW -Argument $env:TUPIL_TASK_ARGUMENTS -WorkingDirectory $env:LOCALAPPDATA
    $trigger = New-ScheduledTaskTrigger -AtLogOn -User $account
    $principal = New-ScheduledTaskPrincipal -UserId $account -LogonType Interactive -RunLevel Limited
    $settings = New-ScheduledTaskSettingsSet -StartWhenAvailable -ExecutionTimeLimit ([TimeSpan]::Zero) `
        -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries -RestartCount 3 `
        -RestartInterval (New-TimeSpan -Minutes 1) -MultipleInstances IgnoreNew
    Register-ScheduledTask -TaskName $name -Action $action -Trigger $trigger -Principal $principal -Settings $settings | Out-Null
} elseif ($task.State -eq 'Disabled') {
    Enable-ScheduledTask -InputObject $task | Out-Null
}
Start-ScheduledTask -TaskName $name
""")
    deadline = time.monotonic() + 51
    while time.monotonic() < deadline:
        invite = active_invite()
        if invite:
            print(f"TUPIL is running in the background and will start at Windows sign-in.\nPrivate invite link: {invite}")
            return
        time.sleep(1)
    raise RuntimeError(f"TUPIL did not finish starting. Check {LOG}")


def quit_website():
    if task_command("if ($task) { 'registered' }") != "registered":
        print("The TUPIL background task has not been set up.")
        return
    STOP.touch()
    task_command(r"""
Disable-ScheduledTask -InputObject $task | Out-Null
$deadline = (Get-Date).AddSeconds(49)
while (Get-CimInstance Win32_Process -Filter "Name='pythonw.exe'" | Where-Object {
    $_.CommandLine -and $_.CommandLine.EndsWith($env:TUPIL_TASK_ARGUMENTS, [StringComparison]::OrdinalIgnoreCase)
}) {
    if ((Get-Date) -ge $deadline) { throw 'The website did not stop in time.' }
    Start-Sleep -Milliseconds 500
}
""")
    if active_invite():
        raise RuntimeError("The website is still responding and may be shutting down. Run Quit again to check.")
    print("TUPIL background website stopped. Windows sign-in startup is disabled.")


def run_background():
    STATUS.unlink(missing_ok=True)
    if STOP.exists():
        return
    with LOG.open("w", encoding="utf-8", buffering=1) as log:
        original = sys.stdout, sys.stderr
        sys.stdout = sys.stderr = log
        try:
            run_server()
        except (OSError, RuntimeError, subprocess.TimeoutExpired) as exc:
            print(exc, file=log)
            raise SystemExit(1)
        finally:
            STATUS.unlink(missing_ok=True)
            sys.stdout, sys.stderr = original


def run_server():
    executable = shutil.which("ngrok")
    if not executable:
        raise RuntimeError("ngrok is missing from this user's PATH.")
    if not serve.DRIVE_FOLDER.is_dir():
        print("Results_Main is unavailable. The figure library will work when the Google Drive mount is ready.")
    args = Namespace(results=serve.DRIVE_FOLDER, notes=NOTES, port=PORT, remote=True, share=False)
    kernel = ctypes.windll.kernel32
    kernel.SetThreadExecutionState(0x80000001)  # ES_CONTINUOUS | ES_SYSTEM_REQUIRED
    try:
        with serve.ReviewServer(args) as server:
            def ready(url):
                temp = STATUS.with_suffix(".tmp")
                temp.write_text(json.dumps({"url": url}), encoding="utf-8")
                temp.replace(STATUS)

                def wait_for_stop():
                    while not STOP.exists():
                        time.sleep(.5)
                    server.shutdown()

                threading.Thread(target=wait_for_stop, daemon=True).start()

            try:
                serve.serve_ngrok(server, executable, None, on_ready=ready)
            except RuntimeError:
                if not STOP.exists():
                    raise
    finally:
        kernel.SetThreadExecutionState(0x80000000)  # ES_CONTINUOUS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("start", "quit", "run"))
    action = parser.parse_args().action
    if os.name != "nt":
        parser.error("The background launcher requires Windows.")
    STATE.mkdir(parents=True, exist_ok=True)
    try:
        {"start": start, "quit": quit_website, "run": run_background}[action]()
    except (OSError, RuntimeError, subprocess.TimeoutExpired) as exc:
        parser.exit(1, f"{exc}\n")


if __name__ == "__main__":
    main()
