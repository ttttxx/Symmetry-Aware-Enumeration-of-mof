#!/usr/bin/env python3
"""Tiny local LaTeX live preview server for PyCharm split-view workflows."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import mimetypes
import queue
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import webbrowser
from dataclasses import dataclass, field
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Iterable
from urllib.parse import urlparse


WATCH_EXTENSIONS = {
    ".tex",
    ".bib",
    ".bst",
    ".cls",
    ".sty",
    ".cfg",
    ".def",
    ".png",
    ".jpg",
    ".jpeg",
    ".pdf",
    ".eps",
}


@dataclass
class PreviewState:
    tex_path: Path
    build_dir: Path
    command: list[str]
    lock: threading.Lock = field(default_factory=threading.Lock)
    event_lock: threading.Lock = field(default_factory=threading.Lock)
    events: list[queue.Queue[str]] = field(default_factory=list)
    version: int = 0
    compiling: bool = False
    pending_compile: bool = False
    success: bool = False
    returncode: int | None = None
    pdf_path: Path | None = None
    log: str = "Waiting for the first compile."
    status: str = "idle"
    started_at: float | None = None
    finished_at: float | None = None

    def snapshot(self) -> dict[str, object]:
        with self.lock:
            return {
                "version": self.version,
                "compiling": self.compiling,
                "success": self.success,
                "returncode": self.returncode,
                "status": self.status,
                "pdfReady": bool(self.pdf_path and self.pdf_path.exists()),
                "texPath": str(self.tex_path),
                "buildDir": str(self.build_dir),
                "command": " ".join(self.command),
                "startedAt": self.started_at,
                "finishedAt": self.finished_at,
            }

    def publish(self) -> None:
        payload = json.dumps(self.snapshot(), ensure_ascii=True)
        with self.event_lock:
            dead: list[queue.Queue[str]] = []
            for event_queue in self.events:
                try:
                    event_queue.put_nowait(payload)
                except queue.Full:
                    dead.append(event_queue)
            for event_queue in dead:
                self.events.remove(event_queue)


def project_default_build_dir(tex_path: Path) -> Path:
    digest = hashlib.sha1(str(tex_path.resolve()).encode("utf-8")).hexdigest()[:10]
    root = Path(tempfile.gettempdir()) / "latex-live-preview"
    return root / f"{tex_path.stem}-{digest}"


def iter_watch_files(root: Path, build_dir: Path) -> Iterable[Path]:
    build_dir = build_dir.resolve()
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        try:
            resolved = path.resolve()
        except OSError:
            continue
        if resolved == build_dir or build_dir in resolved.parents:
            continue
        if path.suffix.lower() in WATCH_EXTENSIONS:
            yield path


def snapshot_files(root: Path, build_dir: Path) -> tuple[tuple[str, int, int], ...]:
    items: list[tuple[str, int, int]] = []
    for path in iter_watch_files(root, build_dir):
        try:
            stat = path.stat()
        except OSError:
            continue
        items.append((str(path), stat.st_mtime_ns, stat.st_size))
    return tuple(sorted(items))


def trim_log(text: str, max_chars: int = 24000) -> str:
    if len(text) <= max_chars:
        return text
    return text[-max_chars:]


def run_compile(state: PreviewState) -> None:
    while True:
        with state.lock:
            if state.compiling:
                state.pending_compile = True
                state.status = "compile queued"
                queued = True
            else:
                queued = False
        if queued:
            state.publish()
            return

        with state.lock:
            state.compiling = True
            state.pending_compile = False
            state.started_at = time.time()
            state.finished_at = None
            state.status = "compiling"
            state.log = "Compiling..."
            state.version += 1
        state.publish()

        state.build_dir.mkdir(parents=True, exist_ok=True)
        command = state.command
        try:
            completed = subprocess.run(
                command,
                cwd=state.tex_path.parent,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            output = completed.stdout or ""
            expected_pdf = state.build_dir / f"{state.tex_path.stem}.pdf"
            success = completed.returncode == 0 and expected_pdf.exists()
            with state.lock:
                state.success = success
                state.returncode = completed.returncode
                state.pdf_path = expected_pdf if expected_pdf.exists() else state.pdf_path
                state.finished_at = time.time()
                state.status = "compiled" if success else "compile failed"
                state.log = trim_log(output)
                state.compiling = False
                state.version += 1
                should_repeat = state.pending_compile
                state.pending_compile = False
        except Exception as exc:  # noqa: BLE001 - show local tool failures in the UI.
            with state.lock:
                state.success = False
                state.returncode = None
                state.finished_at = time.time()
                state.status = "preview tool error"
                state.log = trim_log(f"{type(exc).__name__}: {exc}")
                state.compiling = False
                state.version += 1
                should_repeat = state.pending_compile
                state.pending_compile = False

        state.publish()
        if not should_repeat:
            return


def trigger_compile(state: PreviewState) -> None:
    thread = threading.Thread(target=run_compile, args=(state,), daemon=True)
    thread.start()


def watch_for_changes(state: PreviewState, interval: float) -> None:
    root = state.tex_path.parent
    previous = snapshot_files(root, state.build_dir)
    while True:
        time.sleep(interval)
        current = snapshot_files(root, state.build_dir)
        if current != previous:
            previous = current
            trigger_compile(state)


def response_bytes(text: str) -> bytes:
    return text.encode("utf-8")


def render_index(state: PreviewState) -> str:
    title = f"LaTeX Preview - {state.tex_path.name}"
    escaped_title = html.escape(title)
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escaped_title}</title>
  <style>
    :root {{
      color-scheme: light dark;
      --bg: #f6f7f9;
      --panel: #ffffff;
      --text: #1f2933;
      --muted: #697386;
      --line: #d7dce3;
      --accent: #1f6feb;
      --ok: #167c3b;
      --bad: #b42318;
    }}
    @media (prefers-color-scheme: dark) {{
      :root {{
        --bg: #111418;
        --panel: #181d23;
        --text: #eef2f6;
        --muted: #a6b0bd;
        --line: #303844;
        --accent: #6aa6ff;
        --ok: #4cc38a;
        --bad: #ff8a80;
      }}
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      min-height: 100vh;
      background: var(--bg);
      color: var(--text);
      font: 13px/1.45 -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }}
    .shell {{
      height: 100vh;
      display: grid;
      grid-template-rows: auto minmax(0, 1fr);
    }}
    .toolbar {{
      display: flex;
      align-items: center;
      gap: 10px;
      padding: 9px 10px;
      border-bottom: 1px solid var(--line);
      background: var(--panel);
    }}
    .title {{
      min-width: 0;
      flex: 1;
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
      font-weight: 600;
    }}
    .status {{
      min-width: 110px;
      text-align: center;
      color: var(--muted);
      white-space: nowrap;
    }}
    .status.ok {{ color: var(--ok); }}
    .status.bad {{ color: var(--bad); }}
    button, a.button {{
      appearance: none;
      border: 1px solid var(--line);
      border-radius: 6px;
      background: transparent;
      color: var(--text);
      padding: 5px 9px;
      text-decoration: none;
      cursor: pointer;
      font: inherit;
      white-space: nowrap;
    }}
    button:hover, a.button:hover {{ border-color: var(--accent); color: var(--accent); }}
    .main {{
      min-height: 0;
      display: grid;
      grid-template-rows: minmax(0, 1fr) auto;
    }}
    iframe {{
      width: 100%;
      height: 100%;
      border: 0;
      background: #707780;
    }}
    .log {{
      display: none;
      max-height: 34vh;
      overflow: auto;
      padding: 10px 12px;
      border-top: 1px solid var(--line);
      background: var(--panel);
      white-space: pre-wrap;
      font: 12px/1.45 ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
    }}
    .log.visible {{ display: block; }}
  </style>
</head>
<body>
  <div class="shell">
    <header class="toolbar">
      <div class="title" title="{html.escape(str(state.tex_path))}">{html.escape(state.tex_path.name)}</div>
      <div id="status" class="status">connecting</div>
      <button id="compile" title="Run latexmk now">Compile</button>
      <button id="toggle-log" title="Show or hide the latest latexmk log">Log</button>
      <a id="open-pdf" class="button" href="/pdf" target="_blank" title="Open the current PDF alone">PDF</a>
    </header>
    <main class="main">
      <iframe id="viewer" title="PDF preview"></iframe>
      <pre id="log" class="log"></pre>
    </main>
  </div>
  <script>
    const statusEl = document.getElementById('status');
    const viewer = document.getElementById('viewer');
    const logEl = document.getElementById('log');
    const openPdf = document.getElementById('open-pdf');
    let currentVersion = -1;

    function statusClass(state) {{
      if (state.compiling) return '';
      if (state.success) return 'ok';
      if (state.returncode !== null) return 'bad';
      return '';
    }}

    function update(state) {{
      statusEl.textContent = state.status || 'idle';
      statusEl.className = 'status ' + statusClass(state);
      if (state.pdfReady && state.version !== currentVersion) {{
        currentVersion = state.version;
        const url = '/pdf?v=' + encodeURIComponent(String(state.version));
        viewer.src = url;
        openPdf.href = url;
      }}
      fetch('/log').then(r => r.text()).then(text => {{ logEl.textContent = text; }});
    }}

    document.getElementById('compile').addEventListener('click', () => {{
      fetch('/compile', {{ method: 'POST' }});
    }});

    document.getElementById('toggle-log').addEventListener('click', () => {{
      logEl.classList.toggle('visible');
    }});

    const events = new EventSource('/events');
    events.onmessage = event => update(JSON.parse(event.data));
    events.onerror = () => {{ statusEl.textContent = 'reconnecting'; statusEl.className = 'status bad'; }};
    fetch('/state').then(r => r.json()).then(update);
  </script>
</body>
</html>"""


class PreviewHandler(BaseHTTPRequestHandler):
    server: "PreviewHTTPServer"

    def log_message(self, format: str, *args: object) -> None:
        sys.stderr.write("[%s] %s\n" % (self.log_date_time_string(), format % args))

    def send_body(
        self,
        body: bytes,
        content_type: str,
        status: HTTPStatus = HTTPStatus.OK,
        extra_headers: dict[str, str] | None = None,
    ) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        if extra_headers:
            for name, value in extra_headers.items():
                self.send_header(name, value)
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:  # noqa: N802 - BaseHTTPRequestHandler API.
        parsed = urlparse(self.path)
        path = parsed.path
        state = self.server.preview_state

        if path == "/":
            self.send_body(response_bytes(render_index(state)), "text/html; charset=utf-8")
            return

        if path == "/state":
            self.send_body(
                response_bytes(json.dumps(state.snapshot(), ensure_ascii=True)),
                "application/json; charset=utf-8",
            )
            return

        if path == "/log":
            with state.lock:
                log = state.log
            self.send_body(response_bytes(log), "text/plain; charset=utf-8")
            return

        if path == "/events":
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", "text/event-stream")
            self.send_header("Cache-Control", "no-store")
            self.send_header("Connection", "keep-alive")
            self.end_headers()
            event_queue: queue.Queue[str] = queue.Queue(maxsize=8)
            with state.event_lock:
                state.events.append(event_queue)
            try:
                event_queue.put(json.dumps(state.snapshot(), ensure_ascii=True))
                while True:
                    try:
                        payload = event_queue.get(timeout=15)
                        self.wfile.write(f"data: {payload}\n\n".encode("utf-8"))
                    except queue.Empty:
                        self.wfile.write(b": heartbeat\n\n")
                    self.wfile.flush()
            except (BrokenPipeError, ConnectionResetError):
                pass
            finally:
                with state.event_lock:
                    if event_queue in state.events:
                        state.events.remove(event_queue)
            return

        if path == "/pdf":
            with state.lock:
                pdf_path = state.pdf_path
            if not pdf_path or not pdf_path.exists():
                self.send_body(
                    b"PDF is not ready yet. Check the compile log.",
                    "text/plain; charset=utf-8",
                    HTTPStatus.NOT_FOUND,
                )
                return
            data = pdf_path.read_bytes()
            content_type = mimetypes.guess_type(pdf_path.name)[0] or "application/pdf"
            self.send_body(
                data,
                content_type,
                extra_headers={"Content-Disposition": f'inline; filename="{pdf_path.name}"'},
            )
            return

        self.send_body(b"Not found", "text/plain; charset=utf-8", HTTPStatus.NOT_FOUND)

    def do_HEAD(self) -> None:  # noqa: N802 - BaseHTTPRequestHandler API.
        parsed = urlparse(self.path)
        path = parsed.path
        state = self.server.preview_state

        if path == "/":
            body = response_bytes(render_index(state))
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            return

        if path == "/pdf":
            with state.lock:
                pdf_path = state.pdf_path
            if not pdf_path or not pdf_path.exists():
                self.send_response(HTTPStatus.NOT_FOUND)
                self.send_header("Content-Type", "text/plain; charset=utf-8")
                self.send_header("Content-Length", "0")
                self.send_header("Cache-Control", "no-store")
                self.end_headers()
                return
            content_type = mimetypes.guess_type(pdf_path.name)[0] or "application/pdf"
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(pdf_path.stat().st_size))
            self.send_header("Content-Disposition", f'inline; filename="{pdf_path.name}"')
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            return

        self.send_response(HTTPStatus.NOT_FOUND)
        self.send_header("Content-Length", "0")
        self.send_header("Cache-Control", "no-store")
        self.end_headers()

    def do_POST(self) -> None:  # noqa: N802 - BaseHTTPRequestHandler API.
        if urlparse(self.path).path == "/compile":
            trigger_compile(self.server.preview_state)
            self.send_body(b"queued", "text/plain; charset=utf-8")
            return
        self.send_body(b"Not found", "text/plain; charset=utf-8", HTTPStatus.NOT_FOUND)


class PreviewHTTPServer(ThreadingHTTPServer):
    allow_reuse_address = True

    def __init__(self, server_address: tuple[str, int], state: PreviewState) -> None:
        super().__init__(server_address, PreviewHandler)
        self.preview_state = state


def build_latexmk_command(tex_path: Path, build_dir: Path, engine: str) -> list[str]:
    latexmk = shutil.which("latexmk")
    if not latexmk:
        raise SystemExit("latexmk was not found on PATH. Install MacTeX/TeX Live or add latexmk to PATH.")
    command = [
        latexmk,
        "-pdf",
        "-interaction=nonstopmode",
        "-halt-on-error",
        "-synctex=1",
        f"-outdir={build_dir}",
        tex_path.name,
    ]
    if engine in {"xelatex", "lualatex"}:
        command.insert(2, f"-{engine}")
    return command


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Serve an auto-refreshing PDF preview for a LaTeX file.",
    )
    parser.add_argument("tex", nargs="?", default="paper/main.tex", help="Main .tex file.")
    parser.add_argument("--host", default="127.0.0.1", help="Bind host. Default: 127.0.0.1")
    parser.add_argument("--port", type=int, default=8765, help="Bind port. Default: 8765")
    parser.add_argument("--interval", type=float, default=1.0, help="Watch interval in seconds.")
    parser.add_argument("--build-dir", default=None, help="Output directory. Default: system temp.")
    parser.add_argument(
        "--engine",
        choices=["pdflatex", "xelatex", "lualatex"],
        default="pdflatex",
        help="latexmk engine flag. Default: pdflatex",
    )
    parser.add_argument("--open", action="store_true", help="Open the preview URL in the default browser.")
    parser.add_argument("--no-initial-compile", action="store_true", help="Wait for a file change before compiling.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    tex_path = Path(args.tex).expanduser().resolve()
    if not tex_path.exists():
        print(f"TeX file does not exist: {tex_path}", file=sys.stderr)
        return 2
    if tex_path.suffix.lower() != ".tex":
        print(f"Expected a .tex file: {tex_path}", file=sys.stderr)
        return 2

    build_dir = Path(args.build_dir).expanduser().resolve() if args.build_dir else project_default_build_dir(tex_path)
    command = build_latexmk_command(tex_path, build_dir, args.engine)
    state = PreviewState(tex_path=tex_path, build_dir=build_dir, command=command)
    try:
        server = PreviewHTTPServer((args.host, args.port), state)
    except OSError as exc:
        print(f"Could not start preview server on {args.host}:{args.port}: {exc}", file=sys.stderr)
        print("Try another port, for example: --port 8766", file=sys.stderr)
        return 2

    watcher = threading.Thread(target=watch_for_changes, args=(state, args.interval), daemon=True)
    watcher.start()
    if not args.no_initial_compile:
        trigger_compile(state)

    url = f"http://{args.host}:{args.port}/"
    print(f"LaTeX live preview: {url}", flush=True)
    print(f"TeX source: {tex_path}", flush=True)
    print(f"Build dir:  {build_dir}", flush=True)
    print(
        "Open the URL in PyCharm's browser tab, split it to the right of your .tex editor, and edit normally.",
        flush=True,
    )
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping preview server.")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
