from __future__ import annotations

import sys


def _safe_text(msg: str, stream) -> str:
    """Keep piped/legacy Windows consoles from failing on status symbols."""
    encoding = getattr(stream, "encoding", None) or "utf-8"
    return str(msg).encode(encoding, errors="replace").decode(encoding)


def get_logger():
    """Return a lightweight console logger with optional rich styling."""
    try:
        from rich.console import Console
        from rich.theme import Theme

        theme = Theme(
            {
                "info": "cyan",
                "warn": "yellow",
                "error": "bold red",
                "ok": "green",
            }
        )
        console = Console(theme=theme)

        class _RichLogger:
            def info(self, msg: str):
                console.print(_safe_text(msg, sys.stdout), style="info", markup=False)

            def warn(self, msg: str):
                console.print(_safe_text(msg, sys.stdout), style="warn", markup=False)

            def error(self, msg: str):
                console.print(_safe_text(msg, sys.stdout), style="error", markup=False)

            def ok(self, msg: str):
                console.print(_safe_text(msg, sys.stdout), style="ok", markup=False)

        return _RichLogger()
    except Exception:

        class _PlainLogger:
            def info(self, msg: str):
                print(_safe_text(f"INFO: {msg}", sys.stdout))

            def warn(self, msg: str):
                print(_safe_text(f"WARN: {msg}", sys.stdout))

            def error(self, msg: str):
                print(_safe_text(f"ERROR: {msg}", sys.stderr), file=sys.stderr)

            def ok(self, msg: str):
                print(_safe_text(f"OK: {msg}", sys.stdout))

        return _PlainLogger()
