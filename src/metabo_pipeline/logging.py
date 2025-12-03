from __future__ import annotations

import sys


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
                console.print(f"[info]{msg}")

            def warn(self, msg: str):
                console.print(f"[warn]{msg}")

            def error(self, msg: str):
                console.print(f"[error]{msg}")

            def ok(self, msg: str):
                console.print(f"[ok]{msg}")

        return _RichLogger()
    except Exception:

        class _PlainLogger:
            def info(self, msg: str):
                print(f"INFO: {msg}")

            def warn(self, msg: str):
                print(f"WARN: {msg}")

            def error(self, msg: str):
                print(f"ERROR: {msg}", file=sys.stderr)

            def ok(self, msg: str):
                print(f"OK: {msg}")

        return _PlainLogger()
