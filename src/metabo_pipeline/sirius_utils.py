from __future__ import annotations

import os
import platform
import shutil
from pathlib import Path
from typing import List, Optional, Tuple

def _collect_default_candidates(system: str, exe_name: str) -> List[Tuple[str, str]]:
    candidates: List[Tuple[str, str]] = []
    # Shared environment hints
    sirius_home = os.environ.get("SIRIUS_HOME")
    if sirius_home:
        home_path = Path(sirius_home).expanduser()
        exe = home_path / "bin" / exe_name
        candidates.append((str(exe), "env:SIRIUS_HOME"))
    home_bin = Path.home() / "sirius" / "bin" / exe_name
    candidates.append((str(home_bin), "home:sirius/bin"))
    if system == "Windows":
        program_files = [os.environ.get("ProgramFiles"), os.environ.get("ProgramFiles(x86)"), os.environ.get("LOCALAPPDATA")]
        for root in program_files:
            if not root:
                continue
            root_path = Path(root)
            for parts in (
                ("SIRIUS", exe_name),
                ("sirius", exe_name),
                ("SIRIUS", "bin", exe_name),
                ("sirius", "bin", exe_name),
                ("Programs", "sirius", exe_name),
                ("Programs", "Sirius", exe_name),
            ):
                candidates.append((str(root_path.joinpath(*parts)), f"default:{root_path.name}/{'/'.join(parts)}"))
    elif system == "Darwin":
        apps = [Path("/Applications"), Path.home() / "Applications"]
        for base in apps:
            for bundle in ("sirius.app", "SIRIUS.app"):
                app_path = base / bundle / "Contents" / "MacOS" / "sirius"
                candidates.append((str(app_path), f"bundle:{bundle}"))
        brew_path = Path("/usr/local/bin") / exe_name
        candidates.append((str(brew_path), "usr-local"))
    else:
        for path in (
            Path("/usr/local/bin") / exe_name,
            Path("/usr/bin") / exe_name,
            Path("/opt/sirius/bin") / exe_name,
        ):
            candidates.append((str(path), "linux-default"))
    return candidates

def guess_sirius_executable(preferred: Optional[str] = None) -> Tuple[str, Optional[str]]:
    """Return a likely path to the SIRIUS executable and a note describing the source."""
    system = platform.system()
    exe_name = "sirius.exe" if system == "Windows" else "sirius"
    base = preferred or exe_name

    candidates: List[Tuple[str, str]] = []
    seen: set[str] = set()

    def add(path: Optional[str], note: str) -> None:
        if not path:
            return
        norm = str(Path(path).expanduser())
        if norm in seen:
            return
        seen.add(norm)
        candidates.append((norm, note))

    add(preferred, "cli-option")
    for env_var in ("SIRIUS_EXECUTABLE", "SIRIUS_EXE"):
        add(os.environ.get(env_var), f"env:{env_var}")
    add(shutil.which(base), "which")
    add(shutil.which(exe_name), "which-default")
    for path, note in _collect_default_candidates(system, exe_name):
        add(path, note)

    for cand, note in candidates:
        path = Path(cand)
        try:
            if path.is_file():
                return str(path), note
        except OSError:
            pass
        resolved = shutil.which(str(path))
        if resolved:
            return resolved, note

    fallback = str(Path(base).expanduser())
    return fallback, "fallback"
