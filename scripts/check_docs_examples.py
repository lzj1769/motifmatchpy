"""Run documented CLI and Python examples against the installed package.

Run from any directory with the documentation environment's Python. Examples
execute in a temporary directory so output files never overwrite user data.
"""

from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FENCES = re.compile(r"^```(sh|python|text)\n(.*?)^```", re.MULTILINE | re.DOTALL)


def main() -> None:
    commands = 0
    python_blocks = 0
    previous_cwd = Path.cwd()
    with tempfile.TemporaryDirectory(prefix="motifmatchpy-docs-") as directory:
        work = Path(directory)
        shutil.copytree(ROOT / "docs/examples", work / "docs/examples")
        try:
            os.chdir(work)
            for page in sorted((ROOT / "docs").rglob("*.md")):
                blocks = list(FENCES.finditer(page.read_text()))
                for index, block in enumerate(blocks):
                    language, code = block.groups()
                    label = f"{page.relative_to(ROOT)} block {index + 1}"
                    if language == "python":
                        exec(compile(code, label, "exec"), {"__name__": "__main__"})
                        python_blocks += 1
                    elif language == "sh":
                        last_result = None
                        for line in code.replace("\\\n", " ").splitlines():
                            argv = shlex.split(line, comments=True)
                            if not argv or argv[0] != "motifmatchpy":
                                continue
                            last_result = subprocess.run(
                                [sys.executable, "-m", "motifmatchpy.cli", *argv[1:]],
                                check=True, capture_output=True, text=True,
                            )
                            commands += 1
                        # Compare published stdout examples, including CSV headers.
                        if last_result is not None and index + 1 < len(blocks):
                            following = blocks[index + 1]
                            if following.group(1) == "text" and "-o" not in argv:
                                expected = following.group(2).strip()
                                actual = last_result.stdout.strip()
                                if actual != expected:
                                    raise AssertionError(
                                        f"{label}: expected {expected!r}, got {actual!r}"
                                    )
                        if last_result is not None and "-o" in argv:
                            output = work / argv[argv.index("-o") + 1]
                            if not output.is_file() or not output.read_text().strip():
                                raise AssertionError(f"{label}: missing or empty output {output}")
        finally:
            os.chdir(previous_cwd)
    print(f"Documentation examples passed: {commands} CLI commands, {python_blocks} Python blocks")


if __name__ == "__main__":
    main()
