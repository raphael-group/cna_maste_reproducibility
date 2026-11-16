"""The analysis imports no cnamaste: cnamaste only runs the method, and what the analysis takes from it is ported in
`scripts/cnamaste_port/`.

Referee: every import statement of the analysis' modules, the port and the tests, read from their source.
"""

from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ("scripts/analysis", "scripts/metrics", "scripts/plots", "scripts/cnamaste_port", "tests")
"""The analysis' code; `scripts/baselines/` runs inside the methods' environments and is left out."""


def get_cnamaste_imports(path: Path) -> list[str]:
    """Get the imports of `cnamaste` in one file.

    Args:
        path: a Python file.

    Returns:
        Each import of `cnamaste` or a submodule, as `<line>: <module>`.
    """
    found = []
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Import):
            modules = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            modules = [node.module]
        else:
            continue
        found += [f"{node.lineno}: {m}" for m in modules if m.split(".")[0] == "cnamaste"]
    return found


def test_analysis_imports_no_cnamaste() -> None:
    """No module of the analysis, the port or the tests imports `cnamaste`."""
    found = {
        str(path.relative_to(ROOT)): imports
        for directory in ANALYSIS
        for path in sorted((ROOT / directory).rglob("*.py"))
        if (imports := get_cnamaste_imports(path))
    }
    assert not found, f"cnamaste imported by the analysis: {found}"
