"""The stub's docstrings are the runtime ones (#474).

The docstrings are written twice by hand: as Rust doc comments, which
become each member's runtime `__doc__`, and in `geomcore/__init__.pyi`,
which IDEs show. stubtest compares the signatures only, so this compares
the text. A public class of the stub is checked with the members it
inherits from its stub-only bases (the `_*Methods` method sets).
"""

import ast
import inspect
from pathlib import Path

import geomcore

STUB = Path(__file__).resolve().parents[1] / "geomcore" / "__init__.pyi"


def _base_name(node: ast.expr) -> str | None:
    """The class a base expression names: `Name` or `Name[...]`."""
    if isinstance(node, ast.Subscript):
        node = node.value
    return node.id if isinstance(node, ast.Name) else None


def _classes() -> dict[str, ast.ClassDef]:
    tree = ast.parse(STUB.read_text())
    return {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}


def _member_docstrings(name: str, classes: dict[str, ast.ClassDef]) -> dict[str, str]:
    """Each member's stub docstring, own members first, then those of the
    stub bases in order, as Python's attribute lookup would find them."""
    found: dict[str, str] = {}
    node = classes[name]
    for item in node.body:
        if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
            docstring = ast.get_docstring(item)
            if docstring is not None and item.name not in found:
                found[item.name] = docstring
    for base in node.bases:
        base_name = _base_name(base)
        if base_name in classes:
            for member, docstring in _member_docstrings(base_name, classes).items():
                found.setdefault(member, docstring)
    return found


def _runtime_docstring(value: object) -> str | None:
    docstring = getattr(value, "__doc__", None)
    return inspect.cleandoc(docstring) if docstring else None


def test_stub_docstrings_match_the_runtime_ones():
    classes = _classes()
    public = [name for name in classes if not name.startswith("_")]
    assert len(public) > 50, f"only {len(public)} public classes parsed"
    checked = 0
    drift = []
    for name in public:
        runtime = getattr(geomcore, name, None)
        if runtime is None:
            continue
        class_docstring = ast.get_docstring(classes[name])
        if class_docstring is not None:
            checked += 1
            if _runtime_docstring(runtime) != class_docstring:
                drift.append(name)
        for member, docstring in _member_docstrings(name, classes).items():
            value = inspect.getattr_static(runtime, member, None)
            if value is None:
                continue
            if isinstance(value, (staticmethod, classmethod)):
                value = value.__func__
            checked += 1
            if _runtime_docstring(value) != docstring:
                drift.append(f"{name}.{member}")
    assert checked > 500, f"only {checked} docstrings compared"
    assert drift == [], "stub docstrings that differ from the runtime ones:\n" + "\n".join(drift)
