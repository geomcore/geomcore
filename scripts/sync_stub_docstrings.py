#!/usr/bin/env python3
"""Copies the runtime docstrings of the Python bindings into their stub.

Editors that read only `geomcore-py/geomcore/__init__.pyi` show its
docstrings on hover, so each class, function, method and property of the
stub carries the docstring the built module gives it (#334). Run this
after changing a docstring on the Rust side, with the module built
(`maturin develop`); `test_stub_docstrings_are_the_runtime_ones` fails
until the stub matches.

An overloaded function carries its docstring on its first overload, and
the stub-only method sets (`_CurveMethods`, ...) keep the docstrings
written there, which the parity test compares with each class's.
"""

import ast
import inspect
import sys
from pathlib import Path

import geomcore

STUB = Path(__file__).resolve().parent.parent / "geomcore-py" / "geomcore" / "__init__.pyi"


def runtime_docstring(owner, name):
    """The cleaned docstring of `owner.name` as the module has it, or None."""
    member = inspect.getattr_static(owner, name, None)
    if member is None:
        return None
    if isinstance(member, (staticmethod, classmethod)):
        member = member.__func__
    docstring = getattr(member, "__doc__", None)
    if not isinstance(docstring, str) or not docstring or is_generic(name, member):
        return None
    return inspect.cleandoc(docstring)


def is_generic(name, member):
    """Whether `member` is a special method carrying the interpreter's own
    docstring ("Return repr(self)."), which says nothing of geomcore."""
    if not (name.startswith("__") and name.endswith("__")):
        return False
    if type(member).__name__ == "wrapper_descriptor":
        return True
    reference = getattr(object, name, None)
    return reference is not None and member.__doc__ == reference.__doc__


def formatted(docstring, indent):
    """`docstring` as the stub writes it, opening at `indent` columns."""
    lines = docstring.splitlines()
    if len(lines) == 1:
        return f'"""{lines[0]}"""'
    margin = " " * indent
    body = "\n".join(margin + line if line else "" for line in lines[1:])
    return f'"""{lines[0]}\n{body}\n{margin}"""'


def edits_for(node, docstring, lines):
    """The (start, end, text) edit that gives `node` its docstring, ranges
    as (line, column) with lines from 1, or None when it already has it."""
    first = node.body[0]
    indent = node.col_offset + 4
    text = formatted(docstring, indent)
    if (
        isinstance(first, ast.Expr)
        and isinstance(first.value, ast.Constant)
        and isinstance(first.value.value, str)
    ):
        if ast.get_docstring(node) == docstring:
            return None
        return ((first.lineno, first.col_offset), (first.end_lineno, first.end_col_offset), text)
    # The body opens with `...` or a statement: the docstring goes before it.
    header_end = lines[first.lineno - 1][: first.col_offset].rstrip()
    if header_end.endswith(":"):
        # `def f(): ...` or `class C: ...` on one line.
        start = (first.lineno, len(header_end))
        if len(node.body) == 1 and isinstance(first, ast.Expr) and first.value.value is Ellipsis:
            return (start, (first.end_lineno, first.end_col_offset), f"\n{' ' * indent}{text}")
        return (start, (first.lineno, first.col_offset), f"\n{' ' * indent}{text}\n{' ' * indent}")
    return ((first.lineno, first.col_offset), (first.lineno, first.col_offset), f"{text}\n{' ' * indent}")


def stub_edits(source):
    tree = ast.parse(source)
    lines = source.splitlines()
    edits = []

    def visit(node, owner, name):
        docstring = runtime_docstring(owner, name)
        if docstring is not None:
            edit = edits_for(node, docstring, lines)
            if edit is not None:
                edits.append(edit)

    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and hasattr(geomcore, node.name):
            # The runtime function: its docstring, as for any member.
            docstring = getattr(geomcore, node.name).__doc__
            if docstring:
                edit = edits_for(node, inspect.cleandoc(docstring), lines)
                if edit is not None:
                    edits.append(edit)
        elif isinstance(node, ast.ClassDef) and hasattr(geomcore, node.name):
            owner = getattr(geomcore, node.name)
            docstring = owner.__doc__
            if docstring:
                edit = edits_for(node, inspect.cleandoc(docstring), lines)
                if edit is not None:
                    edits.append(edit)
            seen = set()
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and item.name not in seen:
                    seen.add(item.name)  # the first overload only
                    visit(item, owner, item.name)
    return edits


def apply(source, edits):
    lines = source.splitlines(keepends=True)
    offsets = [0]
    for line in lines:
        offsets.append(offsets[-1] + len(line))

    def position(line_and_column):
        line, column = line_and_column
        return offsets[line - 1] + len(lines[line - 1].encode()[:column].decode())

    for start, end, text in sorted(edits, reverse=True):
        source = source[: position(start)] + text + source[position(end) :]
    return source


def main():
    source = STUB.read_text(encoding="utf-8")
    edits = stub_edits(source)
    STUB.write_text(apply(source, edits), encoding="utf-8")
    print(f"{len(edits)} docstrings written into {STUB.name}", file=sys.stderr)


if __name__ == "__main__":
    main()
