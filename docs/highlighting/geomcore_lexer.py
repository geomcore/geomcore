"""Pygments lexer that colours geomcore API names in docs snippets.

Stock Python highlighting leaves almost every identifier as a plain Name,
so Material paints `Point`, `Circle` and `eval_point` the same as local
variables. This lexer reads the package stubs (no compiled wheel needed)
and retags public types and callables the way build123d's docs lexer does.
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

from pygments.lexers.python import PythonLexer
from pygments.token import Name

_STUB_DIR = Path(__file__).resolve().parents[2] / "geomcore-py" / "geomcore"
_CLASS_RE = re.compile(r"^class\s+([A-Za-z_]\w*)\b", re.M)
_TOP_DEF_RE = re.compile(r"^def\s+([A-Za-z_]\w*)\b", re.M)
_CONSTANT_RE = re.compile(
    r"^\s+([A-Z][A-Z0-9_]*)\s*:\s*ClassVar\b",
    re.M,
)

_installed = False
_original_get_lexer_by_name = None
_original_find_lexer_class_by_name = None


@lru_cache(maxsize=1)
def _api_symbols() -> tuple[frozenset[str], frozenset[str], frozenset[str]]:
    """Classes, module-level functions and ALL_CAPS constants from the stubs."""
    classes: set[str] = set()
    functions: set[str] = set()
    constants: set[str] = set()
    if not _STUB_DIR.is_dir():
        return frozenset(), frozenset(), frozenset()
    for stub in sorted(_STUB_DIR.glob("*.pyi")):
        text = stub.read_text(encoding="utf-8")
        classes.update(_CLASS_RE.findall(text))
        functions.update(_TOP_DEF_RE.findall(text))
        constants.update(_CONSTANT_RE.findall(text))
    functions -= classes
    constants -= classes
    return frozenset(classes), frozenset(functions), frozenset(constants)


class GeomcorePythonLexer(PythonLexer):
    """Python lexer with geomcore types, callables and constants tagged."""

    name = "Geomcore Python"
    aliases = ["python", "py", "python3", "py3", "geomcore"]
    url = "https://geomcore.github.io/geomcore/"
    version_added = ""

    def get_tokens_unprocessed(self, text):
        classes, functions, constants = _api_symbols()
        after_dot = False
        for index, token, value in super().get_tokens_unprocessed(text):
            if value == ".":
                after_dot = True
                yield index, token, value
                continue

            if token in (Name, Name.Other) and value.isidentifier():
                if after_dot:
                    # Attribute access: colour like an IDE (class / const /
                    # method), matching build123d's docs lexer.
                    if value in constants or value.isupper():
                        yield index, Name.Constant, value
                    elif value[0].isupper():
                        yield index, Name.Class, value
                    else:
                        yield index, Name.Function, value
                    after_dot = False
                    continue

                if value in classes:
                    yield index, Name.Class, value
                    continue
                if value in functions:
                    yield index, Name.Function, value
                    continue
                if value in constants:
                    yield index, Name.Constant, value
                    continue

            after_dot = False
            yield index, token, value


def install() -> None:
    """Make ```python fences use GeomcorePythonLexer for this process."""
    global _installed, _original_get_lexer_by_name, _original_find_lexer_class_by_name
    if _installed:
        return

    import pygments.lexers as lexers

    _original_get_lexer_by_name = lexers.get_lexer_by_name
    _original_find_lexer_class_by_name = lexers.find_lexer_class_by_name

    def get_lexer_by_name(alias, **options):
        if alias in GeomcorePythonLexer.aliases:
            return GeomcorePythonLexer(**options)
        return _original_get_lexer_by_name(alias, **options)

    def find_lexer_class_by_name(alias):
        if alias in GeomcorePythonLexer.aliases:
            return GeomcorePythonLexer
        return _original_find_lexer_class_by_name(alias)

    lexers.get_lexer_by_name = get_lexer_by_name
    lexers.find_lexer_class_by_name = find_lexer_class_by_name
    try:
        import pymdownx.highlight as highlight

        highlight.get_lexer_by_name = get_lexer_by_name
    except ImportError:
        pass
    _installed = True
