"""MkDocs hooks for the geomcore docs site."""

from __future__ import annotations

import sys
from pathlib import Path


def on_config(config):
    """Register the geomcore-aware Python lexer before Markdown is rendered."""
    docs = Path(__file__).resolve().parent
    if str(docs) not in sys.path:
        sys.path.insert(0, str(docs))
    from highlighting.geomcore_lexer import install

    install()
    return config
