"""Every Python block of the docs runs as written (#331).

Each page's first Python block is its setup: every other block runs after
it alone, as a reader pasting that block would, and the page runs in
order. A block that imports `geomcore_viewer` is skipped where the viewer
is not installed."""

import importlib.util
import re
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
PAGES = [ROOT / "README.md", ROOT / "geomcore-py" / "README.md"] + sorted(
    (ROOT / "docs").rglob("*.md")
)
FENCE = re.compile(r"^(?P<indent>[ \t]*)```(?P<language>[\w+-]*)[ \t]*$")
HAS_VIEWER = importlib.util.find_spec("geomcore_viewer") is not None


def python_blocks(page):
    """The page's Python blocks, as (line, code), indentation removed."""
    blocks, lines = [], page.read_text(encoding="utf-8").splitlines()
    index = 0
    while index < len(lines):
        opening = FENCE.match(lines[index])
        if not opening or not opening["language"]:
            index += 1
            continue
        indent, start, body = opening["indent"], index + 1, []
        index += 1
        while index < len(lines) and lines[index].rstrip() != indent + "```":
            body.append(lines[index])
            index += 1
        index += 1
        if opening["language"] == "python":
            blocks.append((start, textwrap.dedent("\n".join(body))))
    return blocks


def run(page, blocks, namespace):
    for line, code in blocks:
        if "geomcore_viewer" in code and not HAS_VIEWER:
            pytest.skip("geomcore_viewer is not installed")
        try:
            exec(compile(code, f"{page.relative_to(ROOT)}:{line}", "exec"), namespace)
        except Exception as error:
            raise AssertionError(
                f"the block at {page.relative_to(ROOT)}:{line} raised "
                f"{type(error).__name__}: {error}"
            ) from error


PAGES_WITH_BLOCKS = [page for page in PAGES if python_blocks(page)]
BLOCKS_AFTER_SETUP = [
    pytest.param(page, number, id=f"{page.relative_to(ROOT)}:{line}")
    for page in PAGES_WITH_BLOCKS
    for number, (line, _) in enumerate(python_blocks(page))
    if number > 0
]


@pytest.mark.parametrize(
    "page", PAGES_WITH_BLOCKS, ids=[str(page.relative_to(ROOT)) for page in PAGES_WITH_BLOCKS]
)
def test_each_page_runs_in_order(page):
    run(page, python_blocks(page), {"__name__": "__docs__"})


@pytest.mark.parametrize("page,number", BLOCKS_AFTER_SETUP)
def test_each_block_runs_after_its_page_setup_alone(page, number):
    blocks = python_blocks(page)
    run(page, [blocks[0], blocks[number]], {"__name__": "__docs__"})


def test_the_pages_hold_python_blocks():
    # A guard against the extraction silently finding nothing.
    assert sum(len(python_blocks(page)) for page in PAGES) >= 40
