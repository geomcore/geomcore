"""MkDocs hooks for the geomcore docs site."""

from __future__ import annotations

import html
import re
import sys
from pathlib import Path
from urllib.parse import quote, urljoin


def on_config(config):
    """Register the geomcore-aware Python lexer before Markdown is rendered."""
    here = Path(__file__).resolve().parent
    if str(here) not in sys.path:
        sys.path.insert(0, str(here))
    from highlighting.geomcore_lexer import install

    install()
    return config


# A fenced ``viewer`` block becomes one interactive geomcore-viewer iframe.
#
#     ```viewer
#     scene: assets/scenes/projection-feet.json
#     title: Ellipse with the four projection feet of one query point
#     height: 480
#     ```
#
# ``scene`` is relative to the docs root (or an absolute URL), and ``height``
# is optional. Keeping the scene path, iframe, caption and full-screen link in
# one place means every embed is consistent and authors never hand-write HTML.
#
# The scanner below is fence-aware: a ```viewer shown inside a longer code
# fence (as the viewer page documents the syntax) is left alone.
_OPENING_FENCE = re.compile(r"^(?P<ticks>`{3,})(?P<info>[^\r\n]*)$")


def _closes_fence(line: str, ticks: str) -> bool:
    stripped = line.rstrip("\r\n").strip()
    return (
        stripped.startswith("`" * len(ticks))
        and stripped == "`" * len(stripped)
        and len(stripped) >= len(ticks)
    )


def _parse_viewer_options(body: str) -> dict[str, str]:
    options: dict[str, str] = {}
    for raw_line in body.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        key, separator, value = line.partition(":")
        if not separator:
            raise ValueError(
                f"viewer block: expected 'key: value', got {raw_line!r}"
            )
        options[key.strip()] = value.strip()
    return options


def _render_viewer_block(body: str, config) -> str:
    options = _parse_viewer_options(body)
    allowed = {"scene", "title", "height"}
    unknown = set(options) - allowed
    if unknown:
        names = ", ".join(sorted(unknown))
        raise ValueError(
            f"viewer block: unknown option(s) {names}; use scene, title or height"
        )

    site_url = config.get("site_url") or ""
    viewer_url = config.get("extra", {}).get(
        "viewer_url", "https://geomcore.github.io/geomcore-viewer/"
    )
    scene = options.get("scene", "")
    if scene:
        scene_url = scene if "://" in scene else urljoin(site_url, scene)
        separator = "&" if "?" in viewer_url else "?"
        viewer_url = f"{viewer_url}{separator}scene={quote(scene_url, safe=':/?=&%#')}"

    title = options.get("title", "geomcore-viewer")
    height = options.get("height", "480")
    if not height.isdigit() or not 200 <= int(height) <= 1200:
        raise ValueError(
            f"viewer block: height must be a whole number from 200 to 1200, got {height!r}"
        )

    return (
        f'<figure class="gc-viewer" style="--gc-viewer-height: {height}px">\n'
        f'<iframe src="{html.escape(viewer_url, quote=True)}" '
        f'title="{html.escape(title, quote=True)}" loading="lazy" '
        'allow="fullscreen" allowfullscreen></iframe>\n'
        "<figcaption>"
        f"<span>{html.escape(title)}</span>"
        f'<a href="{html.escape(viewer_url, quote=True)}" target="_blank" '
        'rel="noopener">Open in a new tab ↗</a>'
        "</figcaption>\n"
        "</figure>\n"
    )


def _expand_viewer_blocks(markdown: str, config) -> str:
    lines = markdown.splitlines(keepends=True)
    output: list[str] = []
    index = 0
    while index < len(lines):
        opening = _OPENING_FENCE.match(lines[index].rstrip("\r\n"))
        if opening is None:
            output.append(lines[index])
            index += 1
            continue

        ticks = opening.group("ticks")
        info = opening.group("info").strip()
        if info != "viewer":
            # Copy a non-viewer fence whole, so a ```viewer inside it is text.
            output.append(lines[index])
            index += 1
            while index < len(lines):
                output.append(lines[index])
                closed = _closes_fence(lines[index], ticks)
                index += 1
                if closed:
                    break
            continue

        body_start = index + 1
        body_end = body_start
        while body_end < len(lines) and not _closes_fence(lines[body_end], ticks):
            body_end += 1
        if body_end == len(lines):
            raise ValueError("viewer block: missing closing fence")
        output.append(_render_viewer_block("".join(lines[body_start:body_end]), config))
        index = body_end + 1

    return "".join(output)


def on_page_markdown(markdown: str, *, page, config, files) -> str:
    """Expand interactive ``viewer`` blocks before Markdown is rendered."""
    return _expand_viewer_blocks(markdown, config)
