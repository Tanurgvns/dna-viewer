"""Builds the HTML for the scrolling double-helix view."""

import json
from pathlib import Path

_TEMPLATE = (Path(__file__).parent / "helix_view.html").read_text(encoding="utf-8")


def _js(value) -> str:
    # json.dumps gives a safe JS literal; escaping "</" stops it closing the <script> tag
    return json.dumps(value).replace("</", "<\\/")


def helix_html(sequence: str, file_name: str) -> str:
    return _TEMPLATE.replace("__SEQ__", _js(sequence)).replace("__NAME__", _js(file_name))
