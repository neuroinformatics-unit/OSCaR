from pathlib import Path

from oscar.line import Line


def test_line_statistics():
    Line(Path("data/test-data-1.csv"), "Line-1")
