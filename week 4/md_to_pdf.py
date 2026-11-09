#!/usr/bin/env python3
"""Convert a Markdown file (with LaTeX math and images) to PDF.

Usage:
    python md_to_pdf.py [input.md] [output.pdf]

Defaults to ps_4_solutions.md -> ps_4_solutions.pdf next to this script.
Requires pandoc and a XeLaTeX installation (e.g. TeX Live or MiKTeX).
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

LATEX_HEADER = r"""
\usepackage{fvextra}
\DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,breakanywhere,commandchars=\\\{\},fontsize=\footnotesize}
\usepackage{etoolbox}
\AtBeginEnvironment{longtable}{\small\renewcommand{\arraystretch}{1.5}}
\usepackage{float}
\floatplacement{figure}{H}
"""


def convert(md_path: Path, pdf_path: Path) -> None:
    for tool in ("pandoc", "xelatex"):
        if shutil.which(tool) is None:
            sys.exit(f"Error: '{tool}' was not found on PATH.")
    if not md_path.is_file():
        sys.exit(f"Error: input file not found: {md_path}")

    with tempfile.TemporaryDirectory() as tmp:
        header = Path(tmp) / "header.tex"
        header.write_text(LATEX_HEADER, encoding="utf-8")
        cmd = [
            "pandoc", md_path.name,
            "--from", "markdown+tex_math_dollars+pipe_tables+implicit_figures",
            "--output", str(pdf_path.resolve()),
            "--pdf-engine", "xelatex",
            "--highlight-style", "tango",
            "--include-in-header", str(header),
            "-V", "geometry:a4paper,margin=2.2cm",
            "-V", "fontsize=11pt",
            "-V", "colorlinks=true",
        ]
        # Run inside the Markdown's folder so relative image paths resolve.
        result = subprocess.run(cmd, cwd=md_path.parent, capture_output=True, text=True)

    if result.returncode != 0:
        sys.exit(f"pandoc failed:\n{result.stderr}")
    if result.stderr.strip():
        print(result.stderr.strip())
    print(f"Wrote {pdf_path}")


def main() -> None:
    here = Path(__file__).resolve().parent
    md_path = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else here / "ps_4_solutions.md"
    pdf_path = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else md_path.with_suffix(".pdf")
    convert(md_path, pdf_path)


if __name__ == "__main__":
    main()
