"""
Render one scene with ManimGL (-w and a quality flag), failing loudly on any LaTeX error.
ManimGL reuses latex_cache/working.dvi when a compile fails while an old DVI exists, so a broken
formula can silently show the previous one; here the old DVI is deleted before every compile.
usage (cwd = a run dir containing custom_config.yml):  python strict_render.py <scene file> <class> [-l|--hd]
"""
import sys
from pathlib import Path

scene_file, scene_name = sys.argv[1], sys.argv[2]
quality = sys.argv[3] if len(sys.argv) > 3 else "-l"
sys.argv = ["manimgl", scene_file, scene_name, "-w", quality]

import manimlib.utils.tex_file_writing as tfw  # noqa: E402  (manimlib reads sys.argv on import)
from manimlib.config import manim_config  # noqa: E402

_compile = tfw.full_tex_to_svg


def strict_full_tex_to_svg(full_tex, compiler="latex", message=""):
    work = Path(manim_config.directories.latex_cache)
    for ext in (".dvi", ".xdv"):
        stale = work / f"working{ext}"
        if stale.exists():
            stale.unlink()
    return _compile(full_tex, compiler, message)


tfw.full_tex_to_svg = strict_full_tex_to_svg

from manimlib.__main__ import main  # noqa: E402

main()
