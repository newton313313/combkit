"""
make_api.py —— 用 ast 讀原始碼，重新產生 API.md（不需要 import，也不需要 manim）

    python tools/make_api.py
"""
import ast
from pathlib import Path

DST = Path(__file__).resolve().parents[1]
PKG = DST / "combkit"
DRAW = PKG / "draw"
# 先列這些（依閱讀順序），其餘新加入的模組自動接在後面，不用改這支程式
PURE_FIRST = ["counting", "tableaux", "rep_data", "raney", "tikz"]
DRAW_FIRST = ["styles", "scene_base", "circle_points", "polygon_groups", "root_figure", "lattice_paths",
              "catalan_draw", "mit_draw", "fulton_draw", "rep_draw", "formulas"]


def discover(folder, first):
    found = sorted(f.stem for f in folder.glob("*.py") if f.stem != "__init__")
    return [m for m in first if m in found] + [m for m in found if m not in first]


PUREMODS = discover(PKG, PURE_FIRST)
DRAWMODS = discover(DRAW, DRAW_FIRST)


def first_line(node) -> str:
    doc = ast.get_docstring(node)
    return doc.strip().splitlines()[0].strip() if doc else ""


lines = ["# combkit API 索引", "",
         "由 `tools/make_api.py` 用 ast 從原始碼自動產生；說明取自各函式 docstring 的第一行。手寫的介紹請看 README.md。", ""]
order = [(f"combkit.{n}", PKG / f"{n}.py") for n in PUREMODS] + \
        [(f"combkit.draw.{n}", DRAW / f"{n}.py") for n in DRAWMODS]
for mod, path in order:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    lines += [f"## `{mod}`", ""]
    mdoc = first_line(tree)
    if mdoc:
        lines += [mdoc, ""]
    consts = [t.id for n in tree.body if isinstance(n, ast.Assign) for t in n.targets
              if isinstance(t, ast.Name) and t.id.isupper()]
    if consts:
        lines += ["常數：" + "、".join(f"`{c}`" for c in consts), ""]
    for n in tree.body:
        if isinstance(n, ast.FunctionDef) and not n.name.startswith("_"):
            lines.append(f"- `{n.name}({ast.unparse(n.args)})` — {first_line(n)}")
        elif isinstance(n, ast.ClassDef):
            lines.append(f"- **class `{n.name}`** — {first_line(n)}")
            for m in n.body:
                if isinstance(m, ast.FunctionDef) and not m.name.startswith("_"):
                    args = ast.unparse(m.args).replace("self, ", "").replace("self", "")
                    lines.append(f"  - `.{m.name}({args})` — {first_line(m)}")
    lines.append("")
(DST / "API.md").write_text("\n".join(lines), encoding="utf-8")

n_func = sum(1 for l in lines if l.startswith("- `") or l.startswith("  - `"))
print(f"API.md: {n_func} functions/methods")
