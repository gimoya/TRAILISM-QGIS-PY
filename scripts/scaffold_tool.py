#!/usr/bin/env python3
"""
Scaffold a new TRAILISM tool's documentation in one command.

What it does:
- Creates a sidecar Markdown next to the .py script with reasonable content extracted from the class
  (title, algorithm id, short help, inputs/outputs/parameters when discoverable).
- Creates a docs page under docs/tools/ that includes the sidecar via MkDocs snippets.
- Inserts the new page into mkdocs.yml navigation under the Tools section (if not already present).
- Refreshes the README tools list from all .py files that have a same-named .md (title and description from .md or .py).
- Optionally runs a local mkdocs build for validation.

Usage examples:
  python scripts/scaffold_tool.py --file TRAILISM__New_Tool.py
  python scripts/scaffold_tool.py --file path/to/TRAILISM__New_Tool.py --title "TRAILISM: New Tool" --slug new-tool --validate

This script is intentionally dependency-light. No YAML libs required; nav is updated via safe text surgery.
"""

import argparse
import os
import re
import subprocess
import sys
from typing import Dict, List, Optional, Tuple


def read_text(path: str) -> str:
    with open(path, "r", encoding="utf-8") as f:
        return f.read()


def write_text(path: str, content: str) -> None:
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)


def infer_display_name(py_src: str) -> Optional[str]:
    # Try displayName: return self.tr('...') or return '...'
    m = re.search(r"def\s+displayName\s*\(.*?\):\s*return\s+self\.tr\(\s*(['\"])(.*?)\1\s*\)", py_src, re.DOTALL)
    if m:
        return m.group(2).strip()
    m = re.search(r"def\s+displayName\s*\(.*?\):\s*return\s*(['\"])(.*?)\1", py_src)
    if m:
        return m.group(2).strip()
    return None


def infer_algorithm_id(py_src: str) -> Optional[str]:
    m = re.search(r"def\s+name\s*\(.*?\):\s*return\s*(['\"])(.*?)\1", py_src)
    return m.group(2).strip() if m else None


def infer_short_help(py_src: str) -> Optional[str]:
    # shortHelpString returns self.tr('...') possibly concatenated. Find tr( and match to closing ), skipping string literals.
    m = re.search(r"def\s+shortHelpString\s*\([^)]*\):\s*return\s+self\.tr\s*\(", py_src)
    if not m:
        return None
    start = m.end()
    depth = 1
    i = start
    while i < len(py_src) and depth > 0:
        # Skip string literals so parens inside don't break depth
        if py_src[i] in ("'", '"') and (i == 0 or py_src[i - 1] != "\\"):
            quote = py_src[i]
            i += 1
            while i < len(py_src):
                if py_src[i] == "\\":
                    i += 2
                    continue
                if py_src[i] == quote:
                    i += 1
                    break
                i += 1
            continue
        if py_src[i] == "(":
            depth += 1
        elif py_src[i] == ")":
            depth -= 1
        i += 1
    if depth != 0:
        return None
    body = py_src[start : i - 1]
    # Extract concatenated string literals inside tr(...)
    lits = re.findall(r"(['\"])(.*?)\1", body, re.DOTALL)
    text = "".join(seg for _q, seg in lits)
    # Normalize whitespace
    text = re.sub(r"\r\n|\r", "\n", text)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip()


def infer_parameters(py_src: str) -> List[Dict[str, Optional[str]]]:
    """Extract parameter declarations from addParameter(QgsProcessingParameterXxx(...)).
    Returns list of dicts with keys: key, kind, label, default.
    """
    out: List[Dict[str, Optional[str]]] = []
    # Roughly match blocks: QgsProcessingParameterXxx( self.KEY, self.tr('Label'...) , ... )
    pattern = re.compile(
        r"QgsProcessingParameter(?P<kind>[A-Za-z]+)\(\s*(?:self\.)?(?P<key>[A-Z0-9_]+)\s*,\s*(?:self\.tr\()?(?P<q>['\"])(?P<label>.*?)(?P=q)\)?(?P<rest>[\s\S]*?)\)",
        re.DOTALL,
    )
    for m in pattern.finditer(py_src):
        kind = m.group("kind")
        key = m.group("key")
        label = m.group("label").strip()
        rest = m.group("rest")
        default = None
        dm = re.search(r"defaultValue\s*=\s*([^,\)\n]+)", rest)
        if dm:
            default = dm.group(1).strip()
        out.append({"key": key, "kind": kind, "label": label, "default": default})
    return out


def infer_io(py_src: str) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    """Return (inputs, outputs) from parameter type names.
    input kinds include FeatureSource, RasterLayer, etc.; outputs are FeatureSink.
    """
    params = infer_parameters(py_src)
    inputs: List[Dict[str, str]] = []
    outputs: List[Dict[str, str]] = []
    for p in params:
        kind = (p.get("kind") or "").lower()
        if "sink" in kind:
            outputs.append({"key": p["key"], "label": p.get("label", "")})
        else:
            inputs.append({"key": p["key"], "label": p.get("label", ""), "kind": p.get("kind", "")})
    return inputs, outputs


def derive_slug(title: str) -> str:
    base = title.strip().lower()
    base = re.sub(r"[^a-z0-9]+", "-", base)
    base = re.sub(r"-+", "-", base).strip("-")
    return base or "tool"


def make_sidecar_content(title: str, alg_id: Optional[str], help_text: Optional[str], inputs, outputs, params) -> str:
    lines: List[str] = []
    lines.append(f"## {title}")
    lines.append("")
    if help_text:
        # First paragraph as summary
        summary = help_text.split("\n\n")[0].strip()
        lines.append(summary)
        lines.append("")
    if alg_id:
        lines.append("### Algorithm id")
        lines.append(f"`{alg_id}`")
        lines.append("")
    if inputs:
        lines.append("### Inputs")
        for i in inputs:
            kind = i.get("kind", "")
            label = i.get("label", "")
            lines.append(f"- **{label}** ({kind})")
        lines.append("")
    if params:
        lines.append("### Parameters")
        for p in params:
            label = p.get("label") or p.get("key")
            kind = p.get("kind") or "Param"
            default = p.get("default")
            if default:
                lines.append(f"- **{label}** ({kind}, default {default})")
            else:
                lines.append(f"- **{label}** ({kind})")
        lines.append("")
    if outputs:
        lines.append("### Outputs")
        for o in outputs:
            lines.append(f"- **{o.get('label','')}**")
        lines.append("")
    # Usage block
    if alg_id:
        # Gather parameter keys for a minimal skeleton
        keys = [p.get("key") for p in params if p.get("key")]
        mapping = ",\n    ".join([f"'{k}': ..." for k in keys])
        if mapping:
            mapping += ",\n    'OUTPUT': 'memory:'"
        else:
            mapping = "'OUTPUT': 'memory:'"
        lines.append("### Usage")
        lines.append("- Python:")
        lines.append("```python")
        lines.append("import processing")
        lines.append(f"processing.run('{alg_id}', {{\n    {mapping}\n}})")
        lines.append("```")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def _first_heading_text(md_content: str) -> Optional[str]:
    """First line that starts with ## (with optional leading whitespace); strip the # and whitespace."""
    for line in md_content.splitlines():
        s = line.strip()
        if s.startswith("##"):
            return s.lstrip("#").strip()
    return None


def _first_description_from_md(md_content: str) -> str:
    """First paragraph of body text: after ### What it does, or first non-heading paragraph. Single line, max ~120 chars."""
    paragraphs = re.split(r"\n\s*\n", md_content)
    for p in paragraphs:
        lines = p.strip().splitlines()
        if not lines:
            continue
        first = lines[0].strip()
        if first.startswith("#"):
            if "what it does" in first.lower() and len(lines) > 1:
                desc = " ".join(lines[1:]).strip()
                if len(desc) > 120:
                    desc = desc[:117] + "..."
                return desc
            continue
        desc = p.strip()
        if len(desc) > 120:
            desc = desc[:117] + "..."
        return desc
    return ""


def discover_tools(repo_root: str) -> List[Tuple[str, str, str]]:
    """Find all TRAILISM__*.py that have a same-named .md. Return list of (base_name, title, description) sorted by base_name."""
    tools: List[Tuple[str, str, str]] = []
    if not os.path.isdir(repo_root):
        return tools
    for name in sorted(os.listdir(repo_root)):
        if not name.startswith("TRAILISM__") or not name.endswith(".py"):
            continue
        base = name[:-3]
        md_name = base + ".md"
        md_path = os.path.join(repo_root, md_name)
        py_path = os.path.join(repo_root, name)
        if not os.path.isfile(md_path) or not os.path.isfile(py_path):
            continue
        md_content = read_text(md_path)
        py_src = read_text(py_path)
        title = _first_heading_text(md_content) or infer_display_name(py_src) or base
        description = _first_description_from_md(md_content)
        if not description and infer_short_help(py_src):
            help_text = infer_short_help(py_src)
            first = help_text.split("\n\n")[0].strip()
            if len(first) > 120:
                first = first[:117] + "..."
            description = first
        if not description:
            description = "No description"
        tools.append((base, title, description))
    return tools


def update_readme_tools_list(repo_root: str) -> None:
    readme_path = os.path.join(repo_root, "README.md")
    if not os.path.isfile(readme_path):
        return
    content = read_text(readme_path)
    begin_marker = "<!-- TOOLS LIST BEGIN -->"
    end_marker = "<!-- TOOLS LIST END -->"
    if begin_marker not in content or end_marker not in content:
        return
    tools = discover_tools(repo_root)
    if not tools:
        return
    bullets = "\n".join(f"- {title} — {desc}" for (_base, title, desc) in tools)
    intro_match = re.search(
        r"This repo currently contains (?:eight|\d+) tools:\s*" + re.escape(begin_marker),
        content,
    )
    if not intro_match:
        return
    suffix_start = content.find(end_marker) + len(end_marker)
    new_content = (
        content[: intro_match.start()]
        + f"This repo currently contains {len(tools)} tools:\n\n"
        + begin_marker
        + "\n"
        + bullets
        + "\n"
        + end_marker
        + content[suffix_start:]
    )
    write_text(readme_path, new_content)


def update_mkdocs_nav(mkdocs_path: str, title: str, slug: str) -> None:
    y = read_text(mkdocs_path)
    if f"tools/{slug}.md" in y:
        return
    lines = y.splitlines()
    out_lines: List[str] = []
    i = 0
    inserted = False
    while i < len(lines):
        line = lines[i]
        out_lines.append(line)
        if not inserted and line.strip() == "- Tools:":
            # Insert after the Tools: header; maintain existing 6-space indentation for items
            out_lines.append(f"      - {title}: tools/{slug}.md")
            inserted = True
        i += 1
    if not inserted:
        # Append a nav section if Tools not found
        out_lines.append("nav:")
        out_lines.append("  - Tools:")
        out_lines.append(f"      - {title}: tools/{slug}.md")
    write_text(mkdocs_path, "\n".join(out_lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description="Scaffold a new tool's docs")
    ap.add_argument("--file", required=True, help="Path to the new tool .py script")
    ap.add_argument("--title", default=None, help="Display title (defaults to displayName() or filename)")
    ap.add_argument("--slug", default=None, help="Docs page slug (defaults from title)")
    ap.add_argument("--validate", action="store_true", help="Run mkdocs build --strict after generation")
    args = ap.parse_args()

    script_path = os.path.normpath(args.file)
    if not os.path.isfile(script_path):
        print(f"ERROR: File not found: {script_path}", file=sys.stderr)
        return 2

    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    py_src = read_text(script_path)

    # Derive properties
    display_name = args.title or infer_display_name(py_src) or os.path.splitext(os.path.basename(script_path))[0]
    alg_id = infer_algorithm_id(py_src)
    short_help = infer_short_help(py_src)
    params = infer_parameters(py_src)
    inputs, outputs = infer_io(py_src)
    slug = args.slug or derive_slug(display_name)

    # Sidecar path (next to .py)
    base_name = os.path.splitext(os.path.basename(script_path))[0]
    sidecar_name = base_name + ".md"
    sidecar_path = os.path.join(os.path.dirname(script_path), sidecar_name)

    if os.path.exists(sidecar_path):
        print(f"NOTE: Sidecar already exists: {sidecar_path} (will overwrite)")
    sidecar_md = make_sidecar_content(display_name, alg_id, short_help, inputs, outputs, params)
    write_text(sidecar_path, sidecar_md)

    # Docs page
    docs_page_dir = os.path.join(repo_root, "docs", "tools")
    docs_page_path = os.path.join(docs_page_dir, f"{slug}.md")
    if os.path.exists(docs_page_path):
        print(f"NOTE: Docs page already exists: {docs_page_path} (will overwrite)")
    include_rel = sidecar_name  # MkDocs snippets base_path includes '.' so root-relative sidecar name works
    docs_page = f"# {display_name}\n\n--8<-- \"{include_rel}\"\n\n"
    write_text(docs_page_path, docs_page)

    # mkdocs nav
    mkdocs_path = os.path.join(repo_root, "mkdocs.yml")
    update_mkdocs_nav(mkdocs_path, display_name, slug)

    # README tools list (from all .py + same-named .md in repo root)
    update_readme_tools_list(repo_root)

    print(f"Created sidecar: {os.path.relpath(sidecar_path, repo_root)}")
    print(f"Created docs page: {os.path.relpath(docs_page_path, repo_root)}")
    print("Updated mkdocs.yml nav")
    print("Updated README.md tools list")

    if args.validate:
        try:
            subprocess.check_call([sys.executable, "-m", "mkdocs", "build", "--strict"], cwd=repo_root)
            print("MkDocs build OK")
        except Exception as e:
            print(f"MkDocs build failed: {e}", file=sys.stderr)
            return 3

    return 0


if __name__ == "__main__":
    sys.exit(main())


