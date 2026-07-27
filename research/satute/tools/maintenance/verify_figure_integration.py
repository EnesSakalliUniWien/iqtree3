#!/usr/bin/env python3
"""Verify the active LaTeX figure registry and narrative order.

Adapted from the figure-integration gate in the ``Is it rotation`` manuscript.
This version is intentionally read-only: it does not copy assets or maintain
hash manifests.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Any


PROJECT = Path(__file__).resolve().parents[2]
REGISTRY = PROJECT / "manuscript" / "figures" / "registry.json"
INCLUDE_PATTERN = re.compile(r"\\input\{([^}]+)\}")
GRAPHIC_PATTERN = re.compile(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}")


def document_sources(root: Path) -> list[tuple[Path, str]]:
    """Return active TeX sources in document order."""

    output: list[tuple[Path, str]] = []
    seen: set[Path] = set()

    def visit(path: Path) -> None:
        resolved = path.resolve()
        if resolved in seen:
            return
        seen.add(resolved)
        text = resolved.read_text(encoding="utf-8")
        output.append((resolved, text))
        for raw in INCLUDE_PATTERN.findall(text):
            candidate = PROJECT / "manuscript" / raw
            if candidate.suffix != ".tex":
                candidate = candidate.with_suffix(".tex")
            if not candidate.is_file():
                raise FileNotFoundError(f"active LaTeX input is missing: {candidate}")
            visit(candidate)

    visit(root)
    return output


def active_assets(sources: list[tuple[Path, str]]) -> list[str]:
    assets: list[str] = []
    for _, text in sources:
        for raw in GRAPHIC_PATTERN.findall(text):
            path = Path(raw)
            if not path.is_absolute():
                path = Path("manuscript") / path
            assets.append(path.as_posix())
    return assets


def has_label(text: str, label: str) -> bool:
    escaped = re.escape(label)
    return bool(
        re.search(rf"\\label\{{{escaped}\}}", text)
        or re.search(rf"\\cref@old@label\{{{escaped}\}}", text)
    )


def load_registry() -> dict[str, Any]:
    registry = json.loads(REGISTRY.read_text(encoding="utf-8"))
    if registry.get("schema_version") != 1:
        raise ValueError("figure registry schema_version must be 1")
    if not isinstance(registry.get("figures"), list) or not registry["figures"]:
        raise ValueError("figure registry must contain a non-empty figures list")
    return registry


def main() -> int:
    failures: list[str] = []
    try:
        registry = load_registry()
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"Figure integration audit failed:\n  - {error}", file=sys.stderr)
        return 1

    required = {
        "id",
        "kind",
        "order",
        "label",
        "tex_source",
        "assets",
        "renderer",
        "source_data",
        "takeaway",
    }
    seen_ids: set[str] = set()
    orders: dict[str, list[int]] = {"main": [], "supplement": []}
    registered_assets: set[str] = set()

    documents: dict[str, list[tuple[Path, str]]] = {}
    for kind, relative in registry.get("documents", {}).items():
        root = PROJECT / relative
        try:
            documents[kind] = document_sources(root)
        except (OSError, FileNotFoundError) as error:
            failures.append(str(error))

    for entry in registry["figures"]:
        missing = required - set(entry)
        if missing:
            failures.append(f"{entry.get('id', '<unknown>')}: missing fields {sorted(missing)}")
            continue
        figure_id = str(entry["id"])
        kind = str(entry["kind"])
        if figure_id in seen_ids:
            failures.append(f"duplicate figure id: {figure_id}")
        seen_ids.add(figure_id)
        if kind not in orders:
            failures.append(f"{figure_id}: kind must be main or supplement")
            continue
        orders[kind].append(int(entry["order"]))
        if not str(entry["takeaway"]).strip():
            failures.append(f"{figure_id}: takeaway must not be empty")

        tex_path = PROJECT / entry["tex_source"]
        renderer = PROJECT / entry["renderer"]
        if not tex_path.is_file():
            failures.append(f"{figure_id}: missing TeX source: {entry['tex_source']}")
            tex_text = ""
        else:
            tex_text = tex_path.read_text(encoding="utf-8")
        if not renderer.is_file():
            failures.append(f"{figure_id}: missing renderer: {entry['renderer']}")
        if entry["label"] and not has_label(tex_text, str(entry["label"])):
            failures.append(f"{figure_id}: label {entry['label']!r} is absent from its TeX source")

        for source_data in entry["source_data"]:
            if not (PROJECT / source_data).is_file():
                failures.append(f"{figure_id}: missing source data: {source_data}")
        for asset in entry["assets"]:
            registered_assets.add(str(asset))
            if not (PROJECT / asset).is_file():
                failures.append(f"{figure_id}: missing figure asset: {asset}")
            raw = Path(asset).relative_to("manuscript").as_posix()
            if raw not in GRAPHIC_PATTERN.findall(tex_text):
                failures.append(f"{figure_id}: TeX source does not include {asset}")

    for kind, values in orders.items():
        expected = list(range(1, len(values) + 1))
        if sorted(values) != expected:
            failures.append(f"{kind} figure orders must be consecutive: {values}")

    active_by_document = {
        kind: active_assets(sources) for kind, sources in documents.items()
    }
    active_all = {
        asset for assets in active_by_document.values() for asset in assets
    }
    unregistered = sorted(active_all - registered_assets)
    inactive = sorted(registered_assets - active_all)
    if unregistered:
        failures.append(f"active LaTeX includes unregistered figures: {unregistered}")
    if inactive:
        failures.append(f"registered figures are not active LaTeX inputs: {inactive}")

    for kind in ("main", "supplement"):
        ordered_entries = sorted(
            (entry for entry in registry["figures"] if entry.get("kind") == kind),
            key=lambda entry: int(entry["order"]),
        )
        active = active_by_document.get(kind, [])
        positions: list[int] = []
        for entry in ordered_entries:
            candidate_positions = [
                active.index(asset) for asset in entry["assets"] if asset in active
            ]
            if candidate_positions:
                positions.append(min(candidate_positions))
        if positions != sorted(positions):
            failures.append(f"{kind} registry order disagrees with active LaTeX figure order")

    if failures:
        print("Figure integration audit failed:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 1

    print(
        "Figure integration audit passed: "
        f"{len(registry['figures'])} displays, {len(registered_assets)} unique active assets."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
