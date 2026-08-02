"""
3D structure visualization using MolViewSpec (declarative Mol* view
specs) — renders only the chain(s) you specify, each in its own color.

Requires the `viewer` extra: ``pip install -e ".[viewer]"`` (installs
`molviewspec`). Rendering inline in a notebook additionally needs
`IPython`, which is already present in any Jupyter/Colab environment.
"""

from __future__ import annotations

from typing import Optional, Sequence, Union

from ._client import ISDA_BASE_URL, logger

#: 10-color categorical palette, cycled if there are more chains than colors.
_DEFAULT_CHAIN_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]


def _require_molviewspec():
    try:
        import molviewspec as mvs
        return mvs
    except ImportError as exc:
        raise ImportError(
            "Structure visualization requires molviewspec. Install it with: "
            'pip install -e ".[viewer]"  (or: pip install molviewspec)'
        ) from exc


def _structure_download_url(pdb_id: str, source: str = "isda") -> str:
    pdb_id_lower = pdb_id.lower()
    if source == "isda":
        return f"{ISDA_BASE_URL}/download.{pdb_id_lower}.cif"
    if source == "rcsb":
        return f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    raise ValueError(f"Unknown source '{source}'; expected 'isda' or 'rcsb'.")


def build_chain_view(
    pdb_id: str,
    chain_ids: Union[str, Sequence[str]],
    source: str = "isda",
    colors: Optional[Sequence[str]] = None,
    background_color: Optional[str] = None,
):
    """
    Build a MolViewSpec state showing only the given chain(s) of a PDB
    structure — each in cartoon representation with a distinct color, and
    everything else (other chains, solvent, etc.) left out of the scene.

    Args:
        pdb_id: PDB accession, e.g. "6q0j".
        chain_ids: A chain ID (e.g. "A") or list of chain IDs to display.
            Matched against `auth_asym_id` — the same chain identifiers
            used elsewhere in this package (`pdb_auth_chain` from
            `get_residue_map`, `auth_chain_id` in
            `generate_mutated_structure_script`, `chain_id` in
            `calculate_sasa_for_chain`, etc).
        source: Where to fetch the structure from: "isda" (default) uses
            the IBDC ISDA download endpoint (`{ISDA_BASE_URL}/download.<pdb_id>.cif`);
            "rcsb" fetches from RCSB instead. The download itself happens
            client-side in the browser/Mol*, not in Python.
        colors: Optional list of colors (hex or CSS names), one per
            chain, cycled if there are more chains than colors. Defaults
            to a 10-color categorical palette.
        background_color: Optional canvas background color.

    Returns:
        A `molviewspec` builder (`Root`) instance. Pass it to
        `show_structure_chains` to render inline in a notebook, to
        `save_structure_html` to save a standalone HTML file, or call its
        own `.molstar_html(...)` / `.molstar_notebook(...)` methods
        directly for more control.
    """
    mvs = _require_molviewspec()

    if isinstance(chain_ids, str):
        chain_ids = [chain_ids]
    chain_ids = list(chain_ids)
    if not chain_ids:
        raise ValueError("chain_ids must contain at least one chain identifier.")

    palette = list(colors) if colors else _DEFAULT_CHAIN_COLORS

    builder = mvs.create_builder()
    url = _structure_download_url(pdb_id, source=source)
    model = builder.download(url=url).parse(format="mmcif").model_structure()

    for i, chain_id in enumerate(chain_ids):
        color = palette[i % len(palette)]
        (
            model.component(selector=mvs.ComponentExpression(auth_asym_id=chain_id))
            .representation(type="cartoon")
            .color(color=color)
        )

    if background_color:
        builder.canvas(background_color=background_color)

    return builder


def show_structure_chains(
    pdb_id: str,
    chain_ids: Union[str, Sequence[str]],
    source: str = "isda",
    colors: Optional[Sequence[str]] = None,
    background_color: Optional[str] = None,
    width: Union[int, str] = 950,
    height: Union[int, str] = 600,
) -> None:
    """
    Render only the given chain(s) of a PDB structure inline, in a
    Jupyter notebook or Google Colab cell, via MolViewSpec + Mol*.

    Args:
        pdb_id, chain_ids, source, colors, background_color: see
            `build_chain_view`.
        width, height: Viewer iframe size in pixels (int) or any CSS size
            string (e.g. "100%").

    Returns:
        None — the viewer is displayed as a side effect (via
        `IPython.display`), matching `molviewspec`'s own
        `molstar_notebook` behavior. Requires `IPython` to be installed
        (already the case in any Jupyter/Colab environment).
    """
    builder = build_chain_view(
        pdb_id, chain_ids, source=source, colors=colors, background_color=background_color
    )
    builder.molstar_notebook(width=width, height=height)


def save_structure_html(
    pdb_id: str,
    chain_ids: Union[str, Sequence[str]],
    path: str,
    source: str = "isda",
    colors: Optional[Sequence[str]] = None,
    background_color: Optional[str] = None,
) -> str:
    """
    Save a standalone HTML file showing only the given chain(s) of a PDB
    structure — open it in any browser, no Jupyter/IPython required.

    Args:
        pdb_id, chain_ids, source, colors, background_color: see
            `build_chain_view`.
        path: Output file path, e.g. "structure_view.html".

    Returns:
        `path`, for convenient chaining.
    """
    mvs = _require_molviewspec()
    builder = build_chain_view(
        pdb_id, chain_ids, source=source, colors=colors, background_color=background_color
    )
    html = mvs.molstar_html(builder)
    with open(path, "w") as fh:
        fh.write(html)
    logger.info("Saved structure viewer HTML to %s", path)
    return path
