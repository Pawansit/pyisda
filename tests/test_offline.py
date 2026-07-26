"""
Offline tests: cover the pure-Python logic (sequence mutation, mutation
analysis, ChimeraX script generation, residue-map merging) without hitting
the ISDA API. Run with: pytest tests/
"""

import pandas as pd
import pytest

import pyisda as ibdc


def test_mutate_sequence():
    seq = {"sequence": "MKTAYIAKQRQ"}
    muts = pd.DataFrame({"Alt_AA": ["Gly", "Trp"], "Position": [1, 5]})
    result = ibdc.mutate_sequence(seq, muts)
    assert result[0] == "G"
    assert result[4] == "W"
    assert len(result) == len(seq["sequence"])


def test_mutate_sequence_ignores_unknown_codes():
    seq = {"sequence": "MKTAYIAKQRQ"}
    muts = pd.DataFrame({"Alt_AA": ["Xyz"], "Position": [2]})
    result = ibdc.mutate_sequence(seq, muts)
    assert result == seq["sequence"]  # unmappable mutation skipped


def test_one_letter_three_letter_roundtrip():
    assert ibdc.one_letter_to_three("A") == "Ala"
    assert ibdc.three_letter_to_one("Ala") == "A"


def test_analyze_mutation_properties():
    df = pd.DataFrame({
        "Ref_AA": ["Ala", "Gly"],
        "Alt_AA": ["Val", "Ser"],
    })
    out = ibdc.analyze_mutation_properties(df)
    assert "Mutation_Effect_Type" in out.columns
    assert out.loc[0, "Ref_Prop"] == "Hydrophobic"
    assert out.loc[1, "Ref_Prop"] == "Hydrophilic"


def test_analyze_mutation_properties_empty_returns_none():
    assert ibdc.analyze_mutation_properties(pd.DataFrame()) is None
    assert ibdc.analyze_mutation_properties(None) is None


def test_get_mutation_summary():
    df = pd.DataFrame({
        "Clinical Significance": ["Pathogenic", "Benign", "Pathogenic"],
        "Position": [10, 20, 10],
        "Consequence Type": ["missense", "missense", "nonsense"],
    })
    summary = ibdc.get_mutation_summary(df)
    assert summary["mutation_count"] == 3
    assert summary["significance_status"]["Pathogenic"] == 2
    assert summary["maximum_mutation_spot"][10] == 2


def test_get_mutation_summary_missing_columns():
    assert ibdc.get_mutation_summary(pd.DataFrame({"Position": [1]})) is None


def test_filter_by_significance():
    df = pd.DataFrame({"Clinical Significance": ["Pathogenic", "Benign", "Likely pathogenic"]})
    out = ibdc.filter_by_significance(df, "pathogenic")
    assert len(out) == 2


def test_merge_mutations_with_pdb_map():
    mutation_df = pd.DataFrame({"Position": [10, 20], "Alt_AA": ["Gly", "Trp"]})
    residue_map_df = pd.DataFrame({
        "unp_residue": [10, 20],
        "pdb_residue": [110, 120],
        "pdb_auth_chain": ["A", "A"],
    })
    merged = ibdc.merge_mutations_with_pdb_map(mutation_df, residue_map_df)
    assert list(merged["pdb_residue"]) == [110, 120]


def test_merge_mutations_with_pdb_map_after_boolean_filter():
    # Regression test: residue_map_df sliced via boolean indexing (e.g.
    # df[df['pdb_auth_chain'] == chain_id]) used to raise
    # `KeyError: 'Position'` during the merge, because renaming the first
    # column via `.columns.values[0] = ...` corrupts pandas' internal
    # column-lookup index on a filtered slice even though `.columns`
    # itself looks correct.
    mutation_df = pd.DataFrame({"Position": ["1", "2", "4"], "Alt_AA": ["Gly", "Trp", "Ala"]})
    residue_map_df = pd.DataFrame({
        "unp_residue": ["1", "2", "3", "4", "5"],
        "pdb_residue": ["1", "2", "3", "4", "5"],
        "pdb_auth_chain": ["A", "A", "B", "A", "B"],
        "pdb_chain": ["C", "C", "D", "C", "D"],
    })
    filtered = residue_map_df[residue_map_df["pdb_auth_chain"] == "A"]

    merged = ibdc.merge_mutations_with_pdb_map(mutation_df, filtered)
    assert list(merged["Position"]) == ["1", "2", "4"]
    assert list(merged["Alt_AA"]) == ["Gly", "Trp", "Ala"]


def test_merge_mutations_with_pdb_map_rejects_wrong_dataframe():
    # Regression test: passing a computational-mutations-style DataFrame
    # (which already has its own 'Position' column, no 'unp_residue')
    # instead of get_residue_map() output used to fail deep inside pandas
    # with a cryptic "column label 'Position' is not unique" ValueError,
    # because the function blindly renamed whatever the *first* column
    # was to 'Position'. It should now fail immediately with a clear,
    # actionable message instead.
    annotated_mutations = pd.DataFrame({"Position": ["10"], "Alt_AA": ["Gly"]})
    comp_predictions = pd.DataFrame({
        "protein_variant": ["A10G"], "Ref_AA": ["A"], "Position": ["10"], "Alt_AA": ["G"],
    })

    with pytest.raises(ValueError, match="merge_experimental_with_computational"):
        ibdc.merge_mutations_with_pdb_map(annotated_mutations, comp_predictions)


def test_plot_structural_coverage_runs():
    import matplotlib
    matplotlib.use("Agg")
    records = [
        {"PDB": "1M17", "UNIPROT_START": 696, "UNIPROT_END": 1022, "PDB_CHAIN": "A"},
        {"PDB": "3W2S", "UNIPROT_START": 25, "UNIPROT_END": 646, "PDB_CHAIN": "A"},
    ]
    ax = ibdc.plot_structural_coverage(records, sequence_length=1210)
    assert ax is not None
    assert ax.get_xlim() == (0.0, 1210.0)


def test_plot_structural_coverage_empty():
    import matplotlib
    matplotlib.use("Agg")
    ax = ibdc.plot_structural_coverage([])
    assert ax is not None


def test_plot_mutation_lollipop_runs():
    import matplotlib
    matplotlib.use("Agg")
    df = pd.DataFrame({
        "Position": [10, 10, 20],
        "Clinical Significance": ["Pathogenic", "Benign", "Benign"],
    })
    ax = ibdc.plot_mutation_lollipop(df, sequence_length=100)
    assert ax is not None


def test_plot_sasa_profile_runs():
    import matplotlib
    matplotlib.use("Agg")
    import numpy as np
    result = ibdc.SASAResult(
        residue_ids=np.arange(1, 11),
        residue_names=np.array(["ALA"] * 10),
        sasa_per_residue=np.random.rand(10) * 100,
    )
    ax = ibdc.plot_sasa_profile(result)
    assert ax is not None


def test_plot_protein_overview_runs():
    import matplotlib
    matplotlib.use("Agg")
    import numpy as np
    records = [{"PDB": "1M17", "UNIPROT_START": 696, "UNIPROT_END": 1022, "PDB_CHAIN": "A"}]
    mut_df = pd.DataFrame({"Position": [700], "Clinical Significance": ["Pathogenic"]})
    sasa_result = ibdc.SASAResult(
        residue_ids=np.arange(1, 11),
        residue_names=np.array(["ALA"] * 10),
        sasa_per_residue=np.random.rand(10) * 100,
    )
    fig = ibdc.plot_protein_overview(
        sequence_length=1210, pdb_records=records, mutation_df=mut_df, sasa_result=sasa_result, title="Test"
    )
    assert fig is not None
    assert len(fig.axes) == 3


def test_plot_protein_overview_requires_data():
    with pytest.raises(ValueError):
        ibdc.plot_protein_overview(sequence_length=100)


def test_plot_structural_coverage_html_runs():
    pytest.importorskip("plotly")
    records = [
        {"PDB": "1M17", "UNIPROT_START": 696, "UNIPROT_END": 1022, "PDB_CHAIN": "A"},
        {"PDB": "3W2S", "UNIPROT_START": 25, "UNIPROT_END": 646, "PDB_CHAIN": "A"},
    ]
    fig = ibdc.plot_structural_coverage_html(records, sequence_length=1210)
    assert fig is not None
    assert len(fig.data) == 1  # one Bar trace holding all records


def test_plot_structural_coverage_html_empty():
    pytest.importorskip("plotly")
    fig = ibdc.plot_structural_coverage_html([])
    assert fig is not None


def test_plot_mutation_lollipop_html_runs():
    pytest.importorskip("plotly")
    df = pd.DataFrame({
        "Position": [10, 10, 20],
        "Clinical Significance": ["Pathogenic", "Benign", "Benign"],
    })
    fig = ibdc.plot_mutation_lollipop_html(df, sequence_length=100)
    assert fig is not None


def test_plot_sasa_profile_html_runs():
    pytest.importorskip("plotly")
    import numpy as np
    result = ibdc.SASAResult(
        residue_ids=np.arange(1, 11),
        residue_names=np.array(["ALA"] * 10),
        sasa_per_residue=np.random.rand(10) * 100,
    )
    fig = ibdc.plot_sasa_profile_html(result)
    assert fig is not None


def test_plot_protein_overview_html_runs():
    pytest.importorskip("plotly")
    import numpy as np
    records = [{"PDB": "1M17", "UNIPROT_START": 696, "UNIPROT_END": 1022, "PDB_CHAIN": "A"}]
    mut_df = pd.DataFrame({"Position": [700], "Clinical Significance": ["Pathogenic"]})
    sasa_result = ibdc.SASAResult(
        residue_ids=np.arange(1, 11),
        residue_names=np.array(["ALA"] * 10),
        sasa_per_residue=np.random.rand(10) * 100,
    )
    fig = ibdc.plot_protein_overview_html(
        sequence_length=1210, pdb_records=records, mutation_df=mut_df, sasa_result=sasa_result, title="Test"
    )
    assert fig is not None


def test_save_html_writes_file(tmp_path):
    pytest.importorskip("plotly")
    records = [{"PDB": "1M17", "UNIPROT_START": 696, "UNIPROT_END": 1022, "PDB_CHAIN": "A"}]
    fig = ibdc.plot_structural_coverage_html(records, sequence_length=1210)
    out = tmp_path / "coverage.html"
    ibdc.save_html(fig, str(out))
    assert out.exists()
    assert out.stat().st_size > 0


def test_fetch_structure_cif_builds_correct_url(tmp_path, monkeypatch):
    import pyisda.sasa as sasa_mod

    captured = {}

    def fake_download_file(url, dest_path, timeout=30):
        captured["url"] = url
        captured["dest_path"] = dest_path
        with open(dest_path, "w") as f:
            f.write("fake cif")
        return dest_path

    monkeypatch.setattr(sasa_mod, "download_file", fake_download_file)
    path = sasa_mod.fetch_structure_cif("6GEL", output_dir=str(tmp_path))

    assert captured["url"] == "https://ibdc.dbt.gov.in/isda/api/download.6gel.cif"
    assert path.endswith("6gel.cif")


def test_get_pdb_summary_builds_correct_url(monkeypatch):
    import pyisda.structure as structure_mod

    captured = {}

    def fake_get_json(url, timeout=10, params=None):
        captured["url"] = url
        captured["params"] = params
        return {"pdb_id": "6q0j"}

    monkeypatch.setattr(structure_mod, "get_json", fake_get_json)
    result = structure_mod.get_pdb_summary("6Q0J", additional_outputs="micromolecular_data")

    assert captured["url"] == "https://ibdc.dbt.gov.in/isda/api/pdb_summary/6q0j/"
    assert captured["params"] == {"additional_outputs": "micromolecular_data"}
    assert result == {"pdb_id": "6q0j"}


def test_get_pdb_summary_joins_list_of_additional_outputs(monkeypatch):
    import pyisda.structure as structure_mod

    captured = {}

    def fake_get_json(url, timeout=10, params=None):
        captured["params"] = params
        return {}

    monkeypatch.setattr(structure_mod, "get_json", fake_get_json)
    structure_mod.get_pdb_summary("6q0j", additional_outputs=["micromolecular_data", "ligands"])

    assert captured["params"] == {"additional_outputs": "micromolecular_data,ligands"}


def test_get_pdb_summary_returns_none_on_failure(monkeypatch):
    import pyisda.structure as structure_mod
    from pyisda._client import ISDARequestError

    def fake_get_json(url, timeout=10, params=None):
        raise ISDARequestError("boom")

    monkeypatch.setattr(structure_mod, "get_json", fake_get_json)
    assert structure_mod.get_pdb_summary("6q0j") is None


def test_get_residue_map_filters_by_auth_chain_id(monkeypatch):
    import pyisda.mapper as mapper_mod

    fake_response = {
        "structuralMapping": {
            "pdb_data": [
                {"PDB": "1M17", "UNIPROT_START": 1, "UNIPROT_END": 3, "PDB_START": 1,
                 "AUTH_CHAIN": "A", "PDB_CHAIN": "A"},
                {"PDB": "1M17", "UNIPROT_START": 1, "UNIPROT_END": 3, "PDB_START": 1,
                 "AUTH_CHAIN": "B", "PDB_CHAIN": "B"},
            ]
        }
    }

    monkeypatch.setattr(mapper_mod, "get_json", lambda url, **kw: fake_response)

    full_map = ibdc.get_residue_map("P00533", "1m17")
    assert {r["pdb_auth_chain"] for r in full_map} == {"A", "B"}

    chain_a_only = ibdc.get_residue_map("P00533", "1m17", auth_chain_id="A")
    assert {r["pdb_auth_chain"] for r in chain_a_only} == {"A"}
    assert len(chain_a_only) == 3


def test_merge_experimental_with_computational():
    experimental_df = pd.DataFrame({
        "Position": ["10", "20"],
        "Alt_AA": ["Gly", "Trp"],   # 3-letter, as returned by get_mutation_table
        "Ref_AA": ["Ala", "Ser"],
        "Clinical Significance": ["Pathogenic", "Benign"],
    })
    computational_df = pd.DataFrame({
        "Position": ["10", "20", "30"],
        "Alt_AA": ["G", "W", "A"],   # 1-letter, as returned by get_computational_mutations
        "Ref_AA": ["A", "S", "M"],
        "am_pathogenicity": [0.9, 0.1, 0.5],
    })

    merged = ibdc.merge_experimental_with_computational(experimental_df, computational_df)

    assert list(merged["Position"]) == ["10", "20"]
    assert list(merged["Alt_AA"]) == ["G", "W"]  # canonical 1-letter after normalization
    assert list(merged["am_pathogenicity"]) == [0.9, 0.1]
    assert "Ref_AA_experimental" in merged.columns
    assert "Ref_AA_computational" in merged.columns


def test_merge_experimental_with_computational_inner_join_drops_unmatched():
    experimental_df = pd.DataFrame({"Position": ["10", "99"], "Alt_AA": ["Gly", "Trp"]})
    computational_df = pd.DataFrame({"Position": ["10"], "Alt_AA": ["G"], "score": [0.5]})

    merged = ibdc.merge_experimental_with_computational(experimental_df, computational_df, how="inner")
    assert len(merged) == 1
    assert merged.iloc[0]["Position"] == "10"


def test_merge_experimental_with_computational_missing_column_raises():
    experimental_df = pd.DataFrame({"Position": ["10"]})  # missing Alt_AA
    computational_df = pd.DataFrame({"Position": ["10"], "Alt_AA": ["G"]})
    with pytest.raises(ValueError):
        ibdc.merge_experimental_with_computational(experimental_df, computational_df)


def test_subset_by_positions():
    df = pd.DataFrame({"Position": ["10", "20", "30"], "Alt_AA": ["G", "T", "A"]})
    subset = ibdc.subset_by_positions(df, [10, 30])
    assert list(subset["Position"]) == ["10", "30"]


def test_subset_by_positions_missing_column_raises():
    df = pd.DataFrame({"Alt_AA": ["G"]})
    with pytest.raises(ValueError):
        ibdc.subset_by_positions(df, [10])


def test_generate_mutated_structure_script(tmp_path):
    mt = pd.DataFrame({
        "pdb_auth_chain": ["A", "A", "B"],
        "Alt_AA": ["Gly", "Trp", "Ala"],
        "pdb_residue": [12, 45, 5],
    })
    script_path = ibdc.generate_mutated_structure_script(mt, "1abc", "A", output_dir=tmp_path)
    content = script_path.read_text()
    assert "open 1abc.cif" in content
    assert "swapaa /A:12 Gly" in content
    assert "swapaa /A:45 Trp" in content
    assert "swapaa /B:5 Ala" not in content  # filtered to chain A only
    assert content.strip().endswith("cli quit")
