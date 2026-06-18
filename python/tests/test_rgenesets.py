import pandas as pd

from tackle2 import get_geneset_collections, geneset_membership


def _hallmark_spec():
    return [{"category": "H", "subcategory": "", "collapse": False}]


def test_get_geneset_collections_returns_dataframes():
    collections = get_geneset_collections(_hallmark_spec(), species="Homo sapiens")
    assert "H_" in collections
    hallmark = collections["H_"]
    assert isinstance(hallmark, pd.DataFrame)
    assert not hallmark.empty
    assert {"gs_name", "gs_collection", "gs_subcollection"}.issubset(hallmark.columns)


def test_geneset_membership_collapses_genes():
    collections = get_geneset_collections(_hallmark_spec())
    hallmark = collections["H_"]
    subset = hallmark[hallmark["gs_name"] == hallmark["gs_name"].iloc[0]]
    collapsed = geneset_membership(subset)
    assert collapsed
    genes = next(iter(collapsed.values()))
    assert isinstance(genes, list)
    assert genes
