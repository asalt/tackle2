import importlib


def test_tackle2_exposes_python_package_namespace():
    legacy = importlib.import_module("python")

    import tackle2

    assert tackle2.GenesetAPI is legacy.GenesetAPI
    assert tackle2.get_geneset_collections is legacy.get_geneset_collections
    assert tackle2.geneset_membership is legacy.geneset_membership


def test_tackle2_submodules_alias_python_modules():
    import python.cli as legacy_cli

    import tackle2.cli as new_cli
    import tackle2.python as tackle2_python
    import tackle2.python.cli as nested_cli

    assert new_cli is legacy_cli
    assert nested_cli is legacy_cli
    assert tackle2_python.__package__ == "tackle2.python"
