"""Top-level Python API for tackle2.

Historically the code shipped under the ``python`` package.  This shim keeps
those modules importable while exposing the public API through ``tackle2`` so
``pip install .`` provides a predictable import path:

.. code-block:: python

    from tackle2 import get_geneset_collections
    collections = get_geneset_collections([...])

Every submodule under ``python`` remains available as ``tackle2.<module>`` and
``tackle2.python`` is an alias of the original package for compatibility.
"""

from __future__ import annotations

import importlib
import importlib.abc
import importlib.machinery
import importlib.util
import sys
import types
from pathlib import Path

_python_pkg = importlib.import_module("python")

# Re-export the documented public API.
from python import GenesetAPI, get_geneset_collections, geneset_membership

__all__ = getattr(_python_pkg, "__all__", [])

# Ensure ``import tackle2.<module>`` searches the original package directory.
_pkg_dir = Path(__file__).resolve().parent
__path__ = [str(_pkg_dir)]

# Alias the historical package so ``tackle2.python`` resolves as expected
# without double-importing submodules.


class _PythonPackageShim(types.ModuleType):
    def __init__(self, target):
        super().__init__(name=f"{__name__}.python")
        self._target = target
        self.__doc__ = target.__doc__
        self.__all__ = getattr(target, "__all__", None)
        self.__package__ = f"{__name__}.python"
        self.__path__ = []
        self.__spec__ = importlib.machinery.ModuleSpec(
            name=self.__package__, loader=None, is_package=True
        )
        self.__spec__.submodule_search_locations = []

    def __getattr__(self, attr):
        return getattr(self._target, attr)

    def __dir__(self):
        return dir(self._target)


sys.modules[f"{__name__}.python"] = _PythonPackageShim(_python_pkg)


def _module_exists_locally(qualname: str) -> bool:
    if not qualname:
        return False
    parts = qualname.split(".")
    candidate = _pkg_dir.joinpath(*parts)
    file_candidate = candidate.with_suffix(".py")
    if file_candidate.exists():
        return True
    if candidate.exists() and candidate.joinpath("__init__.py").exists():
        return True
    return False


class _PythonAliasLoader(importlib.abc.Loader):
    def __init__(self, alias_name: str, target_name: str):
        self.alias_name = alias_name
        self.target_name = target_name

    def create_module(self, spec):
        module = importlib.import_module(self.target_name)
        sys.modules[self.alias_name] = module
        return module

    def exec_module(self, module):
        sys.modules[self.alias_name] = module


class _PythonAliasFinder(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        prefix = __name__ + "."
        if not fullname.startswith(prefix):
            return None
        suffix = fullname[len(prefix) :]
        if not suffix:
            return None
        if suffix == "python":
            return None
        if suffix.startswith("python."):
            legacy_suffix = suffix[len("python.") :]
        else:
            legacy_suffix = suffix
        if _module_exists_locally(suffix):
            return None

        legacy_name = f"python.{legacy_suffix}"
        legacy_spec = importlib.util.find_spec(legacy_name)
        if legacy_spec is None:
            return None

        loader = _PythonAliasLoader(fullname, legacy_name)
        alias_spec = importlib.machinery.ModuleSpec(
            fullname,
            loader,
            origin=legacy_spec.origin,
            is_package=legacy_spec.submodule_search_locations is not None,
        )
        if legacy_spec.submodule_search_locations is not None:
            alias_spec.submodule_search_locations = legacy_spec.submodule_search_locations
        return alias_spec
_alias_finder = _PythonAliasFinder()
if not any(isinstance(finder, _PythonAliasFinder) for finder in sys.meta_path):
    sys.meta_path.append(_alias_finder)


def __getattr__(name: str):
    """Fall back to attributes defined on the legacy ``python`` package."""

    if hasattr(_python_pkg, name):
        return getattr(_python_pkg, name)
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


def __dir__():
    return sorted(set(globals()) | set(dir(_python_pkg)))
