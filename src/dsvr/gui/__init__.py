"""GUI package for inspecting DSVR run outputs.

The modules in this package are split into streamlit-free logic
(``inventory``, ``lineage``, ``tables``, ``anomalies``) that can be imported
and tested without the optional ``gui`` dependency, and the Streamlit views
under ``ui`` which import ``streamlit`` lazily.
"""

from dsvr.gui.inventory import RunInventory
from dsvr.gui.lineage import VariantLineage, parse_variant_id

__all__ = ["RunInventory", "VariantLineage", "parse_variant_id"]
