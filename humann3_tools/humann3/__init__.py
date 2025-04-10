# humann3_tools/humann3/__init__.py
"""HUMAnN3 processing functions."""

# Import placeholder for pathway and gene processing
# These need to be implemented
try:
    from humann3_tools.humann3.pathway_processing import process_pathway_abundance
    from humann3_tools.humann3.gene_processing import process_gene_families
    from humann3_tools.humann3.join_unstratify import process_join_unstratify
except ImportError:
    pass
