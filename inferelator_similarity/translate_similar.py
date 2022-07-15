from .utils import SPECIES_COL, GROUP_COL, GENE_COL, group_stats

def similar_translations(
    genes,
    similarity_table,
    translate_from,
    translate_to,
    verbose=False
):

    # Subset to species of interest
    tf = similarity_table[SPECIES_COL] == translate_from
    tt = similarity_table[SPECIES_COL] == translate_to

    similarity_table = similarity_table[tf | tt].copy()

    if verbose:
        print(
            f"Mapping {translate_from} [{tf.sum()}] to "
            f"{translate_to} [{tt.sum()}]"
        )

    # Find 1-1 mappings
    _mappable_groups = group_stats(
        similarity_table
    )

    _mappable_groups = _mappable_groups[
        (_mappable_groups == 2).all(axis=1)
    ]

    if verbose:
        print(
            f"{_mappable_groups.sum()} 1 to 1 matches for mapping"
        )

    similarity_table = similarity_table[
        similarity_table[GROUP_COL].isin(_mappable_groups.index)
    ]

    # Join the perfect matches on GroupID
    # To make a translation table
    return similarity_table[
        similarity_table[SPECIES_COL] == translate_from
    ].set_index(
        GROUP_COL,
        drop=True
    ).join(
        similarity_table[
            similarity_table[SPECIES_COL] == translate_to
        ].set_index(
            GROUP_COL,
            drop=True
        ),
        rsuffix = "_to"
    ).set_index(
        GENE_COL,
        drop=True
    )[GENE_COL + "_to"].reindex(
        genes
    ).dropna()
