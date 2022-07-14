import pandas as pd
import anndata as ad
from typing import Union

from .blast import BLAST_HEADER
from .blast.fasta import fasta_translation_table, translate_blast_table
from .utils import build_joint_sets


BLAST_IDS = BLAST_HEADER[0:2]
MERGE_COL = 'merge_to'
CHROMOSOME_KEY = 'chromosome'

def collapse_genes(
    expression_data: Union[pd.DataFrame, ad.AnnData],
    collapse_table: pd.DataFrame,
    fasta_files: Union[list, None] = None,
    fasta_header_key: Union[str, None] = None,
    agg_func: Union[str, callable] = 'sum',
    axis: int = 0,
    verbose: bool = False
):

    if verbose:
        vprint = print
    else:
        vprint = lambda *x: 0

    # Check and see if it's an expression dataframe
    # Or an AnnData object
    _is_ad = isinstance(expression_data, ad.AnnData)
    if _is_ad:
        _data = expression_data.to_df()
    elif _is_ad and axis != 0:
        raise ValueError("Cannot merge AnnData on axis 1")
    else:
        _data = expression_data.copy()

    # Get a table of chromosome IDs if FASTA files are provided
    # And swap qseqid and sseqid for anything that's got no
    # chromosome
    if fasta_files is not None:

        collapse_table = collapse_table.join(
            fasta_translation_table(
                fasta_files,
                CHROMOSOME_KEY
            ),
            on=BLAST_HEADER[0],
            how='left'
        )

        _swap_idx = collapse_table[CHROMOSOME_KEY].isna()
        _swap_idx |= collapse_table[CHROMOSOME_KEY] == ""

        collapse_table.loc[
            _swap_idx, BLAST_IDS
        ] = collapse_table.loc[_swap_idx, BLAST_IDS[::-1]]

    # Translate based on fasta header
    if fasta_files is not None and fasta_header_key is not None:

        vprint(f"Translating BLAST table to FASTA {fasta_header_key}")
        translate_blast_table(
            collapse_table,
            fasta_translation_table(
                fasta_files,
                fasta_header_key
            ),
            cols_to_translate=BLAST_IDS,
            inplace=True
        )

    # Assign the first qseqid found in a joint set as the to-merge
    # There's probably an easier way to do this

    _collapse_line = build_joint_sets(
        collapse_table,
        include_species=False
    )

    _all_genes = pd.concat((
        collapse_table[BLAST_HEADER[0]],
        collapse_table[BLAST_HEADER[1]]
    )).drop_duplicates()

    _collapse_line = {
        _all_genes[_all_genes.isin(s)].iloc[0]: s
        for s in _collapse_line
        if len(s) > 1
    }

    _collapse_line = pd.DataFrame.from_dict(
        {from_gene: to_gene
        for to_gene, s in _collapse_line.items()
        for from_gene in s},
        orient='index',
        columns=[MERGE_COL]
    )

    if axis == 0 or axis == 'all':

        _data = _agg_data(
            _data,
            agg_func,
            _collapse_line.reindex(
                _data.columns
            )[MERGE_COL].fillna(
                _data.columns.to_series()
            ),
            axis=0
        )

        _data.columns.name = None

    if axis == 1 or axis == 'all':

        _data = _agg_data(
            _data,
            agg_func,
            _collapse_line.reindex(
                _data.index
            )[MERGE_COL].fillna(
                _data.index.to_series()
            ),
            axis=1
        )

        _data.index.name = None

    vprint(f"Merged {expression_data.shape} data to {_data.shape}")

    if isinstance(expression_data, ad.AnnData):

        expression_data = ad.AnnData(
            _data,
            obs=expression_data.obs,
            var=expression_data.var.reindex(_data.columns),
            obsm=expression_data.obsm,
            obsp=expression_data.obsp,
            dtype=expression_data.dtype
        )

    else:

        expression_data = _data

    return expression_data, _collapse_line


def _agg_data(
    data,
    agg_func,
    grouper,
    axis=0
):

    if agg_func == "nonzero" and axis == 1:
        return _nonzero_data(data, grouper)
    elif agg_func == "nonzero" and axis == 0:
        return _nonzero_data(data.T, grouper).T

    elif axis == 1:
        return data.groupby(grouper).agg(agg_func)
    elif axis == 0:
        return data.T.groupby(grouper).agg(agg_func).T


def _nonzero_data(
    data,
    grouper
):

    data = (data != 0).groupby(
        grouper
    ).agg(
        'sum'
    ) > 0

    return data.astype(int)
