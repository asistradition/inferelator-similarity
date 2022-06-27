import pandas as pd
import re

DEFAULT_EXTRACTORS = {
    'id': '>(.*?) ',
    'gene': 'gene:(.*?) ',
    'transcript': 'transcript:(.*?) ',
    'gene_biotype': 'gene_biotype:(.*?) '
}

RAW = 'FASTA'

def parse_fasta_headers(fasta_file, extractors=None):

    with open(fasta_file) as fh:
        df = pd.DataFrame(
            [line.strip() for line in fh
             if len(line) > 0 and line.startswith(">")],
             columns=[RAW]
        )

    extractors = extractors if extractors is not None else DEFAULT_EXTRACTORS

    for col, regex in extractors.items():
        df[col] = df[RAW].str.extract(
            regex,
            expand=False,
            flags=re.IGNORECASE
        )

    return df.drop([RAW], axis=1)

def fasta_translation_table(fasta_files, fasta_key):

    return pd.concat([
        parse_fasta_headers(x) for x in fasta_files
    ]).set_index('id', drop=True)[[fasta_key]]

def translate_blast_table(
    blast_results,
    translation_table,
    cols_to_translate,
    inplace=True
):

    if not inplace:
        blast_results = blast_results.copy()

    _joined_col = translation_table.columns.values[0]
    for c in cols_to_translate:

        col = blast_results[[c]].join(
            translation_table,
            on=c,
            how='left'
        ).fillna("")

        blast_results[c] = col[_joined_col]

    return blast_results
