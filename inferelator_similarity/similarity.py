from typing import Union, Tuple
import pandas as pd
import argparse
import collections

from .blast import pairwise_blasts, BLAST_HEADER, IDENT_COL
from .blast.fasta import fasta_translation_table, translate_blast_table
from .utils import build_joint_sets

def similarity_groups_execute():
    """
    Execute similarity groups using command line arguments
    """
    ap = argparse.ArgumentParser(
        description="Pairwise BLAST of multiple FASTA files"
    )

    ap.add_argument(
        "-f",
        "--fasta",
        dest="fasta",
        help="FASTA File Paths",
        metavar="FILE",
        required=True,
        nargs="+"
    )

    ap.add_argument(
        "-s",
        "--species",
        dest="species",
        help="Species names",
        required=True,
        nargs="+"
    )

    ap.add_argument(
        "--prot",
        dest="prot",
        action='store_const',
        help="FASTA files are Protein",
        const=True,
        default=False
    )

    ap.add_argument(
        "--verbose",
        dest="verbose",
        action='store_const',
        help="Verbose flag",
        const=True,
        default=False
    )

    ap.add_argument(
        "--fasta_header_key",
        dest="header",
        help="Gene ID field in fasta header",
        default=None
    )

    ap.add_argument(
        "--out",
        dest="out",
        help="Output TSV File Name",
        default="./similar_gene_groups.tsv"
    )

    ap.add_argument(
        "--raw_tsv",
        dest="raw",
        help="Raw BLAST Output TSV File Name",
        default=None
    )

    ap.add_argument(
        "--identical_tsv",
        dest="identical",
        help="Identical Gene BLAST Output TSV File Name",
        default=None
    )

    args = ap.parse_args()

    all_gene_results, identicals, process_results = similarity_groups(
        args.species,
        args.fasta,
        fasta_type="prot" if args.prot else "nucl",
        fasta_header_key=args.header,
        verbose=args.verbose,
        identical_feature_file=args.identical
    )

    all_gene_results.to_csv(
        args.out,
        sep="\t",
        index=False
    )

    if args.raw is not None:
        print(
            f"Writing {len(process_results)} BLAST results"
            f"to {args.raw}"
        )

        process_results.to_csv(
            args.raw,
            sep="\t",
            index=False
        )

    if args.identical is not None:

        identicals.to_csv(
            args.identical,
            sep="\t",
            index=False
        )


def similarity_groups(species: list,
    fasta_files: list,
    fasta_type: str = 'prot',
    fasta_header_key: Union[str, None] = None,
    verbose: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Group genes from multiple species by sequence similarity

    :param species: List of species names or IDs
    :type species: list
    :param fasta_files: List of paths to species fasta files
    :type fasta_files: list
    :param fasta_type: Fasta type, 'nucl' or 'prot',
        defaults to 'prot'
    :type fasta_type: str, optional
    :param fasta_header_key: Use different field in FASTA header,
        defaults to None
    :type fasta_header_key: Union[str, None], optional
    :param verbose: Print verbose messages, defaults to False
    :type verbose: bool, optional
    :return: 3-column dataframes with 'Gene', 'Species', and 'GroupID'
        First dataframe returned has all genes
        Second is blast result tables
    :rtype: Tuple[pd.DataFrame, pd.DataFrame]
    """

    if verbose:
        print("Pairwise BLAST:")
        print("\n".join(
            f"{x}:\t{y}" for x, y in zip(species, fasta_files)
        ))

    _processed_results = []

    # Blast every FASTA file against each other
    # Returns (query_species, db_species, result_df)
    blast_results = pairwise_blasts(
        species,
        fasta_files,
        fasta_type=fasta_type
    )

    if fasta_header_key is not None:

        if verbose:
            print(f"Translating to {fasta_header_key} names")

        translator = fasta_translation_table(
            fasta_files,
            fasta_header_key
        )

        _tcols = ['qseqid', 'sseqid'] + species

        for _, _, df in blast_results:
            translate_blast_table(
                df,
                translator,
                cols_to_translate=_tcols,
                inplace=True
            )

    # Get the gene names from within-species blast results
    # And store in a dict
    species_gene_names = {
        q: pd.Index(z['qseqid'].unique())
        for q, d, z in blast_results
        if q == d
    }

    _processed_results.append(
        pd.concat([
            _process_self_blast(z)
            for q, d, z in blast_results
            if q == d
        ])
    )

    identical_features = pd.concat([
        _process_self_blast(z, 100)
        for q, d, z in blast_results
        if q == d
    ])[BLAST_HEADER + [IDENT_COL, 'qspecies', 'sspecies']].copy()

    # Use high-stringency similarity on within-species
    # blast results to group genes initially
    self_results = build_joint_sets(
        _processed_results[-1]
    )

    if verbose:
        _n_merge = sum(len(x) - 1 for x in self_results if len(x) > 1)
        print(f"{_n_merge} genes are >98% identical within species "
              "and are premerged")

    _processed_results.append(
        pd.concat([
            z
            for q, d, z in blast_results
            if q != d
        ])
    )

    # Then use lower stringency for on between-species
    # blast results to group genes
    all_results = build_joint_sets(
        _process_between_blast(
            _processed_results[-1],
            species
        ),
        add_to_list=self_results
    )

    # Turn each group of genes into a separate dataframe
    all_results = [
        pd.DataFrame(list(group), columns=["Species", "Gene"])
        for group in all_results
    ]

    # Add a unique group ID
    for i, r in enumerate(all_results):
        r['GroupID'] = i

    if verbose:
        _group_size = collections.Counter(map(len, all_results))

        print("Similarity grouping complete:")
        for i, j in _group_size.most_common(5):
            print(f"\t{i} Gene Groups: {j} [{i * j} Total Genes]")

    # Concat all the group dataframes into a single dataframe
    all_results = pd.concat(
        all_results
    ).reset_index(
        drop=True
    )

    all_gene_results = []
    group_id_start = all_results['GroupID'].max() + 1

    for sp, sp_gene_names in species_gene_names.items():
        _add_df = pd.DataFrame(
            sp_gene_names[
                ~sp_gene_names.isin(all_results['Gene'])
            ],
            columns=['Gene']
        )
        _add_df['Species'] = sp
        _add_df['GroupID'] = range(
            group_id_start,
            group_id_start + _add_df.shape[0]
        )
        group_id_start += _add_df.shape[0]
        all_gene_results.append(_add_df)

    if verbose:
        _mono_genes = sum(map(len, all_gene_results))
        print(f"{_mono_genes} genes have no similarity")

    all_gene_results = pd.concat(
        [all_results] + all_gene_results
    ).reset_index(
        drop=True
    )

    if verbose:
        _perfect_line = all_gene_results.groupby(
            'GroupID'
        )['Species'].agg(
            ['count', 'nunique']
        ) == len(species)

        _perfect_line = _perfect_line.all(axis=1).sum()

        print(
            f"{_perfect_line} groups "
            f"[{_perfect_line * len(species)} genes] "
            f"are {len(species)}-way perfect matches"
        )

    return all_gene_results, identical_features, pd.concat(_processed_results)

def _process_self_blast(
    blast_results,
    threshold=98
):

    bsr = blast_results.loc[
        blast_results['qseqid'] != blast_results['sseqid'],
        :
    ]

    return bsr.loc[
        bsr['true_ident'] >= threshold,
        :
    ].copy()


def _process_between_blast(
    blast_results: pd.DataFrame,
    species_list: list,
    threshold: Union[int, float] = 50
) -> pd.DataFrame:
    """
    Process a dataframe of blast results to find only
    the highest similarity from one gene to another

    :param blast_results: Blast results
    :type blast_results: pd.DataFrame
    :param species_list: List of species columns
    :type species_list: list
    :param threshold: Minimum threshold for pairwise similarity
        in percent, defaults to 50
    :type threshold: Union[int, float], optional
    :return: Blast results with only the highest hit for each gene
    :rtype: pd.DataFrame
    """

    # Filter all hits below threshold true identity
    # Then take the highest true identity as the most similar
    # for each gene
    bbr = blast_results.loc[
        blast_results['true_ident'] > threshold,
        :
    ].copy().sort_values(
        by=['true_ident'], ascending=False
    ).drop_duplicates(
        subset=['qseqid', 'sseqid']
    )

    _max_result_idx = bbr.groupby(
        species_list
    )['true_ident'].transform(
        max
    ) == bbr['true_ident']

    bbr = bbr[
        _max_result_idx
    ].drop_duplicates()

    return bbr


if __name__ == '__main__':
    similarity_groups_execute()
