from typing import Union, Tuple
import pandas as pd

from .blast import pairwise_blasts
from .blast.fasta import fasta_translation_table, translate_blast_table

def similarity_groups(species: list,
    fasta_files: list,
    fasta_type: str = 'prot',
    fasta_header_key: Union[str, None] = None
) -> Tuple[pd.DataFrame, pd.DataFrame]:

    # Blast every FASTA file against each other
    # Returns (query_species, db_species, result_df)
    blast_results = pairwise_blasts(
        species,
        fasta_files,
        fasta_type=fasta_type
    )

    if fasta_header_key is not None:
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

    # Use high-stringency similarity on within-species
    # blast results to group genes initially
    self_results = _build_joint_sets(
        pd.concat([
            _process_self_blast(z)
            for q, d, z in blast_results
            if q == d
        ])
    )

    # Then use lower stringency for on between-species
    # blast results to group genes
    all_results = _build_joint_sets(
        _process_between_blast(
            pd.concat([
                z
                for q, d, z in blast_results
                if q != d
            ]),
            species
        )#,
        #add_to_list=self_results
    )

    # Turn each group of genes into a separate dataframe
    all_results = [
        pd.DataFrame(list(group), columns=["Species", "Gene"])
        for group in all_results
    ]

    # Add a unique group ID
    for i, r in enumerate(all_results):
        r['GroupID'] = i

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

    all_gene_results = pd.concat(
        [all_results] + all_gene_results
    ).reset_index(
        drop=True
    )

    return all_gene_results, all_results

def _process_self_blast(
    blast_results,
    threshold=98
):

    bsr = blast_results.loc[
        blast_results['qseqid'] != blast_results['sseqid'],
        :
    ]

    return bsr.loc[
        bsr['true_ident'] > threshold,
        :
    ].copy()


def _process_between_blast(
    blast_results,
    species_list,
    threshold=50
):

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

def _build_joint_sets(bd, add_to_list=None):
    """
    Build a list of sets, where each set is highly similar genes,
    from a blast dataframe

    :param bd: Preprocessed blast data
    :type bd: pd.DataFrame
    :param add_to_list: A list of sets to add to, defaults to None
    :type add_to_list: list(set()), optional
    :return: A list of sets, where sets are highly similar genes
    :rtype: list(set())
    """

    valid_genes = {}
    result_dict = {}

    if add_to_list is not None:
        for s in add_to_list:
            for g in s:
                result_dict[g] = s
                valid_genes[g] = True

    # Make dict of sets where each set is group of homologs
    # and the elements of the set all point to the same set
    # Stash the sets in a list for easy access later
    for qgene, sgene, qspecies, sspecies in zip(
        bd['qseqid'],
        bd['sseqid'],
        bd['qspecies'],
        bd['sspecies']
    ):

        i, j = (qspecies, qgene), (sspecies, sgene)

        joint = set([i, j])

        try:
            result_dict[i].update(joint)
        except KeyError:
            result_dict[i] = joint

        try:
            result_dict[j].update(joint)
        except KeyError:
            result_dict[j] = joint

        result_dict[i].update(result_dict[j])
        result_dict[j] = result_dict[i]

        # Make sure result_dict for all the elements in this set are fused
        for k in result_dict[i].copy():
            try:
                if id(result_dict[i]) != id(result_dict[k]):
                    result_dict[i].update(result_dict[k])
            except KeyError:
                pass

            result_dict[k] = result_dict[i]

        valid_genes[i] = True
        valid_genes[j] = True

    # One last unification pass
    for i in list(result_dict.keys()):
        for k in result_dict[i].copy():

            try:
                if id(result_dict[i]) != id(result_dict[k]):
                    result_dict[i].update(result_dict[k])
            except KeyError:
                pass

            result_dict[k] = result_dict[i]

    # Make a list of unique sets from the duplicate references
    # Use a dict with ids to keep track of what we've seen before
    result_list = []
    set_list = {}
    for _, x in result_dict.items():
        try:
            if set_list[id(x)]:
                continue
        except KeyError:
            result_list.append(x)
            set_list[id(x)] = True

    return result_list
