import tempfile
import itertools
import tqdm
import re
import pandas as pd

from ._blast_wrapper import makeblastdb, BLAST_HEADER, blastn, blastp

def pairwise_blasts(
    species,
    fasta_files,
    fasta_type='prot',
    blast_self=False
):

    if fasta_type == 'prot':
        blaster = blastp
    elif fasta_type == 'nucl':
        blaster = blastn
    else:
        raise ValueError(
            "fasta_type must be 'nucl' or 'prot'; "
            f"{fasta_type} passed"
        )

    with tempfile.TemporaryDirectory() as db_holder:

        ### MAKE BLAST DATABASES ###
        dbs = [
            makeblastdb(
                ff,
                sp,
                dbtype=fasta_type,
                out_path=db_holder
            )
            for sp, ff in zip(species, fasta_files)
        ]

        _total = len(species) ** 2

        ### PAIRWISE BLAST ###
        blast_results = [
            (q, d, _process_blast(
                blaster(ff, db),
                q,
                d,
                all_species=species
            ))
            for (q, ff), (d, db) in tqdm.tqdm(
                itertools.product(
                    zip(species, fasta_files),
                    zip(species, dbs)
                ),
                total = _total
            )
        ]

        return blast_results

def _strip_prefs(x):
    rep = x.str.extract(
        '\|(.*?)\|',
        expand=False,
        flags=re.IGNORECASE
    )

    _na_idx = pd.isna(rep)
    rep.loc[_na_idx] = x.loc[_na_idx]

    return rep

def _process_blast(
    data,
    qspecies,
    sspecies,
    strip_prefixes=True,
    all_species=None,
):
    """
    Take a blast dataframe and reprocess it to have species
    information in an easily parsed format

    :param data: Blast dataframe (outfmt 6)
    :type data: pd.DataFrame
    :param qspecies: Query species
    :type qspecies: str
    :param sspecies: Database species
    :type sspecies: str
    :param strip_prefixes: Remove database prefixes, keeping
        only what is between bars (|GENE|), defaults to True
    :type strip_prefixes: bool, optional
    :return: Reprocessed dataframe
    :rtype: pd.DataFrame
    """
    if strip_prefixes:
        data['qseqid'] = _strip_prefs(data['qseqid'])
        data['sseqid'] = _strip_prefs(data['sseqid'])

    data['qspecies'] = qspecies
    data['sspecies'] = sspecies
    data[qspecies] = data['qseqid']
    data[sspecies] = data['sseqid']

    data['true_ident'] = data['pident'] * data['length'] / data['qlen']

    if all_species is not None:
        for sp in all_species:
            if sp != qspecies and sp != sspecies:
                data[sp] = ""

        _process_cols = ['qspecies', 'sspecies', 'true_ident']

        data = data.reindex(
            BLAST_HEADER + _process_cols + all_species,
            axis=1
        )

    return data
