import subprocess
import io
import os
import joblib
import pandas as pd

BLAST_HEADER = [
    "qseqid",
    "sseqid",
    "pident",
    "qlen",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore"
]

BLAST_OUTFMT = " ".join([
    "6",
    "qseqid",
    "sseqid",
    "pident",
    "qlen",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore"
    ]
)

BLASTP_DEFAULTS = {
    "-evalue": "1e-10",
    "-max_target_seqs": "3"
}

def blastp(
    query_file,
    target_db,
    num_threads=-1,
    **kwargs
):
    return _blaster(
        'blastp',
        query_file,
        target_db,
        num_threads,
        **kwargs
    )

def blastn(
    query_file,
    target_db,
    num_threads=-1,
    **kwargs
):
    return _blaster(
        'blastn',
        query_file,
        target_db,
        num_threads,
        **kwargs
    )


def _blaster(
    blast_type,
    query_file,
    target_db,
    num_threads,
    **kwargs
):

    if num_threads == -1:
        num_threads = joblib.cpu_count()

    blast_cmd = [
        blast_type,
        '-query', query_file,
        '-db', target_db,
        '-num_threads', str(num_threads),
        '-outfmt', BLAST_OUTFMT,
    ]

    for cmd, val in kwargs.items():
        blast_cmd.extend([
            _dash(cmd), str(val)
        ])

    for cmd, val in BLASTP_DEFAULTS.items():
        if cmd not in blast_cmd:
            blast_cmd.extend([
                cmd, str(val)
            ])

    proc = subprocess.run(
        blast_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='utf-8',
    )

    if int(proc.returncode) != 0:
        print(proc.stderr)
        raise RuntimeError(f"{blast_type} failed: {proc.args}")

    return pd.read_csv(
        io.StringIO(proc.stdout),
        names=BLAST_HEADER,
        sep="\t"
    )

def makeblastdb(fasta, title, out_path=".", dbtype='prot', **kwargs):

    out_path = os.path.join(out_path, title)

    makeblastdb_cmd = [
        "makeblastdb",
        "-parse_seqids",
        "-dbtype", dbtype,
        "-title", title,
        "-out", out_path,
        "-in", fasta
    ]

    for cmd, val in kwargs.items():
        makeblastdb_cmd.extend([
            _dash(cmd), str(val)
        ])

    proc = subprocess.run(
        makeblastdb_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='utf-8',
    )

    if int(proc.returncode) != 0:
        print(proc.stderr)
        raise RuntimeError(f"makeblastdb failed: {proc.args}")

    return out_path

def _dash(x):
    if len(x) == 0:
        raise ValueError("string is length 0")

    if x.startswith("--"):
        return x[1:]
    elif x.startswith("-"):
        return x
    else:
        return '-' + x
