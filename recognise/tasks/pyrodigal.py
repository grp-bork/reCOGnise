# pylint: disable=E0401
""" Module to add pyrodigal gene calling functionality """

import pyrodigal

from ..utils import read_fasta


def run_pyrodigal(f, genome_id, output_dir):
    """ Call genes with pyrodigal. """
    gf = pyrodigal.GeneFinder(mask=True)

    ids, seqs = zip(*read_fasta(f))
    _ = gf.train(*seqs)

    output_dir.mkdir(exist_ok=True, parents=True,)

    faa_out = open(output_dir / f"{genome_id}.faa", "wt", encoding="UTF-8",)
    ffn_out = open(output_dir / f"{genome_id}.ffn", "wt", encoding="UTF-8",)
    gff_out = open(output_dir / f"{genome_id}.gff", "wt", encoding="UTF-8",)

    with faa_out, ffn_out, gff_out:
        for sid, seq in zip(ids, seqs):
            sid = sid[:sid.find(" ")]
            genes = gf.find_genes(seq)
            genes.write_translations(faa_out, sid)
            genes.write_genes(ffn_out, sid)
            genes.write_gff(gff_out, sid, full_id=False,)
