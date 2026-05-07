""" Utility functions """

import gzip
import itertools as it


GZIP_MAGIC = b"\x1f\x8b\x08"

def read_fasta(f):
    """ Read fasta file. """

    with open(f, "rb") as _in:
        is_gzipped = _in.read(3) == GZIP_MAGIC

    open_f = gzip.open if is_gzipped else open

    with open_f(f, "rt", encoding="utf-8",) as _in:
        gen = it.groupby(_in, key=lambda x:x[0]==">",)
        for is_header, item in gen:
            if is_header:
                seqid = next(item).rstrip()[1:]
                _, seqs = next(gen)
                yield seqid, "".join(s.rstrip() for s in seqs)
