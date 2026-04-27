import itertools as it


def read_fasta(f):
    with open(f) as _in:
        gen = it.groupby(_in, key=lambda x:x[0]==">",)
        for is_header, item in gen:
            if is_header:
                seqid = next(item).rstrip()[1:]
                _, seqs = next(gen)
                yield seqid, "".join(s.rstrip() for s in seqs)
