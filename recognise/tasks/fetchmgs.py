import pathlib
import shutil
import subprocess

try:
    from fetchmgs.fetchmgs import extraction_genes
except ImportError:
    FETCHMGS_LOADED = False
else:
    FETCHMGS_LOADED = True

from ..utils import read_fasta


def call_fetch_mgs(protein_file, gene_file, cog_dir, cpus, cleanup=True,):
    # fetchMG.pl -o \${bin_id}_cogs -t 5 -m extraction -d genecalls/\${bin_id}.extracted.fna genecalls/\${bin_id}.extracted.faa
    fetchmgs_cmd = [
        "fetchMGs.pl",
        "-o", cog_dir,
        "-t", f"{cpus}",
        "-m", "extraction",
        # "-x", "/usr/bin",
        "-d", f"{gene_file}",
        f"{protein_file}",
    ]

    try:
        _ = subprocess.run(
            fetchmgs_cmd,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        raise ValueError(f"fetchMGs: {e.returncode}:\n{e.output}") from e
    else:
        if cleanup:
            shutil.rmtree(cog_dir / "hmmResults")
            shutil.rmtree(cog_dir / "temp")
        
    # if fetchmgs_proc.returncode != 0:
    #     raise ValueError(f"<pre>fetchMGs error\n\n{' '.join(fetchmgs_cmd)}\n\n{fetchmgs_proc.stdout.decode()}</pre>")


def run_fetchmgs(protein_file, gene_file, outdir, cpus, cleanup=True, very_best_hits_only=True,):
    if not FETCHMGS_LOADED:
        call_fetch_mgs(protein_file, gene_file, outdir, cpus, cleanup=cleanup,)
    else:
        extraction_genes(
            [pathlib.Path(protein_file),],
            [pathlib.Path(gene_file),],
            outdir,
            cpus,
            very_best=very_best_hits_only,
        )

        marker_file = outdir / f"{pathlib.Path(protein_file).name}.fetchMGs.fna"
        if not marker_file.is_file():
            raise ValueError(f"{marker_file} does not exist.")
        
        cog_files = {}
        for seqid, seq in read_fasta(marker_file):
            p = seqid.rfind(".")
            seqid, cog = seqid[:p], seqid[p + 1:]
            cfile = cog_files.get(cog)
            if cfile is None:
                cfile = cog_files[cog] = open(outdir / f"{cog}.fna", "wt", encoding="utf-8")
            print(f">{seqid}\n{seq}\n", end="", file=cfile)
        
        for f in cog_files.values():
            f.close()
        



