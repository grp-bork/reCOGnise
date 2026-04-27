import subprocess


def call_prodigal(genome, protein_file, gene_file, gff_file=None):
    # prodigal -i \$(basename ${genome_fna} .gz) -f gff -o ${genome_id}/${genome_id}.gff -a ${genome_id}/${genome_id}.faa -d ${genome_id}/${genome_id}.ffn
    # prodigal -i \$(basename \$genome_file .gz) -f gff -o prodigal/\$genome_id/\$genome_id.gff -a prodigal/\$genome_id/\$genome_id.faa -d prodigal/\$genome_id/\$genome_id.ffn

    gff_params = ["-f", "gff", "-o", gff_file,] if gff_file is not None else []

    try:
        _ = subprocess.run(
            [
                "prodigal",
                "-i",
                genome,
                "-a",
                protein_file,
                "-d",
                gene_file,
            ] + gff_params,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            check=True, capture_output=True,
        )
    except subprocess.CalledProcessError as e:
        raise ValueError(f"prodigal: {e.returncode}:\n{e.output}") from e

    # if prodigal_proc.returncode != 0:
    #     raise ValueError(f"<pre>Prodigal error\n\n{prodigal_proc.stdout}</pre>")
