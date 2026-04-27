import subprocess


def call_fetch_mgs(protein_file, gene_file, cog_dir, cpus):
    # fetchMG.pl -o \${bin_id}_cogs -t 5 -m extraction -d genecalls/\${bin_id}.extracted.fna genecalls/\${bin_id}.extracted.faa
    fetchmgs_cmd = [
        "fetchMGs.pl",
        "-o", cog_dir,
        "-t", f"{cpus}",
        "-m", "extraction",
        "-x", "/usr/bin",
        "-d", f"{gene_file}",
        f"{protein_file}",
    ]		

    try:
        _ = subprocess.run(
            fetchmgs_cmd,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            check=True, capture_output=True,   
        )
    except subprocess.CalledProcessError as e:
        raise ValueError(f"fetchMGs: {e.returncode}:\n{e.output}") from e
        
    # if fetchmgs_proc.returncode != 0:
    #     raise ValueError(f"<pre>fetchMGs error\n\n{' '.join(fetchmgs_cmd)}\n\n{fetchmgs_proc.stdout.decode()}</pre>")
