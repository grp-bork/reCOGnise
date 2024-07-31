# reCOGnise workflow
<table>
  <tr width="100%">
    <td width="150px">
      <a href="https://www.bork.embl.de/"><img src="https://www.bork.embl.de/assets/img/normal_version.png" alt="Bork Group Logo" width="150px" height="auto"></a>
    </td>
    <td width="425px" align="center">
      <b>Developed by the <a href="https://www.bork.embl.de/">Bork Group</a></b><br>
      Raise an <a href="https://github.com/grp-bork/reCOGnise/issues">issue</a> or <a href="mailto:N4M@embl.de">contact us</a><br><br>
      See our <a href="https://www.bork.embl.de/services.html">other Software & Services</a>
    </td>
    <td width="250px">
      Contributors:<br>
      <ul>
        <li>
          <a href="https://github.com/cschu/">Christian Schudoma</a> <a href="https://orcid.org/0000-0003-1157-1354"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID icon" width="20px" height="20px"></a><br>
        </li>
        <li>
          <a href="https://github.com/danielpodlesny/">Daniel Podlesny</a> <a href="https://orcid.org/0000-0002-5685-0915"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID icon" width="20px" height="20px"></a><br>
        </li>
      </ul>
    </td>
    <td width="250px">
      Collaborators:<br>
      <ul>
        <li>
          <a href="https://github.com/fullama/">Anthony Fullam</a> <a href="https://orcid.org/0000-0002-0884-8124"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID icon" width="20px" height="20px"></a><br>
        </li>
        <li>
          <a href="">Askarbek Orakov</a> <a href="https://orcid.org/0000-0001-6823-5269"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID icon" width="20px" height="20px"></a><br>
        </li>
        <li>
          <a href="https://github.com/defleury">Sebastian Schmidt</a> <a href="https://orcid.org/0000-0001-8587-4177"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID icon" width="20px" height="20px"></a><br>
        </li>
      </ul>
    </td>

  </tr>
  <tr>
    <td colspan="4" align="center">The development of this workflow was supported by <a href="https://www.nfdi4microbiota.de/">NFDI4Microbiota <img src="https://github.com/user-attachments/assets/1e78f65e-9828-46c0-834c-0ed12ca9d5ed" alt="NFDI4Microbiota icon" width="20px" height="20px"></a> 
</td>
  </tr>
</table>

---
#### Description
`reCOGnise` is a tool/pipeline for species assignment (specI cluster) of microbial genomes using COG marker genes. `reCOGnise workflow` is a port of a workflow that was used e.g. for species assignment of the [proGenomes database](https://progenomes.embl.de).


---
#### Citation
This workflow: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13143009.svg)](https://doi.org/10.5281/zenodo.13143009)

Also cite:
```
Fullam A, Letunic I, Schmidt TSB, et al. proGenomes3: approaching one million accurately and consistently annotated high-quality prokaryotic genomes. Nucleic Acids Res. 2023;51(D1):D760-D766. doi:10.1093/nar/gkac1078
Coelho LP, Alves R, Del Río ÁR, et al. Towards the biogeography of prokaryotic genes. Nature. 2022;601(7892):252-256. doi:10.1038/s41586-021-04233-4
```

---
# Overview
![reCOGnise Workflow Diagram](https://raw.githubusercontent.com/grp-bork/reCOGnise/main/docs/reCOGnise_workflow.svg)

---
# Requirements

`reCOGnise` requires a docker/singularity installation. All dependencies are contained in the `reCOGnise` docker container.

Dependencies are

* `prodigal`
* `fetchMGS.pl`
* `MAPseq`

---
# Usage
## Cloud-based Workflow Manager (CloWM)
This workflow will be available on the CloWM platform (coming soon).

## Command-Line Interface (CLI)

You can either clone this repository from GitHub and run it as follows
```
git clone https://github.com/grp-bork/reCOGnise.git
nextflow run /path/to/reCOGnise --input_dir /path/to/genome/fastas --output_dir /path/to/output_dir
```

Input genome fasta files have to have one of the following file endings: `{fna,fasta,fa,fna.gz,fasta.gz,fa.gz}`. Alternatively, you can set the pattern with
`params.file_pattern = "**.{<comma-separated-list-of-file-endings>}"`.


Or, you can have nextflow pull it from github and run it from the `$HOME/.nextflow` directory.
```
nextflow run grp-bork/reCOGnise --input_dir /path/to/genome_files --output_dir path/to/output_dir
```