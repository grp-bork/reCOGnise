Output
======

`reCOGnise` produces the following outputs **if** a species can be assigned to an input genome:

### Prodigal Annotation
* `<sample>.faa` -- predicted proteins
* `<sample>.ffn` -- predicted genes
* `<sample>.gff` -- gene annotation

### COG Marker Genes
* `<sample>.cogs.txt` -- marker gene table

### Species Assignment
* `<sample>.specI.status.OK` -- text file signifying if species assignment was succesful
* `<sample>.specI.status` -- text file with status of the species assignment
* `<sample>.specI.txt` -- species (specI cluster) assignment
