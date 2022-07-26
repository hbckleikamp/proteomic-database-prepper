# proteomic-database-prepper
this is a small auxiliary python script that helps with preparing proteomic databases for annotation.
It includes options to:

- Equate I and J to L
- Remove Ambiguous amino acids
- Add decoy sequences (reversed or scrambled)
- Select only bacterial and archaeal sequences 


The parsed NCBI taxonomy file needed to select only bacterial and archaeal sequences is constructed with https://github.com/hbckleikamp/NCBI2Lineage.
