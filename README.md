# ProFlex Protein Flexibility Alphabet
**ProFlex Alphabet for Protein Flexibility Description**

### In order to install the proflex toolkit simply type:
```bash
pip install proflex
```
If problems occur during installation or specifically with structural comparisons this is almost certainly due to pymol2 installation issues. In those cases, please proceed to pymol installation via freely available wheels by following these instructions: https://github.com/cgohlke/pymol-open-source-wheels?tab=readme-ov-file

### Upon installation, the toolkit can be imported and integrated into various workflows. To query a PDB against a proflex database requires only three lines of code:

```python
import proflex as pf
pq = pf.ProFlexQuery("/path/to/database")
pq.query_pdb("input.pdb")
```
As part of this repo we provide a precompile SWISS-PROT proflex database that can be downloaded. To prepare this database for full use, PDB structures need to be downloaded to the PDB subdirectory like so:
```bash
wget https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/swissprot_pdb_v4.tar
```
Simply unpack all structures into the PDB subdirectory of the database. The PDB filenames are referenced by the database and can be retrieved when queries are performed.

### Given a multifasta of proflex sequences, databases can be created like so:
```python
import proflex as pf
new_database = pf.NGramDatabase()
new_database.parse_multifasta("/path/to/multifasta.fasta")
new_database.save_to_directory("/path/to/output/directory")
```
The newly created database can now be queried. In the case of custom databases be sure to provide the mutlifasta sequence file and PDB files in the same directory to allow for full functionality.

Other methods exist for backtranslation etc and are fully available in the installed package. We provide the empirically defined percentiles from our study in the source code which can be easily accessed as follows:

```python
import proflex as pf
print(pf.ProFlex.PERCENTILES)
```





