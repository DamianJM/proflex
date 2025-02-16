# ProFlex Protein Flexibility Alphabet
**ProFlex Alphabet for Protein Flexibility Description**

<img src="https://raw.githubusercontent.com/DamianJM/proflex/main/img/PF_logo.jpg" width="300" height="300"></img> 

# PUBLICATIONS

[Magill, Damian, and Timofey Skvortsov. "Decoding Protein Dynamics: ProFlex as a Linguistic Bridge in Normal Mode Analysis." bioRxiv (2024): 2024-09.](https://www.biorxiv.org/content/10.1101/2024.09.21.614246v1)

# TABLE OF CONTENTS

- [INSTALLATION](#INSTALLATION)
 
- [USAGE](#USAGE)
    - [Querying ProFlex Databases](#Queries)
    - [Creating ProFlex Databases](#Databases)
    - [Precompiled Database Access](#Precompiled-Database-Access)
    - [ProFlex Translation](#Translation)

# INSTALLATION

### In order to install the proflex toolkit simply type:
```bash
pip install proflex
```

Installation should take no longer than a minute.

If problems occur during installation or specifically with structural comparisons this is almost certainly due to pymol2 installation issues. In those cases, please proceed to pymol installation via freely available wheels by following these instructions: https://github.com/cgohlke/pymol-open-source-wheels?tab=readme-ov-file

# USAGE

### Queries

### Upon installation, the toolkit can be imported and integrated into various workflows. To query a PDB against a proflex database requires only three lines of code:

```python
from proflex import ProFlexQuery as pf
pq = pf("/path/to/database")
pq.query_pdb("input.pdb")
```

This provides a report output in HTML format similar to the below:

<img src="https://raw.githubusercontent.com/DamianJM/proflex/main/img/FiguresProFlex.jpg" width="600" height="800"></img> 

Structural and sequence alignments are provided along with top ProFlex query hits and a backtranslated RMSF value comparison with that of the top query hit.

As part of this repo we provide a precompile SWISS-PROT proflex database that can be downloaded. To prepare this database for full use, PDB structures need to be downloaded to the PDB subdirectory like so:
```bash
wget https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/swissprot_pdb_v4.tar
```
Simply unpack all structures into the PDB subdirectory of the database. The PDB filenames are referenced by the database and can be retrieved when queries are performed.

### Databases

### Given a multifasta of proflex sequences, databases can be created like so:
```python
from proflex import NGramDatabase
new_database = NGramDatabase()
new_database.parse_multifasta("/path/to/multifasta.fasta")
new_database.save_to_directory("/path/to/output/directory")
```
The newly created database can now be queried. In the case of custom databases be sure to provide the mutlifasta sequence file and PDB files in the same directory to allow for full functionality.

### Precompiled Database Access

In order to obtain the pre-compiled SWISS-PROT database you will need to use github LFS when cloning this repo. You can install this via the following link: https://git-lfs.com/

Once installed, in order to obtain the database file from the pointer you can do the following:

```bash
git clone https://github.com/DamianJM/proflex.git
cd proflex
git lfs pull
```

### Translation

Other methods exist for backtranslation etc and are fully available in the installed package. We provide the empirically defined percentiles from our study in the source code which can be easily accessed as follows:

```python
from proflex import ProFlex as pf
print(pf.PERCENTILES)
```

Given a set of RMSF values you can derived ProFlex sequences as follows:

```python
from proflex import ProFlex as pf
proflex = pf()
rmsf = [21, 1.2, 10, 1, 7] # input RMSF values
print(proflex.encode_sequence(rmsf_values=rmsf)) # Output ProFlex
```

## A ProFlex online application for direct PDB submission is under development and will be made available soon!
