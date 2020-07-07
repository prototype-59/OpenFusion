### Alice in Wonderland syndrome

<p align="center">
<img src="../img/alice.jpg" height="200"/>
</p>
Alice in Wonderland syndrome (AWS) is a rare condition that causes episodes of distorted perception. You may feel larger or smaller than you actually are.

#### Creating Alice in Wonderland syndrome database

1. Thsi project contains two dictionary files with tab-separated list of synonymous terms: `disease.txt` and `people.txt`. These are terms the OpenFusion will be searching for in the text.
   <br>

2. The project file, e.g. `Alice.yml` contains project details. Please note that text following the # is a comment and ignored by the program.
```yml
---
database: aiwc.db       # name of your SQLite database
email: opefusionx@gmail.com   # email required for PubMed
query: Alice in Wonderland Syndrome   # PubMed query string
dictionary:             # create glossary of terms
- name: People
  file: people.txt
- name: Disease
  file: disease.txt
homographs: False  # True/False for homographs anotations
run:               # what modules to run
- getPubMedCorpus  # get articles online from PubMed
#- getFileCorpus   # get articles from a file
- addDictionary    # create glossary of terms
- annotate         # annotate PubMed articles
```
- `database` section specifies the database that will be created as the output (`aiwc.db`) as well as the name of the text file in MEDLINE format containing articles downloaded from PubMed (`aiwc.db.txt`).
- `email` section is required by PubMed (https://pubmed.ncbi.nlm.nih.gov/). Put your own email address (no PubMed account is needed).
- `query` contains your search terms as keywords or as a PubMed search expression.
- `dictionary` section list your dictionaries. `name` specifies the name you wish to give to your dictionary, and `file` is the file name containing dictionary terms.
- `homographs` section can be set to `True` or `False`. If `True`, terms with the same spelling but with diferent meanings which are listed in different dictionaries, will be annotated. If `False`, only the first term found in the dictionary list  will be annotaded. `True` can slows program significantly.
- `run` section specifies what modules you wish to run. The first time you will run three modules as shown in the listing. This will create the database `aiwc.db` and the file `aiwc.db.txt`. If you wish to add more dictionaries or edit existing ones, you can re-run program as many times you wish (without querying the PubMed) by opting-out the `getPubMedCorpus` option and enalbling the `getFileCorpus`. This option will re-create the database from the file `aiwc.db.txt`.
3. Run:

       ./OpenFusion.py -p Alice.yml

   The program will create SQLite database specified in the project file. You can open this database with the DB Browser for SQLite (https://sqlitebrowser.org). You can use this database in your biomedical projects, machine learning research etc.

    
#### Citation
Get it here:  [![DOI](https://zenodo.org/badge/248162501.svg)](https://zenodo.org/badge/latestdoi/248162501)