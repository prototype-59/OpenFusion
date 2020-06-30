# OpenFusion
## Text Analytics Tools

OpenFusion is a text analytics tool for creating **SQLite** biomedical databases on specific topic by using **PubMed** as a source. This  database can be used locally or on the web.

1. retreive PubMed articles with `getPubmed.py`
2. create dictionaries with `addDictionary.py`
3. annotate text with `annotate.py`
4. see the result in SQLite database `myDatabase.db`

### Processing steps

1. Retreive PubMed articles related to your topic with the `getPubmed.py` tool. For example, to create a database on Alpers Disease, type

        ./getPubmed.py -db my.db -e myemail@email.com -q "Alpers Disease"

    `-db my.db`: database name (if not exists, it will be crearted)\
    `-e myemail@email.com`: your email required by PubMed\
    `-q "Alpers Disease"`: PubMed query string\
    In addition to database, this tool creates a text file `my.db.txt` containing all retreived articles.<br>
    If you already have a file and want to create a new database, type

       ./getPubmed.py-db db_name -f PubMed_filename
    <br>
2. Create files containing dictionary terms, e.g. genes.txt, diseases.txt, etc. Each line must contain a word or phrase, which can be followed by tab-separated synonyms, e.g.

       Alpers Disease       Progressive Cerebral Poliodystrophy
    
    Add these dictionaries, one by one, with the `addDictionary.py` tool to your database created in the first step. For example,

       ./addDictionary.py -db tm.db -name "Human Genes" -color "#d6eaf8"  -f genes.txt
    `-db my.db`: database name\
    `-name "Human Genes"`: dictionary name\
    `-color "#d6eaf8"`: dictionary display color\
    `-f genes.txt`: name of the file with dictionary terms.


    
### Citation
Get it here:  [![DOI](https://zenodo.org/badge/248162501.svg)](https://zenodo.org/badge/latestdoi/248162501)