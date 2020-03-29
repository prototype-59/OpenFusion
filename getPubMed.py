#!/usr/bin/env python
'''
Fetch PubMed articles according to the given query string and create a database tm.db
'''

__author__ = "Aleksandar Radovanovic"
__copyright__ = "Copyright 2020, project covmt"
__license__ = "MIT"
__version__ = "1.0"
__maintainer__ = "Aleksandar Radovanovic"
__email__ = "aleksandar@radovanovic.com"
__status__ = "Development"

import glob, os, sqlite3, sys, getopt, dateutil.parser
from Bio import Entrez, Medline

# get input argument: search string
if len(sys.argv) != 2:
    print( 'program usage: ' + sys.argv[0] +  '"PubMed query string"')
    sys.exit()
query = sys.argv[1]

# get pubmed articles vie entrez utility
Entrez.email = "aleksandar.radovanovic@kaust.edu.sa"
search_results = Entrez.read(
    Entrez.esearch(
        db="pubmed", term=query, datetype="pdat", usehistory="y"
    )
)
count = int(search_results["Count"])
print("Found %i results" % count)

batch_size = 1000
out_handle = open("pubmed_result.txt", "w")
for start in range(0, count, batch_size):
    end = min(count, start + batch_size)
    print("Downloading record %i to %i" % (start + 1, end))
    fetch_handle = Entrez.efetch(
        db="pubmed",
        rettype="medline",
        retmode="text",
        retstart=start,
        retmax=batch_size,
        webenv=search_results["WebEnv"],
        query_key=search_results["QueryKey"],
    )
    data = fetch_handle.read()
    fetch_handle.close()
    out_handle.write(data)
out_handle.close()

# process pubmed articles and write them into database
database = "tm.db"
pubmedTable = """
CREATE TABLE pubmed (
	pmid INTEGER NOT NULL PRIMARY KEY UNIQUE,
    title TEXT,
    abstract TEXT,
    date TEXT
);"""
pubmedInsert = "INSERT INTO pubmed  (pmid, title, abstract, date)  VALUES (?,?,?,?)"

# remove old and create a new database
if os.path.exists(database):
    os.remove(database)
conn = sqlite3.connect(database)
c = conn.cursor()
c.execute(pubmedTable)

with open("pubmed_result.txt") as handle:
    for record in Medline.parse(handle):
        c.execute( pubmedInsert, [record["PMID"],record["TI"],record["AB"],dateutil.parser.parse(record["DP"]).isoformat().split("T")[0]] )

conn.commit()
conn.close()
