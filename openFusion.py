#!/usr/bin/env python
'''
Fetch PubMed articles according to the given query string and create a database
'''

__author__ = "Aleksandar Radovanovic"
__copyright__ = "Copyright 2020, project covmt"
__license__ = "MIT"
__version__ = "1.0"
__maintainer__ = "Aleksandar Radovanovic"
__email__ = "aleksandar@radovanovic.com"
__status__ = "Development"

import argparse, glob, os, sqlite3, sys, getopt, dateutil.parser, re
from Bio import Entrez, Medline

def main():
    args = getArgs()
    if args.r and args.db != None and args.query != None and args.email != None:
        # check email format
        regex = '^\w+([\.-]?\w+)*@\w+([\.-]?\w+)*(\.\w{2,3})+$'
        if not re.search(regex, args.email):
            print("Invalid email address given!")
            parser.print_help()
        else:
            getPubMed(args.db, args.query, args.email)
    return

#------------------------------------------------------------------------------
# get program arguments
#------------------------------------------------------------------------------

def getArgs():
    optional = parser._action_groups.pop() 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-db", dest="db", help="SQLite database name with no file extension.",required=True)
    optional.add_argument("-r", dest="r", help="PubMed retreival required and/or",action='store_true')
    optional.add_argument("-m", dest="m", help="Text mining required",action='store_true')
    optional.add_argument("-q", dest="query", help="PubMed query string.")
    optional.add_argument("-e", dest="email", help="Your email (required by the PubMed).")
    parser._action_groups.append(optional) 
    args = parser.parse_args()
    return args

#------------------------------------------------------------------------------
# PubMed retreival
#------------------------------------------------------------------------------

def getPubMed(db, query, email):
    file = db + ".txt"
    database = db + ".db"

    # get pubmed articles vie entrez utility
    Entrez.email = email
    search_results = Entrez.read(
        Entrez.esearch(
        db="pubmed", term=query, datetype="pdat", usehistory="y"
        )
    )
    count = int(search_results["Count"])
    print("Found %i results" % count)

    batch_size = 1000
    out_handle = open(file, "w")
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

    # process pubmed articles and write them into a database
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

    with open(file) as handle:
        for record in Medline.parse(handle):
            c.execute( pubmedInsert, [record["PMID"],record["TI"],record["AB"],dateutil.parser.parse(record["DP"]).isoformat().split("T")[0]] )

    conn.commit()
    conn.close()
    return


#------------------------------------------------------------------------------
# execute the main function
#------------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='PubMed text mining utility.')
    main()