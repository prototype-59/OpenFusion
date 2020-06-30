#!/usr/bin/env python
'''
Fetch PubMed articles according to the given query string and create a database
'''

__author__ = "Aleksandar Radovanovic"
__copyright__ = "Copyright 2020, project covmt"
__license__ = "MIT"
__version__ = "1.1"
__maintainer__ = "Aleksandar Radovanovic"
__email__ = "aleksandar@radovanovic.com"
__status__ = "Development"

import argparse, glob, os, sqlite3, sys, getopt, dateutil.parser, re
from Bio import Entrez, Medline

def main():
    args = getArgs()
    if args.f != None :
        processPubmed(args.db, args.f)
    elif args.query != None and args.email != None :
        regex = '^\w+([\.-]?\w+)*@\w+([\.-]?\w+)*(\.\w{2,3})+$' # check email format
        if not re.search(regex, args.email):
            print("Invalid email address given!")
            sys.exit()
        else:
            getPubMed(args.db, args.query, args.email)
            processPubmed(args.db, args.db + ".txt")
    else :
        print("Invallid combination of arguments. Use either:\n   -db db_name -f PubMed_filename or \n  -db db_name -e your_email -q \"PubMed query string\"")
    return

#------------------------------------------------------------------------------
# get program arguments
#------------------------------------------------------------------------------

def getArgs():
    optional = parser._action_groups.pop() 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-db", dest="db", help="SQLite database name, e.g. myPubMed.db",required=True)
    optional.add_argument("-f", dest="f", help="File containing PubMed abstracts")
    optional.add_argument("-q", dest="query", help="PubMed query string, e.g. \"Alpers Disease\".")
    optional.add_argument("-e", dest="email", help="Your email (required by PubMed).")
    parser._action_groups.append(optional)
    args = parser.parse_args()
    return args

#------------------------------------------------------------------------------
# convert PubMed dates into ISO standard (e.g. 1969 Sep-Oct  to yyyy-mm-dd)
#------------------------------------------------------------------------------

def date_to_isodate(str) :
    month = ['none','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    date = str.split(" ")
    str = date[0]
    if len(date) > 1 :
        date[1] = date[1].replace("Fall", "Sep").replace("Summer", "Jun").replace("Winter", "Dec").replace("Spring", "Mar")
        str += "-" + f"{month.index(date[1][0:3]):02d}"
        if len(date) > 2 :
            regex = '^[0-9]+$'
            if re.search(regex, date[2]):
                str += "-" + f"{int(date[2]):02d}"
    return str


#------------------------------------------------------------------------------
# PubMed retreival
#------------------------------------------------------------------------------

def getPubMed(db, query, email):
    file = db + ".txt"

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
    return

def processPubmed(db, file) :
    # process pubmed articles and write them into a database
    pubmedTable = """
        CREATE TABLE IF NOT EXISTS corpus (
	        id INTEGER NOT NULL PRIMARY KEY UNIQUE,
            author TEXT,
            title TEXT,
            abstract TEXT,
            body TEXT,
            date TEXT,
            annotation TEXT
    );"""
    pubmedInsert = "INSERT INTO corpus  (id, author, title, abstract, date)  VALUES (?,?,?,?,?)"

    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute(pubmedTable)
    
    with open(file) as handle:
        for record in Medline.parse(handle):
            # csv list of authors's full names (fau)
            fau = ""
            ab = ""
            if "FAU" in record :
                fau = "|".join(map(str, record["FAU"])).replace(',','').replace('|',', ')
            if "AB" in record :
                ab = record["AB"]
            date = date_to_isodate(record["DP"])
            c.execute( pubmedInsert, [record["PMID"], fau, record["TI"],ab,date] )


    conn.commit()
    conn.close()
    return


#------------------------------------------------------------------------------
# execute the main function
#------------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get PubMed articles.')
    main()

