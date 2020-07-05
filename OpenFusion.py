#!/usr/bin/env python

__author__ = "Aleksandar Radovanovic"
__copyright__ = "Copyright 2020, project OpenFusion"
__license__ = "MIT"
__version__ = "1.1"
__maintainer__ = "Aleksandar Radovanovic"
__email__ = "aleksandar@radovanovic.com"
__status__ = "Development"

import os, sys, argparse, sqlite3, json, yaml, re
from Bio import Entrez, Medline

# some constants
WIN_SIZE = 500  # text window size
_end = '_end_'  # trie structure related variables
_id = 0


#------------------------------------------------------------------------------
# convert PubMed dates into ISO standard (e.g. 1969 Sep-Oct  to yyyy-mm-dd)
#------------------------------------------------------------------------------

def date_to_isodate(str) :
    str = str.replace("-", " ").replace("/", " ").replace(".", " ")
    month = {'Fall':'9','Autumn':'9','Aut':'9','Winter':'12','Spring':'3','Summer':'6', 'Jan':'1','Feb':'2','Mar':'3','Apr':'4','May':'5','Jun':'6','Jul':'7','Aug':'8','Sep':'9','Oct':'10','Nov':'11','Dec':'12'}
    date = str.split(" ")

    if date[0].isdigit() :  
        str = date[0]
    else :
        str = date[1]       # sometimes month/season comes before the year
        date[1] = date[0]

    if len(date) > 1 :
        if date[1].isdigit() :
            str += "-" + f"{int(date[1]):02d}"
        else:
            date[1] = re.sub('({})'.format('|'.join(map(re.escape, month.keys()))), lambda m: month[m.group()], date[1])
            str += "-" + f"{int(date[1]):02d}"
        if len(date) > 2 and date[2].isdigit() :
            str += "-" + f"{int(date[2]):02d}"
    return str

#------------------------------------------------------------------------------
# PubMed retreival to a file
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

#------------------------------------------------------------------------------
# create a database from a file retreieved from PubMed
#------------------------------------------------------------------------------

def processPubmedFile(db) :
    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute("DROP TABLE IF EXISTS corpus")
    conn.commit()

    sql = """
        CREATE TABLE IF NOT EXISTS corpus (
	        pmid INTEGER NOT NULL PRIMARY KEY UNIQUE,
            author TEXT,
            article TEXT,
            date TEXT,
            annotation JSON
    );"""
    c.execute(sql)
    conn.commit()
  
    file = file = db + ".txt"
    with open(file) as handle:
        for record in Medline.parse(handle):
            if "TI" not in record : # skip articles with no title
                continue

            fau = ""         # csv list of authors's full names (fau)
            article = record["TI"]
            # some articles are in between [] - remove it
            if record["TI"].startswith('[') and record["TI"].endswith('].') :
                article = article[1:-2] + '.'
            if "FAU" in record :
                fau = "|".join(map(str, record["FAU"])).replace(',','').replace('|',', ')
            if "AB" in record :
                article += '\u00B6' + record["AB"]
            date = date_to_isodate(record["DP"])
            c.execute("INSERT INTO corpus  (pmid, author, article, date)  VALUES (?,?,?,?)", [record["PMID"], fau, article, date])
    conn.commit()
    conn.close()
    return

#------------------------------------------------------------------------------
# drop dictionaries if exists
#------------------------------------------------------------------------------

def dropDictionary(db) :
    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute("DROP TABLE IF EXISTS dictionary")
    c.execute("DROP TABLE IF EXISTS glossary")
    conn.commit()
    return

#------------------------------------------------------------------------------
# drop and add dictionaries
#------------------------------------------------------------------------------

def addDictionary(db, file, name) :
    # 10 colors and shapes for dictionaries
    colors = ["#ffc75f", "#ff9671", "#ff8066", "#fbeaff", "#f9f871", "#e4f0f5", "#1b4f72", "#d5cabd", "#845ec2", "#66095d"]
    shapes = ["ellipse", "star", "round-triangle", "rectangle", "round-pentagon", "vee", "tag", "barrel", "rhomboid", "round-diamond"]

    conn = sqlite3.connect(db)
    c = conn.cursor()

    sql = """
        CREATE TABLE IF NOT EXISTS dictionary (
        did INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
        dictionary TEXT,
        color TEXT,
        shape TEXT
    );
    """
    c.execute(sql)
    sql = """
        CREATE TABLE IF NOT EXISTS glossary (
        tid INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
        did INTEGER,
        term TEXT
    );
    """
    c.execute(sql)
    conn.commit()

    c.execute("INSERT INTO dictionary (dictionary) VALUES (?)", [name])
    conn.commit()
    c.execute("SELECT last_insert_rowid()")
    dictionary = c.fetchone()
    
    # dictionary id (did), color and shape (useful fo graphs e.g. cytoscape.js)
    did = dictionary[0]
    color = colors[did % 10]
    shape = shapes[did % 10]
    c.execute("UPDATE dictionary SET color =?, shape = ? WHERE did = ? ", [color, shape, did])

    # insert dictionary terms
    data = []
    with open(file, mode="r", encoding="utf-8" ) as f:
        for line in f:
            data.append([did,line.strip()])

    # write terms
    c.executemany("INSERT INTO glossary (did, term) VALUES (?,?)", data)
    conn.commit()
    return

#------------------------------------------------------------------------------
# annotation functions
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Build dictionay trie
#------------------------------------------------------------------------------

def make_trie(dictionary):
    root = dict()
    for word in dictionary:
        current_dict = root
        for letter in word[1]:
            current_dict = current_dict.setdefault(letter, {})
        current_dict[_end] = _end
        current_dict[_id] = word[0]
    return root

#------------------------------------------------------------------------------
# Check if a given term is in a dictionary:
# False: if not found, term_id: if found
#------------------------------------------------------------------------------

def in_trie(trie, word):
    current_dict = trie
    for letter in word:
        if letter not in current_dict:
            return False
        current_dict = current_dict[letter]
    #return _end in current_dict
    if _end in current_dict:
        return current_dict[_id]
    else:
        return False

#------------------------------------------------------------------------------
# annotate text in database
# Each term found is recorded as: pmid, did, tid, idx, term_synonym
#------------------------------------------------------------------------------

def annotate(db):
    if not os.path.exists(db):
        print("Database " + db + " does not exists!")
        sys.exit()
    conn = sqlite3.connect(db)
    c = conn.cursor()
    
    # reset previous annotation
    c.execute("DROP TABLE IF EXISTS annotation")
    conn.commit()
    c.execute('UPDATE corpus SET annotation=NULL')
    conn.commit()

    sql = """
        CREATE TABLE annotation (
        id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
        pmid INTEGER,
        did INTEGER, 
        tid INTEGER,
        idx INTEGER,
        term TEXT
    );
    """
    c.execute(sql)
    conn.commit()

    # loop dictionaries 
    c.execute('SELECT did FROM dictionary')
    dictionary_list = [i[0] for i in c.fetchall()]
    
    # annotate for each of the dictionaries
    for did in dictionary_list :

        # get terms in this dictionary
        c.execute('SELECT * FROM glossary WHERE did =?', [did])
        
        # get dictionary terms
        dictionary = []
        for row in c :
            tid = row[0]
            terms = row[2].split('\t') # synonyms
            for term in terms :
                dictionary.append([tid, term])
        
        # contruct the dictionary trie
        trie = make_trie(dictionary)

        # Annotating
        data = []
        punctuation = {'.\u00B6': '  ', '\u00B6': ' ', '? ': '  ', ', ':'  ', '. ':'  ', '; ': '  '}
        c.execute('SELECT pmid, article FROM corpus')
        for record in c:
            # prepare text by replacing special characters and punctuations
            article = record[1] + ' '
            article = re.sub('({})'.format('|'.join(map(re.escape, punctuation.keys()))), lambda m: punctuation[m.group()], article)
            end_idx = len(article) - 1   # index of the last character in article

            # loop the artilce by making text window of 100 (the longest phrase)
            window_start = 0
            window_end = WIN_SIZE + 1 if end_idx > WIN_SIZE else end_idx + 1

            # find the whole word from the right
            while article[window_end:window_end+1] != ' ' and window_end > 2 :
                window_end -= 1

            while window_start < end_idx  :
                while window_end > window_start + 1:
                    text = article[window_start:window_end]
                    term_id = in_trie(trie, text)
                    if term_id :
                        data.append([record[0], did, term_id, window_start, text])
                    window_end -= 1

                # shift the window to the right
                while article[window_start:window_start+1] != ' ' and window_start < end_idx :
                    window_start += 1
                window_start += 1

                window_end = window_start + WIN_SIZE + 1
                if window_end > end_idx : window_end = end_idx + 1
                while article[window_end:window_end+1] != ' ' and window_end > 2 :
                    window_end -= 1
    
        c.executemany("INSERT INTO annotation (pmid, did, tid, idx, term) VALUES (?,?,?,?,?)", data)
        conn.commit()

    conn.close()
    return

#------------------------------------------------------------------------------
# update corpus with articles annotation 
#------------------------------------------------------------------------------

def corpusUpdate(db) :
    conn = sqlite3.connect(db)
    c = conn.cursor()

    # get article list
    c.execute('SELECT pmid FROM corpus')
    pmid_list = [i[0] for i in c.fetchall()]
    
    # update corpus
    for pmid in pmid_list :
        c.execute('SELECT did,tid,idx,term FROM annotation WHERE pmid = ?', [pmid])
        data = c.fetchall()
        if data :
            json_data = json.dumps(data)
            c.execute('UPDATE corpus SET annotation =? WHERE pmid = ?', [json_data, pmid])
    conn.commit()
    conn.close()
    return

# term to pmid_list mappings
def term_to_pmid(db) :
    conn = sqlite3.connect(db)
    c = conn.cursor()

    c.execute("DROP TABLE IF EXISTS term")
    conn.commit()
    sql = """
        CREATE TABLE IF NOT EXISTS term (
        tid INTEGER NOT NULL PRIMARY KEY UNIQUE,
        did INTEGER,
        term TEXT,
        count INTEGER,
        pmid_list TEXT
    );
    """
    c.execute(sql)
    conn.commit()

    sql = """
        INSERT INTO term (tid,did,term,count,pmid_list)
        SELECT T1.tid,T2.did,T2.term, COUNT(T1.pmid) AS count,GROUP_CONCAT(T1.pmid,',') AS pmid_list
        FROM (SELECT DISTINCT tid,pmid FROM annotation) T1
        INNER JOIN glossary T2 ON T1.tid = T2 .tid 
        GROUP BY T1.tid
    """
    c.execute(sql)
    conn.commit()
    return

#------------------------------------------------------------------------------
# term pair to pmid_list mappings
#------------------------------------------------------------------------------

def termpair_to_pmid(db) :
    conn = sqlite3.connect(db)
    c = conn.cursor()

    c.execute("DROP TABLE IF EXISTS termpair")
    conn.commit()
    sql = """
        CREATE TABLE IF NOT EXISTS termpair (
        id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
        tid_1 INTEGER,
        did_1 INTEGER,
        term_1 TEXT,
        tid_2 INTEGER,
        did_2 INTEGER,
        term_2 TEXT,
        count INTEGER,
        pmid_list TEXT
    );
    """
    c.execute(sql)
    conn.commit()

    # loop throught pmid_list and find common pmid(s)
    c.execute("SELECT * FROM term WHERE count > 1")
    termsA = c.fetchall()
    termsB = termsA

    data = []
    for termA in termsA:
        pmid_listA = list(termA[4].split(","))
        for termB in termsB:
            pmid_listB = list(termB[4].split(","))
            common = list(set(pmid_listA).intersection(pmid_listB))
            if common:
                if termA[0] < termB[0] :
                    data.append([termA[0], termA[1], termA[2], termB[0], termB[1], termB[2],len(common), ','.join(common)])
                else:
                    data.append([termB[0], termB[1], termB[2], termA[0], termA[1], termA[2],len(common), ','.join(common)])

    conn.executemany("INSERT INTO termpair (tid_1,did_1,term_1,tid_2,did_2,term_2,count,pmid_list) VALUES (?,?,?,?,?,?,?,?)", data)
    conn.commit()

    # cleanup termpair by deleting self-referencing and duplicate rows
    conn.execute("DELETE FROM termpair WHERE tid_1=tid_2")
    conn.execute("DELETE FROM termpair WHERE rowid NOT IN (SELECT MIN(rowid) FROM termpair GROUP BY tid_1,tid_2)")
    conn.commit()
    return

#------------------------------------------------------------------------------
# the main function
#------------------------------------------------------------------------------

def main():
    required = parser.add_argument_group('required arguments')
    required.add_argument("-p", dest="p", help="YAML configuration file.",required=True)    
    args = parser.parse_args()
    config = yaml.safe_load(open(args.p))
    
    if "getPubMedCorpus" in config['run'] :
        getPubMed(config['database'],config['query'],config['email'])
        processPubmedFile(config['database'])

    if "getFileCorpus" in config['run'] :
        processPubmedFile(config['database'])

    if "addDictionary" in config['run'] :
        dropDictionary(config['database'])
        for dictionary in config['dictionary'] :
            addDictionary(config['database'],dictionary['file'],dictionary['name'])
    
    if "annotate" in config['run'] :
        annotate(config['database'])
        corpusUpdate(config['database'])
        term_to_pmid(config['database'])
        termpair_to_pmid(config['database'])

    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='OpenFusion text mining utility.')
    main()