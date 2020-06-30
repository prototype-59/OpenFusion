#!/usr/bin/env python
'''
Populate a dictionary table in a given database
'''

__author__ = "Aleksandar Radovanovic"
__copyright__ = "Copyright 2020, project covmt"
__license__ = "MIT"
__version__ = "1.1"
__maintainer__ = "Aleksandar Radovanovic"
__email__ = "aleksandar@radovanovic.com"
__status__ = "Development"

import argparse, glob, os, sqlite3, sys, getopt


#------------------------------------------------------------------------------
# get program arguments
#------------------------------------------------------------------------------

def getArgs():
  required = parser.add_argument_group('required arguments')
  required.add_argument("-name", dest="dname", help="dictionary name, e.g. \"human genes\"",required=True)
  required.add_argument("-color", dest="color", help="dictionary hex color, e.g. \"#d6eaf8\"",required=True)
  required.add_argument("-db", dest="db", help="SQLite database name, e.g. myPubMed.db",required=True)
  required.add_argument("-f", dest="f", help="file containing list of terms as tab separated synonyms",required=True)
  args = parser.parse_args()
  return args

#------------------------------------------------------------------------------
# main function
#------------------------------------------------------------------------------
def main():
  args = getArgs()

  # open/create database and create dictionary table if not exits
  conn = sqlite3.connect(args.db)
  c = conn.cursor()
  sql = """
    CREATE TABLE IF NOT EXISTS dictionary (
      did INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
      dictionary TEXT,
      color TEXT
  );
  """
  c.execute(sql)
  sql = """
      CREATE TABLE IF NOT EXISTS term (
      tid INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
      did INTEGER,
      term TEXT
  );
  """
  c.execute(sql)
  conn.commit()

  # check if dictionary already exists and retreive did
  c.execute("SELECT * FROM dictionary WHERE dictionary=?", [args.dname] )
  dictionary = c.fetchone()
  if dictionary == None:
    c.execute("INSERT INTO dictionary (dictionary,color) VALUES (?,?)", [args.dname, args.color])
    conn.commit()
    c.execute("SELECT last_insert_rowid()")
    dictionary = c.fetchone()

  did = dictionary[0]
    
  # insert dictionary terms
  data = []
  with open(args.f, mode="r", encoding="utf-8" ) as f:
    for line in f:
      data.append([did,line.strip()])

  # write terms
  c.executemany("INSERT INTO term (did, term) VALUES (?,?)", data)
  conn.commit()
  return

#------------------------------------------------------------------------------
# execute the main function
#------------------------------------------------------------------------------
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Create a dictionary of terms')
  main()
