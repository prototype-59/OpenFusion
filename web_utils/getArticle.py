#!/usr/bin/python

__description__ = "get the specified article  article from a database in JSON format"
__author__ = "Aleksandar Radovanovic"
__copyright__ = "Copyright 2020, project OpenFusion"
__license__ = "MIT"
__version__ = "1.1"
__maintainer__ = "Aleksandar Radovanovic"
__email__ = "aleksandar@radovanovic.com"
__status__ = "Development"

import cgi, cgitb, sqlite3, json

# get pmid
form = cgi.FieldStorage() 
pmid = form.getvalue('pmid')

# connect to the database  - put your database name
conn = sqlite3.connect("aiwc.db")
c = conn.cursor()

# get dictionary names
c.execute('SELECT did,dictionary FROM dictionary')
dictionary = c.fetchall()

# get article
c.execute('SELECT * FROM corpus WHERE pmid=?',[pmid])
article = c.fetchone()
text = article[2]
annotation = json.loads(article[4])
annotation.sort(reverse = True, key=lambda x: x[3])
term_start = float('inf')   # prevent tags overlaping
for term in annotation :
    if term[2] < term_start :
        text = text[:term[3]+1] + '</span>' + text[term[3]+1:]
        text = text[:term[2]] + '<span title="' + dictionary[term[0]-1][1] + '" class="did-' + str(term[0]) + '" id="tid-' + str(term[1]) + '">' + text[term[2]:]
        term_start = term[2]

full_text = text.split(u'\u00B6')
data = {
    'pmid': article[0],
    'author':  article[1].encode('utf8'),
    'title' : full_text[0].encode('utf8'),
    'abstract' : full_text[1].encode('utf8'),
    'date': article[3]
}

print "Content-type: application/json\r\n\r\n"
print json.dumps(data)
