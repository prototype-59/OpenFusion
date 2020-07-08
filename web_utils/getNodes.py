#!/usr/bin/env python

__description__ = "get all nodes connected to the node with the given id"
__author__ = "Aleksandar Radovanovic"
__copyright__ = "Copyright 2020, project OpenFusion"
__license__ = "MIT"
__version__ = "1.1"
__maintainer__ = "Aleksandar Radovanovic"
__email__ = "aleksandar@radovanovic.com"
__status__ = "Development"


import cgi, cgitb, sqlite3, json

# get the id of origin node
form = cgi.FieldStorage() 
id = form.getvalue('id')

# connect to the database  - put your database name
conn = sqlite3.connect("aiwc.db")
conn.row_factory = sqlite3.Row
c = conn.cursor()

# get dictionary names
c.execute('SELECT * FROM dictionary')
dictionary = c.fetchall()

# get nodes connected to the given one

nodes = []
edges = []
c.execute("""
    SELECT id AS rec_id, tid_2 AS id, term_2 AS term, did_2 AS did, count
    FROM termpair  WHERE tid_1 = ? ORDER BY count DESC LIMIT 50"""
    ,[id])
records = c.fetchall()
for rec in records:
    synonyms = rec['term'].split('\t')
    nodes.append({
        'data': {
            'rec_id': int(rec['rec_id']),
            'id': int(rec['id']),
            'term':synonyms[0],
            'did': int(rec['did']),
            'count': int(rec['count']),
            'color': dictionary[int(rec['did'])-1]['color'],
            'shape': dictionary[int(rec['did'])-1]['shape'],
        }
    })
    edges.append({
        'data': {
            'id': "e" + str(rec['rec_id']),
            'color': dictionary[int(rec['did'])-1]['color'],
            'source': id,
            'target': rec['id'],
            'weight': rec['count']
        }
    })
data = {'nodes': nodes, 'edges': edges}
print "Content-type: application/json\r\n\r\n"
print (json.dumps(data))
