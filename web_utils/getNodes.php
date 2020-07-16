<?php
/*
__description__ = "get all nodes connected to the node with the given id"
__author__ = "Aleksandar Radovanovic"
__copyright__ = "Copyright 2020, project OpenFusion"
__license__ = "MIT"
__version__ = "1.1"
__maintainer__ = "Aleksandar Radovanovic"
__email__ = "aleksandar@radovanovic.com"
__status__ = "Development"
*/

// get the id of origin node
$input = & $_GET;
if ( empty ($input) ) { $input = & $_POST; }
$id = $input['id'];

# open database and get the artilce - put your database name
$dbh = new SQLite3('aiwc.db');

# get dictionaries
$sth = $dbh->prepare('SELECT * FROM dictionary');
$result = $sth->execute();
while ( $dictionary[] = $result->fetchArray(SQLITE3_ASSOC) ) {  }

# get nodes connected to the given one
$sth = $dbh->prepare("
    SELECT id AS rec_id, tid_2 AS id, term_2 AS term, did_2 AS did, pmid_count
    FROM termpair 
    WHERE tid_1 = :id 
    UNION
    SELECT id AS rec_id, tid_1 AS id, term_1 AS term, did_1 AS did, pmid_count
    FROM termpair 
    WHERE tid_2 = :id 
    GROUP BY rec_id
    ORDER BY pmid_count DESC LIMIT 50
");
$sth->bindValue(':id', $id);
$result = $sth->execute();
$data = array();
$nodes = array(); 
$edges = array(); 

while ( $data = $result->fetchArray(SQLITE3_ASSOC) ) 
{
    $data['color'] = $dictionary[$data['did']-1]['color'];
    $data['shape'] = $dictionary[$data['did']-1]['shape'];
    $synonyms = preg_split('/\t/', $data['term']);
    $data['term'] = $synonyms[0];   // only the first term name will be used
    $nodes[] = array('data' => $data);
    $edges[] = array('data' => array('id' => "e" .$data['rec_id'], 'color' => $data['color'],'source' => $id .'', 'target' => $data['id'], 'weight' => $data['pmid_count']));
};
$graph = array(
        'nodes' =>  $nodes,
        'edges' => $edges
);
print json_encode($graph,JSON_PRETTY_PRINT);
?>
