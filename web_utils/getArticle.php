<?php
/*
__description__ = "get the specified article  article from a database in JSON format"
__author__ = "Aleksandar Radovanovic"
__copyright__ = "Copyright 2020, project OpenFusion"
__license__ = "MIT"
__version__ = "1.1"
__maintainer__ = "Aleksandar Radovanovic"
__email__ = "aleksandar@radovanovic.com"
__status__ = "Development"
*/

// get pmid
$input = & $_GET;
if ( empty ($input) ) { $input = & $_POST; }

# open database and get the artilce - put your database name
$dbh = new SQLite3('aiwc.db');

# get dictionary names
$sth = $dbh->prepare('SELECT did,dictionary FROM dictionary');
$result = $sth->execute();
while ( $dictionary[] = $result->fetchArray(SQLITE3_NUM) ) {  }
array_pop( $dictionary );

$sth = $dbh->prepare('SELECT * FROM corpus WHERE pmid=:pmid');
$sth->bindValue(':pmid', $input['pmid']);
$result = $sth->execute();
$record = $result->fetchArray(SQLITE3_ASSOC); 

// tag the article
$text = $record['article'];
$annotation = json_decode($record['annotation'], true);
usort($annotation, 'sort_by_last_index');

$term_start = PHP_INT_MAX;   // prevent tags overlaping
foreach($annotation as $term) {
    if ($term[2] < $term_start) {
        $text = mb_substr_replace($text, '</span>', $term[3]+1, 0); 
        $text = mb_substr_replace($text, '<span title="' . $dictionary[$term[0]-1][1] . '"class="did-' . $term[0] . '" id="tid-' . $term[1] . '">', $term[2], 0);
        $term_start = $term[2];
    }
}
$full_text = mb_split('\u00B6', $text);
$data = array (
    'pmid'      => $record['pmid'],
    'author'    => $record['author'],
    'title'     => $full_text[0],
    'abstract'  => $full_text[1],
    'date'      => $record['date']
);
header("Content-Type: application/json;charset=utf-8");
print json_encode($data,JSON_PRETTY_PRINT);


//-----------------------------------------------------------------------------
// reverse array of arrays sort based on the 3rd index
//-----------------------------------------------------------------------------

function sort_by_last_index($a, $b)
{
  return strnatcmp($b[3], $a[3]);
}

//-----------------------------------------------------------------------------
// utf8 substring replace
// GitHub: https://github.com/fluxbb/utf8/blob/master/functions/substr_replace.php
//-----------------------------------------------------------------------------

function mb_substr_replace($str, $repl, $start, $length = null)
{
    preg_match_all('/./us', $str, $ar);
    preg_match_all('/./us', $repl, $rar);
    $length = is_int($length) ? $length : utf8_strlen($str);
    array_splice($ar[0], $start, $length, $rar[0]);
    return implode($ar[0]);
}

?>