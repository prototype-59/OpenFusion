### Web utilities

This directory contains number of utilities written in Python and PHP that you can use to publish your database online. 

- `getArticle.py` and `getArticle.php` is the same utility that displays selected article on the web page. These utilities are called via GET or POST method by supplying the article pmid. As illustrated in file `example.html`, you can make an asynchronous Ajax call as:

      getArticle.php?pmid=10767914 

  Annotations are highlighted as:
  ```html
  <span title="Disease" class="did-2" id="tid-4">Alice in Wonderland</span>
  ```
  where `did-2` is the dictionary number and `tid-4` is the term number in the database. Colors for each of the dictionaries are defined is `dictionaries.css` stylesheet file. By using jquery you can attach various events to dictionary classes and to term ids.