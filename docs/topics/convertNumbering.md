**convertNumbering** - *convertNumbering: IMGT-Kabat number conversion*

Description
--------------------

Converts numbering systems like Kabat or IMGT using these conventions:
http://www.imgt.org/IMGTScientificChart/Numbering/IMGT-Kabat_part1.html
with Gaps (unoccupied positions) shown by "G" and Asterisks (*) shown by "S": 
arbitrary mappings (multiple possible "to" values) represented with "NA"


Usage
--------------------
```
convertNumbering(locus, from, to, calls)
```

Arguments
-------------------

locus
:   string indicating heavy ("IGH") or light chains ("IGK" or "IGL)

from
:   string indicating numbering system to convert to ("IMGT" or "KABAT")

to
:   string indicating original numbering system ("IMGT" or "KABAT")

calls
:   vector of strings representing original numbering




Value
-------------------

A vector of string indicating the corresponding numbering



Examples
-------------------

```R
convertNumbering("IGH", "IMGT", "KABAT", c("51", "23", "110"))

```


```
[1] "46" "22" "G" 

```


```R
convertNumbering("IGH", "KABAT", "IMGT", c("51", "23", "G"))

```


```
[1] "56" "24" "NA"

```








