# proqueryote

A tool for selecting and downloading prokaryote proteomic and genomic sequence files from NCBI.

## Origin

proqueryote was developed by Stephen Wissow to provide initial data selection and acquisition functionality for a data pipeline used in a group project for an interdisciplinary bioinformatics course at the University of New Hampshire in Spring 2020.

The course was MCBS 913: Bioinformatics, with instructors Kelley Thomas and Toni Westbrook, and the members of the group were Talha Siddique, Ryan Wilmot, Stephen Wissow, and Yibo Xu.

The idea for using NCBI's `prokaryotes.txt` table of contents file came from [Joseph Sevigny](https://github.com/Joseph7e/MDIBL-T3-WGS-Comparative), and was suggested by Ryan Wilmot.

The idea for using NCBI's `taxdump.tar.gz` nodes and names database to add taxonomic rank information to the `prokaryotes.txt` file was suggested by Toni Westbrook.

## Installation

proqueryote was developed using Python 3.7.

## Usage
```
proqueryote.py <query_file> [--fna] [-v[v[v[v]]]]
```
## Options
* `<query_file>` is a plain text file in query format (see below)
* Use the `--fna` flag to download genome instead of proteome data for matching records.
* `[-v[v[v[v]]]]` Four additional levels of verbosity are available. `-v` and `-vv` are useful if you want to know what's happening while you're sitting watching the cursor blink at you. `-vvv` and `-vvvv` may be helpful for debugging.

## Query File
The query file format provides a simple and flexible way to specify search criteria with AND and OR Boolean relationships. The format and syntax are described first, and then the semantics.

The query file is a series of one or more queries:
```
<Query>
[<Query> ...]
```
Each query begins with the `query` keyword on the first line, and is followed on subsequent lines by one or more criteria, one per line:
```
query
<Criterion>
[<Criterion> ...]
```
Each criterion appears in its own line, and is a space delimited string. The first word is the column name, followed by one or more valid value for that column:
```
<column_name> <Value> [<Value> ...]
```
A data record from the prokaryotes database matches a...
* Criterion if it matches ANY of that Criterion's Values (boolean OR, inclusive);
* Query if it matches ALL of that Query's Criterions (boolean AND);
* Query File if it matches ANY of that Query File's Queries (boolean OR, inclusive).

Additional syntax notes:
* Blank lines, and lines beginning with octothorpe character ("#"), are ignored.
* The keyword `query`, all <`column_name>`s, and all `<Value>`s are case sensitive.
* No other words may appear on the "query" keyword line.
* Neither a `<column_name>` nor a `<Value>` may contain a space " " character.
* Every `<column_name>` specified in the query file must exist in the prokaryotes database.

Here is an example query file:
```
query
Genus Actinotelluria Blastococcus Cumulibacter Geodermatophilus Modestobacter

query
Phylum Fusobacteria

query
Family Simpsons
Genus male
```
In the above example query file, rows that match any of the five genera in the first query, or are of phylum Fusobacteria, or are male members of the Simpsons family, will be returned.

## Results
Results will be downloaded to a new folder in the working directory, with the following name format to help distinguish among runs:
```
proqueryote-<TYPE>-<NUM_DATA_TABLE_RECORDS>-<DATE_TIME_STAMP>
```
...where `<TYPE>` is each `faa` or `fna`, and `<NUM_DATA_TABLE_RECORDS>` is the number of records in the prokaryotes database that matched your query. For example,:
```
proqueryote-faa-79-2020-04-24-1242.20/
```
