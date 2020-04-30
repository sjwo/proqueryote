#!/usr/bin/env python

# proqueryote
# UNH MBCBS 913 Spring 2020
# Stephen Wissow


import os
import sys
import datetime
import argparse
import re
import pandas as pd
from typing import List
from pathlib import Path
from collections import namedtuple
Values = List[str]

# CACHE FILE LOCATIONS
HOME = str(Path.home())
CACHE = HOME + '/.proqueryote/'
PROKARYOTES = CACHE + 'prokaryotes.txt'
NODES_FILE = CACHE + 'nodes.dmp'
NAMES_FILE = CACHE + 'names.dmp'
AUGMENTED = CACHE + 'prokaryotes_taxonomic.txt'
TAXDUMP_FILENAME = 'taxdump.tar.gz'
TAXDUMP = CACHE + '/' + TAXDUMP_FILENAME
CACHE_FILES = [PROKARYOTES, NAMES_FILE, NODES_FILE, AUGMENTED]
CACHE_URLS = [
  'ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt',
  'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/' + TAXDUMP_FILENAME
]
FAA_EXTENSION = '/*protein.faa.gz'
FNA_EXTENSION = '/*_cds_from_genomic.fna.gz'
# ORDER LISTING OF RANKS USED BY NODES.DMP FROM TAXDUMP.TAR.GZ FROM NCBI
RANK_ORDER = ['superkingdom', 'kingdom', 'superphylum', 'phylum', 'subphylum', 'superclass', 'class', 'subclass', 'infraclass', 'order', 'superfamily', 'family', 'subfamily', 'genus', 'subgenus', 'species group', 'species subgroup', 'species', 'subspecies']

# RANKS FOR WHICH TO ADD COLUMNS TO PROKARYOTES TABLE
def RANKS():
   return {'species':'','genus':'','family':'','phylum':''}

epilog = """
See query file format specification in associated README.
"""

parser = argparse.ArgumentParser(prog="column_search", usage="column_search <queries_file>", description="Searches prokaryotes_taxanomic.txt file in working directory for all rows matching queries in given <queries_file>.", epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("queries_file", type=str)
parser.add_argument('--verbose', '-v', action='count', default=0,)
parser.add_argument('--fna', action='store_const', const=True, default=False)
args = parser.parse_args()

def log(msg, level=0):
   if (args.verbose >= level): print(msg)

# TEST CACHE PRESENCE
def testCachePresence():
  log('testing cache presence', 2)
  if (os.path.isdir(CACHE)):
    log('testCachePresence: cache directory present', 3)
    for each in CACHE_FILES:
      if not os.path.exists(each):
        log(f'testCachePresence: {each} is missing from cache', 3)
        return False
      else:
        log(f'testCachePresence: {each} file is present', 3)
    log('testCachePresence: cache found to be present. OK.', 3)
    return True
  else:
    log('testCachePresence: cache directory missing', 3)
    return False

def updateCache():
  log('Updating local cache.')
  # remove old cache, create new
  # (will need to make this more fine-grained if we ever
  # store config files in the hidden folder)
  os.system(f'rm -r {CACHE} 2> /dev/null')
  os.mkdir(CACHE)
  
  # download prokaryotes and nodes DB and extract
  working_directory = os.getcwd()
  os.chdir(CACHE)
  log('Downloading resources for local cache', 1)
  for resource in CACHE_URLS:
    log(f'{resource}', 2)
    os.system(f'wget --quiet {resource}')
  os.system(f'tar -z -xf {TAXDUMP}')
  os.chdir(working_directory)

  # created augmented version of table, with new taxonomic rank columns
  taxonomifier = Taxonomifier()
  taxonomifier.produceTaxonomic()


'''
==========
TAXONOMIFY
==========
'''
ParentAndRank = namedtuple('ParentAndRank', 'parent rank')

def getLineCount():
   return sum(1 for record in open(args.proks_file, 'r'))

class Taxonomifier(object):
  def __init__(self):
    log(f'Preparing to augment table with additional taxonomic rank columns', 1)
    log(f'Loading local cache', 2)
    log(f'Loading {PROKARYOTES}', 3)
    self.proks = open(PROKARYOTES, 'r').readlines()
    self.nodes = self.getNodesMap()
    self.names = self.getNamesMap()
    log(f'Read in {len(self.proks)} lines from {PROKARYOTES}', 3)

  @staticmethod
  def getNamesMap():
    log(f'Loading {NAMES_FILE}', 3)
    names = dict()
    with open(NAMES_FILE, 'r') as names_f:
        names_d = names_f.readlines()
        for record in names_d:
          cells = re.split(r'\t\|\t', record)
          taxid = cells[0]
          name = cells[1]
          name_type = cells[3]
          if (name_type.startswith('scientific name')):
              names.setdefault(int(taxid), name)
    return names

  @staticmethod
  def getNodesMap():
    log(f'Loading {NODES_FILE}', 3)
    nodes = dict()
    with open(NODES_FILE, 'r') as nodes_f:
        nodes_d = nodes_f.readlines()
        for record in nodes_d:
          cells = re.split(r'\t\|\t', record)
          this_taxid = cells[0]
          this_parent_taxid = cells[1]
          this_rank = cells[2]
          parent_rank = ParentAndRank(int(this_parent_taxid), this_rank)
          nodes.setdefault(int(this_taxid),parent_rank)
    return nodes

  def getName(self, taxid):
    """
    Retrieves scientific name corresponding to given taxid.
    """
    return self.names[int(taxid)]

  def getParent(self, taxid):
    return self.nodes[int(taxid)].parent

  def getRank(self, taxid):
    return self.nodes[int(taxid)].rank

  @staticmethod
  def getTaxID(proks_record):
    """Parse TaxID field from string of row of prokaryotes.txt
    """
    return proks_record.split('\t')[1]

  def getOneTaxonomy(self, taxid):
    """Searches nodes.dmp for desired ranks, pulling corresponding names from names.dmp.

    :return: A TSV string: "<species>\t<genus>\t<family>\t<phylum>"
    """
    ranks = RANKS()
    searchTaxID = str(taxid)
    log(f'getOneTaxonomy: taxid {taxid}', 4)
    # continue searching until we've found an entry for each desired taxonomic rank
    while ('' in ranks.values()):
        log(f'getOneTaxonomy: searching for {searchTaxID}', 4)
        if (searchTaxID == '1'):
          log(f'getOneTaxonomy: arrived at root, curtailing search with extant results: {ranks}', 4)
          break
        try:
          this_parent_taxid = self.getParent(searchTaxID)
        except KeyError:
          log(f'getOneTaxonomy: caught KeyError raised by getParent: TaxID {searchTaxID} is not listed in {args.nodes_file}.', 4)
          break
        this_rank = self.getRank(searchTaxID)
        log(f'getOneTaxonomy: {searchTaxID} is a {this_rank} with parent {this_parent_taxid}', 4)
        # TODO: may want to alter this to confirm that ranks[this_rank] == '',
        # to prevent over-writing of multiple entries exist in nodes.dmp, not all correct?:
        if (this_rank in ranks):
          name = self.getName(searchTaxID)
          ranks[this_rank] = name
          log(f'getOneTaxonomy: added {searchTaxID} : {name} as {this_rank}. we have thus far found: {ranks}', 4)
        searchTaxID = str(this_parent_taxid)
    log(f'getOneTaxonomy: found {ranks}', 4)
    return(ranks)

  def produceTaxonomic(self):
    log('Building new table', 1)
    # costruct new header row
    old_header = self.proks[0].rstrip()
    new_cols = RANKS().keys()
    log('old header:\n' + old_header, 3)
    log('new columns:\n' + str(new_cols), 4)
    new_header = old_header
    for rank in RANK_ORDER:
      if (rank in new_cols):
        new_header += '\t' + rank.title()
    new_header += '\n'
    log('new header:\n' + new_header, 3)

    # start writing to new file
    with open(AUGMENTED, 'w') as out_f:
      # write new header row
      log(f'About to write new header to new, augmented table:\n{new_header}', 4)
      out_f.write(new_header)
      
      # construct new entries
      for record in self.proks[1:]:
        log('old row:\n' + record, 4)
        taxid = self.getTaxID(record)
        taxonomy = self.getOneTaxonomy(taxid)
        new_record = record.rstrip()
        for rank in RANK_ORDER:
          if (rank in taxonomy.keys()):
              new_record += '\t' + taxonomy[rank]
        new_record += '\n'
        log('new record:\n' + new_record, 4)
        out_f.write(new_record)
    log(f'New table with additional taxonomic columns:\n\t{AUGMENTED}', 2)





'''
=====
QUERY
=====
'''

class Criterion(object):
  '''
  Specifies acceptable values for a given column, for a given query.
  '''
  def __init__(self, column: str, values: Values):
    self.column = column
    self.values = values
  def __repr__(self):
    return f'{self.column}: {self.values}'

class Query(object):
  '''
  Specifies all criteria (in form of Criterion objects) for a given query.
  '''
  def __init__(self):
    self.criteria = list()
  
  def addCriterion(self, criterion: Criterion):
    self.criteria.append(criterion)
  
  def __repr__(self):
    return 'Query:\n' + str(self.criteria)

def parseQueries(queries_file):
  log(f'Parsing query file: {queries_file}', 1)
  queries = list()
  with open(queries_file, 'r') as queries_f:
    for line in queries_f.readlines():
      log("  starting loop", 3)
      line = line.rstrip()
      # skip blank and comment lines
      if (line == '' or line[0] == '#'):
        log("  found empty line or comment; skipping", 3)
        continue
      # start a new query
      elif (line == 'query'):
        log("  found new query", 3)
        queries.append(Query())
      # append criterion to current query
      else:
        log("  adding to criterion to current query", 3)
        delimited = re.split(' ', line)
        queries[len(queries) - 1].addCriterion(Criterion(delimited[0],delimited[1:]))        
  log(f'parsed queries:\n{queries}', 3)
  return queries


class Data(object):
  '''
  API for searching and returning data from a TSV-formatted file with a header row (such as prokaryotes.txt from NCBI, optionally with added taxonomic rank columns).
  '''
  def __init__(self, data_file, queries_file):
    super().__init__()
    log('Preparing to process query using local cache.')
    self.data_file = data_file
    dtype = {
      '#Organism/Name' : 'str',
      'TaxID' : 'int',
      'BioProject Accession' : 'str',
      'BioProject ID' : 'str', # contains '-'
      'Group' : 'str',
      'SubGroup' : 'str',
      'Size (Mb)' : 'float',
      'GC%' : 'str', # contains '-'
      'Replicons' : 'str',
      'WGS' : 'str',
      'Scaffolds' : 'str', # contains '-'
      'Genes' : 'str', # contains '-'
      'Proteins' : 'str', # contains '-'
      'Release Date' : 'str',
      'Modify Date' : 'str',
      'Status' : 'str',
      'Centert' : 'str',
      'BioSample Accession' : 'str',
      'Assembly Accession' : 'str',
      'Reference' : 'str',
      'FTP Path' : 'str',
      'Pubmed ID' : 'str', # must be parsed as CSV ints
      'Strain' : 'str',
      'Phylum' : 'str',
      'Family' : 'str',
      'Genus' : 'str',
      'Species' : 'str',
    }
    parse_dates = ['Release Date', 'Modify Date']
    log('Loading local cache in preparation for querying.', 1)
    self.table = pd.read_csv(
      data_file,
      sep='\t',
      dtype=dtype,
      parse_dates=parse_dates
    )
    self.queries = parseQueries(queries_file)
  
  def columns(self):
    return self.table.columns

  def criterionTruthOf(self, criterion):
    criterion_truth = None
    if (len(criterion.values) > 0):
      try:
        criterion_truth = self.table[criterion.column] == criterion.values[0]
        log('Found first value', 3)
        for additional_value in criterion.values[1:]:
          criterion_truth |= self.table[criterion.column] == additional_value
          log('Found additional value', 3)
      except KeyError:
        log(f'\n ** KeyError: The query file specifies a column name ("{criterion.column}") that does not exist in in the data file ("{self.data_file}"). Aborting. **')
        log('\n    Available columns are:\n',1)    
        {log('        ' + column, 1) for column in self.columns()}
        print()
        sys.exit()
    else:
      log(f'No values found for this criterion: {criterion}. Please fix the query file. Aborting.')
      sys.exit()
    return criterion_truth

  def queryTruthOf(self, query: Query):
    query_truth = None
    if (len(query.criteria) > 0):
      log('Found first criterion', 3)
      query_truth = self.criterionTruthOf(query.criteria[0])
      for additional_criterion in query.criteria[1:]:
        log('Found additional criterion', 3)
        query_truth &= self.criterionTruthOf(additional_criterion)
    else:
      log(f'No criteria found for this query: {query}. Please fix the query file. Aborting.')
      sys.exit()
    return query_truth


  def select_by_queries(self):
    log('Applying parsed queries to data', 3)
    total_truth = None
    if (len(self.queries) > 0):
      log('Found first query', 3)
      total_truth = self.queryTruthOf(self.queries[0])
      for additional_query in self.queries[1:]:
        log('Found additional query', 3)
        total_truth |= self.queryTruthOf(additional_query)
    else:
      log('No queries found. Exiting.')
      sys.exit()
    log('Returning all matching species', 1)
    return self.table[total_truth]
        

def main():
  if (args.fna):
    log(f'Configured to download FNA genomes', 3)
  else:
    log(f'Configured to download FAA proteomes', 3)
  
  def fileExtension():
    if (args.fna):
      return FNA_EXTENSION
    else:
      return FAA_EXTENSION

  if not testCachePresence():
    updateCache()
  
  data = Data(AUGMENTED, args.queries_file)
  selected = data.select_by_queries()
  log(f'Selected these records:\n{selected}', 3)
  count = len(selected)

  def dataType():
    if (args.fna):
      return 'genomes'
    else:
      return 'proteomes'

  print(f'\nFound {count} species. Download available {dataType()}? [y/N]')
  while (True):
    choice = sys.stdin.readline().rstrip().lower()
    if ((choice == 'n') or (choice == '')):
      print('Aborting.')
      return
    elif (choice != 'y'):
      print('Please type "y" or "n" and then press ENTER:')
      continue
    else:
      break
    
  urls = selected['FTP Path']

  stamp = datetime.datetime.now().strftime('%Y-%m-%d-%H%M.%S')
  
  def directoryStub():
    if (args.fna):
      return 'fna'
    else:
      return 'faa'
  
  folder = f'proqueryote-{directoryStub()}-{count}-{stamp}'
  os.system(f'mkdir {folder}')
  os.chdir(folder)
  for url in urls:
    sys.stdout.write('.')
    sys.stdout.flush()
    prot_url = url + fileExtension()
    os.system(f'wget --quiet {prot_url}')
  os.system('gunzip -- *')
  os.chdir('..')
  print()

def test_queries_parser():
  # print(parseQueries(args.queries_file))
  # WARNING: DOES NOT USE CACHE. INSTEAD USES AUGMENTED TABLE IN
  #           WORKING DIRECTORY.
  data = Data('prokaryotes_taxonomic.txt', args.queries_file)
  print(data.queries)

if __name__ == "__main__":
    main()