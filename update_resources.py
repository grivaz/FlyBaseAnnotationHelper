# Copyright 2023 Charlie Grivaz
 #
 # This file is part of Fly Base Annotation Helper
 #
 # Fly Base Annotation Helper is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 #
 # Fly Base Annotation Helper is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License
 # along with Fly Base Annotation Helper. If not, see <http://www.gnu.org/licenses/>.

import configparser
import csv
import os

from tqdm import tqdm
import pickle
import re

config = configparser.ConfigParser()
config.read("config/config.ini")

def getPMIDtoPMCID(path: str):
    """Creates a pmid to pmcid dictionary

    Parameters:
        path, str
            The path to the file containing the pmid to pmcid mapping. Typically, this file is called "PMC-ids.csv"
    """
    pmid_to_pmcid = {}
    with open(path, "r") as f:
        length = len(f.readlines())
    with open(path, "r") as f:
        rd = csv.DictReader(f)
        for row in tqdm(rd, desc=f"Reading {path}", unit=" lines", total=length):
            pmcid = row["PMCID"]
            pmid = row["PMID"]
            if pmid and pmcid:
                pmid_to_pmcid[pmid] = pmcid
    #pickle the dictionary
    # get the directory path from the config parameter
    dir_path = os.path.dirname(config.get('PICKLES', 'PMC_ids_dict'))

    # create the directory if it doesn't exist
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

    # open the file for writing
    with open(config.get('PICKLES','PMC_ids_dict'), "wb") as out:
        pickle.dump(pmid_to_pmcid, out)

def getGenesDict(gene_synonyms_path: str, current_genes_path: str):
    """Creates a dictionary of gene symbols to their flybase id

    Parameters:
        gene_synonyms_path, str
            The path to the file containing the gene synonyms. Typically, this file is called "fb_synonym_fb_[DATE].tsv"
        current_genes_path, str
            The path to the file containing the current genes. Typically, this file is called "currentDmelHsap.txt"
    """
    # tsv column indices
    PRIMARY_FBID = 0
    CURRENT_SYMBOL = 2
    CURRENT_FULLNAME = 3
    FULLNAME_SYNONYMS = 4
    SYMBOL_SYNONYM = 5
    relevant_genes = set()
    with open(current_genes_path, "r") as current_genes_file:
        for line in current_genes_file:
            gene = line.rstrip()
            relevant_genes.add(gene)
    gene_dict = dict()
    fbid_to_symbol = dict()

    with open(gene_synonyms_path, "r") as gene_file:
        length = len(gene_file.readlines())
    with open(gene_synonyms_path, "r") as gene_file:
        rd = csv.reader(filter(lambda row: row[0] != '#', gene_file), delimiter="\t", quotechar='"')
        for row in tqdm(rd, desc="Making genes dictionary", total=length, unit=" lines"):
            if (len(row) == 6 and row[PRIMARY_FBID].startswith("FBgn")):
                fbid = row[PRIMARY_FBID]
                if fbid in relevant_genes:
                    fullname_synonyms = row[FULLNAME_SYNONYMS]
                    if len(fullname_synonyms) > 0:
                        for syn in re.split("\|", fullname_synonyms):
                            gene_dict[syn] = fbid

                    fullname = row[CURRENT_FULLNAME]
                    if len(fullname) > 0:
                        gene_dict[fullname] = fbid  # fullname takes precedence over symbol when ambiguous.
                        # Last ambiguous fullname wins.

                    symbol_synonyms = row[SYMBOL_SYNONYM]
                    if len(symbol_synonyms) > 0:
                        for syn in symbol_synonyms.split("|"):  # here we ignore commas, which is not the same as
                            # before, just because it seems like this one is less comma separated, who knows
                            gene_dict[syn] = fbid

                    symbol = row[CURRENT_SYMBOL]
                    gene_dict[symbol] = fbid  # this should be last to have absolute precedence
                    fbid_to_symbol[fbid] = symbol
        #pickle the dictionaries
        # get the directory paths from the config parameter
        dir_path = os.path.dirname(config.get('PICKLES', 'gene_dict'))
        fbid_to_symbol_dir_path = os.path.dirname(config.get('PICKLES', 'fbid_to_symbol_dict'))

        # create the directories if it doesn't exist
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        if not os.path.exists(fbid_to_symbol_dir_path):
            os.makedirs(fbid_to_symbol_dir_path)

        # output the dictionaries
        with open(config.get('PICKLES','gene_dict'), "wb") as out:
            pickle.dump(gene_dict, out)
        with open(config.get('PICKLES','fbid_to_symbol_dict'), "wb") as out:
            pickle.dump(fbid_to_symbol, out)

if os.path.exists(config.get('PUBMED','PMC_ids')) and os.path.exists(config.get('FLYBASE','gene_synonyms'))\
        and os.path.exists(config.get('FLYBASE','current_genes')):
    getPMIDtoPMCID(config.get('PUBMED','PMC_ids'))
    getGenesDict(config.get('FLYBASE','gene_synonyms'), config.get('FLYBASE','current_genes'))
else:
    if not os.path.exists(config.get('PUBMED','PMC_ids')):
        print(config.get('PUBMED','PMC_ids') + " not found. Please add proper path in config.ini")
    if not os.path.exists(config.get('FLYBASE','gene_synonyms')):
        print(config.get('FLYBASE','gene_synonyms') + " not found. Please add proper path in config.ini")
    if not os.path.exists(config.get('FLYBASE','current_genes')):
        print(config.get('FLYBASE','current_genes') + " not found. Please add proper path in config.ini")