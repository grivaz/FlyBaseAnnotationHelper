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

import argparse
import configparser
import pickle
import sys
import requests
import xmltodict
import time
import os
import get_genes
import csv
import deep_learning
import tqdm

arg_parser = argparse.ArgumentParser(description="Gets gene candidates from xml papers")
arg_parser.add_argument("input", type=argparse.FileType('r'), help="A text file containing one pmid per line")
cmd_args = arg_parser.parse_args()

config_parser = configparser.ConfigParser()
config_parser.read("config/config.ini")

#get pmid to pmcid pickle dictionary
with open(config_parser.get('PICKLES','PMC_ids_dict'), "rb") as f:
    pmid_to_pmcid_dict = pickle.load(f)

def getFtpPath(pmcid: str):
    """Returns the ftp path to a paper given its pmcid

    Parameters:
        pmcid, str
            The pmcid of the paper
    """
    pmc_api_url = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=" + pmcid
    response = requests.get(pmc_api_url)
    dict_data = xmltodict.parse(response.content)
    ftplink = ""
    if "error" in dict_data['OA']:
        print("ERROR: " + dict_data['OA']['error']['@code'] + " " + pmcid, file=sys.stderr)
    else:
        link = dict_data['OA']['records']['record']['link']
        if isinstance(link, list):
            ftplink = link[0]['@href']
        else:
            ftplink = link['@href']
        assert (".tar.gz" in ftplink)
    time.sleep(config_parser.getint('PARAMETERS','sleep_time_between_requests'))
    return ftplink

def download(ftp: str):
    """Downloads a paper given its ftp path

    Parameters:
        ftp, str
            The ftp path to the paper
    """
    wget = f"wget -nc --timeout=10 -P {config_parser.get('PATHS', 'corpus')} {ftp}"
    os.system(wget)

def getXmlFromTar(pmcid: str):
    """Extracts the xml file from the tar.gz file

    Parameters:
        pmcid, str
            The pmcid of the paper
    """
    f = f"{config_parser.get('PATHS', 'corpus')}/{pmcid}.tar.gz"
    untar = "tar -xf " + f + " -C " + config_parser.get('PATHS', 'corpus')
    os.system(untar)
    # make xml directory if it doesn't exist
    if not os.path.exists(config_parser.get('PATHS', 'xml')):
        os.makedirs(config_parser.get('PATHS', 'xml'))
    copy_xml = "cp " + config_parser.get('PATHS', 'corpus') + "/" + pmcid + "/" + "*.nxml " + config_parser.get('PATHS', 'xml') \
               + "/" + pmcid + ".nxml"
    os.system(copy_xml)

def removeFiles(pmcid: str):
    """Removes paper files that were downloaded

    Parameters:
        pmcid, str
            The pmcid of the paper
    """
    f = f"{config_parser.get('PATHS', 'corpus')}/{pmcid}.tar.gz"
    try:
        os.remove(f)
    except FileNotFoundError:
        pass
    f = f"{config_parser.get('PATHS', 'corpus')}/{pmcid}"
    # remove extracted directory
    os.system("rm -r " + f)
    f = f"{config_parser.get('PATHS', 'xml')}/{pmcid}.nxml"
    try:
        os.remove(f)
    except FileNotFoundError:
        pass

#unpickle gene dictionary
with open(config_parser.get('PICKLES','gene_dict'), "rb") as f:
    gene_dict = pickle.load(f)

#unpickle fbid to symbol dictionary
with open(config_parser.get('PICKLES','fbid_to_symbol_dict'), "rb") as f:
    fbid_to_symbol = pickle.load(f)

if config_parser.getboolean('PARAMETERS', 'use_deep_learning'):
    deep_learning.initialize(config_parser.get('PATHS', 'deep_learning_model'))

results = {}
with open(cmd_args.input.name, "r") as f:
    input_list = f.readlines()
for pmid in tqdm.tqdm(input_list, desc="Processing articles", total=len(input_list)):
    if pmid.strip() not in pmid_to_pmcid_dict:
        # print it to the standard error stream
        print("no pmcid for " + pmid.strip(), file=sys.stderr)
    else:
        pmcid = pmid_to_pmcid_dict[pmid.strip()]
        ftp = getFtpPath(pmcid)
        download(ftp)
        getXmlFromTar(pmcid)
        result = None
        if config_parser.getboolean('PARAMETERS', 'use_deep_learning'):
            result = deep_learning.get_genes_with_dl(os.path.join(config_parser.get('PATHS', 'xml'), pmcid + ".nxml"), gene_dict, fbid_to_symbol)
        else:
            result = get_genes.get_genes(os.path.join(config_parser.get('PATHS', 'xml'), pmcid + ".nxml"),
                                    gene_dict, config_parser.get('PARAMETERS', 'snippet_type'),
                                    config_parser.getboolean('PARAMETERS', 'output_gene_occurence'),
                                    config_parser.getboolean('PARAMETERS', 'output_gene_frequency'),
                                    config_parser.getboolean('PARAMETERS', 'output_word_frequency'),
                                    config_parser.getboolean('PARAMETERS', 'output_raw_occurence'))
        results[pmid.strip()] = result
        if config_parser.getboolean('PARAMETERS', 'remove_files'):
            removeFiles(pmcid)


with open(config_parser.get('PATHS', 'output'), 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL, escapechar='\\')
    # write data to file
    # if we're using deep learning
    if config_parser.getboolean('PARAMETERS', 'use_deep_learning'):
        for pmid in results:
            for fbgn in results[pmid]:
                #when using deep learning, we only output the pmid, fbgn, and confidence
                writer.writerow([pmid, fbgn, results[pmid][fbgn]])

    else:
        for pmid in results:
            confidences = results[pmid][0]
            occurrences = results[pmid][1]
            for fbgn in confidences:
                scores = []
                if get_genes.GENES in confidences[fbgn]:
                    scores.append(confidences[fbgn][get_genes.GENES])
                if get_genes.WORD in confidences[fbgn]:
                    scores.append(confidences[fbgn][get_genes.WORD])
                if get_genes.RAW in confidences[fbgn]:
                    scores.append(confidences[fbgn][get_genes.RAW])
                snippet_type = config_parser.get('PARAMETERS', 'snippet_type')
                if snippet_type != 'none' and config_parser.getboolean('PARAMETERS', 'output_gene_occurence'):
                    for i, occurrences_for_gene in enumerate(occurrences[fbgn]):
                        genes_occurrence = occurrences_for_gene[0]
                        snippet = occurrences_for_gene[1]
                        writer.writerow([pmid, fbgn, genes_occurrence, snippet]+scores)
                elif snippet_type != 'none' or config_parser.getboolean('PARAMETERS', 'output_gene_occurence'):
                    for occurrence in occurrences[fbgn]:
                        writer.writerow([pmid, fbgn, occurrence]+scores)
                else:
                    writer.writerow([pmid,fbgn]+scores)


