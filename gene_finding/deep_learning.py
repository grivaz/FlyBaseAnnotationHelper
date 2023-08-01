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
from transformers import AutoTokenizer, AutoModelForSequenceClassification, pipeline
import os
import sys
import pubmed_parser as pp
import typing
import get_genes


model_pipeline = None
tokenizer_kwargs = {'padding': "max_length", 'truncation': True, 'max_length': 512}
def initialize(path_to_model):
    global model_pipeline
    if model_pipeline is None:
        os.environ["TOKENIZERS_PARALLELISM"] = "false"
        scibert = "allenai/scibert_scivocab_uncased"
        tokenizer = AutoTokenizer.from_pretrained(scibert, model_max_length=512)
        model = AutoModelForSequenceClassification.from_pretrained(path_to_model, num_labels=2)
        model_pipeline = pipeline("sentiment-analysis", model=model, tokenizer=tokenizer, device=-1)

def get_gene(fbrf, candidates, fbid_to_symbol):
    occurrences = set()
    for occurrence in candidates:
        occurrences.add(occurrence)
    out = fbid_to_symbol[fbrf] + " " + " ".join(str(candidates.count(cand)) + " " + cand for cand in occurrences)
    return out
 
def get_genes_with_dl(paper_file: str, gene_dict: typing.Dict[str, str], fbid_to_symbol: typing.Dict[str, str]):
    if os.path.isfile(paper_file) and paper_file.endswith("nxml"):
        # get text
        pubmed_dict = pp.parse_pubmed_xml(paper_file)  # dictionary output
        abstract = pubmed_dict["abstract"]  # abstract
        _, candidates = get_genes.get_genes(paper_file, gene_dict, 'none', True, False, False, False)
        results = {}
        for fbrf in candidates:
            gene = get_gene(fbrf, candidates[fbrf], fbid_to_symbol)
            input = gene + ". " + abstract
            results[fbrf] = transform(input)
        return results
    else:
        print("File does not exist or is not an nxml file: " + paper_file, file=sys.stderr)

def transform(input):
    #print error if model_pipeline is not initialized
    if model_pipeline is None:
        raise Exception("model_pipeline not initialized")
    result = model_pipeline(input, **tokenizer_kwargs)
    if result[0]["label"] == "LABEL_1":
        return result[0]["score"]
    else:
        return 1 - result[0]["score"]