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

from xml.etree import ElementTree

import typing

RAW = "raw"
WORD = "word"
GENES = "genes"

exceptions = None
BODY = "body"
SEC = "sec"



def is_exception(gene_canditate: str, exceptions: typing.List[str], exceptions_path: str) -> bool:
    """Indicates whether the given gene candidate matches the given exception list

    A candidate matches an exception if the exception is part of the candidate, and starts and end at the same place as
    the candidate, or in a non-alphanumeric character boundary within the candidate.

    Parameters:
    -----------
    gene_candidate, str
        The candidate as found in the document, i.e. an italized string that matches a gene synonym

    exception, List[str]

        A list of exceptions strings that should not be considered real gene candidates

    Returns:
    --------
    bool
        whether the gene_candidate should be considered an exception
    """

    if exceptions is None:
        exceptions = set()
        with open(exceptions_path, "r") as f:
            for ex in f.readlines():
                ex = ex.strip()
                exceptions.add(ex)

    for ex in exceptions:
        position = gene_canditate.find(ex)
        if position != -1:
            if position == 0 or not gene_canditate[position - 1].isalnum():
                if position + len(ex) == len(gene_canditate) -1 or not gene_canditate[position + len(ex)].isalnum():
                    return True
    return False

def get_genes(paper_file: str, gene_dict: typing.Dict[str, str], snippet_type: str, output_gene_occurrence: bool,
              gene_freq: bool, word_freq: bool, raw_occurrences: bool, exceptions_path: str) -> typing.Dict[str, float]:
    """
    Gets the genes that a paper should be tagged with

    A gene synonym appears at least once in the body of the paper minus the introduction. It is in italics and makes up
    the whole span of the italics except for possible white spaces. The length of its occurrence is more than one letter,
    and it is not part of the exceptions list. Confidence is the number of occurrences of the gene in the paper,
    normalized by the length of the paper.

    :param paper_file: the location of the paper as an xml file
    :param gene_dict: a dictionary of gene synonyms to fbid of the gene
    :param snippet_type: can be 'long', 'short', or 'none'.
    :param output_gene_occurrence: outputs the exact way the gene is spelled in the paper
    :param gene_freq: whether to use gene frequency to compute confidence
    :param word_freq: whether to use word frequency to compute confidence
    :param raw_occurrences: whether to output raw occurrences count
    :return: a dictionary of gene fbid to their confidence for the given paper
    """
    if not (snippet_type == 'long' or snippet_type == 'short' or snippet_type == 'none'):
        raise ValueError("snippet_type must be 'long', 'short', or 'none'")
    cands = set()
    with open(paper_file, "r") as p:
        xml_content = p.read()
        size = len(xml_content.split(" "))
    tree = ElementTree.fromstring(xml_content)
    parent_map = {c: p for p in tree.iter() for c in p}

    def is_in_relevant_section(node):
        in_body = False
        n = node
        sec = None
        while n in parent_map:
            if n.tag == SEC:
                sec = n
            elif n.tag == BODY:
                in_body = True
                break
            n = parent_map[n]
        if n.tag == BODY:
            in_body = True
        if not in_body:
            return False
        assert n.tag == BODY
        if sec == list(n)[0]: # somewhat crude way of guessing if it's in the introduction
            return False
        return True

    tags = dict()
    snippet_dict = dict()
    relevant_gene_mentions = 0

    for node in tree.iter('italic'):
        if node.text:
            gene_canditate = node.text.strip()
            if gene_canditate in gene_dict and len(gene_canditate) > 1 and not is_exception(gene_canditate, 
                                                                                            exceptions, exceptions_path):
                cands.add(gene_canditate)
                if is_in_relevant_section(node):
                    gene = gene_dict[gene_canditate]
                    if not gene in snippet_dict:
                        snippet_dict[gene] = []
                    if snippet_type != 'none':
                        if snippet_type == 'long':
                            pm = parent_map[node]
                            s = ' '.join([t for t in pm.itertext()])
                            if not output_gene_occurrence:
                                snippet_dict[gene].append(s)
                            else:
                                snippet_dict[gene].append((gene_canditate, s))
                        else:
                            if not output_gene_occurrence:
                                snippet_dict[gene].append(gene_canditate + node.tail[:100] if node.tail else "")

                            else:
                                snippet_dict[gene].append((gene_canditate,
                                                           gene_canditate + (node.tail[:100] if node.tail else "")))
                    elif output_gene_occurrence:
                        snippet_dict[gene].append(gene_canditate)
                    if gene in tags:
                        tags[gene] += 1
                    else:
                        tags[gene] = 1
                    relevant_gene_mentions += 1

    output = dict()
    for gene, occurrences in tags.items():
        confidences = {}
        if gene_freq:
            confidences[GENES] = occurrences / relevant_gene_mentions
        if word_freq:
            confidences[WORD] = occurrences / size
        if raw_occurrences:
            confidences[RAW] = occurrences
        output[gene] = confidences
    return output, snippet_dict

