# Fly Base Annotation Helper
The Fly Base Annotation Helper is a Python script that helps annotate Drosophila genes in scientific papers. Given a list 
of PubMed IDs (pmids), the script downloads the corresponding papers, parses their content, and identifies genes that 
appear in italics, in non-italicized text, and in the main sections of the paper (excluding the introduction).

The script outputs a TSV file with each line corresponding to a gene occurrence in the paper, unless the configuration 
states that no gene occurrence or snippet should be output, in which case each line corresponds to one gene in one paper.
The output includes the following columns:

PMID of the paper
FBGNID of the gene
Gene occurrence (the actual spelling of the gene in the text), if desired
Snippet (can be configured for long, short, or none)
Gene frequency (number of gene occurrences in the paper that are the given gene, if desired)
Word frequency (number of words in the paper that are the gene, if desired)
Raw count (number of occurrences of the gene in the paper, if desired)
The last three columns can be used as a confidence measure.

## Getting Started
Before running the Fly Base Annotation Helper, make sure that you have Python 3 installed on your system. You can download
it from the official website. You will also need to install the required dependencies listed in the requirements.txt file
by running the following command in your terminal:

```pip install -r requirements.txt```

## Usage
To use the Fly Base Annotation Helper, you first need to generate the necessary resources by running the 
update_resources.py script. This script looks up the necessary FlyBase files from config.ini and generates two pickles: 
a pmid to pmcid dictionary and a gene dictionary that maps different spellings of the genes to their FBGNID.

Once the resources are generated, you can run the main script fly_base_annotation_helper.py. The script can be 
configured via config.ini. You will also need to pass it a file containing one pmid per line, for the papers for which you
want the genes to be found. The configuration file includes options to specify paths to FlyBase files, paths to pickles,
output paths, and various parameters, such as the politeness parameter for requests to the NCBI server and whether to 
output gene occurrence, snippet, gene frequency, word frequency, and raw count. If you use the machine learning algorithm, the output will be paper, gene and confidence.

## Output
The output of the Fly Base Annotation Helper is a TSV file with the columns described above. The output can be used help
human annotators tag the papers with the genes that they are about.

## Contributions and Issues
If you have any questions or issues with the Fly Base Annotation Helper, please feel free to open an issue on the [GitHub 
repository](https://github.com/grivaz/FlyBaseAnnotationHelper). Contributions are also welcome via pull requests.

## License
The Fly Base Annotation Helper is released under the GPL 3.0 License.
