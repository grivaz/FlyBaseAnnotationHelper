from setuptools import setup, find_packages

setup(
    name='gene_finding',
    version='0.1',
    author='Charlie Grivaz',
    description='Finds drosophilia genes in papers using a fine-tuned BERT model',
    packages=find_packages(),
    install_requires=[
        'transformers',
        'pubmed_parser',
        'typing;python_version<"3.5"',
    ],
)

