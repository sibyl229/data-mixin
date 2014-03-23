A Pipeline for Multiple Data Extraction & Integration
========

This pipeline is extracted from a previous project, where I needed to compute various features for human and yeast proteins, and look for correlation between protein half-life and these features.

The pipeline aims to **organize** scripts that handle different types of data, **automate** data extraction and integration, and improve **extensibility**.

I don't know if this would be useful for anyone, but I present this experimental pipeline here.



## Usage

- Clone the repo

- At the project root, run `python code/main.py` and follow the prompt.



## How to extend it to include a new feature and/or a new species

### Add a new species (skip if only adding a new feature)

- create a subdirectory under raw_inputs/ for data of the new species

- add a module (a Python file) under code/species/, say `sp1.py`

- edit this species in the `species_choices` variable in `code/main.py` (more details there)

- include a MySPC class (in the module you create), that extends SpeciesConnectivityMixin of species.base
(This enables the pipeline to dynamically determine what Data objects are available that are associated with the species)

- follow the instructions below to add features



### Add a new feature
- Create a new class in the species module, say `code/species/sp1`, by extending **MySPC** and additionally

  	 - if the feature score already is in csv/tsv file, where each row includes the quantitative measurements and the identifier that you intend to keep, extend the AbstractData class (refer to data.py)

	 - otherwise, extend the AbstractComputedData class (refer to compute_data.py), AND edit main() in compute_data.py to instantiate the new class (which computes scores for the feature)

- run code/main.py 


## Project Organization

`raw_inputs/species_name/`
	   	 
- ALL data gathered, AND feature scores computed from those data

`input/species_name/`

- automatically generated from the pipeline that serves as input for further analysis
- contains preprocessed data (feature scores for **a select list of proteins**)	      	


`code/`

- `main.py` a command-line interface

- `data.py` and `compute_data.py` deal with extracting and computing various features from raw data. They contain AbstractData, AbstractComputedData, and various subclasses of them, which are extended and customized in species-specific modules under `code/species/`.

- `DataRemap.py` takes the master list of proteins as template and fill in feature scores for them (producing files with .derived_features.tsv suffix under input/species_name/)

`code/species/`
- customization for species-specific data handling 

`learning_setup.py` 
- combine various features for learning, normalize the features

`linear.py` 
- tries to play with some learning algorithms		 	 



## Dependencies
### Python Dependencies
** required
numpy
scipy
biopython	# if computing features from Fasta files
matplotlib	# for plotting

** likely needed
scikit-learn 	# for machine learning

if you have pip setup, all packages can be installed with:
pip install package_name




## Citations

1. Belle A, Tanay A, Bitincka L, Shamir R, O’Shea EK. 2006. Quantification of protein half-lives in the budding yeast proteome. Proc Natl Acad Sci U S A 103:13004–13009.

2. Cherry JM, Hong EL, Amundsen C, Balakrishnan R, Binkley G, Chan ET, Christie KR, Costanzo MC, Dwight SS, Engel SR, Fisk DG, Hirschman JE, Hitz BC, Karra K, Krieger CJ, Miyasato SR, Nash RS, Park J, Skrzypek MS, Simison M, Weng S, Wong ED. 2011. Saccharomyces Genome Database: the genomics resource of budding yeast. Nucl. Acids Res. gkr1029.

3. The UniProt Consortium. 2014. Activities at the Universal Protein Resource (UniProt). Nucleic Acids Research 42:D191–D198.

4. Di Domenico T, Walsh I, Martin AJM, Tosatto SCE. 2012. MobiDB: a comprehensive database of intrinsic protein disorder annotations. Bioinformatics 28:2080–2081.

5. Venne AS, Vögtle F-N, Meisinger C, Sickmann A, Zahedi RP. 2013. Novel Highly Sensitive, Specific, and Straightforward Strategy for Comprehensive N-Terminal Proteomics Reveals Unknown Substrates of the Mitochondrial Peptidase Icp55. J. Proteome Res. 12:3823–3830.

6. Burton JL, Solomon MJ. 2001. D box and KEN box motifs in budding yeast Hsl1p are required for APC-mediated degradation and direct binding to Cdc20p and Cdh1p. Genes Dev 15:2381–2395.
