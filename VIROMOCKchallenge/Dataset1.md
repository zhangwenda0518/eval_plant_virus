## Dataset 1: Viral concentration at the strain level (CTV on Citrus)

### Introduction to the dataset

The real dataset is composed of mixed infections of *Citrus tristeza virus* (CTV), *Citrus vein enation virus* (CVEV), *Citrus exocortis viroid* (CEVd), *Citrus viroid III* (CVd-III) and *Hop stunt viroid (HSVd)* on citrus. Three other viruses/viroids (*Citrus endogenous pararetrovirus* (CitPRV), *Citrus dwarfing viroid* (CDVd) and *Apple hammerhead viroid-like RNA* (AHVd RNA)) have been detected by NGS, but not confirmed with other methods (possibly the result of cross-contamination). **The challenge addressed is the ability to detect several strains showing different concentrations**, using CTV. The goal was therefore to add artificial CTV strains at different frequencies in the real datasets.

First, all the CTV genomes available on NCBI have been downloaded and compared to the CTV strain already presents in the real dataset. Then, a phylogenetic tree has been built in [Geneious](https://www.geneious.com/) using [PHYML](https://github.com/stephaneguindon/phyml) (Figure 1). The CTV strain already present is surrounded in blue in the tree, while the 3 new selected strains are surrounded in red. All these CTV strains belong to different phylogenetic groups, so that it should still be possible to identify each strain.

![Figure1](images/Dataset1_Figure1_TreeCTVstrains_2.png)
*Figure 1: Phylogenetic tree of CTV. The strain indicated in blue is the consensus sequence derived from the real dataset. The strains indicated in red are the strains which were artificially added to the dataset.*

The percentage of identity between the four strains at the nucleotide level ranges between 79 and 87%. Figure 2 shows the percentage of similarity between the strains along the genome, using the CTV strain already present in the sample as reference sequence to which the other strains were compared.

![Figure2](images/Dataset1_Figure2_SimilarityPlotCTVstrains_2.png)
*Figure 2: Similarity plot of the full-length sequence of the CTV strains.*

### Artificial reads added to the dataset

An overview of the expected composition of the dataset can be found in Table 1.

*Table 1: Composition of Dataset 1. The total number of reads in this dataset is 2x21,703,43 (2x150 bp).*

|     Virus/Viroid                                     	|     Reads type                          	|     Number of artificial reads added    	|     Expected proportion of artificial CTV reads (%)    	|     Expected average number of reads per position    	|
|------------------------------------------------------	|-----------------------------------------	|-----------------------------------------	|--------------------------------------------------------	|------------------------------------------------------	|
|     *Citrus endogenous pararetrovirus* (CitPRV)      	|     Real but unconfirmed<sup>1</sup>    	|                                         	|                                                        	|                                                      	|
|     *Citrus dwarfing viroid* (CDVd)                  	|     Real but unconfirmed<sup>1</sup>    	|                                         	|                                                        	|                                                      	|
|     *Apple hammerhead viroid-like RNA* (AHVd RNA)    	|     Real but unconfirmed<sup>1</sup>    	|                                         	|                                                        	|                                                      	|
|     *Citrus vein enation virus* (CVEV)               	|     Real                                	|                                         	|                                                        	|                                                      	|
|     *Citrus exocortis viroid* (CEVd)                 	|     Real                                	|                                         	|                                                        	|                                                      	|
|     *Citrus viroid III* (CVd-III)                     |     Real                                	|                                         	|                                                        	|                                                      	|
|     *Hop stunt viroid* (HSVd)                        	|     Real                                	|                                         	|                                                        	|                                                      	|
|     *Citrus tristeza virus* (CTV)                    	|     Real                                	|                                         	|                                                        	|                                                      	|
|     *Citrus tristeza virus* (CTV) JQ911663           	|     Artificial                          	|     67838                               	|     69.75                                              	|     530.3                                            	|
|     *Citrus tristeza virus* (CTV) KU883267           	|     Artificial                          	|     28418                               	|     29.22                                              	|     221.7                                            	|
|     *Citrus tristeza virus* (CTV) MH323442           	|     Artificial                          	|     1002                               	|     1.03                                               	|     7.8                                              	|

<sup>1</sup> There are a few reads of this virus/viroid present, but these could not be confirmed by qPCR and are most probably the result of cross-contamination.

### Results of different labs analyzing the dataset

The different labs that analyzed the dataset are listed in Table 2.

*Table 2: Participants to the VIROMOCK challenge of Dataset 1.*

| Participant                      	| Institute 	| Country 	| email                                  	|
|----------------------------------	|-----------	|---------	|----------------------------------------	|
| Lucie Tamisier                   	| ULg       	| Belgium 	| <lucie.tamisier@uliege.be>             	|
| Annelies Haegeman, Yoika Foucart 	| ILVO      	| Belgium 	| <annelies.haegeman@ilvo.vlaanderen.be> 	|
| Alex Hu | USDA-APHIS | USA | <Xiaojun.hu@usda.gov> |


The observed values of the different participants are listed in Table 3.

*Table 3: Observed composition of Dataset 1 after analysis by different labs.*

| Institute/Lab | Observed closest   NCBI accession | Observed   proportion of real + artificial CTV reads (%) | Observed   proportion of artificial CTV reads (%) | Were you able to   detect all strains? |
|---------------|-----------------------------------|----------------------------------------------------------|---------------------------------------------------|----------------------------------------|
| ILVO          | HF679486.1                        | NA                                                       | NA                                                | yes                                    |
| ILVO          | S67442.1                          | NA                                                       | NA                                                |                                        |
| ILVO          | AB054623.1                        | NA                                                       | NA                                                |                                        |
| ILVO          | KJ810553.1                        | NA                                                       | NA                                                |                                        |
| ILVO          | AF260651                          | 9.31                                                     | NA                                                |                                        |
| ILVO          | JQ911663                          | 63.29                                                    | 69.79                                             |                                        |
| ILVO          | KU883267                          | 26.46                                                    | 29.18                                             |                                        |
| ILVO          | MH323442                          | 0.94                                                     | 1.03                                              |                                        |
| USDA-APHIS | HF679486.1 | NA   | NA    | yes |
| USDA-APHIS | FJ751933.1 | NA   | NA    |     |
| USDA-APHIS | S76452.1   | NA   | NA    |     |
| USDA-APHIS   | AB054611.1 | NA   | NA    |     |
| USDA-APHIS   | KC748391   | 4.22 | NA    |     |
| USDA-APHIS   | JQ911663   | 66.7 | 69.64 |     |
| USDA-APHIS   | KU883267   |   28 | 29.24 |     |
| USDA-APHIS   | MH323442   | 1.07 |  1.12 |     |

### Comments of different labs while analyzing the dataset

Here you can read some comments the participants had while analyzing the dataset.

| Institute/Lab     | Comments                                                                                                                                                        |
|-------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| USDA-APHIS | Citrus   exocortis Yucatan viroid (FJ751933.1) was detected. It's differerent from   Citrus exocortis viroid. Several partial CEVd genomes were also assembled. |

### How to participate to the VIROMOCK challenge

The dataset can be downloaded [here](https://doi.org/10.5061/dryad.0zpc866z8) (click the arrow next to "November 2, 2021" to download each dataset separately).

If you finish your analysis, we encourage you to submit your results through [this Google Sheet](https://docs.google.com/spreadsheets/d/1ooAtkGojO0P0cfZC4YObk5S7wUH82xFAC-A8fRdjRFk/edit#gid=0).

The Google Sheet will allow you to share your results in detail. Only the green columns are required. However, we encourage you to give as much information as possible.

After submission of the Google Sheet, your results will be processed and added to Table 3.
