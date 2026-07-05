## Dataset 4: Close viroids (GYSVd-1 and GYSVd-2 on Grapevine)

### Introduction to the dataset

The real dataset is composed of mixed infections of *Grapevine red blotch virus* (GRBV), *Grapevine rupestris stem pitting-associated virus* (GRSPaV), *Hop stunt viroid* (HSVd) and *Grapevine yellow speckle viroid 1* (GSYVd-1) on grapevine. **The challenge addressed is the ability to detect closely related viroids, using Grapevine yellow speckle viroid 1 (GYSVd-1) and Grapevine yellow speckle viroid 2 (GYSVd-2)**. 117 sequences of GYSVd-2 have been downloaded from NCBI and aligned against the consensus sequence of GYSVd-1 from the real dataset. The GYSVd-2 reference DQ377131 has been selected to be added in the dataset 4. This reference shows a pairwise nucleotide identity of 73.9% with the consensus sequence of the naturally present GYSVd-1, a portion of the two genomes being very similar while the other part shows more differences (Figure 1).

![Figure1](images/Dataset4_Figure1_SimilarityPlotViroids_2.png)
*Figure 1: Similarity plot of the full-length sequence of GYSVd-1 and GYSVd-2 – DQ377131 reference, GYSVd-1 consensus sequence being the query sequence.*

### Artificial reads added to the dataset

An overview of the expected composition of the dataset can be found in Table 1.

*Table 1: Composition of Dataset 4. The total number of reads in this dataset is 2x10,054,658 (2x75 bp).*


|                         Virus/Viroid                        | Reads type | Number of   artificial reads added | Expected average number of reads per position |
|:-----------------------------------------------------------:|:----------:|:----------------------------------:|:---------------------------------------------:|
|              *Grapevine red blotch virus* (GRBV)              |    Real    |                                    |                                               |
| *Grapevine rupestris stem pitting-associated virus* (GRSPaV)  |    Real    |                                    |                                               |
|                   *Hop stunt viroid* (HSVd)                   |    Real    |                                    |                                               |
|          *Grapevine yellow speckle viroid 1* (GSYVd-1)          |    Real    |                                    |                                               |
|          *Grapevine yellow speckle viroid 2* (GSYVd-2)          | Artificial |                2306                |                     530.3                     |


The different labs that analyzed the dataset are listed in Table 2.

*Table 2: Participants to the VIROMOCK challenge of Dataset 4.*

| Participant                      	| Institute 	| Country 	| email                                  	|
|----------------------------------	|-----------	|---------	|----------------------------------------	|
| Lucie Tamisier                   	| ULg       	| Belgium 	| <lucie.tamisier@uliege.be>             	|
| Annelies Haegeman, Yoika Foucart 	| ILVO      	| Belgium 	| <annelies.haegeman@ilvo.vlaanderen.be> 	|
| Alex Hu | USDA-APHIS | USA | <xiaojun.hu@usda.gov> |
| Miguel Morard| ValGenetics SL | Spain | mmorard@valgenetics.com |

The observed values of the different participants are listed in Table 3.

*Table 3: Observed composition of Dataset 4 after analysis by different labs.*

| Institute/Lab | Observed closest   NCBI accession | Observed   proportion of GSYVd reads (%)<sup>1</sup> | Were   you able to detect the two viroids? |
|:-------------:|-----------------------------------|------------------------------------------|--------------------------------------------|
|      ILVO     | KY426921                          | NA                                       | yes                                        |
|      ILVO     | MG938315                          | NA                                       |                                            |
|      ILVO     | GQ995465                          | NA                                       |                                            |
|      ILVO     | KP010010                          | 41.84                                    |                                            |
|      ILVO     | MG780425                          | 58.16                                    |                                            |
| USDA-APHIS | MF276895.1 | NA   | yes |
| USDA-APHIS  | MG938313   | NA   |     |
|  USDA-APHIS  | X06873.1   | NA   |     |
|  USDA-APHIS  | KU880716.1 | 54.1 |     |
|  USDA-APHIS  | DQ377129   | 45.9 |     |
| ValGenetics | MF276895.1 | NA   | yes |
| ValGenetics | MG938313   | NA   |     |
| ValGenetics | X06873.1   | NA   |     |
| ValGenetics | KU880712.1 | 50.20 |     |
| ValGenetics | DQ377129.1   | 49.80 |     |

<sup>1</sup> This number represents the number of GYSVd filtered reads mapped for each accession against the total number of GYSVd filtered reads.

### Comments of different labs while analyzing the dataset

Here you can read some comments the participants had while analyzing the dataset.

| Institute/Lab     | Comments                                              |
|-------------------|-------------------------------------------------------|
| USDA, APHIS, PGQP | DQ377129   is the same as GYSVd-2 reference DQ377131. |

### How to participate to the VIROMOCK challenge

The dataset can be downloaded [here](https://doi.org/10.5061/dryad.0zpc866z8) (click the arrow next to "November 2, 2021" to download each dataset separately).

If you finish your analysis, we encourage you to submit your results through [this Google Sheet](https://docs.google.com/spreadsheets/d/1Bfx_34APkOxX23bYP8T299P737gXjmqI6zuHgfJw6aA/edit#gid=0).

The Google Sheet will allow you to share your results in detail. Only the green columns are required. However, we encourage you to give as much information as possible.

After submission of the Google Sheet, your results will be processed and added to Table 3.
