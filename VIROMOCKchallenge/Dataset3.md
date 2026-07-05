## Dataset 3: Viral concentration at the species level and non-complete genome (mixed infections on Grapevine)

### Introduction to the dataset

The real dataset is composed of mixed infections of *Grapevine rupestris vein feathering virus* (GRVFV), *Grapevine rupestris stem pitting-associated virus* (GRSPaV), *Grapevine leafroll-associated virus 2* (GLRaV2), *Hop stunt viroid* (HSVd) and *Grapevine yellow speckle viroid 1* (GYSVd1) on grapevine. **The challenge addressed is the ability to detect several viral/viroid species showing different frequencies and incomplete genome coverages**. After the analyse of the real dataset with [Geneious](https://www.geneious.com/), 35,729 reads were mapped on the 3 viral genomes (GRVFV, GRSPaV and GLRaV2), the number of mapped reads ranging between 6,486 for GRVFV and 21,570 for GLRaV2. In order to obtain incomplete and/or low coverages, 4,000 viral reads among the 35,729 mapped have been randomly selected, using the “Randomly Sample Sequences” tool in [Geneious](https://www.geneious.com/). The old mapped reads have been removed from the dataset, and the new 4,000 selected reads have been added. The final coverages obtained for the 3 viruses (GRSPaV, GRVFV and GLRaV2) are shown in Figure 1.

![Figure1](images/Dataset_3_Figure1_Coverages.png)
*Figure 1: Mapping of filtered reads of dataset 3 against (a) GRSPaV, (b) GRVFV and (c) GLRaV2 references.*

### Artificial reads added to the dataset

An overview of the expected composition of the dataset can be found in Table 1.

*Table 1: Composition of Dataset 3. The total number of reads in this dataset is 2x24,526,416 (2x150 bp).*

|                         Virus/Viroid                         |      Reads type      |
|:------------------------------------------------------------:|:--------------------:|
| *Grapevine rupestris stem pitting-associated virus* (GRSPaV) | Real, but rarified |
|      *Grapevine rupestris vein feathering virus* (GRVFV)     |  Real, but rarified  |
|      *Grapevine leafroll-associated virus 2* (GLRaV2)        |  Real, but rarified  |
|                    *Hop stunt viroid* (HSVd)                   |         Real         |
|          *Grapevine yellow speckle viroid 1* (GYSVd1)          |         Real         |


The different labs that analyzed the dataset are listed in Table 2.

*Table 2: Participants to the VIROMOCK challenge of Dataset 3.*

| Participant                      	| Institute 	| Country 	| email                                  	|
|----------------------------------	|-----------	|---------	|----------------------------------------	|
| Lucie Tamisier                   	| ULg       	| Belgium 	| <lucie.tamisier@uliege.be>             	|
| Annelies Haegeman, Yoika Foucart 	| ILVO      	| Belgium 	| <annelies.haegeman@ilvo.vlaanderen.be> 	|
| Alex Hu | USDA-APHIS | USA | <xiaojun.hu@usda.gov> |
| Miguel Morard| ValGenetics SL | Spain | mmorard@valgenetics.com |


The observed values of the different participants are listed in Table 3.

*Table 3: Observed composition of Dataset 3 after analysis by different labs.*



| Institute/Lab | Observed closest   NCBI accession | Observed   proportion of VIRAL reads (%)<sup>1</sup> | Observed   proportion of VIROID reads (%)<sup>1</sup> | Were you able to   detect all viruses/viroids? |
|:-------------:|-----------------------------------|------------------------------------------|-------------------------------------------|------------------------------------------------|
|      ILVO     | KX274274.1                        | 19                                       | NA                                        | yes                                            |
|      ILVO     | MF000326.1                        | 32.5                                     | NA                                        |                                                |
|      ILVO     | JX513891.1                        | 48.5                                     | NA                                        |                                                |
|      ILVO     | MF576417.1                        | NA                                       | 69.35                                     |                                                |
|      ILVO     | KP010010.1                        | NA                                       | 30.65                                     |                                                |
| USDA-APHIS | KX274274.1 | 10.65 | NA    | yes |
| USDA-APHIS | KY513702.1 | 16.21 | NA    |     |
| USDA-APHIS | MN548394.1 | 73.1  | NA    |     |
| USDA-APHIS | KR909028.1 | NA    | 71.66 |     |
| USDA-APHIS | EU682454.1 | NA    | 28.34 |     |
| ValGenetics| KX274274.1 | 14.13 | NA    | yes |
| ValGenetics | AY706994.1 | 36.96 | NA    |     |
| ValGenetics | MT899925.1 | 48.91  | NA    |     |
| ValGenetics | KJ466332.1 | NA    | 71.97 |     |
| ValGenetics | KU880715.1 | NA    | 28.03 |     |

<sup>1</sup> This number represents the number of virus (viroid) filtered reads mapped for each accession against the total number of virus (viroid) filtered reads.


### Comments of different labs while analyzing the dataset

Here you can read some comments the participants had while analyzing the dataset.

No comments yet.

### How to participate to the VIROMOCK challenge

The dataset can be downloaded [here](https://doi.org/10.5061/dryad.0zpc866z8) (click the arrow next to "November 2, 2021" to download each dataset separately).

If you finish your analysis, we encourage you to submit your results through [this Google Sheet](https://docs.google.com/spreadsheets/d/1KxO4pvvJFbbRgyprVyCIFX9mygraMcNzCi376cjsdr8/edit?usp=drive_web&ouid=101686715531089518912).

The Google Sheet will allow you to share your results in detail. Only the green columns are required. However, we encourage you to give as much information as possible.

After submission of the Google Sheet, your results will be processed and added to Table 3.
