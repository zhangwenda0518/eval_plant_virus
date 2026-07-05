## Dataset 7: Non-complete genome (TSWV on Tobacco)

### Introduction to the dataset

The real dataset is composed of two variants of *Tomato spotted wilt virus* (TSWV) from tobacco. The genome of TSWV consists of 3 negative ssRNA named S, M and L. The variants diverge only for the L genomic segment, one being of normal length (8,913 bp) and the other being a shorter defective (2,612 bp). The defective variant lacks the genomic part from 760 to 7,060 bp. There is no evidence of variants for the segments M (4,827 bp) and S (2,923 bp).
The real dataset shows already a challenging composition, and has therefore not been spiked with artificial viruses. **The challenge addressed is therefore the ability to detect both a defective and a normal length variant**. Figure 1 shows the reads coverage distribution along the genome of the TSWV.


![Figure1](images/Dataset7_Figure1_CoveragesTSWV2_2.png)
*Figure 1: Reads coverage distribution along the genome of the TSWV, using the NCBI references (a) MF159047 for the L segment, (b) KT717692 for the M segment and (c) KU179543 for the S segment.*


The different labs that analyzed the dataset are listed in Table 1.

*Table 1: Participants to the VIROMOCK challenge of Dataset 7.*

| Participant                      	| Institute 	| Country 	| email                                  	|
|----------------------------------	|-----------	|---------	|----------------------------------------	|
| Lucie Tamisier                   	| ULg       	| Belgium 	| <lucie.tamisier@uliege.be>             	|
| Annelies Haegeman, Yoika Foucart 	| ILVO      	| Belgium 	| <annelies.haegeman@ilvo.vlaanderen.be> 	|
| Alex Hu | USDA-APHIS | USA | <xiaojun.hu@usda.gov> |   |
| Steven Sewe | NRI | UK | <ss0291e@gre.ac.uk> |   |

The observed values of the different participants are listed in Table 2.

*Table 2: Observed composition of Dataset 7 after analysis by different labs.*

| Institute/Lab | Virus/viroid                  | Observed   closest NCBI accession | Observed   number of reads mapped | Were   you able to identify the deleted region from the defective variant? |
|:-------------:|-------------------------------|-----------------------------------|-----------------------------------|----------------------------------------------------------------------------|
|      ILVO     | TSWV-S                        | KU179543                          | 4752                              | yes                                                                        |
|      ILVO     | TSVW-M                        | KT717692                          | 2121                              |                                                                            |
|      ILVO     | TSWV-L   (normal + defective) | MF159047                          | 10107                             |                                                                            |
| USDA-APHIS | TSWV-S                        | DQ915948 | 30408  | yes |
|  USDA-APHIS  | TSVW-M                        | AY744486 | 6290   |     |
|  USDA-APHIS  | TSWV-L   (normal + defective) | D10066   | 151801 |     |
|         NRI         | TSWV-S                        |          | 28082  |   yes  |
|         NRI         | TSVW-M                        |          | 5924   |     |
|         NRI         | TSWV-L   (normal + defective) | D10066   | 70542  |  |

### Comments of different labs while analyzing the dataset

Here you can read some comments the participants had while analyzing the dataset.

No comments yet.

### How to participate to the VIROMOCK challenge

The dataset can be downloaded [here](https://doi.org/10.5061/dryad.0zpc866z8) (click the arrow next to "November 2, 2021" to download each dataset separately).

If you finish your analysis, we encourage you to submit your results through [this Google Sheet](https://docs.google.com/spreadsheets/d/1xvA5AoRrt-yCBTJ1PxlrHt3INtNdP3Mh6L3H7tIleyo/edit#gid=0).

The Google Sheet will allow you to share your results in detail. Only the green columns are required. However, we encourage you to give as much information as possible.

After submission of the Google Sheet, your results will be processed and added to Table 2.
