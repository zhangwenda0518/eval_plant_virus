## Dataset 5: Recombination (PVY on Potato)

### Introduction to the dataset

The real dataset is composed of two *Potato virus Y* (PVY) strains, an NTN strain and the N605 isolate on potato. **The challenge addressed is the ability to detect one recombinant strain and one parental one**. We have chosen to use only one and not two parental strains because the presence of both parental strains plus the recombinant one could have made the challenge too difficult, especially for the identification of the recombinant strain. The parental and recombinant strains have been selected based on the publication of [Hu et al. (2009)](https://www.microbiologyresearch.org/docserver/fulltext/jgv/90/12/3033.pdf?expires=1584980064&id=id&accname=guest&checksum=344B6D652F4BE5D43173AF58746D9423). The parental isolate is AY884983.1, which is an N strain, and the recombinant isolate is EF026076.1, which is a recombinant between an N and an O strains. Both isolates show a pairwise identity of 88.2% at the nucleotide level. The beginning of the genomes (first 2000 nucleotides) are almost the same, while the rest of the genomes show more differences between both sequences (Figure 1). Another challenge of this dataset is the presence of two close PVY isolates. Indeed, the already present N605 isolate and the artificially added AY884983.1 isolate show a pairwise identity around 99% at the nucleotide level.

![Figure1](images/Dataset5_Figure1_SimilarityPlotPVY_2.png)
*Figure 1: Similarity plot of the full-length sequence of PVY AY884983 and EF026076 references, AY884983 being the query sequence.*

### Artificial reads added to the dataset

An overview of the expected composition of the dataset can be found in Table 1.

*Table 1: Composition of Dataset 5. The total number of reads in this dataset is 2x31,277,475 (1x50 bp).*

| Virus/Viroid | Reads type | Number of   artificial reads added | Expected average number of reads per position  |
|:------------:|:----------:|:----------------------------------:|:----------------------------------------------:|
|    *Potato virus Y* (PVY) NTN  |    Real    |                                    |                                                |
|   PVY N605   |    Real    |                                    |                                                |
|      PVY (EF026076.1)    | Artificial |                119874               |                      413.7                     |
|      PVY (AY884983.1)    | Artificial |                29942               |                       103                      |

The sequence of EF026076 PVY strain can be downloaded [here](https://gitlab.com/ahaegeman/PHBN-WP3-VIROMOCKchallenge/-/blob/master/Datasets/fasta/PVY_EF026076_1.fasta).

The sequence of AY884983 PVY strain can be downloaded [here](https://gitlab.com/ahaegeman/PHBN-WP3-VIROMOCKchallenge/-/blob/master/Datasets/fasta/PVY_AY884983_1.fasta).

The different labs that analyzed the dataset are listed in Table 2.

*Table 2: Participants to the VIROMOCK challenge of Dataset 5.*

| Participant                      	| Institute 	| Country 	| email                                  	|
|----------------------------------	|-----------	|---------	|----------------------------------------	|
| Lucie Tamisier                   	| ULg       	| Belgium 	| <lucie.tamisier@uliege.be>             	|
| Annelies Haegeman, Yoika Foucart 	| ILVO      	| Belgium 	| <annelies.haegeman@ilvo.vlaanderen.be> 	|
| Alex Hu | USDA-APHIS | USA | <xiaojun.hu@usda.gov> |   |


The observed values of the different participants are listed in Table 3.

*Table 3: Observed composition of Dataset 5 after analysis by different labs.*

| Institute/Lab | Observed closest   NCBI accession | Observed   length coverage (%) | Observed   proportion of PVY reads (%)<sup>1</sup> | Were   you able to notice that we are dealing with a recombinant strain? |
|:-------------:|-----------------------------------|--------------------------------|----------------------------------------|--------------------------------------------------------------------------|
|      ILVO     | KM396648                          | 47.9                           | 11.6                                   | yes                                                                      |
|      ILVO     | X97895.1                            | 38.1                           | 10.9                                   |                                                                          |
|      ILVO     | EF026076.1                          | 70.2                           | 58.4                                   |                                                                          |
|      ILVO     | AY884983.1                          | 44.6                           | 19.1                                   |                                                                          |
| USDA-APHIS | KM396648 | 99.66 | 20.21 | yes |
|  USDA-APHIS  | X97895   | 99.98 | 35.58 |     |
|  USDA-APHIS  | EF026076 | 99.84 | 16.98 |     |
|  USDA-APHIS  | AY884983 | 100   | 27.23 |     |

<sup>1</sup> This number represents the number of PVY filtered reads mapped for each accession against the total number of PVY filtered reads.


### Comments of different labs while analyzing the dataset

Here you can read some comments the participants had while analyzing the dataset.

| Institute/Lab     | Comments                                                                                                     |
|-------------------|--------------------------------------------------------------------------------------------------------------|
| USDA, APHIS, PGQP | The identity between X97895 and AY884983 is 99%. It's hard to separate them. Difficult to assemble AY884983. |

### How to participate to the VIROMOCK challenge

The dataset can be downloaded [here](https://doi.org/10.5061/dryad.0zpc866z8) (click the arrow next to "November 2, 2021" to download each dataset separately).

If you finish your analysis, we encourage you to submit your results through [this Google Sheet](https://docs.google.com/spreadsheets/d/1WKzS5q2NFT8FNZDlvJQRuqaGjJUswcTAKyvbuP6V1Eg/edit#gid=0).

The Google Sheet will allow you to share your results in detail. Only the green columns are required. However, we encourage you to give as much information as possible.

After submission of the Google Sheet, your results will be processed and added to Table 3.
