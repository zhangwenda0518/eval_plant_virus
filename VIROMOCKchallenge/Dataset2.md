## Dataset 2: mutations (CTV on Citrus)

### Introduction to the dataset

The real dataset is the same as the dataset 1, and is therefore composed of mixed infections of *Citrus tristeza virus* (CTV), *Citrus vein enation virus* (CVEV), *Citrus exocortis viroid* (CEVd), *Citrus viroid III* (CVd-III) and *Hop stunt viroid* (HSVd) on citrus. **The challenge addressed is the ability to identify different types of mutations showing different frequencies**. The creation of this dataset has been inspired by the second dataset of the Proficiency Test No 2.
The CTV - MH323442 strain has been selected and 5 identical sequences of this strain have been used. On each sequence, 1 substitution, 1 bp deletion and 1 bp insertion have been added. The MH323442 strain without mutations has also been added to the final dataset.

### Artificial reads added to the dataset

Table 1 summarizes the different haplotypes created and an overview of the expected composition of the dataset can be found in Table 2.

*Table 1: Description of the 5 haplotypes created.*

|     Haplotype    |     Substitution    |     Insertion    |     Deletion    |     ORF            |     Relative frequency (%)    |
|------------------|---------------------|------------------|-----------------|--------------------|-------------------------------|
|     1            |     A,2164,C        |     +C,2374      |     -A,2696     |     Polyprotein    |     24.39                     |
|     2            |     G,9549,T        |     +A,9867      |     -C,10051    |     RdRp           |     10.24                     |
|     3            |     G,12359,A       |     +T,12434     |     -G,12629    |     P65            |     5.37                      |
|     4            |     T,16242,G       |     +G,16292     |     -T,16619    |     P25            |     1.07                      |
|     5            |     A,17769,T       |     +A,17854     |     -G,17982    |     P20            |     0.39                      |


*Table 2: Composition of Dataset 2. The total number of reads in this dataset is 2x21,756,961 (2x150 bp).*

|            Virus/Viroid            | Reads type | Number of   artificial reads added | Expected average number of reads per position |
|:----------------------------------:|:----------:|:----------------------------------:|:---------------------------------------------:|
| *Citrus vein enation virus* (CVEV) |    Real    |                  -                 |                       -                       |
|  *Citrus exocortis viroid* (CEVd)  |    Real    |                  -                 |                       -                       |
|    *Citrus viroid III* (CVdIII)    |    Real    |                  -                 |                       -                       |
|       *Hop stunt viroid* (HSVd)      |    Real    |                  -                 |                       -                       |
|    *Citrus tristeza virus* (CTV)   |    Real    |                  -                 |                       -                       |
|           CTV - MH323442           | Artificial |               204312                |                     1596.1                    |

### Results of different labs analyzing the dataset

The different labs that analyzed the dataset are listed in Table 3.

*Table 3: Participants to the VIROMOCK challenge of Dataset 2.*

| Participant                      	| Institute 	| Country 	| email                                  	|
|----------------------------------	|-----------	|---------	|----------------------------------------	|
| Lucie Tamisier                   	| ULg       	| Belgium 	| <lucie.tamisier@uliege.be>             	|
| Annelies Haegeman, Yoika Foucart 	| ILVO      	| Belgium 	| <annelies.haegeman@ilvo.vlaanderen.be> 	|
| Alex Hu | USDA-PHIS | USA | <xiaojun.hu@usda.gov> |   |


The observed values of the different participants are listed in Table 4.

*Table 4: Observed mutations present in Dataset 2 for CTV after analysis by different labs.*

| Institute/Lab | Variant type | Variant   | Observed position | Observed   frequency | Were   you able to detect all SNPs? |
|:-------------:|--------------|-----------|-------------------|----------------------|-------------------------------------|
|      ILVO     | substitution | A,2163,C  | 2163              | 26.26                | no                                  |
|      ILVO     | substitution | G,9549,T  | 9549              | 9.64                 |                                     |
|      ILVO     | substitution | G,12359,A | 12359             | 6.97                 |                                     |
|      ILVO     | substitution | T,16242,G | 16242             | 0.78                 |                                     |
|      ILVO     | substitution | A,17769,T | 17769             | 0.77                 |                                     |
|      ILVO     | insertion    | C,2374    | 2373              | 22.68                |                                     |
|      ILVO     | insertion    | A,9867    | 9866              | 9.23                 |                                     |
|      ILVO     | insertion    | T,12434   | 12433             | 5.9                  |                                     |
|      ILVO     | insertion    | G,16292   | 16291             | 1.32                 |                                     |
|      ILVO     | insertion    | A,17854   | 1754              | 0.27                 |                                     |
|      ILVO     | deletion     | -A,2696   | 2695              | 21.01                |                                     |
|      ILVO     | deletion     | -C,10051  | 10050             | 10.88                |                                     |
|      ILVO     | deletion     | -G,12629  | 12628             | 4.74                 |                                     |
|      ILVO     | deletion     | -T,16619  | 17981             | 1.05                 |                                     |
|      ILVO     | deletion     | -G,17982  | 16618             | 0.77                 |                                     |
| USDA-APHIS | substitution |  A,2163,C |  2163 |   25  | yes |
| USDA-APHIS | substitution |  G,9549,T |  9549 |  9.63 |     |
| USDA-APHIS | substitution | G,12359,A | 12359 |  6.49 |     |
| USDA-APHIS | substitution | T,16242,G | 16242 |  0.97 |     |
| USDA-APHIS | substitution | A,17769,T | 17769 |  0.57 |     |
| USDA-APHIS |   insertion  |   C,2374  |  2373 | 23.76 |     |
| USDA-APHIS |   insertion  |   A,9867  |  9866 |  9.57 |     |
| USDA-APHIS |   insertion  |  T,12434  | 12433 |  5.78 |     |
| USDA-APHIS |   insertion  |  G,16292  | 16291 |  1.2  |     |
| USDA-APHIS |   insertion  |  A,17854  | 17853 |  0.19 |     |
| USDA-APHIS |   deletion   |  -A,2696  |  2694 |  21.9 |     |
| USDA-APHIS |   deletion   |  -C,10051 | 10049 | 10.61 |     |
| USDA-APHIS |   deletion   |  -G,12629 | 12627 |  5.31 |     |
| USDA-APHIS |   deletion   |  -T,16619 | 16617 |  0.94 |     |
| USDA-APHIS |   deletion   |  -G,17982 | 17980 |  0.68 |     |

### Comments of different labs while analyzing the dataset

Here you can read some comments the participants had while analyzing the dataset.

| Institute/Lab     | Comments                                                                                                                                                |
|-------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------|
| ILVO              | For SNP T,16242,G; A,17769,T; A,17854; -T,16619 and -G,17982: Not easy to identify, lots of other variants with similar frequencies from the real virus. |
| USDA-APHIS        | For SNP T,16242,G and A,17769,T: Detection at the noise level.                                                                                          |

### How to participate to the VIROMOCK challenge

The dataset can be downloaded [here](https://doi.org/10.5061/dryad.0zpc866z8) (click the arrow next to "November 2, 2021" to download each dataset separately).

If you finish your analysis, we encourage you to submit your results through [this Google Sheet](https://docs.google.com/spreadsheets/d/11oizMUYOkpiO5Zgds-HM7uB-lbN5FqtVLCph4K99gdY/edit#gid=0).

The Google Sheet will allow you to share your results in detail. Only the green columns are required. However, we encourage you to give as much information as possible.

After submission of the Google Sheet, your results will be processed and added to Table 4.
