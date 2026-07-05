## Dataset 6: New strain (PVY on Potato)

### Introduction to the dataset

This dataset is the same as [Dataset 5](Datasets/Dataset5.md), and is therefore composed of one PVY NTN isolate and the PVY N605 strain on potato. **The challenge addressed is the ability to detect a new PVY strain that does not exist in the database, within a dataset already composed of PVY strains**.
The creation of the new PVY strain has been performed as follow:

- First, 397 PVY references belonging to different clades have been downloaded from NCBI. They have been aligned using the [MAFFT](https://www.ebi.ac.uk/Tools/msa/mafft/) alignment program on [Geneious](https://www.geneious.com/), and the FJ214726 strain has been selected for the next steps. This strain has been chosen because it was one of the strains showing the highest number of differences at the nucleotide level (never more than 84% of identity with the other strains).

- PVY encodes a polyprotein. The nucleotide sequence encoding the polyprotein of the FJ214726 strain has first been translated in amino acids. Then, this amino acid sequence has been reverse translated into nucleotides, using the online [EMBOSS Backtranseq tool](https://www.ebi.ac.uk/Tools/st/emboss_backtranseq/). Thanks to the degeneracy of the genetic code, the nucleotide sequence obtained was different from the one of FJ214726.

- The artificial stain obtained was still identical to FJ214726 at the amino acid level. Therefore, a third strain, MH795855 has been used. This strain shows 83.4% identity with FJ214726 at the nucleotide level. A multiple alignment has been performed using [Muscle](https://www.ebi.ac.uk/Tools/msa/muscle/) on [Geneious](https://www.geneious.com/) between the three strains. The 5’ and 3’ UTR regions of MH795855 have been added to the artificial PVY strain. Then, non-synonymous substitutions have been manually added to the new artificial strain. Thanks to the alignment, we could verify that these mutations were non-synonymous and did not generated stop codons. If three different nucleotides were present at a position between the three strains, we sometimes change the nucleotide of the artificial strain into the one of MH795855, so that the new strain differs from FJ214726, but still shows a PVY-like sequence.

The PVY artificial strain obtained shows a percentage of identity with the 397 other PVY strains ranging between 71.862 and 68.062% at the nucleotide level. It is less divergent at the amino acid level, showing a percentage of identity ranging between 98.889 and 83.318%. Figure 1 shows the phylogenetic trees built on [Geneious](https://www.geneious.com/) using [PHYML](https://github.com/stephaneguindon/phyml), with the nucleotide and amino acid alignments of 57 PVY references representative of PVY diversity and the new artificial PVY strain.

![Figure1](images/Dataset6_Figure1_PhylogeneticTreesPVY_2.png)
*Figure 1: Phylogenetic tree of PVY. (a) Tree generated using nucleotide sequences, (b) tree generated using amino acid sequences. The new artificial strain is surrounded in red.*

### Artificial reads added to the dataset

An overview of the expected composition of the dataset can be found in Table 1.

*Table 1: Composition of Dataset 6. The total number of reads in this dataset is 31,327,327 (1x50 bp).*

|        Virus/Viroid       | Reads type | Number of artificial reads added | Expected proportion of artificial   PVY reads (%) | Expected average number of reads   per position |
|:-------------------------:|:----------:|:--------------------------------:|:-------------------------------------------------:|:-----------------------------------------------:|
|          PVY NTN          |    Real    |                                  |                                                   |                                                 |
|          PVY N605         |    Real    |                                  |                                                   |                                                 |
| PVY new artificial strain | Artificial |              199668              |                         27                        |                      1034.6                     |

The sequence of the artificial PVY strain can be downloaded [here](https://gitlab.com/ahaegeman/PHBN-WP3-VIROMOCKchallenge/-/blob/master/Datasets/fasta/Potato_virus_Y_artificial_strain.fasta).

The different labs that analyzed the dataset are listed in Table 2.

*Table 2: Participants to the VIROMOCK challenge of Dataset 6.*

| Participant                      	| Institute 	| Country 	| email                                  	|
|----------------------------------	|-----------	|---------	|----------------------------------------	|
| Lucie Tamisier                   	| ULg       	| Belgium 	| <lucie.tamisier@uliege.be>             	|
| Annelies Haegeman, Yoika Foucart 	| ILVO      	| Belgium 	| <annelies.haegeman@ilvo.vlaanderen.be> 	|
| Alex Hu | USDA-APHIS | USA | <xiaojun.hu@usda.gov> |
| Miguel Morard| ValGenetics SL | Spain | mmorard@valgenetics.com |

The observed values of the different participants are listed in Table 3.

*Table 3: Observed composition of Dataset 6 after analysis by different labs.*


| Institute/Lab | Virus/viroid     | Length   of assembled new strain(bp) | %identity   of assembled sequence to artificially sequence of the new strain | Were you able to   notice that we are dealing with a new strain? |
|:-------------:|------------------|--------------------------------------|------------------------------------------------------------------------------|------------------------------------------------------------------|
|      ILVO     | PVY   NTN        | NA                                   | NA                                                                           | no                                                               |
|      ILVO     | PVY   N605       | NA                                   | NA                                                                           |                                                                  |
|      ILVO     | PVY   artificial |                                      |                                                                              |                                                                  |
| USDA-APHIS | PVY NTN          | NA   | NA    | yes |
|  USDA-APHIS  | PVY   N605       | NA   | NA    |     |
|  USDA-APHIS  | PVY   artificial | 9666 | 73.57 |     |
| ValGenetics | PVY NTN          | NA   | NA    | yes |
|  ValGenetics  | PVY   N605       | NA   | NA    |     |
|  ValGenetics  | PVY   artificial | 9276 | 73.56 |     |

### Comments of different labs while analyzing the dataset

Here you can read some comments the participants had while analyzing the dataset.

| Institute/Lab | Comments                                                                                            |
|---------------|-----------------------------------------------------------------------------------------------------|
| ILVO          | We   did not notice the new strain, hence we did not assign assembled contigs to   each PVY strain. |

### How to participate to the VIROMOCK challenge

The dataset can be downloaded [here](https://doi.org/10.5061/dryad.0zpc866z8) (click the arrow next to "November 2, 2021" to download each dataset separately).

If you finish your analysis, we encourage you to submit your results through [this Google Sheet](https://docs.google.com/spreadsheets/d/1ILZWDYJaEIiBWrJLrMGhQZIOROVKUD4fx_2vG8fVpyc/edit#gid=0).

The Google Sheet will allow you to share your results in detail. Only the green columns are required. However, we encourage you to give as much information as possible.

After submission of the Google Sheet, your results will be processed and added to Table 3.
