## Dataset 10: New strain (PBNSPaV + PPV on Prunus)

### Introduction to the dataset

The real dataset is composed of *Plum bark necrosis stem pitting-associated virus* (PBNSPaV) from prunus. **The limit tested is the ability to detect a new viral strain that does not exist in the database, adding a virus that is not already present in the dataset** (contrary to [Dataset 6](Datasets/Dataset6.md)).

We have chosen to create and add a new strain of *Plum pox virus* (PPV), because it infects prunus and codes for a polyprotein, so we can use the same method as for the creation of the new PVY strain in [Dataset 6](Datasets/Dataset6.md). The creation of the new PPV strain has been performed as follows:

- A sub-sample of 111 PPV references have been downloaded from NCBI. We could not run the analyses on all the PPV sequences, the number being too large. They have been aligned using the [MAFFT](https://www.ebi.ac.uk/Tools/msa/mafft/) alignment program on [Geneious](https://www.geneious.com/), and the MK387313 strain has been selected for the next steps. This strain has been chosen because it was one of the strains showing the highest number of differences at the nucleotide level (never more than 82% of identity with the other strains).

- The nucleotide sequence encoding the polyprotein of the MK387313 strain has first been translated in amino acids. Then, this amino acid sequence has been reverse translated into nucleotides, using the online [EMBOSS Backtranseq tool](https://www.ebi.ac.uk/Tools/st/emboss_backtranseq/). Thanks to the degeneracy of the genetic code, the nucleotide sequence obtained was different from the one of MK387313.

- The artificial stain obtained was still identical to MK387313 at the amino acid level. Therefore, a third strain, AB545926 has been used. This strain shows 78.6% identity with MK387313 at the nucleotide level. A multiple alignment has been performed using [Muscle](https://www.ebi.ac.uk/Tools/msa/muscle/) on [Geneious](https://www.geneious.com/) between the three strains (MK387313, the artificial strain and AB545926). The 5’ and 3’ UTR regions of AB545926 have been added to the artificial PPV strain. Then, non-synonymous substitutions have been manually added to the new artificial strain. Thanks to the alignment, we could verify that these mutations were non-synonymous and did not generated stop codons. If three different nucleotides were present at a position between the three strains, we sometimes changed the nucleotide of the artificial strain into the one of AB545926, so that the new strain differs from MK387313, but still shows a PPV-like sequence.

The PPV artificial strain obtained shows a percentage of identity with the 111 other PPV strains ranging between 73.532 and 69.128% at the nucleotide level. It is less divergent at the amino acid level, showing a percentage of identity ranging between 98.834 and 88.864%. Figure 10 shows the phylogenetic trees built on [Geneious](https://www.geneious.com/) using [PHYML](https://github.com/stephaneguindon/phyml), with the nucleotide and amino acid alignments of the 111 PPV references and the new artificial PPV strain

![Figure1](images/Dataset10_Figure1_PhylogeneticTreesPPV_2.png)
*Figure 1: Phylogenetic tree of PPV. (a) Tree generated using nucleotide sequences, (b) tree generated using amino acid sequences. The new PPV artificial strain is surrounded in red.*

### Artificial reads added to the dataset

An overview of the expected composition of the dataset can be found in Table 1.

*Table 1: Composition of Dataset 10. The total number of reads in this dataset is 2x24,573,681 (1x75 bp).*

| Virus/Viroid | Reads type | Number of artificial reads added | Expected   proportion of viral reads (%) 1 | Expected average number of reads   per position |
|:------------:|:----------:|:--------------------------------:|:------------------------------------------:|:-----------------------------------------------:|
|    PBNSPaV   |    Real    |                                  |                                            |                                                 |
|      PPV     | Artificial |               6002               |                    49.2                    |                        46                       |

The sequence of the artificial PPV strain can be downloaded [here](https://gitlab.com/ahaegeman/PHBN-WP3-VIROMOCKchallenge/-/blob/master/Datasets/fasta/Plum_pox_virus_artificial_strain.fasta).

The different labs that analyzed the dataset are listed in Table 2.

*Table 2: Participants to the VIROMOCK challenge of Dataset 10.*

| Participant                      	| Institute 	| Country 	| email                                  	|
|----------------------------------	|-----------	|---------	|----------------------------------------	|
| Lucie Tamisier                   	| ULg       	| Belgium 	| <lucie.tamisier@uliege.be>             	|
| Annelies Haegeman, Yoika Foucart 	| ILVO      	| Belgium 	| <annelies.haegeman@ilvo.vlaanderen.be> 	|
| Alex Hu 	| USDA-APHIS      	| USA 	| <xiaojun.hu@usda.gov> 	|
| Miguel Morard 	| ValGenetics SL     	| Spain 	| <mmorard@valgenetics.com> 	|


The observed values of the different participants are listed in Table 3.

*Table 3: Observed composition of Dataset 10 after analysis by different labs.*

| Institute/Lab | Virus/viroid | Length   of assembled new virus (bp) | %identity   of assembled sequence to artificially sequence of the new virus | Observed proportion of viral reads (%) | Were   you able to notice that there is a new virus present? |
|:-------------:|--------------|--------------------------------------|-----------------------------------------------------------------------------|--------------------------------------------------------------|------|
|      ILVO     | PBNSPaV      | NA                                   | NA   | 52                                                                       | yes                                                          |
|      ILVO     | new   virus  | 9785                                 | to   be completed       | 48                                                    |                                                              |
|      USDA-APHIS     | PBNSPaV      | NA                                   | NA          | 90.15                                                                | yes                                                          |
|      USDA-APHIS     | new   virus  | 9785                                 | 100                                                           |  9.85 |
|      ValGenetics    | PBNSPaV      | NA                                   | NA          | 89.92                                                                | yes                                                          |
|      ValGenetics    | new   virus  | 9785                                 | 100                                                           |  10.08 |

### Comments of different labs while analyzing the dataset

Here you can read some comments the participants had while analyzing the dataset.

No comments yet.

### How to participate to the VIROMOCK challenge

The dataset can be downloaded [here](https://doi.org/10.5061/dryad.0zpc866z8) (click the arrow next to "November 2, 2021" to download each dataset separately).

If you finish your analysis, we encourage you to submit your results through  [this Google Sheet](https://docs.google.com/spreadsheets/d/1eBBisPW0yqw_548EK-9BYhLjaGHqFaoeoSNt-Lha_AQ/edit#gid=0).

The Google Sheet will allow you to share your results in detail. Only the green columns are required. However, we encourage you to give as much information as possible.

After submission of the Google Sheet, your results will be processed and added to Table 3.
