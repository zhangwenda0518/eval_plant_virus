## PHBN-WP3-VIROMOCKchallenge

### Introduction

The goal of the WP3 of the **Plant Health Bioinformatics Network** (PHBN) project is to help researchers to compare and validate their bioinformatics pipelines for virus detection. Semi-artificial datasets have been created for this purpose. They are composed of a “real” HTS dataset spiked with artificial viral reads. It will allow researchers to adjust their pipeline/parameters as good as possible to approximate the actual viral composition of the semi-artificial datasets. This will happen in the frame of the **VIROMOCK challenge**. Each **(semi-)artificial** dataset addresses one or several challenges that could prevent virus detection or a correct virus identification from HTS data (i.e. low viral concentration, new viral species, non-complete genome). This document describes the general pipeline that has been used to create the semi-artificial datasets and presents each one of them. Another challenge in the plant virology community is the ability to reconstruct viral haplotypes. Completely artificial datasets only composed of viral reads have been created for this purpose, and are described at the end of this document.

### Creation of the (semi-)artificial datasets

The (semi-)artificial datasets of this challenge are based on real datasets which are spiked with artificial data. More detailed info on how these datasets were created can be found [here](Dataset_creation).

### Participation to the VIROMOCK challenge

The (semi-)artificial datasets are **open access** and hence you can use them easily to validate your way of working, test new software, etc. The exact composition of the artificial reads is always shared together with the dataset, hence you can see for yourself whether or not your pipeline is behaving as expected, and what aspects of the analysis are difficult to tackle. All information is found on the dataset page, which can be accessed by clicking the dataset link in Table 1. All datasets can be downloaded in batch [here](https://doi.org/10.5061/dryad.0zpc866z8). We encourage you to also **give us some feedback** on what your results were for each dataset you analyze, and what difficulties you experienced by completing a **Google Form** is added to each dataset page. **General comments, suggestions or questions** can be discussed on our [![Gitter](https://badges.gitter.im/ahaegeman-PHBN-WP3-VIROMOCKchallenge/community.svg)](https://gitter.im/ahaegeman-PHBN-WP3-VIROMOCKchallenge/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge).

### Semi-artificial and real datasets of the VIROMOCK challenge

A total of **10 semi-artificial and real datasets** have been created or selected for the challenge. Seven are semi-artificial datasets and three are real datasets. Table 1 summarizes the main characteristics of each dataset.
Each dataset adresses at least **one challenge** using one or several viruses/viroids. Since we used real data, other viruses/viroids can be part of the background of the datasets and, in most cases, their presence was confirmed by molecular methods. All the viruses/viroids present are listed in the tables describing the composition of the datasets. However, these “background” viruses have not been investigated at the strain or SNP level, because it was not the goal of the datasets to analyse these viruses.

Each dataset name is clickable and will lead to you a page with more information on that dataset. There you can read more details about the construction and/or composition of the dataset, analysis results of different labs, and how to participate to the challenge yourself. By participating, your analysis results and/or comments will appear on the dedicated dataset page for comparison to results from other labs.

*Table 1: Characteristics of each semi-artificial and real dataset.*

|     Dataset                               	|     Researcher, institute, country     	|     Dataset type       	|     Plant species    	|     Virus/Viroids already present           	|     Modification                    	|     Reads (bp)                     	|     Total number of reads                 	|     Challenge                                                              	|
|-------------------------------------------	|----------------------------------------	|------------------------	|----------------------	|---------------------------------------------	|-------------------------------------	|------------------------------------	|-------------------------------------------	|-------------------------------------------------------------------------	|
|     [Dataset1](Datasets/Dataset1.md)      	|     Kris De Jonghe, ILVO, BE           	|     Semi-artificial    	|     Citrus           	|     CTV, CVEV, CEVd, CVd-III, HSVd          	|     Addition of CTV                 	|     2 x 150                        	|     21,703,434 (R1)    21,703,434 (R2)    	|     Different viral concentration (CTV strains) |
|     [Dataset2](Datasets/Dataset2.md)      	|     Kris De Jonghe, ILVO, BE           	|     Semi-artificial    	|     Citrus           	|     CTV, CVEV, CEVd, CVd-III, HSVd          	|     Addition of CTV                 	|     2 x 150                        	|     21,756,961 (R1)    21,756,961 (R2)    	|     Mutation present in different frequencies (CTV haplotypes)                                                            	|
|     [Dataset3](Datasets/Dataset3.md)      	|     Marie Lefebvre, INRA, FR           	|     Semi-artificial    	|     Grapevine        	|     GRSPaV, GLRaV2, GRVFV, HSVd, GYSVd1     	|     Removing of real viral reads    	|     2 x 150                        	|     24,526,416 (R1)    24,526,416 (R2)    	|     Different viral concentration (at the species level) + Non complete virus genome coverage (GRSPaV, GLRaV2 and GRVFV)|
|     [Dataset4](Datasets/Dataset4.md)      	|     Jean-Sébastien Reynard, AGS, CH    	|     Semi-artificial    	|     Grapevine        	|     GRBV, GRSPaV, HSVd, GYSVd-1             	|     Addition of GYSVd-2             	|     2 x 75                         	|     10,054,658 (R1)    10,054,658 (R2)    	|     Viroids with very similar sequence (GYSVd1 and GYSVd2)                                                       	|
|     [Dataset5](Datasets/Dataset5.md)      	|     Denis Kutnjak, NIB, SI             	|     Semi-artificial    	|     Potato           	|     PVY                                     	|     Addition of PVY                 	|     1 x 50                         	|     31,277,475                            	|     Mix of recombinant and parental viral PVY strains                                                       	|
|     [Dataset6](Datasets/Dataset6.md)      	|     Denis Kutnjak, NIB, SI             	|     Semi-artificial    	|     Potato           	|     PVY                                     	|     Addition of PVY                 	|     1 x 50                         	|     31,327,327                            	|     New PVY strain                                                          	|
|     [Dataset7](Datasets/Dataset7.md)      	|     Paolo Margaria, DSMZ, DE           	|     Real               	|     Tobacco          	|     TSWV                                    	|     -                               	|     2 x 301                        	|     1,904,369 (R1)    1,904,369 (R2)      	|     Complete genome + defective form of TSWV                                                 	|
|     [Dataset8](Datasets/Dataset8.md)      	|     Paolo Margaria, DSMZ, DE           	|     Real               	|     Chenopodium      	|     PFBV + mitovirus                        	|     -                               	|     2 x 301                        	|     65,177 (R1)    65,177 (R2)            	|     Cryptic mitovirus virus + low mitovirus concentration                                   	|
|     [Dataset9](Datasets/Dataset9.md)      	|     Nihal Buzkan, UCDAVIS, USA         	|     Real               	|     Pistacio         	|     PiVB                                    	|     -                               	|     1 x 151 (R1)    1 x 84 (R2)    	|     5,259,903 (R1)    5,259,903 (R2)      	|     Concentration of different PiVB genomic segments                   	|
|     [Dataset10](Datasets/Dataset10.md)    	|     Kristian Stevens, UCDAVIS, USA     	|     Semi-artificial    	|     Prunus           	|     PBNSPaV                                 	|     Addition of PPV                 	|     1 x 75                         	|     24,573,681                            	|     New PBNSPaV strain                                                          	|

### Completely artificial datasets

Finally, eight **artificial datasets** only composed of viral reads (no background data) were created. Each dataset consists of a mix of several isolates from the same viral species showing different frequencies. The viral species were selected to be as divergent as possible. Therefore, the selected viruses have: (i) a genome composed of DNA or RNA, (ii) a single or double-stranded genome, (iii) a linear, circular and/or segmented genome, and (iv) show a genome length ranging from 2.8 to 17.1 kb. Viral reads of 150 bp have been synthesized using the ART software, as described [here](Dataset_creation). The reads have been directly derived from the NCBI references and **no SNPs were added**. These datasets can be used for example to test haplotype reconstruction software, the goal being to reconstruct all the isolates present in a dataset. The datasets are presented in Table 2, more information on the exact composition of each dataset can be found by clicking the dataset links in the table.

*Table 2: Characteristics of each artificial dataset.*

| Dataset                            	| Virus                                 	| Number of isolates added 	|
|------------------------------------	|---------------------------------------	|-------------------------	|
| [Dataset11](Datasets/Dataset11.md) 	| *Pepino mosaic virus* (PepMV)         	| 6                       	|
| [Dataset12](Datasets/Dataset12.md) 	| *Cassava mosaic virus*                	| 4                       	|
| [Dataset13](Datasets/Dataset13.md) 	| *Banana streak virus* (BSV)           	| 6                       	|
| [Dataset14](Datasets/Dataset14.md) 	| *Potato virus Y* (PVY)                	| 5                       	|
| [Dataset15](Datasets/Dataset15.md) 	| *Eggplant mottled dwarf virus* (EMDV) 	| 3                       	|
| [Dataset16](Datasets/Dataset16.md) 	| *Bell pepper endornavirus* (BPEV)     	| 4                       	|
| [Dataset17](Datasets/Dataset17.md) 	| *Little cherry virus 1* (LChV1)       	| 5                       	|
| [Dataset18](Datasets/Dataset18.md) 	| *Barley yellow dwarf virus* (BYDV)    	| 6                       	|

### Overview of the challenges

Figure 1 provides an overview of the 18 semi-artificial, real and completely artificial datasets as well as the challenges they addresse.

![Figure1](Datasets/images/ReadMe_Figure1_Overview.png)
*Figure 1: Schematic representation of the bioinformatics challenges presented in this study that could prevent detection of, e.g., viruses, viral strains, viral isolates, SNPs. Each challenge is addressed by at least one dataset. The datasets are either real (blue), semi-artificial (orange) or completely artificial (grey).*

### Publication

<img align="left" width="128" height="90" src="Datasets/images/PCI_Genomics.png">
<br>

This work has been the subject of a publication recommended by [PCI Genomics](https://genomics.peercommunityin.org/), a free recommendation process of scientific preprints based on peer-reviews. The publication can be found [here](https://zenodo.org/record/4584718#.YHRROj869PY) and the recommendation [here](https://genomics.peercommunityin.org/articles/rec?id=14).

<br>


### Other interesting datasets

In 2019, [Massart *et al.*](https://apsjournals.apsnet.org/doi/10.1094/PHYTO-02-18-0067-R) have compared the ability of 21 plant virology laboratories to detect 12 plant viruses through a double-blind large-scale performance test. They have used 3 real datasets of 21- to 24-nucleotide small RNA (sRNA) sequences from three different infected plants (grapevine, potato and apple). From the fastq files of the 3 initial sRNA datasets, the total reads have been rarefied to obtained 3 different sequencing depth: 50K, 250K and 2.5M. Additionally, a second random subsampling has been carried out at the 250K-read depth for the grapevine dataset to create two technical pseudoreplicates. The resulting 10 datasets are publicly available [here](https://github.com/plantvirology/COST_Action_PT/releases).

### Contact

For any question or remark, you can email <lucie.tamisier@inrae.fr>, <sebastien.massart@uliege.be> or <annelies.haegeman@ilvo.vlaanderen.be>.

### Acknowledgements

We wish to thank all researchers that contributed their datasets as starting point (see Table 1).

This work is distributed under a [CC BY 4.0 license](http://creativecommons.org/licenses/by/4.0/).

This research is partially funded by the Belgian Federal Public Service of Health, Food Chain Safety and Environment ([FPS Health](https://www.health.belgium.be/en)) through the contract "RI 18_A-289" and by the [Euphresco](http://www.euphresco.net) project "Plant Health Bioinformatics Network" (PHBN) (2018-A-289). The work was coordinated by the University of Liège (ULG), Belgium ([Lucie Tamisier](mailto:lucie.tamisier@inrae.fr) & [Sébastien Massart](mailto:sebastien.massart@uliege.be)) with Flanders Research Institute for Agriculture and Fisheries (ILVO), Belgium as main partner ([Annelies Haegeman](mailto:annelies.haegeman@ilvo.vlaanderen.be)).
