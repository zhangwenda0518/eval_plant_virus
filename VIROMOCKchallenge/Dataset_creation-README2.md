### General pipeline used to create semi-artificial datasets

#### 1. Analysis of the real datasets with Geneious

Before the creation of the artificial reads, each real dataset has been analysed in order to check that we were able to detect every viruses/viroids and to assess their frequencies. Note that for each dataset, the presence of the viruses/viroids has also been confirmed by PCR and/or ELISA. The datasets have been analysed in [Geneious](https://www.geneious.com/). The main steps performed were:
- Pairing of the reads for paired-end datasets. For each paired-end dataset, all the reads have always been paired (i.e. no unpaired reads have been used).
- Merging / assembling of the paired reads to obtain longer reads for paired-end datasets.
- Trimming of the merged and unmerged reads ([BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/), minimum quality = Q25, minimum length = 20 bp).
- Removing duplicates of the trimmed, merged and unmerged reads. We obtain the final “filtered reads”.
- De novo assembly of the filtered reads using [SPAdes](http://cab.spbu.ru/software/spades/) or [Velvet](https://github.com/dzerbino/velvet/tree/master) for the short reads (50 or 75 bp).
- Blasting of the contigs against the [RefSeq virus database](https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/) (using [tblastx] from the blast+ suite (https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)).
- Mapping of the filtered reads against either the blast references or the contigs (if we have a complete genome, which depends on the dataset). If a read could be mapped with equal match to two or more strains, it was randomly assigned to one of them.
- Assessment of the proportion of each virus/viroid of the dataset.

In the end, 8 real datasets have been analysed. Three were already interesting and showed limits that we wanted to test. Therefore, these 3 datasets (No. 7, 8 and 9) have not been spiked with artificial viral reads or modified. The other 5 real datasets have been used to created 7 semi-artificial datasets.


#### 2. Creation of the semi-artificial datasets with ART and Linux tools
[ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) is a set of simulation tools allowing the generation of synthetic next-generation sequencing reads. It mimics the real sequencing process and allows the creation of FASTQ files with single or paired-end reads of different sizes.

The **first challenge** was to create artificial reads showing the same quality score as the reads from the real datasets. ART allows the creation of our own quality profile files. Therefore, for each dataset, quality profiles have been created using the real reads. In the case of paired-end datasets, quality profiles have been generated both for the reads 1 and reads 2. Then, artificial reads showing the same size and quality score as the one of the real datasets have been created. The number of reads generated depends on the dataset (see following section for details). Figure 1 shows the results of a FastQC analysis on both the real and artificial reads generated for  dataset 1. The quality scores are extremely similar between both type of reads, and we have obtained the same kind of results for the other semi-artificial datasets.

![Figure1](images/Figure1_Quality_scores.png)
*Figure 1: Quality scores of real and artificial reads of the paired-end dataset 1.*

The **second challenge** was to merge artificial and real reads in the same file and with similar headers. Indeed, ART creates its own header for each read and does not allow to change it. The different steps to overcome this issue are summarized in Figure 2. These steps have been performed using command-lines in Linux and are the following:
- Both artificial and real reads have been merged in the same file.
- The old headers have been removed and new headers have been added. These new headers are all the same, they look like a “real” header, but they have a number in the end. If the dataset contains x reads, x different numbers have been generated and put randomly at the end of the headers. Here is an example of headers for a hypothetical paired-end dataset with 3 reads, the number added being the last number before the space:

```
@E00526:48:HT7CTCXY:3:1101:16295:993:1 1:N:0:NGAGCTAC+NTAGCCTT   Read 1
@E00526:48:HT7CTCXY:3:1101:16295:993:1 2:N:0:NGAGCTAC+NTAGCCTT   Read 2
@E00526:48:HT7CTCXY:3:1101:16295:993:3 1:N:0:NGAGCTAC+NTAGCCTT   Read 1
@E00526:48:HT7CTCXY:3:1101:16295:993:3 2:N:0:NGAGCTAC+NTAGCCTT   Read 2
@E00526:48:HT7CTCXY:3:1101:16295:993:2 1:N:0:NGAGCTAC+NTAGCCTT   Read 1
@E00526:48:HT7CTCXY:3:1101:16295:993:2 2:N:0:NGAGCTAC+NTAGCCTT   Read 2
```

- The reads have been sorted according to number of the header. It allows the mixing of the real and artificial reads.

![Figure2](images/Figure2_Dataset_creation.png)
*Figure 2: Final steps of the semi-artificial datasets creation.*
