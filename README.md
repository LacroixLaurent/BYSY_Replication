# BYSY_Replication
### Laurent Lacroix (*laurent.lacroix@inserm.fr*)
***
## About this repository  
To see this text in browser, go there [https://lacroixlaurent.github.io/BYSY_Replication/](https://lacroixlaurent.github.io/BYSY_Replication/), otherwise this repository can be viewed in RStudio.  
This repository contains the scripts used for all the analysis in the one-chromosome yeast replication manuscript [Pellet et al., 2025](https://doi.org/10.1038/XXXXXX).  

Raw data are available from ENA repository under accession number [PRJEB89792](https://www.ebi.ac.uk/ena/).  

This repository documents the bioinformatics workflows used to compare the replication program of a yeast strain with all 16 chromosomes fused (strain SY14, [Shao et al](https://doi.org/10.1038/s41586-018-0382-x)) to the standard 16-chromosome reference strain BY4742.  
Genome sequences for SY14 and BY4742 have been downloaded from public database using the accession number from the original SY14 strain publication [PRJNA429985](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA429985).  

[NFS](https://github.com/LacroixLaurent/NanoForkSpeed) and [nanoTiming](https://github.com/LacroixLaurent/NanoTiming) experiments have been conducted as in [Theulot et al, 2022](https://doi.org/10.1038/s41467-022-31012-0) and [Theulot et al, 2024](https://doi.org/10.1038/s41467-024-55520-3) respectively (see BYSY publication for more details). NFS and nanotiming data have been processed with matching reference genome. Please note that these reference genomes do not contain *pBL-hsvTKCO-hENT1CO* sequence, the plasmid used to enable BrdU usage.  

The procedure are explained in separate notebooks corresponding to the following processes:  

1. [Data transposition](./01_BYSY_Data_Transposition.nb.html)

2. [Genomes annotations](./02_BYSY_Genome_Annotation.nb.html)

3. [Generating RFD profiles](./03_BYSY_RFD_Profiles.nb.html)

4. [Defining active replication origins in SY and/or BY strain](./04_BYSY_Active_Ori.nb.html)

5. [Importing bibliographic data](./05_BYSY_Biblio.nb.html)

6. [Computing OEM](./06_BYSY_OEM.nb.html)

7. [Testing RFD data](./07_BYSY_RFD_Test.nb.html)

8. [Testing Initiation data](./08_BYSY_Init_Test.nb.html)

9. [Testing nanoTiming data](./09_BYSY_nanoT_Test.nb.html)

10. [Generating master table](./10_BYSY_Master_Table.nb.html)

11. [Generating and testing speed map](./11_BYSY_Speed_Map.nb.html)

The scripts starting with *BYSY_FiguresArticle* are the one used to generate the figures.  
The file *Helper_function.R* contains accessory R functions used during the analysis.

* the *Data_raw* folder contains the raw data used to generate the figures and analysis.  
* the *Data* folder contains the data produced by the analysis scripts.  
* the *BigWig* folder contain the genomic track generated during the analysis.  
* the *Genome_annotations* folder contains the genomic annotation  tracks generated during the analysis.  
* the *FigArticle* folder contains the figures generated during the analysis. The figures with *mod* in their name have been manually modified with Inkscape to adjust ARS name positions for readability.  
* the *BigFiles* folder contains *BigWig*, *Data_raw* and *Data* files too big for github and that have to be downloaded from the Zenodo repository [https://doi.org/10.5281/zenodo.15462552](https://doi.org/10.5281/zenodo.15462552).  

***
Code to generate the session info:  
library(sessioninfo)  
session_info() %>% capture.output(file="BYSY_Replication_processing_session_info.txt")  
