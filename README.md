This README describes the scripts used for the analyses in:   
[Diversity of functionally permissive sequences in the receptor-binding site of influenza hemagglutinin](http://www.cell.com/cell-host-microbe/abstract/S1931-3128(17)30204-4)

## ANALYSIS FOR HA RBS DEEP MUTATIONAL SCANNING EXPERIMENTS
This study aims to examine the functional sequence diversity and epistasis of influenza A virus hemagglutinin (HA) receptor-binding site (RBS). Deep mutational scanning experiment was performed for the HA RBS of two strains, namely A/WSN/33 (WSN; H1N1) and A/Hong Kong/1/1968 (HK68; H3N2). The experiment probed for the fitness effect of mutants that contain up to three amino-acid substitutions.

### INPUT FILE
* All sequencing raw reads, which can be downloaded from NIH SRA database [PRJNA353496](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA353496), should be placed in fastq/folder:
  * WSN single mutant input library: fastq/WSN\_HARBS-1\_S12\_L001\_R1\_001.fastq and fastq/WSN\_HARBS-1\_S12\_L001\_R2\_001.fastq
  * WSN single mutant passaged library: fastq/WSN\_HARBS-4\_S15\_L001\_R1\_001.fastq and fastq/WSN\_HARBS-4\_S15\_L001\_R2\_001.fastq
  * WSN double mutant input library: fastq/WSN\_HARBS-5\_S16\_L001\_R1\_001.fastq and fastq/WSN\_HARBS-5\_S16\_L001\_R2\_001.fastq
  * WSN double mutant passaged library: fastq/WSN\_HARBS-2\_S13\_L001\_R1\_001.fastq and fastq/WSN\_HARBS-2\_S13\_L001\_R2\_001.fastq
  * WSN triple mutant input library: fastq/WSN\_HARBS-3\_S14\_L001\_R1\_001.fastq and fastq/WSN\_HARBS-3\_S14\_L001\_R2\_001.fastq
  * WSN triple mutant passaged library: fastq/WSN\_HARBS-6\_S17\_L001\_R1\_001.fastq and fastq/WSN\_HARBS-6\_S17\_L001\_R2\_001.fastq
  * HK68 triple mutant input library: fastq/HK68-Tlib-1\_S1\_L001\_R1\_001.fastq and fastq/HK68-Tlib-1\_S1\_L001\_R2\_001.fastq
  * HK68 triple mutant passaged library (round 1): fastq/HK68-Tlib-2\_S2\_L001\_R1\_001.fastq and fastq/HK68-Tlib-2\_S2\_L001\_R2\_001.fastq
  * HK68 triple mutant passaged library (round 2): fastq/HK68-Tlib-3\_S3\_L001\_R1\_001.fastq and fastq/HK68-Tlib-3\_S3\_L001\_R2\_001.fastq

### ANALYSIS PEPLINE
#### FOR WSN
1. script/WSN\_HARBS\_read2RFindex.py: Converting raw reads to counts and RF index
  * Input file: fastq/WSN\*.fastq
  * Output file: result/WSN\*.count
2. script/WSN\_HARBS\_compileMutFit.py: Compile \*.count files into a single file
  * Input file: result/WSN\*.count
  * Output file: result/WSN\_MutFitTable.tsv
3. script/WSN\_HARBS\_CalMaxFit.py: Calculate the maximum RF index among all mutants that carried a specified substitution
  * Input file: result/WSN\_MutFitTable.tsv
  * Output file: result/WSN\_MaxFitMut\_\*.tsv
4. script/WSN\_HARBS\_CrypticBen.py: Search for the maximum beneficial effect for a given substitution of interest (max. RF increase) among all genetic backgrounds
  * Input file: result/WSN\_MutFitTable.tsv
  * Output file: result/WSN\_crypticben.tsv

#### FOR HK68
1. script/HK68\_HARBS\_read2RFindex.py: Converting raw reads to counts and RF index
  * Input file: fastq/HK68\*.fastq
  * Output file: result/HK68\_Tlib.count
2. script/HK68\_HARBS\_compileMutFit.py: Compile \*.count files into a single file and also generate a combined file with WSN
  * Input file:
    * result/HK68\_Tlib.count
    * result/WSN\_MutFitTable.tsv
  * Output file:
    * result/HK68\_MutFitTable.tsv
    * result/MutCompareTable.tsv
3. script/HK68\_HARBS\_CalMaxFit.py: Calculate the maximum RF index among all mutants that carried a specified substitution
  * Input file: result/HK68\_MutFitTable.tsv
  * Output file: result/HK68\_MaxFitMut\_\*.tsv
4. script/HK68\_HARBS\_CrypticBen.py: Search for the maximum beneficial effect for a given substitution of interest (max. RF increase) among all genetic backgrounds
  * Input file: result/HK68\_MutFitTable.tsv
  * Output file: result/HK68\_crypticben.tsv

#### COMPARISON ANALYSIS
1. script/EpisCount.py: Count reciprocal sign epistasis
  * Input file: 
    * result/WSN\_MutFitTable.tsv
    * result/HK68\_MutFitTable.tsv
  * Output file: result/EpisCountAroundWT.tsv
  
### PLOTTING
#### FOR WSN
* script/WSN\_HARBS\_FitQC.R: Plot the RF index distribution of silent, missense, and nonsense mutations
  * Input file: result/WSN\*.count
  * Output file: graph/WSN\_FitQC\_stripchart.png
* script/WSN\_HARBS\_PlotMaxFit.R: Plot the maximum RF index among all mutants that carried a specified substitution as heatmap and also plot their distribution
  * Input file: result/WSN\_MaxFitMut\_\*.tsv
  * Output file: graph/WSN\_MaxFitMut\_\*.png
* script/WSN\_HARBS\_CrypticBen.R: Plot the maximum beneficial effect for a given substitution of interest (max. RF increase) among all genetic backgrounds
  * Input file: result/WSN\_crypticben.tsv
  * Output file: WSN\_CrypticBen\*.png

#### FOR HK68
* script/HK68\_HARBS\_FitQC.R: Plot the RF index distribution of silent, missense, and nonsense mutations
  * Input file:
    * result/HK68\_Tlib.count
    * result/HK68\_MutFitTable.tsv
  * Output file: graph/HK68\_FitQC\_stripchart.png
* script/HK68\_HARBS\_PlotMaxFit.R: lot the maximum RF index among all mutants that carried a specified substitution as heatmap and also plot their distribution
  * Input file: result/HK68\_MaxFitMut\_\*.tsv
  * Output file: graph/HK68\_MaxFitMut\_\*.png
* script/HK68\_HARBS\_CrypticBen.R: Plot the maximum beneficial effect for a given substitution of interest (max. RF increase) among all genetic backgrounds
  * Input file: result/HK68\_crypticben.tsv
  * Output file: HK68\_CrypticBen\*.png

#### COMPARISON ANALYSIS
* script/EpisPlot.R: Compare reciprocal sign epistasis by barplot and perform Fisher's exact test
  * Input file: result/EpisCountAroundWT.tsv
  * Output file: 
    * graph/WSN\_EpisCount.png
    * graph/HK68\_EpisCount.png

* script/HARBS\_compare.R: Compare RF index of a given variant in WSN and HK68 backgrounds
  * Input file: result/MutCompareTable.tsv
  * Output file: graph/MutFitCompare.png

* script/EpiMap.R: Plot the pairwise position of reciprocal sign epistasis
  * Input file: 
    * result/WSN\_EpisMap.tsv
    * result/HK68\_EpisMap.tsv
  * Output file: 
    * graph/WSN\_EpiMap.png
    * graph/HK68\_EpiMap.png
