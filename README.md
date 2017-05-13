## ANALYSIS FOR HA RBS DEEP MUTATIONAL SCANNING EXPERIMENTS
This study aims to examine the functional sequence diversity and epistasis of influenza A virus hemagglutinin (HA) receptor-binding site (RBS). Deep mutational scanning experiment was performed for the HA RBS of two strains, namely A/WSN/33 (H1N1) and A/Hong Kong/1/1968 (H3N2). The experiment probed for the fitness effect of mutants that contain up to three amino-acid substitutions.

### FILE
* All sequencing raw reads, which can be downloaded from NIH SRA database [PRJNA353496](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA353496), should be placed in fastq/folder:
  * WSN single mutant input library: fastq/WSN\_HARBS-1\_S12\_L001\_R1\_001.fastq
  * WSN single mutant passaged library: fastq/WSN\_HARBS-4\_S15\_L001\_R1\_001.fastq
  * WSN double mutant input library: fastq/WSN\_HARBS-5\_S16\_L001\_R1\_001.fastq
  * WSN double mutant passaged library: fastq/WSN\_HARBS-2\_S13\_L001\_R1\_001.fastq
  * WSN triple mutant input library: fastq/WSN\_HARBS-3\_S14\_L001\_R1\_001.fastq
  * WSN triple mutant passaged library: fastq/WSN\_HARBS-6\_S17\_L001\_R1\_001.fastq
* result/\*.count: 

### ANALYSIS PEPLINE
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

### PLOTTING
* script/WSN\_HARBS\_FitQC.R: Plot the RF index distribution of silent, missense, and nonsense mutations
  * Input file: result/WSN\*.count
  * Output file: graph/WSN\_FitQC\_stripchart.png
* script/WSN\_HARBS\_PlotMaxFit.R: Plot the maximum RF index among all mutants that carried a specified substitution as heatmap and also plot their distribution
  * Input file: result/WSN\_MaxFitMut\_\*.tsv
  * Output file: graph/WSN\_MaxFitMut\_\*.png
* script/WSN\_HARBS\_CrypticBen.R: Plot the maximum beneficial effect for a given substitution of interest (max. RF increase) among all genetic backgrounds
  * Input file: result/WSN\_crypticben.tsv
  * Output file: WSN\_CrypticBen\*.png
