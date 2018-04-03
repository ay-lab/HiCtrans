Requirement:
R "changepoint" package. This can be installed by typing install.packages("changepoint") from R interactive enviornment.
_____________________________________________________________________

To run, please go inside the "Translocation_calling" folder and type :
perl run_HiCtrans.pl T47D.chr7_chr15.input.txt

This will create a bash file with required commands to call translocations from Hi-C data. 
Currently HiCtrans can call translocations from 40Kb resolution file but can call from 5Kb resolution also. Only requirement is to have the 5Kb Genomic_feature file (present in the data folder).
We use this file to filter out low GC content, mappability and blacklisted regions from the inter-chromosomal Hi-C data.
Please change the  "../data/hg19.fa.40000.genome_feature.txt" line in "scripts/transNormScript.pl" file in order to use for other resolutions.
_____________________________________________________________________

To generate the genome feature file for other resolutions, genomes or restriction fragments use the our "./data/utility/create_F_GC_MAP_file.pl" script. For more details check the "readme.txt" under "data"
_____________________________________________________________________

To run individual chromosomes separately, create a *.input.txt file like the following. E.g. The T47D.chr7_chr15.input.txt file has all 
the information to call a translocation from a Hi-C matrix file. The file descriptions is as follows,

T47D.chr7_15.bed #This is the bed file with index informatio. Users can provide the full genome bed file.

T47D.chr7_15.matrix #Hi-C matrix file in sparse matrix format. Unser can provide the full genome Hi-C matrix also. The format of the 
matrix file should be same as that of the example file.
chr7 #Chromosome A
chr15 #Chromosome B
chr7_chr15_Folder #All the results will be placed under this folder once the calling part is finished

All chromosome scanning:

To scan all the chromosomes please create a *.input.txt file inside "Translocation_calling" folder.
In the input.txt write the following lines

Genome.bed
Genome.matrix
all
all

Genome.bed is the full genome bed file with index information similar to T47D.chr7_15.bed file (default 40Kb resolution).
Genome.matrix is the full genome raw Hi-C data in sparse matrix format similar to T47D.chr7_15.matrix file (default 40Kb resolution).

Copy your bed and matrix file inside "Translocation_calling" folder.

Run HiCtrans with the following command

perl run_HiCtrans.pl input.txt (This will create the bash file to run the job)

Note: HiCtrans requires a restriction enzyme specific genome feature file with GC content, mappability information for translocation detection. By default hg19 HindIII restriction enzyme specific file @40kb resolution is provided.
_____________________________________________________________________
