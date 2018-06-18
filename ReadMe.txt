Requirement:
R "changepoint and optparse" package. This can be installed by typing install.packages(c("changepoint","optparse")) from R interactive enviornment.
_____________________________________________________________________
Quick start:

To run, please go inside the "Translocation_calling" folder and type :
perl run_HiCtrans.pl T47D.chr7_chr15.input.txt

This will create a bash file with required commands to call translocations from Hi-C data. 
Currently HiCtrans can call translocations from 40Kb resolution file but can call from 5Kb resolution also. Only requirement is to have the 5Kb Genomic_feature file (present in the data folder).
We use this file to filter out low GC content, mappability and blacklisted regions from the inter-chromosomal Hi-C data.
Please change the  "../data/hg19.fa.40000.genome_feature.txt" line in "scripts/transNormScript.pl" file in order to use for other resolutions.
_____________________________________________________________________
Note:
To generate the genome feature file for other resolutions, genomes or restriction fragments use the our "./data/utility/create_F_GC_MAP_file.pl" script. For more details check the "readme.txt" under "data"
_____________________________________________________________________
Specific chromosome wise translocation detection:

To run individual chromosomes separately, create a *.input.txt file like the following. E.g. The T47D.chr7_chr15.input.txt file has all 
the information to call a translocation from a Hi-C matrix file. The file descriptions is as follows,

T47D.chr7_15.bed  #This is the bed file with index information. Users can provide the full genome bed file.

T47D.chr7_15.matrix #Hi-C matrix file in sparse matrix format. Unser can provide the full genome Hi-C matrix also. The format of the 
matrix file should be same as that of the example file.
chr7              #Chromosome A
chr15             #Chromosome B
chr7_chr15_Folder #All the results will be placed under this folder once the calling part is finished

_____________________________________________________________________
Genome wide translocation detection:

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
Output:

After completion "All.chromosome.Translocation.EntropyFiltered.Final.result" file will be generated. This file will contain the translocation breakpoint coordinates. Following is an example output:

chrA  BoundaryAS  BreakPointA BoundaryAE  chrB  BoundaryBS  BreakPointB BoundaryBE  count box.entropy random.entropy.99uCI  ratio
chr7  81900000  87260000  87340000  chr15 29980000  30100000  34140000  23  0.146 0.028 5.159

chrA  & chrB: Translocated chromosome pairs.

BoundaryAS, BoundaryAE, BoundaryBS and BoundaryBE: chrA breakpoint boundary start, chrA breakpoint boundary end, chrB breakpoint boundary start, and chrB breakpoint boundary end.

BreakPointA, BreakPointB: chrA and chrB breakpoint coordinate.

count: Contact count at BreakPointA and BreakPointB coordinate.

box.entropy: Normalized entropy of counts within BoundaryAS, BoundaryAE, BoundaryBS and BoundaryBE.

random.entropy.99uCI: Normalized entropy + 99% confidence interval of counts of random boxes (Similar area defined by BoundaryAS, BoundaryAE, BoundaryBS and BoundaryBE but excluding all the breakpoint boundaries).

ratio: box.entropy/random.entropy.99uCI (Translocated region will be enriched in heterogeneous mixture of different count values [high entropy] compared to a random region with homogeneous count values [Low entropy])  

_____________________________________________________________________
Juicebox output to HiCtrans matrix format conversion

If you have an output file from Juicebox (e.g. by following way java -jar juicebox_tools.7.5.jar dump observed NONE K562.combined.hic 1 10 BP 50000 chr1_chr10.txt) then you can convert the chr1_chr10.txt file to HiCtrans input files by using "juiceboxToHiCtrans_MatrixFormat.r" program present under scripts folder.

For example to convert the chr1_chr10.txt file use the script in the following way

Rscript juiceboxToHiCtrans_MatrixFormat.r --bedtools bedtools --genomesize hg19.sizes --resolution 50000 --column1 chr1 --column2 chr10 --juicebox_file chr1_chr10.txt

This script will create two files, input to the HiCtrans. 

chr1_chr10.hictrans.index.bed and chr1_chr10.hictrans.matrix file. 

for help just do 

Rscript juiceboxToHiCtrans_MatrixFormat.r -h

Usage: juiceboxToHiCtrans_MatrixFormat.r [options]


Options:
        --bedtools=BEDTOOLS
                Path to bedtools program

        --genomesize=GENOMESIZE
                chromosome wise genome size file. eg. a file like hg19.txt with followign information
        chr1    249250621
        chr2    243199373
        chr3    198022430
        ....


        --resolution=RESOLUTION
                Hi-C resolution

        --column1=COLUMN1
                Column1 chromosome name of juicebox output

        --column2=COLUMN2
                Column2 chromosome name of juicebox output

        --juicebox_file=JUICEBOX_FILE
                Juicebox output file name

        -h, --help
                Show this help message and exit
_____________________________________________________________________
Contact

abhijit@lji.org (Abhijit Chakraborty)
