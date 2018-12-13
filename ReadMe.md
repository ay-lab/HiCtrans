Requirement:
	
	Within R :
	1. changepoint
	2. optparse 
	3. data.table
	4. hashmap
	5. caTools
	6. depmixS4
	7. DEoptimR

	Within Perl :
  	1. Getopt::Long
	2. Parallel::ForkManager

Quick start:

	To run, please go inside the "Translocation_calling" folder and type :

	$ perl HiCtransScript.pl -bed T47D.chr7_15.bed -mat T47D.chr7_15.matrix -chromA chr7 -chrB chr15 -prefix T47D

	This will call and run the HiCtrans program which will generate the tranlocation borders. 
	HiCtrans can call translocations from any resolution file (The example provided is of 40Kb resolution). 
	But for each resolution there needs to a Genomic_feature file (present in the data folder).
	HiCtrans uses this file to filter out low GC content, mappability and blacklisted regions from the inter-chromosomal Hi-C data.
	Please change the  "../data/hg19.fa.40000.genome_feature.txt" line in "scripts/transNormScript.pl" file in order to use for other resolutions.

	Note:
	To generate the genome feature file for other resolutions, genomes or restriction fragments use the our "./data/utility/create_F_GC_MAP_file.pl" script. For more details check the "readme.txt" under "data"


	Starting from this current release, HiCtrans can call translocation breakpoints at the restriction site resolution.
	For that the users needs to provide two separate files:
		
		1. A whitespace/tab separated file that contains, on each line
		
		<readname> <chr1> <pos1> <strand1> <chr2> <pos2> <strand2>
		
		e.g. 
   		SRR6213722.1   chr11   124331538       -       chr11   124345246       -
        	SRR6213722.2   chr1    198436365       -       chr1    199923196       +

		2. Restriction fragments bed file. e.g. 
		
		chr1    0       16007   HIC_chr1_1      0       +
	        chr1    16007   24571   HIC_chr1_2      0       +
        	chr1    24571   27981   HIC_chr1_3      0       +

		

	For details please do, 
	
	$ perl HiCtransScript.pl -h 

	Usage: perl HiCtransScript.pl -bed T47D.chr7_15.bed -mat T47D.chr7_15.matrix -chromA chr7 -chrB chr15 -prefix T47D
               perl HiCtransScript.pl -bed Genome.bed -mat Genome.matrix -chromA chr7 -chrB chr15 -prefix T47D -vp T47D_validpairs.txt -refrags HindIII_REfrags.bed

        -bed: Genome wide binned bed file with index information

        chr1    0       40000   1
        chr1    40000   80000   2
        chr1    80000   120000  3
        chr1    120000  160000  4
        .....

        -mat: Hi-C interaction file

        1       2       234
        1       3       100
        2       3       150
        3       4       110
        .....

        -chromA: chromosome A to be searched in the Hi-C data (either chromosome name or all)

        -chromB: chromosome B to be searched in the Hi-C data (either chromosome name or all)

        NOTE: For all possible pairs provide '-chromA all -chromB all' as argument

        -prefix: Output folder or file name

        -size: chromsome names and their sizes (NOTE: Ordering of the chromosome should be same as that of the binned bed file)

        chr1    249250621
        chr2    243199373
        chr3    198022430
        .....

        -cf: Filter out intergenic-interactions below this contact count (Can be changed based on sample, default is 10)

        -rf: This is the box.entropy by random.entropy ratio (This is will control the noise effect, default is 2)

        -vp: Hi-C validpair fragment file

        SRR6213722.1   chr11   124331538       -       chr11   124345246       -
        SRR6213722.2   chr1    198436365       -       chr1    199923196       +
        .....

        -refrags: Restriction fragment file

        chr1    0       16007   HIC_chr1_1      0       +
        chr1    16007   24571   HIC_chr1_2      0       +
        chr1    24571   27981   HIC_chr1_3      0       +

        NOTE: Breakpoint finding is only possible if validpair and restriction fragment files are not provided. Otherwise only translocation boundaries are provided.

        -cpu: Number of cpu to use (default 1)

        -h: Help


			
	Note:
	T47D.chr7_15.bed is the bed file with index information. Users can provide the full genome bed file.
	T47D.chr7_15.matrix is the Hi-C matrix file in sparse matrix format. Unser can provide the full genome Hi-C matrix also.


Output:

	After completion a *_Folder will be generated which will have all the raw files and a *.BreakPoints.txt file will be generated. 
	This file will contain the translocation breakpoint coordinates. Following is an example output:

	chrA    BoundaryAS      BoundaryAE      chrB    BoundaryBS      BoundaryBE      MaxCount        box.entropy     random.entropy.99uCI    ratio   chromA_BP_start chromA_BP_end   chromB_BP_start chromB_BP_end
	chr7    83060000        87340000        chr15   29980000        31500000        21      0.2095       0.0320      6.5384        87334408        87337551        29962616        29970578

	chrA  & chrB: Translocated chromosome pairs.

	BoundaryAS, BoundaryAE, BoundaryBS and BoundaryBE: chrA/chrB translocation boundary start/end.

	MaxCount: Max xontact count within the translocation boundary.

	box.entropy: Normalized entropy of counts within BoundaryAS, BoundaryAE, BoundaryBS and BoundaryBE.

	random.entropy.99uCI: Normalized entropy + 99% confidence interval of counts of random boxes (Similar area defined by BoundaryAS, BoundaryAE, BoundaryBS and BoundaryBE but excluding all the breakpoint boundaries).

	ratio: box.entropy/random.entropy.99uCI (Translocated region will be enriched in heterogeneous mixture of different count values [high entropy] compared to a random region with homogeneous count values [Low entropy]). 
	
	chromA_BP_start, chromA_BP_end: chrA breakpoint at restriction site resolution.
  
 	chromB_BP_start, chromB_BP_end: chrB breakpoint at restriction site resolution.



Juicebox output to HiCtrans matrix format conversion:

	If you have an output file from Juicebox (e.g. by following way java -jar juicebox_tools.7.5.jar dump observed NONE K562.combined.hic 1 10 BP 50000 chr1_chr10.txt) then you can convert the chr1_chr10.txt file to HiCtrans input files by using "juiceboxToHiCtrans_MatrixFormat.r" program present under scripts folder. For example to convert the chr1_chr10.txt file use the script in the following way,

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
        	        chromosome wise genome size file. eg. a file like hg19.txt with following information
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

Contact

	abhijit@lji.org (Abhijit Chakraborty)
