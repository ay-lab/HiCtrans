chromosome	start	end	index	chromosome_id	start	end	effective_fragment_length	GC_content	mappability	blacklisted_or_not (I:include or not blacklisted)
chr1    0       40000   1       1       0       40000   1000    0.488   0.190   I
chr1    40000   80000   2       1       40000   80000   1000    0.350   0.305   I
chr1    80000   120000  3       1       80000   120000  1000    0.379   0.124   I
chr1    120000  160000  4       1       120000  160000  1000    0.416   0.057   I
chr1    160000  200000  5       1       160000  200000  1000    0.385   0.068   I
chr1    200000  240000  6       1       200000  240000  1000    0.433   0.101   I
chr1    240000  280000  7       1       240000  280000  1000    0.358   0.095   I
chr1    280000  320000  8       1       280000  320000  1000    0.394   0.051   I
chr1    320000  360000  9       1       320000  360000  1000    0.407   0.066   I

To create this genomic feature file use the scripts are under "utility" folder. Open the "create_F_GC_MAP_file.pl" file and change the required variables inside the program as per your requirement. Type "perl create_F_GC_MAP_file.pl" to generate the genomic feature file.

__________________________________________________________________________
Parameters to change in the "create_F_GC_MAP_file.pl" files are:

$frag = "HindIII_resfrag_hg19.bed"; #Change the restriction fragment file as per your experiment. Some restriction fragment files for hg/mm genomes are provided in the "HindIII_resfrag_files.zip" file 

$read_length = 50; #Hi-C fastq read length

$fasta_file = "hg19.fa"; #Needs hg19.fa or equivalent genome sequence file

$map_bedgraph = "hg19.MapabilityAlign50mer.bedGraph"; #Need to download the mappability file separately 

$resolution = 40000; #Change the HiC resolution as per your requirement
 
$chrom_size_file = "chrom_hg19.sizes"; #Chromosome size file

$bedtools = "/mnt/BioApps/bedtools/bin/bedtools"; #bedtools path

$blacklisted_region = "EncodeExcludableRegions.hg19.bed"; # Keep this variable as "No" if you don't have the black listed region information

__________________________________________________________________________

This file can also be created by following the HiCNorm protocol (http://www.people.fas.harvard.edu/~junliu/HiCNorm/estimate_mapp_sub.rar).

Effective fragment length is not nesseccary. So one can change it a constant value instead.
