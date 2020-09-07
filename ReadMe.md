# HiCtrans 

This is an updated version of HiCtrans program. HiCtrans can scan inter-chromosomal Hi-C matrix and report translocations, their breakpoints at restriction site or at any lower resolution.
Check the paper https://doi.org/10.1093/bioinformatics/btx664 for details about the method.

Changes made from previous version:

1. No requirement of genome feature file. 
2. Given a starting Hi-C data, HiCtrans now can scan for translocations at different resolutions and report a translocation observed at multiple resolutions.
3. Users can also check potential translocations at different resolutions.
4. No requirement of perl and its associted libraries.
5. Error handling.
6. Faster.

Result description:

A successfull HiCtrans run will generate the following result files and folders

```bash

<prefix>_hictrans
       .
	\<prefix>_hictrans_<chrA>_<chrB>_<resolution>
               .
		\Lower_Resolution_HiC_Data
			<prefix>_<resolution>.<chrA>_<chrB>.matrix	
			<prefix>_<resolution>.<chrA>_<chrB>_abs.bed
			....
  	       .
		\Translocations
                       .
			\Details	
				<prefix>_<resolution>.<chrA>_<chrB>.Details.txt
		       .	....
			\juicebox_files	
				<prefix>_<resolution>.<chrA>_<chrB>.Translocations_jcbx.txt
				....				
			<prefix>_<resolution>.<chrA>_<chrB>.preCluster.txt
	       .
		\MultiResolution_supported_Translocations
			<prefix>_<resolution>.<chrA>_<chrB>.MultiResolution_Filtered.Translocation.txt

        <prefix>.matrix	
	<prefix>_abs.bed
	<prefix>.<chrA>_<chrB>.mat.txt 
	<prefix>.<chrA>_<chrB>.log.txt
```

NOTE: MultiResolution_supported_Translocations folder is only created when there are such cases.

Once finished check the <prefix>_<resolution>.<chrA>_<chrB>.MultiResolution_Filtered.Translocation.txt file for possible translocations.
\<prefix>\_\<resolution>.\<chrA>\_\<chrB>.MultiResolution_Filtered.Translocation.txt provides strong support for the translocation with any anomaly in the inter Hi-C data.
If there is no multi-resolution supported translocations, users can check the \<prefix>\_\<resolution>.\<chrA>\_\<chrB>.preCluster.txt file 'Translocations' folder.
This file will have all the translocations (BreakPoints and Translocation boxes) found in the chromosomal pair data at different resolutions. 
Users can check the highest resolution in the 'resolution column' (lower the value higher the resolution) for further investigation.
The zscore column repersents the enrichment of counts within the box associated with the translocation. The 'count' column is simply the hic count of the
breakpoint detected within the enriched box. Users can ignore the 'id' column.

For detailed help use the following

```bash

Rscript hictrans.v3.R --help

Usage: hictrans.v3.R [options]

Options:
        --mat=MAT
                An upper triangular Hi-C sparse matrix
                It should have the following columns

                <indexA> <indexB> <count>

                1 1 300
                1 2 30
                1 3 10
                2 2 200
                2 3 20
                3 3 200
                ....


        --bed=BED
                Bed file with index information
                It should have the following columns


                <chr> <start> <end> <index>

                chr1 1 40000 1
                chr1 40000 80000 2
                chr1 80000 120000 3
                ....


        --chrA=CHRA
                Chromosome A name. It will represent the rows in the inter-chromosomal matrix. It should be the <indexA> chromosome.


        --chrB=CHRB
                Chromosome B name. It will represent the columns in the inter-chromosomal matrix. It should be the <indexB> chromosome.


        --prefix=PREFIX
                Prefix of the output file <prefix.chrA_chrB>. All the output files and folders will be generated with this prefix.


        --covq=COVQ
                Quantile value to be subtracted from one dimensional trans-coverage profile [trans.coverage - quantile(trans.coverage, covq)] [default 0.10].

                Bins with very low coverage values are removed with this filter.
                Increasing <covq> value will keep only the most stringent bins in the two chromosome.


        --minzscore=MINZSCORE
                Minimum Zscore of a possible translocation box to be retained [default is 1].

                HiCtrans will find enriched boxes within the inter-chromosomal matrix as potential translocation box.
                The enrichment is calculated as Z-score against a background with all possible similar sized boxes in the inter-chromosomal matrix.
                Increasing <minzscore> value will keep the most enriched trans interacting boxes.

	--minboxsize=MINBOXSIZE
                Minimum size of a possible translocation box relative to its Hi-C resolution [default is 0 i.e. no filtering. If set to non-zero value, then (Breakpoint.start - Breakpoint.end)/HiC.resolution > minboxsize filtering will be applied].

                HiCtrans will find enriched boxes within the inter-chromosomal matrix as potential translocation box.
                The minimum box size threshold will filter out small false positive multi-resolution supported potential translocation boxes.

                Increasing <minboxsize> value will keep the most enriched and larger trans interacting boxes.


        --boxzscore=BOXZSCORE
                Minimum Zscore of a possible translocation box to be retained [default is 1].

                HiCtrans will keep boxes enriched above boxzscore threshold to find translocations among them.
                Increasing <boxzscore> value will keep the most enriched trans interacting boxes.


        --locq=LOCQ
                Top percentile to be reported as possible breakpoints within a translocation box [default top 0.1%]

                For each enriched translocated boxes, HiCtrans will report top <locq>% interacting pairs (Weighted by the frequency of total interaction).
                Decreasing <locq> will reduce the number of reported breakpoints within an enriched trans interacting box. A <locq> value of 0 will report only the top interacting pair


        --mincount=MINCOUNT
                Minimum count of a possible breakpoint to be retained when compared to all possible chrA-chrB interaction [default cutoff is 10].

                This is an absolute minimum count cutoff to filter out breakpoints detected at any resolution.
                Increasing <mincount> value will keep only the most stringent interacting pair.


        --glbq=GLBQ
                Percentile value for minimum count cutoff at each resolution [default cutoff is at top 0.1% of the count distribution].

                This is a relative count cutoff based on the inter-chromosomal count distribution determined for each resolution independently.
                Increasing <glbq> value keep only the most stringent interacting pair.


        --resolutions=RESOLUTIONS
                Comma separated list of integers to be multiplied with the starting Hi-C resolution to get the lower resolutions [default 2,4,5,10].

                HiCtrans will search for enriched trans-interacting boxes and breakpoint finding at different resolutions.
                Provide only integer values in a comma separated list.


        --multires=MULTIRES
                Number of Hi-C resolutions at which the breakpoint should be supported with [default is at least 2 different resolutions].

                HiCtrans wiil find enriched trans-interacting boxes and subsequent breakpoints at different resolutions. The ultimate goal is to find a true translocation supported by multiple resolutions.
                Increasing <multires> value will keep only the enriched boxes and breakpoints supported by at least <multires> number of different resolution.


        --maxres=MAXRES
                Maximum resolution upto which the breakpoint is kept after multi-resolution filtering [default is 3 X user provided Hi-C resolution].

                HiCtrans wiil find enriched trans-interacting boxes and subsequent breakpoints at different resolutions. The ultimate goal is to find a true translocation supported by the highest resolutions.
                Increasing <maxres> value will keep only the enriched boxes and breakpoints supported by upto <maxres> X starting Hi-C resolutions.


        --relevel=RELEVEL
                Should the breakpoints be refined upto restriction-level resolution [default is 'No'; If 'Yes', the following parameters are MUST]


        --fragsFile=FRAGSFILE
                Restriction Fragment file [MUST].

                chr1    0       16007   HIC_chr1_1      0       +
                chr1    16007   24571   HIC_chr1_2      0       +
                chr1    24571   27981   HIC_chr1_3      0       +
                ......


        --chromsize=CHROMSIZE
                Chromosome size file [MUST].

                chr1    249250621
                chr2    243199373
                chr3    198022430
                chr4    191154276
                .....


        --validpair=VALIDPAIR
                Valid pair file of the HiC data [MUST].

                SRR6213722.1    chr11   124331538       -       chr11   124345246       -
                SRR6213722.2    chr1    198436365       -       chr1    199923196       +
                .....


        --clusdist=CLUSDIST
                Distance threshold in basepairs to cluster the nearby breakpoints obtained from multi-resolution filtered (MultiResolution_Filtered.Translocation.txt) or individual Translocations_jcbx.txt files [Default 1Mb]

        --ssA=SSA
                Extend -(ve) bp of the 5' HMM segment border of chromosome A for breakpoint identification. Default 100Kb.

        --seA=SEA
                Extend +(ve) bp of the 3' HMM segment border of chromosome A for breakpoint identification. Default 100Kb.

        --ssB=SSB
                Extend -(ve) bp of the 5' HMM segment border of chromosome B for breakpoint identification. Default 100Kb.

        --seB=SEB
                Extend +(ve) bp of the 3' HMM segment border of chromosome B for breakpoint identification. Default 100Kb.


        --precheck=PRECHECK
                Precheck option will help to restrict HiCtrans search only to chromosome combinations with significant max interaction compared to mean non-zero count value [Default 1. Lower value will increase stringency]

        -h, --help
                Show this help message and exit

```

Users need to run each chromosome pair independently. This is a helper function to generate all the combination of chromosomal pairs and run hictrans.R
HiCtrans can be run with or without the restriction level validpair information file. 

If you don't have the validpair file, please use the following command

```bash
perl -e '@F=`cat $ARGV[0]`; for($i=0; $i<$#F; $i++){chomp $F[$i]; for($j=$i+1; $j<=$#F; $j++){chomp $F[$j]; print "Rscript hictrans.v3.R --mat $ARGV[1] --bed $ARGV[2] --chrA $F[$i] --chrB $F[$j] --prefix $ARGV[3] --resolutions 2,3,4,5,6,8,10 --covq 0.1 --chromsize chrom_hg19.sizes\n";}}' <chrom_name.file> <matrix.file> <bed.file> <prefix>
```

Example

```bash
Rscript ../hictrans.v3.R --mat T47D_20Kb_chr10_chr20.matrix --bed T47D_20Kb_chr10_chr20_abs.bed --chrA chr10 --chrB chr20 --prefix T47D_20Kb_chr10_chr20 --resolutions 2,3,4,5,6,8,10 --covq 0.1 --chromsize chrom_hg19.sizes
```

If you have the validpair file, please use the following

```bash
perl -e '@F=`cat $ARGV[0]`; for($i=0; $i<$#F; $i++){chomp $F[$i]; for($j=$i+1; $j<=$#F; $j++){chomp $F[$j]; print "Rscript ../HiCtrans/hictrans.v3.R --mat $ARGV[1] --bed $ARGV[2] --chrA $F[$i] --chrB $F[$j] --prefix $ARGV[3] --resolutions 2,3,4,5,6,8,10 --covq 0.1 --relevel YES --fragsFile Resfrag_hg19.bed --validpair $ARGV[4] --chromsize chrom_hg19.sizes --precheck 1e-5\n";}}' <chrom_name.file> <matrix.file> <bed.file> <prefix> <validpair.file>
```

Example (Validpair files are not provided in the example folder, but we assume the data by default is processed with HiC-pro and validpair file is available)

```bash
Rscript ../hictrans.v3.R --mat T47D_20Kb_chr10_chr20.matrix --bed T47D_20Kb_chr10_chr20_abs.bed --chrA chr10 --chrB chr20 --prefix T47D_20Kb_chr10_chr20 --resolutions 2,3,4,5,6,8,10 --covq 0.1 --relevel YES --fragsFile HindIII_resfrag_hg19.bed --validpair T47D_validpair.txt --chromsize chrom_hg19.sizes --precheck 1e-5
```


Here, chrom_name.file is a signle column file with chromsome names; matrix and bed files are names of the Hi-C sparse matrix and the associated bed files.
To generate the sparse matrix use the 'build_matrix.cpp' file (compile this program by running 'g++ build_matrix.cpp -o build_matrix' in your command prompt). 
For details of the program check the https://github.com/nservant/HiC-Pro repository.
The input to the build_matrix program is a validpair file described in the help section. 

If you are staring with HiCUP, then use hicup_filter to create valid Hi-C read pairs (generally ends with a name filt.bam or filt.sam).
Then use the following command to generate a validpair file from the filt.bam file

```bash
samtools view filt.bam| awk -v OFS='\t' '{print $1,$3,$4,"+"}' |paste - - |awk -v OFS='\t' '{print $1,$2,$3,$4,$6,$7,$8}' > hictrans.validpair
```

### Note
1. HiCtrans requires 10Kb resolution Hi-C matrix or multiple of 10Kb matrix to begin with
2. The bam file should be sorted based on read name.

# Filtering the result
HiCtrans recommends filtering the black listed regions (+/- 100Kb) from the output result.
These regions tends to produce higher interactions and can artificially appear as translocations.
Cheack the black list here for different genomes 
https://sites.google.com/site/anshulkundaje/projects/blacklists

# R library requirements:
data.table, hashmap, changepoint, hashmap, 
optparse, Rcpp, caTools, depmixS4 DEoptimR

# For troubleshoot:
Abhijit Chakraborty (abhijit@lji.org)
