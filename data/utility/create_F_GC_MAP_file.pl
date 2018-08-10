#Download the file wget ftp://hgdownload.cse.ucsc.edu/gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign50mer.bw
#bigWigToBedGraph wgEncodeCrgMapabilityAlign50mer.bw hg19.MapabilityAlign50mer.bedGraph
#
####### User specific parameters #######
$re_name = "HindIII"; #Change the restriction enzyme name

$frag = "HindIII_resfrag_hg19.bed"; #Change the restriction fragment file as per your experiment. Some restriction fragment files for hg/mm genomes are provided in the "HindIII_resfrag_files.zip" file 

$frag_size = 500;

$read_length = 50; #Hi-C fastq read length

$frag_file_name = "HindIII_hg19.$frag_size.$read_length";

$fasta_file = "hg19.fa"; #Needs hg19.fa or equivalent genome sequence file

$map_bedgraph = "hg19.MapabilityAlign50mer.bedGraph"; #Need to download the mappability file separately 

$resolution = 40000; #Change the HiC resolution as per your requirement
 
$chrom_size_file = "chrom_hg19.sizes"; #Chromosome size file

$bedtools = "/mnt/BioApps/bedtools/bin/bedtools"; #bedtools path

$blacklisted_region = "EncodeExcludableRegions.hg19.bed"; # Keep this variable as "No" if you don't have the black listed region information
#########################################
$frag_half = $frag_size/2;

$name = "$re_name.$read_length"."mer.$resolution";

open (out,">F_GC_MAP.$name.sh");
print out "awk '{print \$1\"\\t\"\$2\"\\t\"\$2+$frag_half\"\\t\"\$4\"\\t\"\$2\"\\t\"\$3\"\\t\"\$3-\$2\"\\n\"\$1\"\\t\"\$3-$frag_half\"\\t\"\$3\"\\t\"\$4\"\\t\"\$2\"\\t\"\$3\"\\t\"\$3-\$2}' $frag|awk '{if(\$2 >= 0){print}}'|sortBed > $frag_file_name.bed\n";
print out "bedtools nuc -fi $fasta_file -bed $frag_file_name.bed > $frag_file_name.GC.bed\n";
print out "bedtools map -a $frag_file_name.GC.bed -b $map_bedgraph -c 4 -o mean > $frag_file_name.GC_Map.bed\n";
print out "perl gc_map_per_fragment.pl $frag_file_name.GC_Map.bed $frag > $frag_file_name.F_GC_MAP.bed\n";
print out "$bedtools makewindows -g $chrom_size_file -w $resolution -i winnum > $fasta_file.$resolution.bed\n";
print out "$bedtools sort -i $fasta_file.$resolution.bed > sorted.$fasta_file.$resolution.bed\n";
print out "$bedtools sort -i $frag_file_name.F_GC_MAP.bed > sorted.$frag_file_name.F_GC_MAP.bed\n";
print out "$bedtools map -a sorted.$fasta_file.$resolution.bed -b sorted.$frag_file_name.F_GC_MAP.bed -c 4 -o mean -null 0 > sorted.$fasta_file.$resolution.GC.bed\n";
if ($blacklisted_region eq "No"){
	print out "$bedtools map -a sorted.$fasta_file.$resolution.bed -b sorted.$frag_file_name.F_GC_MAP.bed -c 5 -o mean -null 0|awk '{print \$0\"\\tI\"}' > sorted.$fasta_file.$resolution.MAP.BL.bed\n";
}
else {
	print out "$bedtools map -a sorted.$fasta_file.$resolution.bed -b sorted.$frag_file_name.F_GC_MAP.bed -c 5 -o mean -null 0 > sorted.$fasta_file.$resolution.MAP.bed\n";
	print out "$bedtools map -a sorted.$fasta_file.$resolution.MAP.bed -b $blacklisted_region -c 4 -o collapse -null I|cut -d, -f 1 > sorted.$fasta_file.$resolution.MAP.BL.bed\n";
}
print out "perl create_genome_feature.pl $fasta_file.$resolution.bed sorted.$fasta_file.$resolution.GC.bed sorted.$fasta_file.$resolution.MAP.BL.bed > $fasta_file.$name.genome_feature.txt\n";
print out "rm sorted.*.bed";
close out;
`chmod 755 F_GC_MAP.$name.sh`;
