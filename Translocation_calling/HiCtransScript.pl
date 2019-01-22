#HiCtrans script for translocation breakpoint finding from Hi-C data
#Author: Abhijit Chakraborty & Ferhat Ay
#Email: abhijit@lji.org & ferhatay@lji.org
#Run perl HiCtransScript.pl for details
#
use Getopt::Long;

sub HELP{
	print "\n\tUsage: perl HiCtransScript.pl -bed T47D.chr7_15.bed -mat T47D.chr7_15.matrix -chromA chr7 -chromB chr15 -prefix T47D -size chrom_hg19.sizes
               perl HiCtransScript.pl -bed Genome.bed -mat Genome.matrix -chromA chr7 -chrB chr15 -prefix T47D -vp T47D_validpairs.txt -refrags HindIII_REfrags.bed -size chrom_hg19.sizes\n\n";
	print "\t-bed: Genome wide binned bed file with index information\n
	chr1    0       40000   1
	chr1    40000   80000   2
	chr1    80000   120000  3
	chr1    120000  160000  4
	.....
	\n"; 
	print "\t-mat: Hi-C interaction file\n
	1	2	234
	1	3	100
	2	3	150
	3	4	110
	.....
	\n";
	print "\t-chromA: chromosome A to be searched in the Hi-C data (either chromosome name or all)\n\n";
	print "\t-chromB: chromosome B to be searched in the Hi-C data (either chromosome name or all)\n\n\tNOTE: For all possible pairs provide '-chromA all -chromB all' as argument\n\n";
	print "\t-prefix: Output folder or file name\n\n";
	print "\t-size: chromsome names and their sizes (NOTE: Ordering of the chromosome should be same as that of the binned bed file)\n
	chr1    249250621
	chr2    243199373
	chr3    198022430
	.....
	\n";
	print "\t-cf: Filter out intergenic-interactions below this contact count (Can be changed based on sample, default is 10)\n\n";
	print "\t-rf: This is the box.entropy by random.entropy ratio (This is will control the noise effect, default is 2)\n\n";
	print "\t-vp: Hi-C validpair fragment file\n
	SRR6213722.1   chr11   124331538       -       chr11   124345246       -
	SRR6213722.2   chr1    198436365       -       chr1    199923196       +
	.....
 	\n";
	print "\t-refrags: Restriction fragment file\n
	chr1    0       16007   HIC_chr1_1      0       +
	chr1    16007   24571   HIC_chr1_2      0       +
	chr1    24571   27981   HIC_chr1_3      0       +
	\n\tNOTE: Breakpoint finding is only possible if validpair and restriction fragment files are not provided. Otherwise only translocation boundaries are provided.\n\n";
	print "\t-cpu: Number of cpu to use (default 1)\n\n";
	print "\t-h: Help\n\n";
	exit;
}


GetOptions(
        'bed=s'     => \my $bed_file,
        'mat=s'     => \my $mat_file,
	'chromA=s'  => \my $chromA,
	'chromB=s'  => \my $chromB,
	'prefix=s'  => \my $prefix,
	'size=s'    => \my $size,
	'cf=s' => \my $count_filter,
	'rf=s' => \my $ratio_filter,
	'vp=s' => \my $valid_pair,
	'refrags=s' => \my $refrags_bed,
	'cpu=s'	=> \my $cpu,
	'h=s'   => \my $h,
)  or HELP;

chomp ($bed_file, $mat_file, $chromA, $chromB, $prefix, $size, $count_filter, $ratio_filter, $valid_pair, $refrags_bed, $h);

if ($cpu eq ""){$cpu = 1;}
if ($count_filter eq ""){$count_filter = 10;}
if ($ratio_filter eq ""){$ratio_filter = 2;}
if ($bed_file eq "" || $mat_file eq "" || $chromA eq "" || $chromB eq "" || $prefix eq "" || $size eq ""){HELP;}

sub run_hictrans{
	my $bed_file  = $_[0];
	my $mat_file  = $_[1];
   	my $chromA    = $_[2];
	my $chromB    = $_[3];
	my $prefix    = $_[4];
	my $size      = $_[5];
	my $cf	      = $_[6];
	my $rf	      = $_[7];
	my $valid_pair  = $_[8];
	my $refrags_bed = $_[9];
	
	print "Running step1: perl ../scripts/transNormScript.pl $bed_file $mat_file $chromA $chromB\n";
        system("perl ../scripts/transNormScript.pl $bed_file $mat_file $chromA $chromB") == 0 or die "Error in step1: $chromA-$chromB\n";

        print "Running step2: Rscript ../scripts/translocationFind.r $chromA-$chromB/$chromA-$chromB.count.matrix $chromA-$chromB\n";
        system("Rscript ../scripts/translocationFind.r $chromA-$chromB/$chromA-$chromB.count.matrix $chromA-$chromB") == 0 or die "Error in step2: $chromA-$chromB\n";

        print "Running step3: Rscript ../scripts/maxCountFilter.r $chromA-$chromB/$chromA-$chromB.fragsfile $chromA-$chromB/$chromA-$chromB.intersfile $chromA $chromB $chromA-$chromB\n";
        system("Rscript ../scripts/maxCountFilter.r $chromA-$chromB/$chromA-$chromB.fragsfile $chromA-$chromB/$chromA-$chromB.intersfile $chromA $chromB $chromA-$chromB") == 0 or die "Error in step3: $chromA-$chromB\n";

        print "Running step4: perl ../scripts/mapBack.pl $chromA $chromB $chromA-$chromB/$chromA-$chromB.fragsfile $chromA-$chromB.tmp.result > $chromA-$chromB.Translocation.result\n";
        system("perl ../scripts/mapBack.pl $chromA $chromB $chromA-$chromB/$chromA-$chromB.fragsfile $chromA-$chromB.tmp.result > $chromA-$chromB.Translocation.result") == 0 or die "Error in step4: $chromA-$chromB\n";

        print "Running step5: Rscript ../scripts/noiseReduction.r -f $chromA-$chromB/$chromA-$chromB.intersfile -t $chromA-$chromB.Translocation.result -s 10 -c $count_filter -n $ratio_filter -o $chromA-$chromB.Translocation.EntropyFiltered.result\n";
        system("Rscript ../scripts/noiseReduction.r -f $chromA-$chromB/$chromA-$chromB.intersfile -t $chromA-$chromB.Translocation.result -s 10 -c $count_filter -n $ratio_filter -o $chromA-$chromB.Translocation.EntropyFiltered.result") == 0 or die "Error in step5: $chromA-$chromB\n";

        if ($valid_pair ne "" && $refrags_bed ne ""){
		if (stat("$chromA-$chromB.Translocation.EntropyFiltered.result") ne ""){
                	print "For details check Rscript ../scripts/translocationREBreakPoint.r --help\n";
                	print "Running step6: Rscript ../scripts/translocationREBreakPoint.r --fragsFile $refrags_bed --chromsize $size --validpair $valid_pair --prefix $prefix --hictrans $chromA-$chromB.Translocation.EntropyFiltered.result\n";
                	system("Rscript ../scripts/translocationREBreakPoint.r --fragsFile $refrags_bed --chromsize $size --validpair $valid_pair --prefix $prefix --hictrans $chromA-$chromB.Translocation.EntropyFiltered.result") == 0 or die "Error in step6: $chromA-$chromB\n";
		}
        }
	$folder = $prefix."_".$chromA."_".$chromB."_Folder";
        if (stat($folder) eq ""){system("mkdir $folder");}
	system("mv $chromA-$chromB/ $folder/");
        system("mv $chromA-$chromB.tmp.result $folder/");
        system("mv $chromA-$chromB.Translocation.result $folder/");
	if (stat("$chromA-$chromB.Translocation.EntropyFiltered.result") ne ""){
        	system("mv $chromA-$chromB.Translocation.EntropyFiltered.result $folder/");
		system("mv $prefix.$chromA* $folder/");
		## Commented on 01/22/2019 ##
		#system("mv $prefix\_$chromA\_$chromB\_* $folder/");
		#system("mv $prefix.$chromB* $folder/");
		## END ##
	}
}

if ($chromA ne "all" && $chromB ne "all"){
	run_hictrans(
		$bed_file,
		$mat_file,		
		$chromA,
		$chromB,
		$prefix,
		$size,
		$cf,
		$rf,
		$valid_pair,
		$refrags_bed);

} elsif ($chromA eq "all" && $chromB eq "all"){

	use Parallel::ForkManager;		
	$MAX_PROCESSES = $cpu;
	$pm = new Parallel::ForkManager($MAX_PROCESSES);

	@chr = `awk '{print \$1}' $size`;	
	$i = 0;
	while ($i < $#chr){
		chomp $chr[$i];
		$j = $i+1;
                while ($j <= $#chr){
                        chomp $chr[$j];
			push @pairs,"$chr[$i]\t$chr[$j]\n";
			$j++;
		}
		$i++;
	}
	foreach $pair (@pairs){
		chomp $pair;
		@data = split(/\s+/,$pair);
		$chrA = @data[0];
		$chrB = @data[1];
		my $pid = $pm->start and next;
		print "Started: $chrA-$chrB\n";
		run_hictrans(
               		$bed_file,
                	$mat_file,
                	$chrA,
                	$chrB,
                	$prefix,
                	$size,
                	$cf,
                	$rf,
                	$valid_pair,
                	$refrags_bed);
		undef @data;
		$pm->finish;
	}
	undef @pairs;
	undef @chr;
} 
