$input = @ARGV[0];
@input = `cat $input`;
chomp (@input[0], @input[1], @input[2], @input[3], @input[4]);
$chrA = @input[2];
$chrB = @input[3];
$folder = @input[4];

$stamp = `date|awk '{print \$3_\$2_\$6}'`; 
$job   = $$;
chomp ($stamp,$job);
open (out,">HiCtrans_$stamp.$job.sh");
if ($chrA ne "all" && $chrB ne "all"){
	print out "perl ../scripts/transNormScript.pl $input\n";
	print out "Rscript ../scripts/translocationFind.r $chrA-$chrB/$chrA-$chrB.count.matrix $chrA-$chrB\n";
	print out "Rscript ../scripts/maxCountFilter.r $chrA-$chrB/$chrA-$chrB.fragsfile $chrA-$chrB/$chrA-$chrB.intersfile $chrA $chrB $chrA-$chrB\n";
	print out "perl ../scripts/mapBack.pl $chrA $chrB $chrA-$chrB/$chrA-$chrB.fragsfile $chrA-$chrB.tmp.result > $chrA-$chrB.Translocation.result\n";
	print out "mkdir $folder\n";
	print out "mv $input $folder/\n";
	print out "mv $chrA-$chrB $folder/\n";
	print out "mv $chrA-$chrB.tmp.result $folder/\n";
	print out "mv $chrA-$chrB.Translocation.result $folder/";
	close out;
}
elsif ($chrA eq "all" && $chrB eq "all"){
	@chr = qw(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX);
	open (out,">HiCtrans_$stamp.$job.sh");
	$i = 0;
	while ($i < $#chr){
		chomp $chr[$i];
		print $chr[$i],"\n";
		$j = $i+1;
		while ($j <= $#chr){
			chomp $chr[$j];
			open (input_write,">$chr[$i]_$chr[$j].input.txt");
			print input_write "@input[0]\n";
			print input_write "@input[1]\n";	
			print input_write "$chr[$i]\n";
			print input_write "$chr[$j]\n";
			print input_write "$chr[$i]_$chr[$j]_Folder";
			close input_write;
			
			print out "perl ../scripts/transNormScript.pl $chr[$i]_$chr[$j].input.txt\n";
        		print out "Rscript ../scripts/translocationFind.r $chr[$i]-$chr[$j]/$chr[$i]-$chr[$j].count.matrix $chr[$i]-$chr[$j]\n";
        		print out "Rscript ../scripts/maxCountFilter.r $chr[$i]-$chr[$j]/$chr[$i]-$chr[$j].fragsfile $chr[$i]-$chr[$j]/$chr[$i]-$chr[$j].intersfile $chr[$i] $chr[$j] $chr[$i]-$chr[$j]\n";
        		print out "perl ../scripts/mapBack.pl $chr[$i] $chr[$j] $chr[$i]-$chr[$j]/$chr[$i]-$chr[$j].fragsfile $chr[$i]-$chr[$j].tmp.result > $chr[$i]-$chr[$j].Translocation.result\n";
        		print out "mkdir $chr[$i]_$chr[$j]_Folder\n";
        		print out "mv $chr[$i]_$chr[$j].input.txt $chr[$i]_$chr[$j]_Folder/\n";
        		print out "mv $chr[$i]-$chr[$j] $chr[$i]_$chr[$j]_Folder/\n";
        		print out "mv $chr[$i]-$chr[$j].tmp.result $chr[$i]_$chr[$j]_Folder/\n";
        		print out "mv $chr[$i]-$chr[$j].Translocation.result $chr[$i]_$chr[$j]_Folder/\n";
			$j++;
		}
		$i++;
	}
	print out "echo \"chrA\t5'-Boundary\tBreakPoint\t3'-Boundary\tchrB\t5'-Boundary\tBreakPoint\t3'-Boundary\tCount\" > All.chromosome.Translocation.result\n";
	print out "cat chr*_chr*_Folder/*.Translocation.result|grep -v \"Count\"|awk '{if(\$9 > 5){print}}' >> All.chromosome.Translocation.result\n";
	close out;
}
`chmod 755 HiCtrans_$stamp.$job.sh`;
