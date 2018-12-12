$chrA = @ARGV[0];
$chrB = @ARGV[1];
@fragsFile = `cat @ARGV[2]`;
@binomFile = `cat @ARGV[3]`;
@result = `grep -v "score" ./$chrA-$chrB/data.txt|grep -v "Inf"|sed -e 's/"//g'|sed 's/:/ /g'|awk '{print \$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10\"\\t\"\$11}'`;

$i = 0;
while ($i <= $#fragsFile){
	chomp $fragsFile[$i];
	@fragsData = split(/\s+/,$fragsFile[$i]);
	$index{@fragsData[3]}{@fragsData[0]} = @fragsData[2];
	undef @fragsData;
	$i++;
}
########### Abhijit: 12/10/2018 ###########
#print "chrA\t5'-Boundary\tBreakPoint\t3'-Boundary\tchrB\t5'-Boundary\tBreakPoint\t3'-Boundary\tCount\n";
print "chrA\t5'-Boundary\t3'-Boundary\tchrB\t5'-Boundary\t3'-Boundary\tMaxCount\n";
###########
$i = 0;
while ($i <= $#result){
	$check = 1;
	chomp $result[$i];
	@resultData = split(/\s+/,$result[$i]);
	$chrA_lb = $index{@resultData[0]}{$chrA};
	$chrA_ub = $index{@resultData[1]}{$chrA};
	$chrB_lb = $index{@resultData[2]}{$chrB};
        $chrB_ub = $index{@resultData[3]}{$chrB}; 
 
        ############ Abhijit: 12/10/2018 ############
        
	$j = 0;
	while ($j <= $#binomFile){
        	chomp $binomFile[$j];
		@binomData = split(/\s+/,$binomFile[$j]);
		$chrA_bp = @binomData[1];
		$chrB_bp = @binomData[3];
		if ($chrA_bp >= $chrA_lb && $chrA_bp <= $chrA_ub && $chrB_bp >= $chrB_lb && $chrB_bp <= $chrB_ub && @binomData[4] > 6){	
			if ($check > @binomData[5]){
				$check = @binomData[5];
   				########### Abhijit: 12/10/2018 ###########
				#$breakPoint{$i} = "@binomData[0]\t$chrA_lb\t@binomData[1]\t$chrA_ub\t\t@binomData[2]\t$chrB_lb\t@binomData[3]\t$chrB_ub\t@binomData[4]";
				$breakPoint{$i} = "@binomData[0]\t$chrA_lb\t$chrA_ub\t@binomData[2]\t$chrB_lb\t$chrB_ub\t@binomData[4]";
				###########
			}
			push @pval,"@binomData[5] #@binomData[0]\t$chrA_lb\t@binomData[1]\t$chrA_ub\t\t@binomData[2]\t$chrB_lb\t@binomData[3]\t$chrB_ub\t@binomData[4]\n";
		}
		undef @binomData;
        	$j++;
	}
	if ($breakPoint{$i} ne ""){
		print "$breakPoint{$i}\n";	
	}
	########################
	$i++;
}
undef %breakPoint;
undef @fragsFile;
undef @binomFile;
undef @result;
#`rm @ARGV[3]`;
