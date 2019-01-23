#author abhijit
#created in 8th Sept, 2016
#

###########################
	## Commented: 12/11/2018 ##
#@file = `cat @ARGV[0]`;
#$bed  = @file[0];
#$mat  = @file[1];
#$chrA = @file[2];
#$chrB = @file[3];

	## Added: 12/11/2018 ##
$bed  = @ARGV[0];
$mat  = @ARGV[1];
$chrA = @ARGV[2];
$chrB = @ARGV[3];
###########################

#Change the $half_resolution as per the experiment
$half_resolution = 2e4

chomp ($bed,$mat,$chrA,$chrB);
if (stat("$chrA-$chrB/") eq ""){
	`mkdir $chrA-$chrB`;
}
open(bed,"$bed");
while (<bed>)
{
	chomp $_;
	@bedData = split(/\s+/,$_);
	if (@bedData[0] eq $chrA)
	{
		if ($bedCountAStart{@bedData[0]} eq "")
		{
			$bedCountAStart{@bedData[0]} = @bedData[3];
			$chrACount = @bedData[3]-1;
		}
		$bedChr{@bedData[3]}   = @bedData[0];
		$bedStart{@bedData[3]} = @bedData[1];
		$bedEnd{@bedData[3]}   = @bedData[2];
		$bedMid{@bedData[3]}   = @bedData[1]+$half_resolution;
		$chrACount++; 
	}
	if (@bedData[0] eq $chrB)
        {
		if ($bedCountBStart{@bedData[0]} eq "")
                {
                        $bedCountBStart{@bedData[0]} = @bedData[3];
                        $chrBCount = @bedData[3]-1;
                }
                $bedChr{@bedData[3]}   = @bedData[0];
                $bedStart{@bedData[3]} = @bedData[1];
                $bedEnd{@bedData[3]}   = @bedData[2];
		$bedMid{@bedData[3]}   = @bedData[1]+$half_resolution;
		$chrBCount++;
        }
	undef @bedData;	
} 
$bedCountAEnd{$chrA} = $chrACount;
$bedCountBEnd{$chrB} = $chrBCount;
close bed;

@hgbed = `cat ../data/hg19.fa.40000.genome_feature.txt`;  #Change here, if you using different genome/resolution and restriction enzyme
$i = 0;
while ($i <= $#hgbed)
{
        $c++;
        chomp $hgbed[$i];
        @bedData = split(/\s+/,$hgbed[$i]);
        $chrInfo{@bedData[3]}  = @bedData[0];
        $chrNum{@bedData[3]}   = @bedData[4];
        $chrStartBin{@bedData[3]} = @bedData[5];
        $chrEndBin{@bedData[3]}   = @bedData[6];
        $chrLen{@bedData[3]}   = @bedData[7];
        $chrGcc{@bedData[3]}   = @bedData[8];
        $chrMap{@bedData[3]}   = @bedData[9];
        $chrRegion{@bedData[3]}   = @bedData[10];
        if ($chrStart{@bedData[0]} eq "")
        {
                $chrStart{@bedData[0]} = @bedData[3];
        }
        if ($chrStart{@bedData[0]} < $c)
        {
                $chrEnd{@bedData[0]}   = @bedData[3];
        }
        undef @bedData;
        $i++;
}

open(chrA,">./$chrA-$chrB/$chrA.original.matrix");
open(chrB,">./$chrA-$chrB/$chrB.original.matrix");
open(chrAB,">./$chrA-$chrB/$chrA-$chrB.original.matrix");
open(mat,"$mat");
while (<mat>)
{
	chomp $_;
	@matData = split(/\s+/,$_);
	if ($bedChr{@matData[0]} eq $chrA && $bedChr{@matData[1]} eq $chrA)
	{
		print chrA "@matData[0]\t@matData[1]\t@matData[2]\n";
	}
	if ($bedChr{@matData[0]} eq $chrB && $bedChr{@matData[1]} eq $chrB)
        {
		print chrB "@matData[0]\t@matData[1]\t@matData[2]\n";
	}
	if ($bedChr{@matData[0]} eq $chrA && $bedChr{@matData[1]} eq $chrB)
        {
		print chrAB "@matData[0]\t@matData[1]\t@matData[2]\n";
	}
	if ($bedChr{@matData[0]} eq $chrB && $bedChr{@matData[1]} eq $chrA)
	{
		print chrAB "@matData[1]\t@matData[0]\t@matData[2]\n";
	}
	undef @matData;	
}
close chrA;
close chrB;
close chrAB;
close mat;

sub generateCisNormFiles{
	my $chr = $_[0];
	chomp $chr;
	open(matrix,"./$chrA-$chrB/$chr.original.matrix");
	while (<matrix>)
	{
	        chomp $_;
	        @matData = split(/\s+/,$_);
	        $matInfo{@matData[0]}{@matData[1]} = @matData[2];
	        $matInfo{@matData[1]}{@matData[0]} = @matData[2];
	        undef @matData;
	}
	close matrix;
	open(feature,">./$chrA-$chrB/$chr.feature");
	open(matrix,">./$chrA-$chrB/$chr.count.matrix");
	open(number,">./$chrA-$chrB/$chr.num.matrix");
	print feature "chr\tstart\tend\tlen\tgcc\tmap\n";
	if ($chr eq $chrA)
	{
		$start = $bedCountAStart{$chr};
		$end   = $bedCountAEnd{$chr};
	}
	elsif ($chr eq $chrB)
	{
		$start = $bedCountBStart{$chr};
                $end   = $bedCountBEnd{$chr};
	}
	print "$chr\t$start\t$end\n";
	my $i = $start;
	while ($i <= $end)
	{
	        if ($chrLen{$i} > 0 && $chrGcc{$i} >= 0.2 && $chrMap{$i} >= 0.5 && $chrRegion{$i} eq "I")
	        {
	                my $j = $start;
	                while ($j <= $end)
	                {
	                        if ($chrLen{$j} > 0 && $chrGcc{$j} >= 0.2 && $chrMap{$j} >= 0.5 && $chrRegion{$j} eq "I")
	                        {
	                                if ($matInfo{$i}{$j} ne "")
	                                {
	                                        if ($i == $j)
	                                        {
	                                                $line .= "0\t";
	                                                $numb .= "$j\t";
	                                        }
	                                        else
	                                        {
	                                                $line .= "$matInfo{$i}{$j}\t";
	                                                $numb .= "$j\t";
	                                        }
	                                }
	                                elsif ($matInfo{$i}{$j} eq "")
	                                {
	                                        $line .= "0\t";
       		                                $numb .= "$j\t";
                                	}
                        	}
                        	$j++;
                	}
                	print feature "$chrNum{$i}\t$chrStartBin{$i}\t$chrEndBin{$i}\t$chrLen{$i}\t$chrGcc{$i}\t$chrMap{$i}\n";
                	chop $line;
                	print matrix "$line\n";
                	print number "$numb\n";
                	$line = undef;
                	$numb = undef;
        	}
        	$i++;
	}
	close feature;
	close matrix;
	close number;
	undef %matInfo;
}


generateCisNormFiles($chrA);
generateCisNormFiles($chrB);

sub generateTransNormFiles{
        my $chra = $_[0];
	my $chrb = $_[1];
        chomp ($chra,$chrb);
        open(matrix,"./$chrA-$chrB/$chra-$chrb.original.matrix");
 	open(inters,">./$chrA-$chrB/$chra-$chrb.intersfile");
        while (<matrix>)
        {
                chomp $_;
                @matData = split(/\s+/,$_);
                $matInfo{@matData[0]}{@matData[1]} = @matData[2];
                $matInfo{@matData[1]}{@matData[0]} = @matData[2];
                undef @matData;
        }
        close matrix;
        open(matrix,">./$chrA-$chrB/$chra-$chrb.count.matrix");
        open(number,">./$chrA-$chrB/$chra-$chrb.num.matrix");
	print "$chra\t$chrb\t$bedCountAStart{$chra}\t$bedCountBEnd{$chrb}\n";	
	#Make sure chromosomes comes as similar as listed in the bed file. Otherwise start and end positions will be problematic.
        my $i = $bedCountAStart{$chra};
        while ($i <= $bedCountAEnd{$chra})
        {
               	if ($chrLen{$i} > 0 && $chrGcc{$i} >= 0.2 && $chrMap{$i} >= 0.5 && $chrRegion{$i} eq "I")
               	{
                       	my $j = $bedCountBStart{$chrb};
                       	while ($j <= $bedCountBEnd{$chrb})
                       	{
                       	        if ($chrLen{$j} > 0 && $chrGcc{$j} >= 0.2 && $chrMap{$j} >= 0.5 && $chrRegion{$j} eq "I")
                       	        {
                       	                if ($matInfo{$i}{$j} ne "")
                       	                {
                       	       	                if ($i == $j)
                       	        	      	{
                        	        		$line .= "0\t";
                        	        	      	$numb .= "$j\t";
                        	        	}
                        	      		else
                        	    		{
                        	       			$line .= "$matInfo{$i}{$j}\t";
                        	        	     	$numb .= "$j\t";
							 print inters "$bedChr{$i}\t$bedMid{$i}\t$bedChr{$j}\t$bedMid{$j}\t$matInfo{$i}{$j}\n";
                        	        	}
                        		}
                        	 	elsif ($matInfo{$i}{$j} eq "")
                        	   	{
                        	       		$line .= "0\t";
                        	         	$numb .= "$j\t";
                               	       	}
                                }
                                $j++;
                        }
                        chop $line;
                        print matrix "$line\n";
                        print number "$numb\n";
                        $line = undef;
                        $numb = undef;
		}
                $i++;
        }
        close matrix;
        close number;
	close inters;
        undef %matInfo;
}

generateTransNormFiles($chrA,$chrB);

@updatedBedA = `grep -v "start" $chrA-$chrB/$chrA.feature |awk '{c++;print "$chrA\\t1\\t"\$2+$half_resolution"\\t"c"\\t1"}'`;
@updatedBedB = `grep -v "start" $chrA-$chrB/$chrB.feature |awk '{c++;print "$chrB\\t1\\t"\$2+$half_resolution"\\t"c"\\t1"}'`;

open (updatedBed,">./$chrA-$chrB/$chrA-$chrB.fragsfile");
print updatedBed @updatedBedA;
print updatedBed @updatedBedB;
close updatedBed;
