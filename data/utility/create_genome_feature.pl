$chr{"chr1"} = 1;
$chr{"chr2"} = 2;
$chr{"chr3"} = 3;
$chr{"chr4"} = 4;
$chr{"chr5"} = 5;
$chr{"chr6"} = 6;
$chr{"chr7"} = 7;
$chr{"chr8"} = 8;
$chr{"chr9"} = 9;
$chr{"chr10"} = 10;
$chr{"chr11"} = 11;
$chr{"chr12"} = 12;
$chr{"chr13"} = 13;
$chr{"chr14"} = 14;
$chr{"chr15"} = 15;
$chr{"chr16"} = 16;
$chr{"chr17"} = 17;
$chr{"chr18"} = 18;
$chr{"chr19"} = 19;
$chr{"chr20"} = 20;
$chr{"chr21"} = 21;
$chr{"chr22"} = 22;
$chr{"chrX"} = 23;
$chr{"chrY"} = 24;
$chr{"chrM"} = 25;

$bed_file = @ARGV[0];
$gc_file  = @ARGV[1];
$map_file = @ARGV[2];
chomp ($bed_file,$gc_file,$map_file);
open (gc_in,"$gc_file");
while (<gc_in>){
	chomp $_;
	@gcData = split(/\s+/,$_);
	$gc_info{@gcData[0]}{@gcData[1]}{@gcData[2]} = sprintf("%0.3f",@gcData[4]);
	undef @gcData;
}
close gc_in; 
open (map_in,"$map_file");
while (<map_in>){
        chomp $_;
        @mapData = split(/\s+/,$_);
        $map_info{@mapData[0]}{@mapData[1]}{@mapData[2]} = sprintf("%0.3f",@mapData[4]);
	$bl_info{@mapData[0]}{@mapData[1]}{@mapData[2]}  = @mapData[5];
        undef @mapData;
}
close map_in;
open (bed_in,"$bed_file");
$c = 0;
while (<bed_in>){
	$c++;
	chomp $_;
	@bedData = split(/\s+/,$_);
	print "@bedData[0]\t@bedData[1]\t@bedData[2]\t$c\t$chr{@bedData[0]}\t@bedData[1]\t@bedData[2]\t1000\t$gc_info{@bedData[0]}{@bedData[1]}{@bedData[2]}\t$map_info{@bedData[0]}{@bedData[1]}{@bedData[2]}\t$bl_info{@bedData[0]}{@bedData[1]}{@bedData[2]}\n";
	undef @bedData;
	$i++;
} 
close bed_in;
