###############################################################################
# Author: Alex Di Genova 
# Laboratory: IARC/SCG/RCG
# Copyright (c)
# year: 2020
###############################################################################
use Data::Dumper;
use Getopt::Std;
use strict;


sub usage {
   print "$0 usage : -a <rf_model_snv_pred.txt> -b <rf_model_indel_pred.txt> -c <cutoff:0.5>  -d <annot> -s <sample> \n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:c:d:s:", \%opts );
if ( !defined $opts{a} or !defined $opts{b} or !defined $opts{d} or !defined $opts{s} ) {
   usage;
}

my $cuttof=0.5;
if(defined $opts{c}){
	$cuttof=$opts{c};
}

#we get the decision
my ($hsnp)=load_vars($opts{a},$cuttof,$opts{s});
my ($hindel)=load_vars($opts{b},$cuttof,$opts{s});



print STDERR "Total somatic SNVs loaded ($opts{s})  :". scalar(keys %{$hsnp}) ."\n";
print STDERR "Total somatic INDELs loaded ($opts{s}):". scalar(keys %{$hindel}) ."\n";
#print Dumper($hash);
open(FILE, $opts{d}) or die "cannot open file $opts{d}\n";

my $somatic=0;
my $bases={A=>1,C=>1,T=>1,G=>1};
my $total_lines=0;
while(my $line=<FILE>){
        chomp $line;
	$total_lines++;
        my @data=split /\t/,$line;
	if($line=~m/^Chr/){
		print $line."\n";
	}else{
		#we check if variants is SNP
		if(defined $bases->{$data[3]} and defined $bases->{$data[4]}){
			if(defined $hsnp->{join("__",$data[0],$data[1])}){
				#print $line." $somatic\n";
				print $line."\n";
				$somatic++;
			}else{
				#print STDERR "NOT FOUND\n".$line."\n";
			}
              }else{
			#we try an indel 
			if(defined $hindel->{join("__",$data[0],$data[1])}){
				print $line."\n";
				$somatic++;
			#annovar encode alt deletions increasing pos+1, example below:
			#chr1    941119  941119  A       -       intronic
			#MESO_084_T      chr1    941118  NN      MODIFIER#	
			#for ref deletions is not the case
			#chr1    3086979 3086979 -       A	
			#MESO_084_T      chr1    3086979 NN      MODIFIER	
			#INDEL case
			}elsif(defined $hindel->{join("__",$data[0],$data[1]-1)} and length($data[3]) > length($data[4])){
				#print $line." $somatic\n";
				print $line."\n";
				$somatic++;
			# handle the MNP case, such as CC>TT
			}elsif(length($data[3]) == length($data[4]) and $data[4] eq "-" and defined $hindel->{join("__",$data[0],$data[1]-1)}){
				#print STDERR "NOT FOUND\n".$line."\n";
				print $line."\n";
				$somatic++;
			}
	      }
	
	}
}


#SOMATIC
print STDERR "A total of $somatic variants were annotated for $opts{s}\n";
#print STDERR "Total lines  $total_lines in $opts{s}\n";

sub load_vars{
	my ($file,$cutoff, $sample)=@_;
	open(FILE,$file) or die "cannot open file $file\n";
	my $hash=();
	while(my $line=<FILE>){
	    	chomp $line;
        	my @data=split /\t/,$line;
		next if($data[1] ne $sample);
		#next if not somatic
		if($data[$#data-3] > $cuttof ){
			$hash->{join("__",$data[2],$data[3])}=$data[$#data-3];
		}
	}
	close(FILE);
	return $hash;
}


