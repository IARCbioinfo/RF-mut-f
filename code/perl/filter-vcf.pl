
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
   print "$0 usage : -a <rf_model_snv_pred.txt> -b <rf_model_indel_pred.txt> -c <cutoff:0.5>  -d <vcf> -s <sample> \n";
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


#print Dumper($hash);
print STDERR "Total somatic SNVs loaded ($opts{s})  :". scalar(keys %{$hsnp}) ."\n";
print STDERR "Total somatic INDELs loaded ($opts{s}):". scalar(keys %{$hindel}) ."\n";
#print Dumper($hash);
open(FILE,"gzip -dc $opts{d}|") or die "cannot open file $opts{d}\n";

my $somatic=0;

my $bases={A=>1,C=>1,T=>1,G=>1};

while(my $line=<FILE>){
	chomp $line;
        my @data=split /\t/,$line;
	if($line =~m/^#/){
		if($line =~m/#CHROM/){
			print "##FORMAT=<ID=RF_SP,Number=A,Type=Float,Description=\"Random Forest probability for somatic class\">\n";
			print $line."\n";
		}else{
		 	print $line."\n";
		}
		
	}else{
		#we check if variants is SNP
		if(defined $bases->{$data[3]} and defined $bases->{$data[4]}){
			if(defined $hsnp->{join("__",$data[0],$data[1])}){
				my $prob=$hsnp->{join("__",$data[0],$data[1])};
				$data[7]="RF_SP=$prob;".$data[7];
				print join("\t",@data)."\n";
				$somatic++;
			}
              }else{
			#we try an indel 
			if(defined $hindel->{join("__",$data[0],$data[1])}){
				my $prob=$hindel->{join("__",$data[0],$data[1])};
				$data[7]="RF_SP=$prob;".$data[7];
				print join("\t",@data)."\n";
				$somatic++;
			}
	      }
	}
}

#SOMATIC
print STDERR "A total of $somatic variants were annotated for $opts{s}\n";




sub load_vars{
	my ($file,$cutoff, $sample)=@_;
	open(FILE,$file) or die "cannot open file $file\n";
	my $hash=();
	while(my $line=<FILE>){
	    	chomp $line;
        	my @data=split /\t/,$line;
		next if($data[1] ne $sample);
		#print $data[$#data-4]."\n";
		#print $data[$#data-3]."\n";
		#next if not somatic
		if($data[$#data-3] > $cuttof ){
			#print join("\t",@data)."\n";
			$hash->{join("__",$data[2],$data[3])}=$data[$#data-3];
		}
	}
	close(FILE);
	return $hash;

}
