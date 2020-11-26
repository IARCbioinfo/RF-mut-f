
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
   print "$0 usage : -a  -b  -c\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:c:", \%opts );
if ( !defined $opts{a}  ) {
   usage;
}

my ($hash)=load_vars($opts{a},$opts{b});
#print Dumper($hash);
open(FILE,"gzip -dc $opts{c}|") or die "cannot open file $opts{c}\n";

while(my $line=<FILE>){
	chomp $line;
        my @data=split /\t/,$line;
	if($line =~m/^#/){
		print $line."\n";
	}else{
		if(defined $hash->{$data[0]}->{$data[1]}){
			print $line."\n";
		}
	}
}







sub load_vars{
	my ($file,$sample)=@_;
	open(FILE,$file) or die "cannot open file $file\n";
	my $hash=();
	while(my $line=<FILE>){
	    	chomp $line;
        	my @data=split /\t/,$line;
		next if ($data[1] ne $sample);
		#next if not somatic
		next if($data[$#data] != 1);
		if($data[2] =~m/chr/){
					
		}else{
			$data[2]="chr".$data[2];
		}
		if($data[2] =~m/MT/){
			#print STDERR $line."\n";
			$data[2]="chrM";
		}
		#print join("\t",@data)."\n";
		$hash->{$data[2]}->{$data[3]}=$data[4];
	}
	close(FILE);
	return $hash;

}
