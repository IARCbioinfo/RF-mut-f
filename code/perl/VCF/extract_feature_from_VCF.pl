###############################################################################
# Author: Alex Di Genova
# Laboratory: GCS/INRIA
# Copyright (c)
# year: 2019
###############################################################################
use Data::Dumper;
use Getopt::Std;
use FindBin;
use lib "$FindBin::Bin";
use VCF;


use strict;

sub usage {
   print "$0 usage : -a  -b  -c\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:c:r:", \%opts );
if ( !defined $opts{a}  ) {
   usage;
}

my $vcf = new VCF($opts{a});
#returns an array with a tag => desc values availables in the VCF file
my $feats=$vcf->get_features_tags();
my $skyp=0; #remove some features from matrix
if(defined $opts{b}){
	$skyp=1;
}
#release
my $release=0;
if(defined $opts{r}){
	$release=1;
}
$vcf->build_feature_table($feats,$skyp,$release);
