
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

my $cuttof=0;
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
open(NPVCF,">".$opts{s}.".non_pass_rf_filter.vcf") or die "cannot open file $opts{s}\n";
open(PFVCF,">".$opts{s}.".pass_rf_filter.vcf") or die "cannot open file $opts{s}\n";

my $somatic=0;
my $total=0;
my $bases={A=>1,C=>1,T=>1,G=>1};
my $is_cod=();
my $coding=0;
my $non_coding=0;
#coding is exonic + splicing
my $reg_types={"intergenic"=>0,"intronic"=>0,"exonic"=>1,"downstream"=>0,"upstream"=>0,
               "splicing"=>1,"UTR3"=>0,"UTR5"=>0,"ncRNA_UTR5"=>0,
               "ncRNA_exonic"=>0,"ncRNA_intronic"=>0,"ncRNA_splicing"=>0};


while(my $line=<FILE>){
	chomp $line;
        my @data=split /\t/,$line;
	if($line =~m/^#/){
		if($line =~m/#CHROM/){
			print PFVCF "##FORMAT=<ID=RF_SP,Number=A,Type=Float,Description=\"Random Forest probability for somatic class\">\n";
			print PFVCF "##FORMAT=<ID=IS_CODING,Number=A,Type=Integer,Description=\"Indicate if the variant is coding or not (Annovar classification)\">\n";
			print PFVCF $line."\n";

			print NPVCF "##FORMAT=<ID=RF_SP,Number=A,Type=Float,Description=\"Random Forest probability for somatic class\">\n";
			print NPVCF "##FORMAT=<ID=IS_CODING,Number=A,Type=Integer,Description=\"Indicate if the variant is coding or not (Annovar classification)\">\n";
			print NPVCF $line."\n";
		}else{
		 	print PFVCF $line."\n";
		 	print NPVCF $line."\n";
		}
		
	}else{
		$total++;
		#we check if the variant is coding or not
		my ($is_coding)=is_coding($data[7],$reg_types);
		#$is_cod->{$is_coding}++;
		#we check if variants is SNP
		if(defined $bases->{$data[3]} and defined $bases->{$data[4]}){
			if(defined $hsnp->{join("__",$data[0],$data[1])}){
				my $prob=$hsnp->{join("__",$data[0],$data[1])};
				$data[7]="IS_CODING=$is_coding;RF_SP=$prob;".$data[7];
				if($is_coding and $prob > 0.5){
					print PFVCF join("\t",@data)."\n";
					$coding++;
					$somatic++;
				}elsif($prob > 0.75){
				 	print PFVCF join("\t",@data)."\n";
					$non_coding++;
					$somatic++;
				}else{
				 	print NPVCF  join("\t",@data)."\n";
				}	
				
			}else{
				 print NPVCF  join("\t",@data)."\n";
			}
              }else{
			#we try an indel 
			if(defined $hindel->{join("__",$data[0],$data[1])}){
				my $prob=$hindel->{join("__",$data[0],$data[1])};
				$data[7]="IS_CODING=$is_coding;RF_SP=$prob;".$data[7];
				if($is_coding and $prob > 0.5){
					print PFVCF join("\t",@data)."\n";
					$coding++;
					$somatic++;
				}elsif($prob > 0.75){
				 	print PFVCF join("\t",@data)."\n";
					$non_coding++;
					$somatic++;
				}else{
				 	print NPVCF join("\t",@data)."\n";
				}
				#print join("\t",@data)."\n";
			}else{
				#var not classified
				 print NPVCF join("\t",@data)."\n";
			}
	      }
	}
}

#SOMATIC
print STDERR "A total of $somatic [coding=$coding, non_coding=$non_coding ] variants were classified as somatic for $opts{s} (Total=$total)\n";
#print Dumper($is_cod);

sub is_coding{
	my ($line,$types)=@_;
	my @tags=split(";",$line);
	my $hash=();
	foreach my $tt(@tags){
		my ($t,$v)=split("=",$tt);
		#next if($v eq ".");
		if($v=~m/x3/){
			my $p=index($v,"x3");
			#print STDERR $v." match"." at pos ".$p." sub ".substr($v,0,$p-1)."\n";
			my $a=substr($v,0,$p-1);
			$v=$a;
		}
		$hash->{$t}=$v;
	}	
	#print Dumper($hash);
	#print Dumper($hash->{"Func.ensGene"});
	#print $hash->{"Func.ensGene"}."\n";
	if(defined $types->{$hash->{"Func.ensGene"}}){
			return  $types->{$hash->{"Func.ensGene"}};
	}else{
		print STDERR $hash->{"Func.ensGene"}." not found in types\n";
		print STDERR Dumper($hash->{"Func.ensGene"});
		print STDERR Dumper($types);
	}
	#return $hash->{"Func.ensGene"};
		
	#return $hash->{"ExonicFunc.ensGene"};
}


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
		#if($data[$#data-3] > $cuttof ){
			#print join("\t",@data)."\n";
			$hash->{join("__",$data[2],$data[3])}=$data[$#data-3];
		#}
	}
	close(FILE);
	return $hash;

}
