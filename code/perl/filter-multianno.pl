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

my $cuttof=0;
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
open(NPVCF,">".$opts{s}.".non_pass_rf_filter.hg38_multianno.txt") or die "cannot open file $opts{s}\n";
open(PFVCF,">".$opts{s}.".pass_rf_filter.hg38_multianno.txt") or die "cannot open file $opts{s}\n";


my $somatic=0;
my $bases={A=>1,C=>1,T=>1,G=>1};
my $total_lines=0;
my $total=0;
#coding is exonic + splicing
my $reg_types={"intergenic"=>0,"intronic"=>0,"exonic"=>1,"downstream"=>0,"upstream"=>0,
               "splicing"=>1,"UTR3"=>0,"UTR5"=>0,"ncRNA_UTR5"=>0,
               "ncRNA_exonic"=>0,"ncRNA_intronic"=>0,"ncRNA_splicing"=>0};


               my $coding=0;
               my $non_coding=0;
while(my $line=<FILE>){
        chomp $line;
	$total_lines++;
        my @data=split /\t/,$line;
	if($line=~m/^Chr/){
    my @h=split("\t",$line);
    #we add the two values in the array
    #splice(@h,@h,4,"RF_SP","IS_CODING");
		print NPVCF join("\t",@h[0 .. 4],"RF_SP","IS_CODING",@h[5 .. $#h])."\n";
    print PFVCF join("\t",@h[0 .. 4],"RF_SP","IS_CODING",@h[5 .. $#h])."\n";
	}else{
    $total++;
		#we check if variants is SNP
    my ($is_coding)=is_coding($data[5],$reg_types);
    #we try a SNV variant
		if(defined $bases->{$data[3]} and defined $bases->{$data[4]}){
			if(defined $hsnp->{join("__",$data[0],$data[1])}){
			  my $prob=$hsnp->{join("__",$data[0],$data[1])};
        if($is_coding and $prob > 0.5){
          #splice(@data,@data,4,$prob,$is_coding);
          print PFVCF join("\t",@data[0 .. 4],$prob,$is_coding,@data[5 .. $#data] )."\n";
          $coding++;
          $somatic++;
        }elsif($prob > 0.75){
          print PFVCF join("\t",@data[0 .. 4],$prob,$is_coding,@data[5 .. $#data] )."\n";
          $non_coding++;
          $somatic++;
        }else{
          print NPVCF join("\t",@data[0 .. 4],$prob,$is_coding,@data[5 .. $#data] )."\n";
        }
			}else{
				#print STDERR "NOT FOUND\n".$line."\n";
        #variant not classified 4 in MESO_084_T
         print NPVCF join("\t",@data[0 .. 4],"NA","NA",@data[5 .. $#data] )."\n";
			}
     }else{
			#we try an indel
      #annovar encode alt deletions increasing pos+1, example below:
			#chr1    941119  941119  A       -       intronic
			#MESO_084_T      chr1    941118  NN      MODIFIER#
			#for ref deletions is not the case
			#chr1    3086979 3086979 -       A
			#MESO_084_T      chr1    3086979 NN      MODIFIER

			if(defined $hindel->{join("__",$data[0],$data[1])}){
				#print $line."\n";
        my $prob=$hindel->{join("__",$data[0],$data[1])};
        if($is_coding and $prob > 0.5){
          #splice(@data,@data,4,$prob,$is_coding);
          print PFVCF join("\t",@data[0 .. 4],$prob,$is_coding,@data[5 .. $#data] )."\n";
          $coding++;
          $somatic++;
        }elsif($prob > 0.75){
          print PFVCF join("\t",@data[0 .. 4],$prob,$is_coding,@data[5 .. $#data] )."\n";
          $non_coding++;
          $somatic++;
        }else{
          print NPVCF join("\t",@data[0 .. 4],$prob,$is_coding,@data[5 .. $#data] )."\n";
        }
				#$somatic++;
      	#INDEL case
			}elsif(defined $hindel->{join("__",$data[0],$data[1]-1)} and length($data[3]) > length($data[4])){
				#print $line." $somatic\n";
        my $prob=$hindel->{join("__",$data[0],$data[1]-1)};
				#print $line."\n";
				#$somatic++;
        if($is_coding and $prob > 0.5){
          #splice(@data,@data,4,$prob,$is_coding);
          print PFVCF join("\t",@data[0 .. 4],$prob,$is_coding,@data[5 .. $#data] )."\n";
          $coding++;
          $somatic++;
        }elsif($prob > 0.75){
          print PFVCF join("\t",@data[0 .. 4],$prob,$is_coding,@data[5 .. $#data] )."\n";
          $non_coding++;
          $somatic++;
        }else{
          print NPVCF join("\t",@data[0 .. 4],$prob,$is_coding,@data[5 .. $#data] )."\n";
        }

			# handle the MNP case, such as CC>TT
			}elsif(length($data[3]) == length($data[4]) and $data[4] eq "-" and defined $hindel->{join("__",$data[0],$data[1]-1)}){
				#print STDERR "NOT FOUND\n".$line."\n";
				#print $line."\n";
        my $prob=$hindel->{join("__",$data[0],$data[1]-1)};
        if($is_coding and $prob > 0.5){
          #splice(@data,@data,4,$prob,$is_coding);
          print PFVCF join("\t",@data[0 .. 4],$prob,$is_coding,@data[5 .. $#data] )."\n";
          $coding++;
          $somatic++;
        }elsif($prob > 0.75){
          print PFVCF join("\t",@data[0 .. 4],$prob,$is_coding,@data[5 .. $#data] )."\n";
          $non_coding++;
          $somatic++;
        }else{
          print NPVCF join("\t",@data[0 .. 4],$prob,$is_coding,@data[5 .. $#data] )."\n";
        }
				#$somatic++;
			}
	   }

	}
}


#SOMATIC
#print STDERR "A total of $somatic variants were annotated for $opts{s}\n";
print STDERR "A total of $somatic [coding=$coding, non_coding=$non_coding ] variants were classified as somatic for $opts{s} (Total=$total)\n";

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
		#if($data[$#data-3] > $cuttof ){
			$hash->{join("__",$data[2],$data[3])}=$data[$#data-3];
		#}
	}
	close(FILE);
	return $hash;
}

#ask if var is coding using the types
sub is_coding{
	my ($line,$types)=@_;
	my @tags=split(";",$line);
	if(defined $types->{$tags[0]}){
			return  $types->{$tags[0]};
	}else{
		print STDERR $tags[0]." not found in types\n";
		print STDERR Dumper($tags[0]);
		print STDERR Dumper($types);
	}
}
