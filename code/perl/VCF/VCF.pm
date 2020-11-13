package VCF;


=head1 NAME

VCF

=head1 DESCRIPTION

This object perform several operation on VCFs files

=head2 Available methods


=cut

use strict;
use Data::Dumper;

sub new{
  my ($packagename, $vcf) = @_;
  my $self = {vcffile => $vcf};
  bless ($self, $packagename);
  return ($self);
}


sub _get_autosomes_X_numbers{
  my @chrs=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X");
  return @chrs;
}

sub _get_autosomes_X_chr{
  my @chrs=("chr1","chr2","chr3","chr4","chr5","chr6","chr7",
  "chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
  "chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX");
  return @chrs;
}

sub _get_hash_autosomes_X_chr{
    my $self= shift;
    my $hchr=();
    foreach my $c($self->_get_autosomes_X_chr()){
        $hchr->{$c}=1;
    }
    return $hchr;
}

sub _parse_ctg_entry{
      my $self=shift;
      my $entry=shift;
      chomp $entry;
      my $s=index($entry,"<");
      my $e=index($entry,">");
      #print join(" ",$s,$e,$entry,substr($entry,$s+1,abs($s+1-$e)))."\n";
      my $ht=();
      foreach my $en(split(",",substr($entry,$s+1,abs($s+1-$e)))){
        my ($t,$v)=split("=",$en);
        #print join(" ",$t,$v);
          $ht->{$t}=$v;
      }
      #print join(" ",$s,$e,$entry,substr($entry,$s+3,abs($s+3-$e)))."\n";
      #print $e;
      #print Dumper($ht);
      return $ht;
}


sub filter_by_region{
    my $self=shift;
    my $prefix=shift;
    my $hchr=$self->_get_hash_autosomes_X_chr();
    #print Dumper($hchr);
    #filter in header
    ##contig=<ID=chr1,length=248956422>
    open(VCF, $self->{vcffile}) or die "cannot open VCF file\n";
    my $p =substr($self->{vcffile},0,length($self->{vcffile})-3)."autosomes_plus_X.vcf";
    open(OUT, ">".$p) or die "cannot create the output file\n";

    while (my $line=<VCF>) {
      if ($line =~m/^#/){
        #print $line;
          if($line=~m/#contig=/){
              my ($e)=$self->_parse_ctg_entry($line);
              print OUT $line if(defined $hchr->{$e->{ID}});
          }else{
            print OUT $line;
          }
      }else{
        #is an entry
        my @tmp=split("\t",$line);
        print OUT $line if(defined $hchr->{$tmp[0]});
      }
    }
    close(OUT);
}
#for the moment insertion and deletions only
sub add_type_entry{
  my $self=shift;
  my $entry=shift;
  my $m=shift;

  chomp $entry;
  my @d=split("\t",$entry);

 #chrX    138269888       .       C       CAATGTCATCAGTTAAGGCAGGAACAGGCCATTTTCACTTCTTTTGTGGTGG    60      .       .       GT      1/1
 #PRECISE;SVMETHOD=PAFTS;SVTYPE=INS;SVLEN=325
 my $info="";
 #deletion
 if(length($d[3]) > length($d[4])){
     $info="PRECISE;SVMETHOD=$m;SVTYPE=DEL;SVLEN=-".length($d[3]);
 }else{
   #insertion
   $info="PRECISE;SVMETHOD=$m;SVTYPE=INS;SVLEN=".length($d[4]);
 }
return join("\t",@d[0 .. 6],$info,@d[8 .. $#d]);

}


sub add_sv_type{
  my $self=shift;
  my $m=shift;
  if(length($m) == 0 ){
    $m="Denovo";
  }

  open(VCF, $self->{vcffile}) or die "cannot open VCF file\n";
  my $p =substr($self->{vcffile},0,length($self->{vcffile})-3)."svtypes.vcf";
  open(OUT, ">".$p) or die "cannot create the output file\n";
  my $types=<<EOF;
##fileformat=VCFv4.3
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=INVDUP,Description="InvertedDUP with unknown boundaries">
##ALT=<ID=TRA,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
##INFO=<ID=UNRESOLVED,Number=0,Type=Flag,Description="An insertion that is longer than the read and thus we cannot predict the full size.">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
EOF
my $info=1;
  while (my $line=<VCF>) {
    if ($line =~m/^#/){
        if($line=~m/^##INFO/ and $info == 1){
           print OUT $types;
           print OUT $line;
           $info++;
        }else{
          print OUT $line;
        }
    }else{
        #we are in an entry
        print OUT $self->add_type_entry($line,$m)."\n";
    }
  }
    close(OUT);

}

#parse tags values from VCF files
sub _get_tags{
    my $self=shift;
    my $entry=shift;
    my $tags=();
    chomp $entry;
    my @d=split("\t",$entry);
    my @vals=split(";",$d[7]);
    foreach my $t(@vals){
      my ($tag,$v)=split("=",$t);
          if(length($v) > 0){
          $tags->{$tag}=$v;
          }else{
          #for flag tags present in VCF files
          $tags->{$tag}=1;
         }
    }
    #we parse the genotype tags
    my @tvals=split(":",$d[8]);
     @vals=split(":",$d[9]);
    #print Dumper(@tvals,@vals);
    my $j=0;
    foreach my $v(@vals){
      my @vv=split(",",$v);
        my $i=1;
        foreach my $val(@vv){
          $tags->{$tvals[$j]."_".$i}=$val;
          $i++;
        }
        $j++;
    }
    #print Dumper($tags);
    return $tags;
}


sub keep_indels_svs{
      my $self=shift;
      open(VCF, $self->{vcffile}) or die "cannot open VCF file\n";
      my $p =substr($self->{vcffile},0,length($self->{vcffile})-3)."sv_insdel.vcf";
      open(OUT, ">".$p) or die "cannot create the output file\n";

      while (my $line=<VCF>) {
        if ($line =~m/^#/){
              print OUT $line;
        }else{
            #we are in an entry
            my $tags=$self->_get_tags($line);
            # print Dumper($tags);
            if($tags->{SVTYPE} eq "INS" or $tags->{SVTYPE} eq "DEL"){
                if(abs($tags->{SVLEN}) >= 50 ){
                  print OUT $line;
                }
            }
        }
      }
        close(OUT);
}

#get a list of features from VCF files
sub get_features_tags{
  my $self=shift;
  open(VCF, "bcftools  view $self->{vcffile} |") or die "cannot open VCF file\n";
  my $types=(); # array of values
  my $htypes=(); # has with counts for each feature stored in VCF file
  while (my $line=<VCF>) {
   if($line=~m/^##INFO|^##FORMAT/){
      $line=~s/##FORMAT=<//;
      $line=~s/##INFO=<//;
      my $tmp=();
      #we catch the description
      if($line =~m/Description="(.*)"/){
        #print $1."\n";
        $tmp->{desc}=$1;
      }
      my ($id,undef,$type)=split(",",$line);
      #ID=F1R2 Type=Integer
      $id=~s/ID=//;
      $type=~s/Type=//;
      $tmp->{id}=$id;
      $tmp->{type}=$type;
      #print join(" ",$tmp->{id},$tmp->{type},$tmp->{desc})."\n";
      push(@{$types},$tmp);
      #we init a counter for each type
      $htypes->{$tmp->{id}}=0;
   }elsif($line !~m/^#/){
     #we parse the tags
     my $tags=$self->_get_tags($line);
     $htypes->{$_}++ foreach(keys %{$tags});
     #print Dumper($tags);
     #if(defined $tags->{BCSQ}){
      # print $tags->{BCSQ}."\n";
     #}
    }
  }
  close(VCF);
  #print Dumper($htypes);
  my $ftypes=();
  foreach my $t(@{$types}){
     #feat not defined in the vcf of tumor-only
     next if($t->{type} eq "NALOD" or $t->{type} eq "NLOD");
     if($htypes->{$t->{id}} > 0){
       print join("\t","#",$self->{vcffile},$t->{id},$t->{type},$htypes->{$t->{id}},$t->{desc})."\n";
       $t->{count}=$htypes->{$t->{id}};
       push(@{$ftypes},$t);
     }
  }
  return $ftypes;
}

#we build the matrix of features
sub build_feature_table{
   my $self=shift;
   my $feats=shift;
   my $skyp=shift;	
   my $rel=shift;
	
   #not common feats NALOD or NALOD
   #my @f= sort @{$feats};
   my %f = map { $_->{id} => 1 } @{$feats};
   my @fs = sort keys %f;
   #print Dumper(@fs);
   my $prefix=$self->{vcffile};
    $prefix=~s/.vcf.bgz//;
   open(VCF, "bcftools  view $self->{vcffile} |") or die "cannot open VCF file\n";
   open(SNV,">".$prefix.".snv.matrix") or die "cannnot create matrix file\n";
   open(INDEL,">".$prefix.".indel.matrix") or die "cannnot create matrix file\n";

   #feat to check
#$VAR1 = 'BCSQ'; #MODIFIER, LOW, MODERATE, HIGH
#$VAR2 = 'CENTROMER'; # 1, 0
#$VAR3 = 'CGENE'; # Cosmic consensus gene
#$VAR4 = 'CNAME'; # CENTROMER name (ignored)
#$VAR5 = 'CONTQ'; # Mutect var
#$VAR6 = 'COSMIC'; # site in COSMIC
#$VAR7 = 'COSMIC_CENSUS_GENE'; # Sites listed in COSMIC_CENSUS_GENE
#$VAR8 = 'DP';
#$VAR9 = 'ECNT';
#$VAR10 = 'GERMQ';
#$VAR11 = 'GNOMAD_AC'; #AC
#$VAR12 = 'GNOMAD_AN'; #Total number of alleles in samples
#$VAR13 = 'MBQ';
#$VAR14 = 'MFRL';
#$VAR15 = 'MMQ';
#$VAR16 = 'MPOS';
#$VAR17 = 'NO_GNOMAD'; # 0
#$VAR18 = 'POPAF';
#$VAR19 = 'RPA';
#$VAR20 = 'RU';
#$VAR21 = 'SEQQ';
#$VAR22 = 'STR';
#$VAR23 = 'STRANDQ';
#$VAR24 = 'STRQ';
#$VAR25 = 'TLOD';
 #my $manuals={"CENTROMER"=>1,"CGENE"=>1,"COSMIC"=>1,"COSMIC_CENSUS_GENE"=>1,"GNOMAD_AC"=>1,
  #            "GNOMAD_AN"=>1,"NO_GNOMAD"=>1,"CNAME"=>1;};
  my $hcsq  = get_consequence_values();
  my $cd_nc = get_coding_noncoding_values();

  #print SNV join(" ","FILE","CHROM","POS","SNVS","BCSQ","CODING","COSMIC_CENSUS_GENE","COSMIC","GNOMAD","CENTROMER","MPOS","DP","GERMQ","SEQQ","STRANDQ","TLOD","AF","ADR","ADA","OR1","OR2","OA1","OA2")."\n";
  if($rel){
   if($skyp){	
  	print SNV join(" ","FILE","CHROM","POS","SNVS","BCSQ","CODING","COSMIC_CENSUS_GENE","COSMIC","GNOMAD","CENTROMER","MPOS","DP","TLOD","AF","ADR","ADA","OR1","OR2","OA1","OA2")."\n";
  	print INDEL join(" ","FILE","CHROM","POS","SNVS","BCSQ","CODING","COSMIC_CENSUS_GENE","COSMIC","GNOMAD","CENTROMER","MPOS","DP","TLOD","AF","ADR","ADA","OR1","OR2","OA1","OA2")."\n";
  }else{
  	print SNV join(" ","FILE","CHROM","POS","SNVS","BCSQ","CODING","COSMIC_CENSUS_GENE","COSMIC","GNOMAD","CENTROMER","MPOS","DP","GERMQ","TLOD","AF","ADR","ADA","OR1","OR2","OA1","OA2")."\n";
  	print INDEL join(" ","FILE","CHROM","POS","SNVS","BCSQ","CODING","COSMIC_CENSUS_GENE","COSMIC","GNOMAD","CENTROMER","MPOS","DP","GERMQ","TLOD","AF","ADR","ADA","OR1","OR2","OA1","OA2")."\n";
  }	
  }else{
  #rel2 t-only files
  if($skyp){	
  print SNV join(" ","FILE","CHROM","POS","SNVS","BCSQ","CODING","COSMIC_CENSUS_GENE","COSMIC","GNOMAD","CENTROMER","MPOS","DP","SEQQ","STRANDQ","TLOD","AF","ADR","ADA","OR1","OR2","OA1","OA2")."\n";
  print INDEL join(" ","FILE","CHROM","POS","SNVS","BCSQ","CODING","COSMIC_CENSUS_GENE","COSMIC","GNOMAD","CENTROMER","MPOS","DP","SEQQ","STRANDQ","TLOD","AF","ADR","ADA","OR1","OR2","OA1","OA2")."\n";
  }else{
  print SNV join(" ","FILE","CHROM","POS","SNVS","BCSQ","CODING","COSMIC_CENSUS_GENE","COSMIC","GNOMAD","CENTROMER","MPOS","DP","GERMQ","SEQQ","STRANDQ","TLOD","AF","ADR","ADA","OR1","OR2","OA1","OA2")."\n";
  print INDEL join(" ","FILE","CHROM","POS","SNVS","BCSQ","CODING","COSMIC_CENSUS_GENE","COSMIC","GNOMAD","CENTROMER","MPOS","DP","GERMQ","SEQQ","STRANDQ","TLOD","AF","ADR","ADA","OR1","OR2","OA1","OA2")."\n";
  }
 }

    #forward
    #AC AG AT CA CG CT GA GC GT TA TC TG
    #revcomp
    #TG TC TA GT GC GA CT CG CA AT AG AC
    #pairs
    #AC TG => AC 1
    #AG TC => AG 2
    #AT TA => AT 3
    #CA GT => CA 4
    #CG GC => CG 5
    #CT GA => CT 6
    #mirror pairs
    #GA CT 6
    #GC CG 5
    #GT CA 4
    #TA AT 3
    #TC AG 2
    #TG AC 1
    #signatures for SNPs
  my $sig={"AC"=>"AC", "TG"=>"AC",
           "AG"=>"AG","TC"=>"AG",
           "AT"=>"AT", "TA"=>"AT",
           "CA"=>"CA", "GT"=>"CA",
           "CG"=>"CG", "GC"=>"CG",
           "CT"=>"CT", "GA"=>"CT"};

   while (my $line=<VCF>){
          next if($line=~m/^#/);
          my $tags=$self->_get_tags($line);
          my @cols=();
          my @data=split("\t",$line);
          #variant effect
          if(defined $tags->{BCSQ}){
              push(@cols,annot_consequence($hcsq,$tags->{BCSQ}));
              #annot if the var is coding or not
              push(@cols,annot_coding($cd_nc,$tags->{BCSQ}));

          }else{
            push(@cols,"MODIFIER");
            #non coding variant
            push(@cols,0);
          }

          #variant in cosmic gene
          if(defined $tags->{COSMIC_CENSUS_GENE}){
                push(@cols,1);
          }else{
            push(@cols,0);
          }

          #variant annotated in cosmic
          if(defined $tags->{COSMIC}){
                push(@cols,1);
          }else{
            push(@cols,0);
          }

          #variant in gnomad
          if(defined $tags->{NO_GNOMAD}){
            push(@cols,0);
          }else{
            push(@cols,$tags->{GNOMAD_AC});
          }

          #variant in centromer
          if(defined $tags->{CENTROMER}){
            push(@cols,1);
          }else{
            push(@cols,0);
          }
          #quality metrics mutec2
          #MESO_093_filtered_PASS_norm.tonly.csq.vcf.bgz   MBQ     Integer 24645   median base quality
          #MESO_093_filtered_PASS_norm.tonly.csq.vcf.bgz   MFRL    Integer 24645   median fragment length
          #MESO_093_filtered_PASS_norm.tonly.csq.vcf.bgz   MMQ     Integer 24645   median mapping quality
          #MESO_093_filtered_PASS_norm.tonly.csq.vcf.bgz   MPOS    Integer 24645   median distance from end of read
          #push(@cols,$tags->{MBQ});
          #push(@cols,$tags->{MFRL});
          #push(@cols,$tags->{MMQ});
          push(@cols,$tags->{MPOS});
          push(@cols,$tags->{DP});
	  if($skyp){
          #push(@cols,$tags->{GERMQ});
		#tumor raw dont have GERMQ 
	  }else{
          push(@cols,$tags->{GERMQ});
	  }
	  if($rel){
		#rel3 dont have this values
	  }else{
          push(@cols,$tags->{SEQQ});
          push(@cols,$tags->{STRANDQ});
	  }
          push(@cols,int($tags->{TLOD}));
          #genotype variables Allelic fraction,  ref and alt depths
          push(@cols,$tags->{AF_1});
          push(@cols,$tags->{AD_1});
          push(@cols,$tags->{AD_2});
          #read orientation supporting the calling
          #Count of read pairs in the F1R2 and F2R1 configurations supporting REF and ALT alleles (F1R2, F2R1)
          #reference
          push(@cols,$tags->{F1R2_1});
          push(@cols,$tags->{F1R2_2});
          #alt
          push(@cols,$tags->{F2R1_1});
          push(@cols,$tags->{F2R1_2});
          #foreach my $f (@fs){
          #  if(defined $tags->{BCSQ}){
          #      annot_consequence($hcsq,$tags->{BCSQ});
          #  }
          #}
          #print Dumper($tags);
          #SNP variant
          if(length($data[3]) == length($data[4]) and length($data[4]) == 1 ){
              #print SNV join(" ",$self->{vcffile},$data[0],$data[1],$data[3].$data[4],@cols)."\n";
              print SNV join(" ",$self->{vcffile},$data[0],$data[1],$sig->{$data[3].$data[4]},@cols)."\n";
        }else{
          #indel variant
               print INDEL join(" ",$self->{vcffile},$data[0],$data[1],"NN",@cols)."\n";
        }
   }

   close (SNV);
   close (INDEL);
   close(VCF);
}

sub annot_coding{
     my $hvals=shift;
     my $vars=shift;
     my $rank={"NC"=>0,"CD"=>1};
     #print Dumper($hvals);
     #print Dumper($vars);
     my @ip=split(",",$vars);
     my $init="NC";
     foreach my $ii (@ip){
         my ($v)= split /\|/,$ii;
           if($v =~m/&/){
            ($v)=split("&",$v);
          }
          $v =~ s/\W//g;
         if(defined $hvals->{$v}){
           if($rank->{$hvals->{$v}} > $rank->{$init}){
              $init=$hvals->{$v};
           }
         }else{
           print STDERR "$v not defined in impact hash\n";
           print STDERR "$v is equal to BCSQ=\@6415 is a bug in MT\n";
           #print STDERR Dumper($vars);
         }
     }
     #print $init." ".join(" ",@ip)."\n";
     return $init eq "NC" ? 0:1;
}


sub annot_consequence{
     my $hvals=shift;
     my $vars=shift;
     my $rank={"MODIFIER"=>0,"LOW"=>1,"MODERATE"=>2,"HIGH"=>3};
     #print Dumper($hvals);
     #print Dumper($vars);
     my @ip=split(",",$vars);
     my $init="MODIFIER";
     foreach my $ii (@ip){
         my ($v)= split /\|/,$ii;
           if($v =~m/&/){
            ($v)=split("&",$v);
          }
          $v =~ s/\W//g;
         if(defined $hvals->{$v}){
           if($rank->{$hvals->{$v}} > $rank->{$init}){
              $init=$hvals->{$v};
           }
         }else{
           print STDERR "$v not defined in impact hash\n";
           print STDERR "$v is equal to BCSQ=\@6415 is a bug in MT\n";
           #print STDERR Dumper($vars);
         }
     }
     #print $init." ".join(" ",@ip)."\n";
     return $init;
}


sub get_consequence_values{
  my $vals={"intergenic"=>"MODIFIER",
            "downstream"=>"MODIFIER",
            "upstream"=>"MODIFIER",
            "intron"=>"MODIFIER",
            "non_coding"=>"MODIFIER",
            "regulatory"=>"MODIFIER",
            "5_prime_utr"=>"MODIFIER",
            "3_prime_utr"=>"MODIFIER",
            "stop_retained"=>"LOW",
            "start_retained"=>"LOW",
            "synonymous"=>"LOW",
            "splice_region"=>"LOW",
            "coding_sequence"=>"MODIFIER",
            "missense"=>"MODERATE",
            "inframe"=>"MODERATE",
            "exon_loss"=>"HIGH",
            "disruptive"=>"HIGH",
            "splice_acceptor"=>"HIGH",
            "splice_donor"=>"HIGH",
            "start_lost"=>"HIGH",
            "stop_lost"=>"HIGH",
            "stop_gained"=>"HIGH",
            "frameshift"=>"HIGH",
            "inframe_insertion"=>"MODERATE",
            "inframe_deletion"=>"MODERATE",
          };
          return $vals;
}

sub get_coding_noncoding_values{
  #CD = coding variant
  #NC = non coding variant
  my $vals={"intergenic"=>"NC",
            "downstream"=>"NC",
            "upstream"=>"NC",
            "intron"=>"NC",
            "non_coding"=>"NC",
            "regulatory"=>"NC",
            "5_prime_utr"=>"NC",
            "3_prime_utr"=>"NC",
            "stop_retained"=>"NC",
            "start_retained"=>"NC",
            "synonymous"=>"CD",
            "splice_region"=>"CD",
            "coding_sequence"=>"CD",
            "missense"=>"CD",
            "inframe"=>"CD",
            "exon_loss"=>"CD",
            "disruptive"=>"CD",
            "splice_acceptor"=>"CD",
            "splice_donor"=>"CD",
            "start_lost"=>"CD",
            "stop_lost"=>"CD",
            "stop_gained"=>"CD",
            "frameshift"=>"CD",
            "inframe_insertion"=>"CD",
            "inframe_deletion"=>"CD",
          };
          return $vals;
}


1; #EOM
