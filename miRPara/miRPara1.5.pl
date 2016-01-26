#! /usr/bin/perl -w
#mirpara.pl
#################################USE############################################
use warnings;
use strict;
use Getopt::Long;
use Algorithm::SVM;
use Algorithm::SVM::DataSet;
use File::chdir;#to change the perl environment
##################################GOCOBLE#######################################
our $abname="noname";
our $species="overall";
our $level=15;
our $share="/home/weasteam/biology/WIV/Simon_Rayner/RNAi/datas/program/miRPara/test/";#share data,,add"/" at last
our $prilength=60;
our %rstcount;
our %mirbase;
our $pmtonly="";
our $svmmodel;
our $pri="";
my %options;
my $path;
my $infile="";
my @splitname;
my @basedata;
my @base;
my $m;
################################PROGROM START###################################
################################################################################
################################################################################
GetOptions \%options, 'version|v' => sub { version('miRPara.pl') },
						'help|h' => sub { usage() },
						'name|n=s' => \$abname,
						'species|s=s' => \$species,
						'level|l=i'=>\$level,
						'pri'=>\$pri,
						'pmt'=>\$pmtonly,
						'prilengthlimit|p=i' => \$prilength or die("Type miRPara.pl -h for help!\n");
#################################data prepare###################################
##infile
$infile=shift or die("Error: data not specified\nRun 'miRPara.pl -h' for help\n");
##enviroment
$path=$infile;
($infile)= ($infile =~ m!([^/]+)$!);#to generate the final file name
$path=~ s/$infile$//g;
if ($path ne ""){
   $CWD = $path;#change the envirioment directory
}
##species
if (index("overallanimalplantvirus",$species) eq "-1"){#to check the species
   print "Please provide right species name!!\n";
   die;
}
##count default
$rstcount{'pri'}=0;
$rstcount{'area'}=0;
$rstcount{'candidate'}=0;
$rstcount{'preserial'}=0;
$rstcount{'prepriid'}="";
if ($pmtonly ne 1){
   #svm data
   $svmmodel = new Algorithm::SVM(Model => $share."models/$species\_$level.model");
   ##blast data
   open (IN,$share."blast.dat");#read the blast files
   @basedata=<IN>;
   close IN;
   foreach (@basedata){
	  $_=~ s/\n//g;
	  if ($_ ne ""){
		 @base=split("\t",$_);
		 $mirbase{lc($base[0])}="$base[1]";
	  }
   }
   @basedata=();
   @base=();
}
##abname
if ($abname eq "noname"){#reset the abname
	$abname=$infile;
	$abname=~s/.\w+$//;
}
##extra parameter
if ($pri eq "1"){
   system "cp $infile $abname\_hairpin.fasta";
   goto step4;
}
##file format
if (index($infile,".pmt") ne -1){
   pmtpredict($infile);
   goto step5;
}
##################################steps#########################################
##################STEP1##################
step1:
@splitname=seqsplit($infile,6000,500);
##################STEP2##################
step2:
callunafd('long');#call unafold
##################STEP3##################
step3:
readlongfile(@splitname);#readlongfile
##################step4##################
step4:
if ($rstcount{'pri'} ne 0 or $pri eq "1"){
	candidate();#followed by parameter and pmt2svmdata
}
##################step5##################
step5:
if ($pmtonly ne 1){
   head();
}
if (-e "tmp*"){
   system "rm tmp*";#delete all tmp files
}
print "\nAll Done.\n";
##################STOP HERE##################
sub usage () {
    print <<EOF;
Usage: miRPara.pl [Options] file [File]

Options:
-v, --version
-h, --help
-n, --name=<abbreviated name of the species>
-s, --species=<overall, animal, plant or virus> (defaults as overall)
-l, --Level=<1..10>(defaults to 7)
-p, --prilengthlimit=<limit to the pri-miRNA length> (defaults to 60)
--pri, --Consider all input sequences as pri-miRNAs
--pmt, --Only calculate the parameters without further prediction

File:
*.fasta, --Only fasta format sequences were allowed
*.pmt, --Predict from the parameter files

EOF
print 'Report bugs to weasteam@gmail.com or raynere@wh.iov.cn', "\n";
    exit;
}
sub version ($) {
	print "\n";
    print "$_[0] (miRPare) 1.5\n";
    print "By Yonggan Wu and Simon Rayner\n";
    print "Copyright (C) 2009\n";
    print "Wuhan Institute of Virology, Chinese Academic of Science.\n";
    print "Wuhan, 430071, China\n\n";
    exit;
}
sub seqsplit{#split sequence into small sequences
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The seqsplit was used to split the given sequence into small ones
   #===========Usage: seqsplit($infile,$abname,splitlength,overlap)
   #===========The $abname will be return to the main program.
   #===========The *_splitseq.nm will be created with the created seq (fasta format)
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my ($infile,$slen,$overlap)=@_;#input the file
   my @data;
   my $inseq;
   my $splitseq;
   my $olseqs;#start
   my $olseqe;#end
   my $midseq;
   my $splitname;
   my @splitname;
   my $inlen;
   my $start=1;
   my $end;
   my @seq;
   my @tmp;
   @seq=readfas($infile);
   open (OUT,">$abname\_splitseq.fasta");
   print OUT "";
   close OUT;
   foreach (@seq){
	  @tmp=split("\t",$_);
	  $inseq=$tmp[1];
	  $inlen=length($inseq);
	  while ($inlen>$slen){
		 $end=$start+$slen-1;
		 $splitname="$tmp[0]\_$start-$end";
		 $olseqs=substr($inseq,0,$overlap);
		 $midseq=substr($inseq,$overlap,$slen-2*$overlap);
		 $olseqe=substr($inseq,$slen-$overlap,$overlap);
		 $inseq=substr($inseq,$slen-$overlap);
		 open (OUT,">>$abname\_splitseq.fasta");
		 print OUT ">$splitname\n";
		 print OUT "$olseqs$midseq$olseqe\n";
		 close OUT;
		 $start+=length($olseqs)+length($midseq);
		 $inlen=length($inseq);
		 push(@splitname,$splitname);
	  }
	  $end=$start+$inlen-1;
	  $splitname="$tmp[0]\_$start-$end";
	  open (OUT,">>$abname\_splitseq.fasta");
	  print OUT ">$splitname\n";
	  print OUT "$inseq\n";
	  close OUT;
	  push(@splitname,$splitname);
   }
   return @splitname;
}
sub callunafd{
	#>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	#===========The callunafd was used to call unafold.pl to generate the second
	#			structure and minimal free energy.
	#===========Usage: seqsplit($abname,parameter)
	#			-long for the long file *_splitseq.nm
	#			-tmp for tmp files
	#===========The the ct2out will be called at the same time to generate the *.out
	#===========The *.str file will be output and other files from unafold.pl will be deleted
	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my $parameter;
	my @seq;
	my @tmp;
    ($parameter)=@_;#-l for the long file, and -s for the small file
    if ($parameter eq "long"){
	  @seq=readfas("$abname\_splitseq.fasta");
	  foreach (@seq){
		 @tmp=split("\t",$_);
		 print "RNA folding with UNAFold to $tmp[0]\n";
		 open (OUT,">tmp.fasta");
		 print OUT ">$tmp[0]\n";
		 print OUT $tmp[1];
		 close OUT;
		 system "UNAFold.pl tmp.fasta > error.log";
		 system "ct2out <tmp.fasta_1.ct> $tmp[0]\.str";#run ct2out
		 system "rm tmp.fas*";
	  }
    }
    if ($parameter eq "tmp"){
	  if (-e "tmp.str"){
		 system "rm tmp.str";
	  }
	  system "UNAFold.pl tmp.fasta > error.log";
	  system "ct2out <tmp.fasta_1.ct> tmp.str";#run ct2out
	  system "rm tmp.fas*";
    }
}
sub readlongfile{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The readlongfile was used to get all hairpin seq from *.out
   #===========Usage: readlongfile($abname,$limit)
   #			--limit the length limit of prilength
   #===========The *_splitseq.nm is need for right information
   #===========The sequence will be output as fasta format
   #			and a *_hairpin.nm will also be output
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #format readlongfile(length)
   my @splitname;
   my $serial;#the number of hairpins
   my @data;#the data of the *.out file
   my $m;#tmp number
   my $n;#tmp number
   my $string;#tmp string
   my %hairpin;#hairpin structure
   my $rp=0;#whether a sequence from repeat or overlap reagion
   my $seq;#each small sequence
   my @sequence;#collection of all seq, to avoid replicate
   my $do=1;
   my $mfe;
   my $splitname;
   (@splitname)=@_;
   open (NAME,">$abname\_hairpin.fasta");#clean the fasta file
   print NAME "";
   close NAME;
   foreach (@splitname){
	  $serial=1;
	  $splitname=$_;
	  print "Read the secondary structure file from $_\n";
	  open (DATA,"$_\.str") or open (DATA,"$_\.out") or open (DATA,"$_\.fas.out");#read the data of *.out file
	  @data=<DATA>;
	  close DATA;
	  $do=0;
	  while ($do <= 1){
		($m,$mfe,@data)=strtrim(0,@data);
		$do=$do+$m;
		if ($m ne 1){
			$hairpin{1}=lc($data[0]);#add data to each hairpin segment
			$hairpin{2}=lc($data[1]);
			$hairpin{3}=lc($data[2]);
			$hairpin{4}=lc($data[3]);
			$hairpin{1}=~s/t/u/g;
			$hairpin{2}=~s/t/u/g;
			$hairpin{3}=~s/t/u/g;
			$hairpin{4}=~s/t/u/g;
			$hairpin{4}=~s/^\s+//;
			($hairpin{4})= ($hairpin{4} =~ m!([^\\]+)$!);#to get the sequence without \
			$n=length($hairpin{3})-length($hairpin{4});
			$hairpin{1}=substr($hairpin{1},$n);
			$hairpin{2}=substr($hairpin{2},$n);
			$hairpin{3}=substr($hairpin{3},$n);
			#length decide
			$string=$hairpin{1}.$hairpin{2}.$hairpin{3}.$hairpin{4};
			$string=~ s/\s+//g;#get rid of black
			$string=~ s/\\+//g;#get rid of \
			$string=~ s/-+//g;#get rid of -
			$string=~ s/\n+//g;#get rid of -
			if (length($string)>=$prilength){#output the sequences
				$rstcount{'pri'}+=1;
			   $seq=hairpin2seq($hairpin{1},$hairpin{2},$hairpin{3},$hairpin{4});
			   foreach (@sequence){#whether it is repeat to provious seq
				  if ($seq eq $_){
					 $rp=1
				  }
			   }
			   if ($rp eq 0){
				  open (HAIRPIN,">>$abname\_hairpin.fasta");
				  print HAIRPIN ">$splitname-mir-$serial","\n";
				  print HAIRPIN $seq,"\n";
				  close HAIRPIN;
				  push (@sequence,$seq);#add the seq to the pool for replication detecting
				  $serial+=1;
			   }
			   $rp=0;
			}
			shift(@data);
			shift(@data);
			shift(@data);
			shift(@data);
		}
	  }
   }
}
sub hairpin2seq {#conver the second structure into a line
	#>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	#===========The hairpin2seq was used to generate a seq from second str
	#===========Usage: readlongfile(seq1,seq2,seq3,seq4)
	#===========A line seq will be resurn.
	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	#formathairpin2seq(seq1,seq2,seq3,seq4)
    my $s1;
    my $s2;
    my $s3;
    my $s4;
    my $seq;#for the sequence
    my $upperseq;#for the upper strand
    my $lowerseq;#for the lower strand
    ($s1,$s2,$s3,$s4)=@_;
	$s1=~ s/\n//g;
	$s2=~ s/\n//g;
	$s3=~ s/\n//g;
	$s4=~ s/\n//g;
    my $m=0;
    while ($m<=(length($s2)-1)){
		$upperseq.=substr($s1,$m,1).substr($s2,$m,1);#get the seq
		$lowerseq.=substr($s3,-($m+1),1).substr($s4,-($m+1),1);
		$m+=1;
    }
    $seq=$upperseq.$lowerseq;
    $seq=~ s/\\//;#get rid of \\
    $seq=~ s/\s//g;#get rid of black
    $seq=~ s/-//g;#get rid of -
    $seq=~ s/\n//g;#get rid of \n,if any
	$upperseq="";
	$lowerseq="";
    return $seq;
}
sub candidate{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The candidate was used to generate the miRNA candidate of a pri-miRNA
   #===========Usage: candidate($abname,$species)
   #===========*_hairpin.nm and *.fas will be used
   #===========Some primary parameter were generated and then put into <parameter>
   #			to calculate more other parameters.
   #			*_candi.nm will be output too.
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my %para;#hash for the parameters of a miRNA
#   my $id="";#the id of the sequence
   my @strdata;#the data of *.out files
   my $seq;
   my $s;#tmp string
   my $m;#tmp number
   my $n;#tmp number
   my $seq5;#the line seq of 5' strnad
   my $seq3;#the line seq of 3' strnad
   my %pstr5;#the position of a nt in its second structure
   my %pstr3;#the position of a nt in its second structure
   my $serial=1;#the serial of each miRNA
   my $position;#the position of each nt in liner sequence
   my @pmt;
   my $tmp;
   my @tmp;
   my @seq;
   open (DATAOUT, ">$abname\_candidate.pmt");#tab separate files
   $m=1;
   @pmt=pmt("all");
   foreach (@pmt){
	  print DATAOUT "[$m]$_\t";
	  $m+=1;
   	}
   print DATAOUT "\n";
   close DATAOUT;
   @pmt="";
   if ($pmtonly ne 1){
	  open (OT,">$abname\_mir_overview_$level.out");
	  print OT "";
	  close OT;
	  open (OT,">$abname\_mir_parameter_$level.out");
	  print OT "";
	  close OT;
   }
   $rstcount{'pri'}=0;
   @seq=readfas("$abname\_hairpin.fasta");
   foreach (@seq){
	  @tmp=split("\t",$_);
	  $para{'priid'}=$tmp[0];
	  $seq=$tmp[1];
#	  $rstcount{'pri'}+=1;
	  open (OUT,">tmp.fasta");
	  print OUT ">tmp\n";
	  print OUT $seq;
	  close OUT;
	  #unafold
	  system "UNAFold.pl tmp.fasta > error.log";
	  system "ct2out <tmp.fasta_1.ct> tmp.str";#run ct2out
	  system "rm tmp.fas*";
#	  callunafd("tmp");
	  open (DATA,"tmp.str");#read the name of the long sequences
	  @strdata=<DATA>;
	  close DATA;
	  ($m,$para{'primfe'},@strdata)=strtrim(1,@strdata);#read the minimal free energy
	  system "rm tmp.str";
	  $para{'pristr_1'}=lc($strdata[0]);#read the second structure
	  $para{'pristr_2'}=lc($strdata[1]);
	  $para{'pristr_3'}=lc($strdata[2]);
	  $para{'pristr_4'}=lc($strdata[3]);
	  $para{'pristr_1'}=~ s/\n//g;
	  $para{'pristr_2'}=~ s/\n//g;
	  $para{'pristr_3'}=~ s/\n//g;
	  $para{'pristr_4'}=~ s/\n//g;
	  $para{'priseq'} = lc($seq);#get the sequence
	  #judge if there is any budding stem
	  $s=join("",$para{'pristr_1'},$para{'pristr_2'},$para{'pristr_3'},$para{'pristr_4'});
	  $s=~ s/[\s\\\.-]+//g;
	  if (length($s) ne length($para{'priseq'})){
		 goto line1;#budding stem
	  }
	  else {#no budding stem
		 $para{'buddingstem'}="NO";#budding stem
		 #generate the line seq of each strand
		 $s=join("",$para{'pristr_1'},$para{'pristr_2'});
		 $s=~ s/[\s\\\.-]+//g;
		 $seq5=substr($para{'priseq'},0,length($s));#get the sequence of 5'strand
		 $seq3=substr($para{'priseq'},length($s));#get the sequence of 3'strand
		 #generate the position of each nt in the second str.
		 $m=0;
		 $n=1;
		 while ($n<=length($seq5)){#get position of each nuclitide in second structure
			if ((substr($para{'pristr_1'},$m,1) ne "-") and (substr($para{'pristr_2'},$m,1) ne "-")){
			   $pstr5{$n}=$m+1;
			   $n +=1;
			}
			$m +=1;
		 }
		 $pstr5{0}=0;
		 $m=0;
		 $n=1;
		 while ($n<=length($seq3)){#get position of each nuclitide in second structure
			if ((substr($para{'pristr_3'},$m,1) ne "-") and (substr($para{'pristr_4'},$m,1) ne "-")){
			   $pstr3{$n}=$m+1;
			   $n +=1;
			}
			$m +=1;
		 }
		 $pstr3{0}=0;
		 #generate the candidate with 20-24 in length
		 #goto line;
		 $position=1;
		 $tmp=$para{'pristr_1'};
		 $tmp=~s/\s+/ /g;
		 $tmp=~s/ $//g;
		 @tmp=split (" ",$tmp);
		 $tmp=pop(@tmp);
		 my $loop=length($tmp)+1;
		 while ((length($seq5)-$position-$loop)>=19){
			$m=20;
			while (($m<=24) and ($position+$m)<=(length($seq5)-$loop)){
			   $para{'miid'}="$abname-MIR-$serial\_$m";
			   #upper position
			   $para{'upperstart'}=$pstr5{$position};#the start position of a miRNA in second str
			   $para{'upperend'}=$pstr5{$position+$m-1};#the start position of a miRNA in second str
			   #mi position
			   $para{'mistart'}=$position;
			   $para{'miend'}=$position+$m-1;
			   #the second structure of the pri-miRNA
			   $para{'pristr_3'}=lc($para{'pristr_3'});
			   $para{'pristr_4'}=lc($para{'pristr_4'});
			   $para{'pristr_1'}=lc(substr($para{'pristr_1'},0,$para{'upperstart'}-1)).uc(substr($para{'pristr_1'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1)).lc(substr($para{'pristr_1'},$para{'upperend'}));
			   $para{'pristr_2'}=lc(substr($para{'pristr_2'},0,$para{'upperstart'}-1)).uc(substr($para{'pristr_2'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1)).lc(substr($para{'pristr_2'},$para{'upperend'}));
			   #the second structure of pre-miRNA
			   $para{'prestr_3'}=lc(substr($para{'pristr_3'},$para{'upperstart'}-1));
			   $para{'prestr_4'}=lc(substr($para{'pristr_4'},$para{'upperstart'}-1));
			   $para{'prestr_1'}=uc(substr($para{'pristr_1'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1)).lc(substr($para{'pristr_1'},$para{'upperend'}));
			   $para{'prestr_2'}=uc(substr($para{'pristr_2'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1)).lc(substr($para{'pristr_2'},$para{'upperend'}));
			   #the second structure of miRNA
			   $para{'mistr_3'}=lc(substr($para{'pristr_3'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1));
			   $para{'mistr_4'}=lc(substr($para{'pristr_4'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1));
			   $para{'mistr_1'}=uc(substr($para{'pristr_1'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1));
			   $para{'mistr_2'}=uc(substr($para{'pristr_2'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1));
			   #length
			   $s=join("",$para{'prestr_1'},$para{'prestr_2'},$para{'prestr_3'},$para{'prestr_4'});
			   $s=~ s/[\s\\-]+//g;
			   $para{'prelength'}=length($s);
			   #seq
			   $para{'preseq'}=substr($para{'priseq'},$position-1,$para{'prelength'});
			   $para{'miseq'}=substr($seq5,$position-1,$m);
			   #strand
			   $para{'strand'}='5';
			   print "Calculating miRNA candidate: $para{'miid'} at $para{'priid'}";
			   parameter($serial,%para);#sent the information to calculate the parameters
			   $m +=1;
			}#20-24
			$serial +=1;
			$position +=1;
		 }#position
		 $serial -=1;
		 #line:
		 $position=length($seq3)-$loop;
		 while ($position>=20){
			$m=20;
			while (($m<=24) and ($position >=$m)){
			   $para{'miid'}="$abname-MIR-$serial\_$m";
			   #upper position
			   $para{'upperstart'}=$pstr3{$position-$m+1};#the start position of a miRNA in second str
			   $para{'upperend'}=$pstr3{$position};#the start position of a miRNA in second str
			   #mi position
			   $para{'mistart'}=length($seq5)+length($seq3)-$position+1;
			   $para{'miend'}=length($seq5)+length($seq3)-$position+1+$m-1;
			   #the second structure of the pri-miRNA
			   $para{'pristr_1'}=lc($para{'pristr_1'});
			   $para{'pristr_2'}=lc($para{'pristr_2'});
			   $para{'pristr_3'}=lc(substr($para{'pristr_3'},0,$para{'upperstart'}-1)).uc(substr($para{'pristr_3'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1)).lc(substr($para{'pristr_3'},$para{'upperend'}));
			   $para{'pristr_4'}=lc(substr($para{'pristr_4'},0,$para{'upperstart'}-1)).uc(substr($para{'pristr_4'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1)).lc(substr($para{'pristr_4'},$para{'upperend'}));
			   #the second structure of pre-miRNA
			   $para{'prestr_1'}=lc(substr($para{'pristr_1'},$para{'upperstart'}-1));
			   $para{'prestr_2'}=lc(substr($para{'pristr_2'},$para{'upperstart'}-1));
			   $para{'prestr_3'}=uc(substr($para{'pristr_3'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1)).lc(substr($para{'pristr_3'},$para{'upperend'}));
			   $para{'prestr_4'}=uc(substr($para{'pristr_4'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1)).lc(substr($para{'pristr_4'},$para{'upperend'}));
			   #the second structure of miRNA
			   $para{'mistr_1'}=lc(substr($para{'pristr_1'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1));
			   $para{'mistr_2'}=lc(substr($para{'pristr_2'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1));
			   $para{'mistr_3'}=uc(substr($para{'pristr_3'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1));
			   $para{'mistr_4'}=uc(substr($para{'pristr_4'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1));
			   #length
			   $s=join("",$para{'prestr_1'},$para{'prestr_2'},$para{'prestr_3'},$para{'prestr_4'});
			   $s=~ s/[\s\\-]+//g;
			   $para{'prelength'}=length($s);
			   #seq
			   $para{'preseq'}=substr($para{'priseq'},$m-$position-$para{'prelength'},$para{'prelength'});
			   $para{'miseq'}=substr($seq3,-$position,$m);
			   #strand
			   $para{'strand'}='3';
			   print "calulating miRNA candidate: $para{'miid'} at $para{'priid'}";
			   parameter($serial,%para);#sent the information to calculate the parameters
			   $m +=1;
		   }
		   $serial +=1;
		   $position-=1;
		 }
	  }#no budding stem
	  line1:
   }#each seq
}#sub
sub parameter{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The parameter was used to generate the miRNA parameter values
   #===========Usage: parameter($abname,$species,%para)
   #===========Some other sub programs will be called:
   #			-gc
   #			-pairs
   #			-gu
   #			-ntcontent
   #			-internalloop
   #			-internalloopnumber
   #			-unpairedbases
   #			-unpairedrate
   #===========The *.pmt file will be created.
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my %para;#a hash for all parameters
   my @strdata;#the data of *.out file
   my $s;#tmp string
   my $s1;#tmp string
   my $s2;#tmp string
   my $m;#tmp number
   my $n;#tmp number
   my @pmt;
   my $pmt;
   my $rst;#the predict result;
   my $serial;
   my $tmp;
   my @tmp;
   ($serial,%para)=@_;#receive the parameter
   #length
   $para{'prilength'}=length($para{'priseq'});
   $para{'prelength'}=length($para{'preseq'});
   $para{'milength'}=length($para{'miseq'});
   #premfe
#   open (OUTPUT,">tmp.fasta");#output the pre seq
#   print OUTPUT ">tmp\n";
#   print OUTPUT $para{'preseq'};
#   close OUTPUT;
#   if (-e "tmp.str"){
#	  system "rm tmp.str";
#   }
#   system "UNAFold.pl tmp.fasta > error.log";
#   system "ct2out <tmp.fasta_1.ct> tmp.str";#run ct2out
#   system "rm tmp.fas*";
#   callunafd('tmp');#call unafold to generate the mfe
#   open (DATA,"tmp.str");
#   @strdata=<DATA>;
#   close DATA;
#   ($m,$para{'premfe'},@strdata)=strtrim(1,@strdata);
   #gc
   $para{'premfe'}=sprintf '%4.2f',-1*(($para{'primfe'}*(-1))*($para{'prelength'}/$para{'prilength'}));
   $para{'prigc'}=gc($para{'priseq'});#prigc
   $para{'pregc'}=gc($para{'preseq'});#pregc
   $para{'migc'}=gc($para{'miseq'});#migc
   #pairs
   $para{'pripairs'}=pairs($para{'pristr_2'},'yes');
   $para{'prepairs'}=pairs($para{'prestr_2'},'yes');
   $para{'mipairs'}=pairs($para{'mistr_2'});
   #gu
   $para{'prigu'}=gu($para{'pristr_2'},$para{'pristr_3'},'yes');
   $para{'pregu'}=gu($para{'prestr_2'},$para{'prestr_3'},'yes');
   $para{'migu'}=gu($para{'mistr_2'},$para{'mistr_3'});
   #ntcontent
   $para{'printcontent_a'}=ntcontent($para{'priseq'},'a');
   $para{'printcontent_u'}=ntcontent($para{'priseq'},'u');
   $para{'printcontent_c'}=ntcontent($para{'priseq'},'c');
   $para{'printcontent_g'}=ntcontent($para{'priseq'},'g');
   $para{'prentcontent_a'}=ntcontent($para{'preseq'},'a');
   $para{'prentcontent_u'}=ntcontent($para{'preseq'},'u');
   $para{'prentcontent_c'}=ntcontent($para{'preseq'},'c');
   $para{'prentcontent_g'}=ntcontent($para{'preseq'},'g');
   $para{'mintcontent_a'}=ntcontent($para{'miseq'},'a');
   $para{'mintcontent_u'}=ntcontent($para{'miseq'},'u');
   $para{'mintcontent_c'}=ntcontent($para{'miseq'},'c');
   $para{'mintcontent_g'}=ntcontent($para{'miseq'},'g');
   #firstbase
   $para{'firstbase'}=uc(substr($para{'miseq'},0,1));
   #BasalSegment
   $s1=$para{'pristr_1'};
   $s1=~ s/[-|\\\^]+//g;#get rid of unnecessary charaters
   $s2=$para{'pristr_4'};
   $s2=~ s/[-|\\\^]+//g;#get rid of unnecessary charaters
   my $basalen1=0;
   my $basalen2=0;
   $basalen1=index($s1," ",0);
   $basalen2=index($s2," ",0);
   if ($basalen1>=$basalen2) {
	  $para{'length_basalsegment'}=$basalen1;
	  $para{'basalend'}=$basalen1;#the end position of basal segment
   }
   else {
	  $para{'length_basalsegment'}=$basalen2;
	  $para{'basalend'}=$basalen2;
   }
	  #if the mature miRNA is locate in the reagin
   if (($para{'length_basalsegment'}+length($para{'prestr_1'}))>length($para{'pristr_1'})){
	  $para{'length_basalsegment'}=0;
   }
   #terminal loop
   
   $tmp=$para{'pristr_1'};
   $tmp=~s/\s+/ /g;
   $tmp=~s/ $//g;
   @tmp=split (" ",$tmp);
   $tmp=pop(@tmp);
   my $terminalloopstem=length($tmp);
   $para{'loopstart'}=length($para{'pristr_1'})-$terminalloopstem;#the start position of the loop
   if (index($para{'pristr_2'},"\\") eq -1) {#no \
	  $para{'length_terminalloop'}=2*$terminalloopstem+2;
   }
   else {
	  $para{'length_terminalloop'}=2*$terminalloopstem+1;
   }
   #lowerstem
   $para{'length_lowerstem'}=$para{'upperstart'}-$para{'basalend'}-1;
   if ($para{'length_lowerstem'}<=0){
	  $para{'length_lowerstem'}="NULL";
   }
   #upperstem
   $para{'length_upperstem'}=$para{'upperend'}-$para{'upperstart'}+1;
   #topstem
   $para{'length_topstem'}=$para{'loopstart'}-$para{'upperend'}-1;
   if ($para{'length_topstem'}<=0){
	  $para{'length_topstem'}="NULL";
   }
   #internalloop
   $para{'priinternalloop'}=internalloop(
	  substr($para{'pristr_1'},$para{'basalend'},$para{'loopstart'}-$para{'basalend'}-1),
	  substr($para{'pristr_4'},$para{'basalend'},$para{'loopstart'}-$para{'basalend'}-1));
   $para{'preinternalloop'}=internalloop(
	  substr($para{'pristr_1'},$para{'upperstart'}-1,$para{'loopstart'}-$para{'upperstart'}),
	  substr($para{'pristr_4'},$para{'upperstart'}-1,$para{'loopstart'}-$para{'upperstart'}));
   $para{'miinternalloop'}=internalloop($para{'mistr_1'},$para{'mistr_4'});
   if ($para{'length_lowerstem'} ne "NULL"){
	  $para{'internalloop_lowerstem'}=internalloop(
		 substr($para{'pristr_1'},$para{'basalend'},$para{'length_lowerstem'}),
		 substr($para{'pristr_4'},$para{'basalend'},$para{'length_lowerstem'}));
   }
   else{
	  $para{'internalloop_lowerstem'}="NULL"
   }
   if ($para{'length_topstem'} ne "NULL"){
	  $para{'internalloop_topstem'}=internalloop(
		 substr($para{'pristr_1'},$para{'upperend'},$para{'length_topstem'}),
		 substr($para{'pristr_4'},$para{'upperend'},$para{'length_topstem'}));
   }
   else{
	  $para{'internalloop_topstem'}="NULL"
   }
   #internalloop number
   $para{'priinternalloopnumber'}=internalloopnumber($para{'pristr_2'})-1;
   $para{'preinternalloopnumber'}=internalloopnumber($para{'prestr_2'})-1;
   $para{'miinternalloopnumber'}=internalloopnumber($para{'mistr_2'});
   if ($para{'length_lowerstem'} ne "NULL"){
	  $para{'internalloopnumber_lowerstem'}=internalloopnumber(
		 substr($para{'pristr_2'},$para{'basalend'},$para{'length_lowerstem'}));
   }
   else{
	  $para{'internalloopnumber_lowerstem'}="NULL"
   }
   if ($para{'length_topstem'} ne "NULL"){
	  $para{'internalloopnumber_topstem'}=internalloopnumber(
		 substr($para{'pristr_2'},$para{'upperend'},$para{'length_topstem'}));
   }
   else{
	  $para{'internalloopnumber_topstem'}="NULL"
   }
   if ($para{'length_basalsegment'}>0){
	  $para{'priinternalloopnumber'}-=1;
	  $para{'preinternalloopnumber'}-=1;
  }
   #unpairedbases
   $para{'priunpairedbases'}=unpairedbases(
	  substr($para{'pristr_1'},$para{'basalend'},$para{'loopstart'}-$para{'basalend'}-1),
	  substr($para{'pristr_4'},$para{'basalend'},$para{'loopstart'}-$para{'basalend'}-1));
  $para{'preunpairedbases'}=unpairedbases(
	  substr($para{'pristr_1'},$para{'upperstart'}-1,$para{'loopstart'}-$para{'upperstart'}),
	  substr($para{'pristr_4'},$para{'upperstart'}-1,$para{'loopstart'}-$para{'upperstart'}));
   $para{'miunpairedbases'}=unpairedbases($para{'mistr_1'},$para{'mistr_4'});
   if ($para{'length_lowerstem'} ne "NULL"){
	  $para{'unpairedbases_lowerstem'}=unpairedbases(
		 substr($para{'pristr_1'},$para{'basalend'},$para{'length_lowerstem'}),
		 substr($para{'pristr_4'},$para{'basalend'},$para{'length_lowerstem'}));
   }
   else{
	  $para{'unpairedbases_lowerstem'}="NULL"
   }
   if ($para{'length_topstem'} ne "NULL"){
	  $para{'unpairedbases_topstem'}=unpairedbases(
		 substr($para{'pristr_1'},$para{'upperend'},$para{'length_topstem'}),
		 substr($para{'pristr_4'},$para{'upperend'},$para{'length_topstem'}));
   }
   else{
	  $para{'unpairedbases_topstem'}="NULL"
   }
   #unpairedrate
   $para{'priunpairedrate'}=unpairedrate($para{'priunpairedbases'},
	  substr($para{'pristr_2'},$para{'basalend'},$para{'loopstart'}-$para{'basalend'}-1),
	  substr($para{'pristr_3'},$para{'basalend'},$para{'loopstart'}-$para{'basalend'}-1));
   $para{'preunpairedrate'}=unpairedrate($para{'preunpairedbases'},
	  substr($para{'pristr_2'},$para{'upperstart'}-1,$para{'loopstart'}-$para{'upperstart'}),
	  substr($para{'pristr_3'},$para{'upperstart'}-1,$para{'loopstart'}-$para{'upperstart'}));
   $para{'miunpairedrate'}=unpairedrate($para{'miunpairedbases'},$para{'mistr_2'},$para{'mistr_3'});
   if ($para{'length_lowerstem'} ne "NULL"){
	  $para{'unpairedrate_lowerstem'}=unpairedrate($para{'unpairedbases_lowerstem'},
		 substr($para{'pristr_2'},$para{'basalend'},$para{'length_lowerstem'}),
		 substr($para{'pristr_3'},$para{'basalend'},$para{'length_lowerstem'}));
   }
   else{
	   $para{'unpairedrate_lowerstem'}="NULL"
   }
   if ($para{'length_topstem'} ne "NULL"){
	   $para{'unpairedrate_topstem'}=unpairedrate($para{'unpairedbases_topstem'},
		 substr($para{'pristr_2'},$para{'upperend'},$para{'length_topstem'}),
		 substr($para{'pristr_3'},$para{'upperend'},$para{'length_topstem'}));
   }
   else{
	   $para{'unpairedrate_topstem'}="NULL"
   }
   #stability
   $s1=lc(join("",substr($para{'mistr_2'},0,4),substr($para{'mistr_3'},0,4)));
   $s2=lc(join("",substr($para{'mistr_2'},-4,4),substr($para{'mistr_3'},-4,4)));
   $s1=~ s/\s+//g;
   $s1=~ s/c/zz/g;
   $s2=~ s/\s+//g;
   $s2=~ s/c/zz/g;
   if ($para{'strand'} eq "5"){
	  if (length($s2) ne 0){
		 $para{'stability'}=sprintf '%.2f',length($s1)/length($s2);
	  }
	  else{
		 $para{'stability'}=-1;
	  }
   }
   else{
	  if (length($s1) ne 0){
		 $para{'stability'}=sprintf '%.2f',length($s2)/length($s1);
	  }
	  else{
		 $para{'stability'}=-1;
	  }
   }
   #overhang
   $s="";
   $m=1;
   while ((length($s) < 2)  and ($para{'upperend'}+$m)<=length($para{'pristr_3'})){
	  if ($para{'strand'} eq "5"){#generate the two nt
		 $s=join("",substr($para{'pristr_3'},$para{'upperend'},$m),
				 substr($para{'pristr_4'},$para{'upperend'},$m));
	  }
	  else{
		 $s=join("",substr($para{'pristr_1'},$para{'upperstart'}-$m-1,$m),
				 substr($para{'pristr_2'},$para{'upperstart'}-$m-1,$m));
	  }
	  $s=~ s/[\s\\-]+//g;
	  $m++;
   }
   if ($s eq ""){#no overhang
	  $para{'penultimateposition'}="NULL";
	  $para{'terminalnucleotide'}="NULL";
   }
   if (length($s) eq 1){#one overhang
	  $para{'penultimateposition'}=$s;
	  $para{'terminalnucleotide'}="NULL";
   }
   if (length($s) eq 2){#two overhang
	  if ($para{'strand'} eq "5"){
		 $para{'penultimateposition'}=substr($s,0,1);
		 $para{'terminalnucleotide'}=substr($s,1,1);
	  }
	  else{
		 $para{'penultimateposition'}=substr($s,1,1);
		 $para{'terminalnucleotide'}=substr($s,0,1);
	  }
   }
   $para{'penultimateposition'}=uc($para{'penultimateposition'});
   $para{'terminalnucleotide'}=uc($para{'terminalnucleotide'});
   if ($pmtonly ne 1){
	  $rst=predict(%para);
	  if ($rst eq "1"){
		 print "....TRUE\n";
		 $rstcount{'candidate'}+=1;
		 if (($serial-$rstcount{'preserial'})>=2){
			$rstcount{'area'}+=1;
		 }
		 $rstcount{'preserial'}=$serial;
		 out(%para);
	  }
	  else{
		 print "....FALSE\n";
	  }
   }
   else{
	  print "\n";
   }
   #output
   @pmt=pmt("all");
   foreach (@pmt){
	  $pmt.= "$para{$_}\t";
   }
   chop($pmt);
   open (DATAOUT, ">>$abname\_candidate.pmt");
   print DATAOUT "$pmt\n";
   close DATAOUT;
}
sub gc{#parameter
	#>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	#===========The gc was used to generate gc content of given seq
	#===========Usage: gc($seq)
	#===========The GC content will be returned
	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	my $seq;#receive the sequences
	my $s;#tmp string
	my $gc;
	($seq)=@_;
	$s=lc($seq);
	$s=~ s/[aut]+//g;
	$gc=sprintf '%4.4f', length($s)/length($seq);#prigc
	return $gc;
}
sub pairs{#parameter
	#>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	#===========The pairs was used to generate the number of pairs of second str
	#===========Usage: pairs($seq,$loop)
	#			$loop is to test whether there is any nt in the loop
	#			yes/no to decide the $loop value
	#===========The number of pair bands will be returned
	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	my $seq;
	my $loop;
	my $pairs;
	($seq,$loop)=@_;
	if (lc($loop) eq "yes"){
		chop($seq);
	}
	$seq=~ s/\s//g;
	$pairs=length($seq);
	return $pairs;
}
sub gu{#parameter
	#>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	#===========The gu was used to generate the number of GU wobbles
	#===========Usage: gu(seq1,seq2,yes/no)
	#			$loop is to test whether there is any nt in the loop
	#			yes/no to decide the $loop value
	#===========The number of GU wobbles will be returned
	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	my $seq1;
	my $seq2;
	my $s;#tmp string
	my $loop;
	my $gu;
	my $u;
	my $a;
	($seq1,$seq2,$loop)=@_;
	if (lc($loop) eq "yes"){
		chop($seq1);
		chop($seq2);
	}
	$s=lc(join("",$seq1,$seq2));
	$s=~ s/\s+//g;
	$u=$s;
	$u=~ s/[ut]+//g;
	$a=$s;
	$a=~ s/[a]+//g;
	$gu=length($a)-length($u);
	return $gu;
}
sub ntcontent{#parameter
	#>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	#===========The ntcontent was used to generate the nt content of four nts
	#===========Usage: ntcontent($seq,$nt)
	#			$nt reprecent A,U,G OR C
	#===========The nt content of given nt will be returned
	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	my $seq;
	my $s;
	my $nt;
	my $ntcontent;
	($seq,$nt)=@_;
	$s=lc($seq);
	if (lc($nt) eq "a"){
		$s=~ s/[a]+//g;
	}
	if ((lc($nt) eq "t") or (lc($nt) eq "u")){
		$s=~ s/[tu]+//g;
	}
	if (lc($nt) eq "c"){
		$s=~ s/[c]+//g;
	}
	if (lc($nt) eq "g"){
		$s=~ s/[g]+//g;
	}
   $ntcontent=sprintf '%4.4f', (length($seq)-length($s))/length($seq);
	return $ntcontent;
}
sub internalloop{#parameter
	#>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	#===========The internalloop was used to generate the size of the biggest internal loop
	#===========Usage: internalloop($seq1,$seq2)
	#===========The size of the biggest internal loop will be returned
	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	#format internalloop(seq1,seq2)
	my $seq1;
	my $seq2;
	($seq1,$seq2)=@_;
	$seq1=~ s/[-]+//g;
	$seq1=~ s/\s+/>/g;
	$seq1=~ s/>$//;
	$seq2=~ s/[-]+//g;
	$seq2=~ s/\s+/>/g;
	$seq2=~ s/>$//;
	my @internalloop=split(">",join(">",$seq1,$seq2));
	my $lengthinternalloop=0;
	foreach (@internalloop) {
		if (length($_)>$lengthinternalloop) {
			$lengthinternalloop=length($_)
		}
	}
	return $lengthinternalloop;
}
sub internalloopnumber{#parameter
	#>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	#===========The internalloopnumber was used to generate the number of internal loop
	#===========Usage: internalloopnumber($seq),the second or the third one
	#===========The number of internal loop will be returned
	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	my $seq;
	my $tmp;
	my $lengthinternalloopnumber;
	($seq)=@_;
	$seq=~ s/\s+/-/g;
	$tmp=$seq;
	$tmp=~s/-//g;
	$lengthinternalloopnumber=length($seq)-length($tmp);
	return $lengthinternalloopnumber;
}
sub unpairedbases{#parameter
	#>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	#===========The unpairedbases was used to generate the size of unpaired bases
	#===========Usage: unpairedbases($seq1,$seq2)
	#===========The size of unpaired bases will be returned
	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	my $seq1;
	my $seq2;
	my $s;
	($seq1,$seq2)=@_;
	$s=join("",$seq1,$seq2);
	$s=~ s/[\s-]+//g;
	return length($s);
}
sub unpairedrate{#parameter
	#>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	#===========The unpairedrate was used to generate the unpaired rate of seq
	#===========Usage: unpairedrate(unpairedbases, seq1,seq2)
	#===========The unpaired rate of seq will be returned
	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	my $unpairedbases;
	my $seq1;
	my $seq2;
	my $s;
	my $rate;
	($unpairedbases,$seq1,$seq2)=@_;
	$s=join("",$seq1,$seq2);
	$s=~ s/\s+//g;
	if (length($s) ne 0) {
		$rate=sprintf '%4.4f',($unpairedbases)/($unpairedbases+length($s));
	}
	else{
		$rate="NULL";
	}
	return $rate;
}
sub predict{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The predict was used to predict the data with svm
   #===========Usage: predict($abname,$species)
   #===========The collect seq will be put in *._mir.nm
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #the species and parameter should be used
   my $rst=0;
   my $label=2;
   my @row;
   my (%para)=@_;
   my @svmpmt;
   my @rangeparameter;
   my %rangeup;
   my %rangedown;
   my $rangeresult;
   my @svmone;
   my $value;
   my $dataset;
   @svmpmt=pmt('svm');
   @rangeparameter=pmt("range");
   %rangeup=range("up");
   %rangedown=range("down");
   $rangeresult=0;
   #range filter
   foreach (@rangeparameter){
	  $rangeresult =compare($para{$_},$rangeup{$_},$rangedown{$_});
	  if ($rangeresult eq 0){
		 $rst=0;
		 goto endline;
	  }
   }
   #prepare the data
   @svmone=();
   foreach (@svmpmt){
	  if (lc($para{$_}) eq "null" or $_ eq "firstbase" or $_ eq "penultimateposition" or $_ eq "terminalnucleotide"){
		 $value=nt2number($para{$_});
	  }
	   else{
		 $value=$para{$_};
	   }
	   push(@svmone,$value);
   }
   $dataset = new Algorithm::SVM::DataSet(Label => $label,Data  => [@svmone]);
   $rst = $svmmodel->predict($dataset);
   endline:
   return $rst;
}
sub nt2number{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The nt2number was used to translate the nt to ACC number
   #===========Usage: nt2number(nt)
   #===========The number of acc will be return
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my $value;
   ($value)=@_;
   $value=lc($value);
   if ($value eq "null"){
	  $value="-1";
   }
   elsif ($value eq ""){
	  $value="";#if nothing
   }
   elsif ($value eq "a"){
	  $value="1";
   }
   elsif ($value eq "c"){
	  $value="2";
   }
   elsif ($value eq "g"){
	  $value="3";
   }
   elsif (($value eq "u") or ($value eq "t")){
	  $value="4";
   }
   elsif ($value eq "n"){
	  $value="0";
   }
   return $value;
}
sub out{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The out was used to generate the result
   #===========Usage: out($abname)
   #===========The *_mir_overview.out and *_mir_parameter.out will be created
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my (%para)=@_;
   my $blast="";
   my $id;
   my $miseq;
   my @pmt;
   my $l;
   my $data;
   if (exists($mirbase{lc($para{'miseq'})})){
	  $blast=$mirbase{lc($para{'miseq'})};
   }
   if (length($para{'miid'})>=20){
	  $id="$para{'miid'}\t";
   }
   else{
	  $id=$para{'miid'}.(" " x (20-length($para{'miid'})));
   }
   $miseq=uc($para{'miseq'}).(" " x (30-length($para{'miseq'})));
   open (OT,">>$abname\_mir_overview_$level.out");
   print OT $id.$miseq.$blast."\n";
   close OT;
   @pmt=pmt('all');
   $l=length($para{'miid'});
   $data=("=" x ((84-$l)/2)).">$para{'miid'}<".("=" x ((84-$l)/2))."\n";
   foreach (@pmt){
	   $data.=uc($_).":".(" " x (30-length($_))).$para{$_}."\n";
   }
   $data.=("-" x 86)."\n\n";
   open (OT,">>$abname\_mir_parameter_$level.out");
   print OT "$data";
   close OT;
}
sub compare{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The compare was used to decide whether a value in in certain range
   #===========Usage: compare(n1,n2,n3,)
   #===========the value 1 for true or 0 for flase will be return
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my $n1="";
   my $n2="";
   my $n3="";
   my $result;
   ($n1,$n2,$n3)=@_;
   if (lc($n1) eq "null" or lc($n1) eq "a" or lc($n1) eq "u" or lc($n1) eq "c" or lc($n1) eq "g" or lc($n1) eq ""){
	  $n1=nt2number($n1);
   }
   if ($n1 eq "-1"){
	  $result=1;
   }
   else{
	  if ($n1>=$n2 and $n1<=$n3){
		 $result=1;
	  }
	  else{
		 $result=0;
	  }
   }
   return $result;
}
sub pmt{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The pmt was used to provide pmt files
   #===========Usage: pmt(parameter)
   #			--all all parameters
   #			--range the range parameters
   #			--svm the svm parameters
   #			--display the display parameters
   #===========the parameters will be return
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my $class;
   my @parameter;
   ($class)=@_;
   if (lc($class) eq "range"){
	  @parameter=("milength","prelength","prilength","length_basalsegment",
				  "length_lowerstem","length_upperstem","length_topstem",
				  "length_terminalloop","mipairs","prepairs","pripairs",
				  "premfe","primfe","migc","prigc","mintcontent_a","mintcontent_c",
				  "mintcontent_g","mintcontent_u","prentcontent_a","prentcontent_c",
				  "prentcontent_g","prentcontent_u","printcontent_a","printcontent_c",
				  "printcontent_g","printcontent_u","miinternalloop","preinternalloop",
				  "internalloop_lowerstem","internalloop_topstem","miinternalloopnumber",
				  "priinternalloopnumber","internalloopnumber_lowerstem",
				  "internalloopnumber_topstem","miunpairedbases","preunpairedbases",
				  "priunpairedbases","unpairedbases_lowerstem","unpairedbases_topstem",
				  "miunpairedrate","preunpairedrate","priunpairedrate","unpairedrate_lowerstem",
				  "unpairedrate_topstem","migu","mistart","miend","upperstart",
				  "upperend","stability","preinternalloopnumber","pregc",
				  "priinternalloop","pregu","prigu");#56
   }
   elsif (lc($class) eq "svm"){
        @parameter=("unpairedrate_lowerstem","prelength","internalloopnumber_lowerstem","length_upperstem",
        "miinternalloop","firstbase","mintcontent_a","migc","pregc","prentcontent_u",
        "prentcontent_a","internalloop_topstem","preunpairedrate","mipairs","prepairs",
        "internalloopnumber_topstem","unpairedrate_topstem","mintcontent_c","mistart",
        "miunpairedrate","mintcontent_g","terminalnucleotide","prentcontent_c",
        "prentcontent_g","mintcontent_u");#25pmt
   }
   elsif (lc($class) eq "all"){
	  @parameter=("miid","priid",
	  "miseq","mistr_1","mistr_2","mistr_3","mistr_4",
	  "preseq","prestr_1","prestr_2","prestr_3","prestr_4",
	  "priseq","pristr_1","pristr_2","pristr_3","pristr_4",
	  "milength","prelength","prilength",
	  "length_basalsegment","length_lowerstem","length_upperstem","length_topstem","length_terminalloop",
	  "mipairs","prepairs","pripairs",
	  "premfe","primfe",
	  "migc","pregc","prigc",
	  "mintcontent_a","mintcontent_c","mintcontent_g","mintcontent_u",
	  "prentcontent_a","prentcontent_c","prentcontent_g","prentcontent_u",
	  "printcontent_a","printcontent_c","printcontent_g","printcontent_u",
	  "miinternalloop","preinternalloop","priinternalloop",
	  "internalloop_lowerstem","internalloop_topstem",
	  "miinternalloopnumber","preinternalloopnumber","priinternalloopnumber",
	  "internalloopnumber_lowerstem","internalloopnumber_topstem",
	  "miunpairedbases","preunpairedbases","priunpairedbases",
	  "unpairedbases_lowerstem","unpairedbases_topstem",
	  "miunpairedrate","preunpairedrate","priunpairedrate",
	  "unpairedrate_lowerstem","unpairedrate_topstem",
	  "migu","pregu","prigu",
	  "strand","firstbase","penultimateposition","terminalnucleotide",
	  "mistart","miend","upperstart","upperend","stability");#77
   }
   else{
	  @parameter=();
   }
   return @parameter;
}
sub range{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The range was used to provide range information
   #===========Usage: range(parameter)
   #			--up all parameters
   #			--down the range parameters
   #===========the corespond value will be return
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my ($id)=@_;
   my @range;
   my %range;
   my @rangeparameter=pmt("range");#the range parameter
   my $m=0;
   if (lc($id) eq "up"){
	  @range=("17","41","49","0","-1","17","-1","3","10","12","15","-172.2",
			  "-307.7","0.16","0.21","0","0","0","0","0.05","0.06","0.08",
			  "0.07","0.05","0.07","0.1","0.08","0","0","-1","-1","0","0",
			  "-1","-1","0","0","0","-1","-1","0","0","0","-1","-1","0",
			  "1","20","1","19","-1","-1","0.2","0","0","0");
   }
   elsif (lc($id) eq "down"){
#	  @range=("27","250.96","379","24","116.42","35","103.85","27.42","23",
#			  "100.42","150","-6.16","-9.83","0.86","0.83","0.57","0.59",
#			  "0.7","0.6","0.43","0.44","0.49","0.47","0.42","0.43","0.5",
#			  "0.45","12","17.42","16","17.42","6","20.85","13.42","12","26",
#			  "57.42","80.92","56.54","52.85","0.51","0.47","0.46","0.64",
#			  "0.75","6","258.38","279.65","117.42","136.42","5.21","14",
#			  "0.85","23.85","9","12");
	  @range=("27","250.96","379","24","116.42","35","103.85","27.42","23",
			  "100.42","150","-6.16","-20","0.86","0.83","0.57","0.59",
			  "0.7","0.6","0.43","0.44","0.49","0.47","0.42","0.43","0.5",
			  "0.45","12","17.42","16","17.42","6","20.85","13.42","12","26",
			  "57.42","80.92","56.54","52.85","0.51","0.47","0.46","0.64",
			  "0.75","6","258.38","279.65","117.42","136.42","5.21","14",
			  "0.85","23.85","9","12");
   }
   for ($m=0;$m<@rangeparameter;$m+=1){
	  $range{$rangeparameter[$m]}=$range[$m];
   }
   return %range;
}
sub strtrim{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The strtrim was used to provide the mfe and first four nts from *.str
   #===========Usage: strtrimstr.dat)
   #===========the corespond value will be return
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my ($onetime,@data)=@_;
   my $tmp;
   my @tmp;
   my $mfe="0";
   my $do=1;
   my $stop=0;
   #first loop
   while ($do eq 1){
	  blacktest:
	  $tmp=$data[0];
	  $tmp=~ s/\n//g;
	  $tmp=~ s/\.//g;
	  $tmp=~ s/\s+//g;
	  if (length($tmp) eq 0){
		 shift(@data);
		 goto blacktest;
	  }
	  if (index($data[0],"dG") ne "-1"){
		 $data[0]=~ s/\s+/ /g;
		 @tmp=split(" ",$data[0]);
		 while ($tmp[0] ne "="){
			shift(@tmp);
		 }
		 $mfe=$tmp[1];
		 shift(@data);
		 goto blacktest;
	  }
	  if (index(lc($data[0]),"structure") ne "-1"){
		 $stop=1;
		 shift(@data);
		 goto blacktest;
		 if ($onetime ne 1){
			goto trimend;
		 }
	  }
	  if(index(lc($data[0]),"________") ne "-1"){
		 $stop=1;
		 goto trimend;
	  }
	  $tmp=$data[0];
	  $tmp=~ s/\d//g;
	  if (length($tmp) ne length($data[0])){
		 shift(@data);
		 goto blacktest;
	  }
	  else{
		 $do=0;
	  }
   }
   trimend:
   return ($stop,$mfe,@data);
}
sub head{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The head was used to generate head file for out file
   #===========Usage: strtrimstr.dat)
   #===========the corespond value will be return
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my $head;
   my @data;
   $head.= "miRNA predicting result by miRPara1.5\n";
   $head.= "By Yonggan Wu and Simon Rayner\n";
   $head.= 'Report bugs to weasteam@gmail.com or raynere@wh.iov.cn'."\n";
   $head.= "Wuhan Institute of Virology, Chinese Academic of Science.\n";
   $head.= "Wuhan, 430071, China\n\n";
   $head.= ("-" x 86)."\n";
   $head.= "Your data was predicted at level: $level (Total 20 levels)\n";
   $head.= "The number of pri-miRNAs, miRNA regions and miRNA candidates are:\n";
   $head.= "pri-miRNA           $rstcount{'pri'}\n";
   $head.= "miRNA candidates    $rstcount{'candidate'}\n";
   $head.= "miRNA region        $rstcount{'area'}\n";
   open (IN,"$abname\_mir_overview_$level.out");
   @data=<IN>;
   close IN;
   unshift (@data,$head,("-" x 86)."\n",
			"Name                miRNA sequences               blast in miRBase13.0\n",
			("-" x 86)."\n");
   open (OUT,">$abname\_mir_overview_$level.out");
   print OUT @data;
   if ($rstcount{'area'} eq 0){
	  print OUT "No miRNA avaiable in your sequence.\n";
	  print OUT "You can try again with lower strict level.\n";
   }
   close OUT;
   open (IN,"$abname\_mir_parameter_$level.out");
   @data=<IN>;
   close IN;
   unshift (@data,$head);
   open (OUT,">$abname\_mir_parameter_$level.out");
   print OUT @data;
   if ($rstcount{'area'} eq 0){
	  print OUT "\nNo miRNA avaiable in your sequence.\n";
	  print OUT "You can try again with lower strict level.\n";
   }
   close OUT;
}
sub readfas{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The readfas was used to read the fasta files;
   #===========Usage: readfas
   #===========the the id and seq will be return
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my ($infile)=@_;
   my @data;
   my $id;
   my $seq="";
   my @seq;
   open (IN,$infile);
   @data=<IN>;
   close IN;
   while (@data>=1){
	  if (index($data[0],">") ne "-1"){
		 if ($seq ne ""){
			push(@seq,"$id\t$seq");
		 }
		 $id=$data[0];
		 $id=~s/\n//g;
		 $id=~s/^>//;
		 if (index($id," ") ne "-1"){
			$id=substr($id,0,index($id," "));
		 }
		 $id=~s/\|/\_/g;
		 $id=~s/^\_//;
		 $id=~s/\_$//;
		 shift(@data);
		 $seq="";
	  }
	  else{
		 $data[0]=~s/\n//g;
		 $data[0]=~s/\s+//g;
		 $seq.=lc($data[0]);
		 shift(@data);
	  }
	  
   }
   $seq=~s/t/u/g;
   $seq=~s/\n//g;
   $seq=~s/\r//g;
   push(@seq,"$id\t$seq");
   return @seq;
}
sub pmtpredict{
   my @data;
   my $m;
   my @pmt=pmt('all');
   my ($infile)=@_;
   my %para;
   my @tmp;
   my $tmp;
   my $rst;
   my $serial=0;
   my $tmp2;
   #clear the files;
   open (OT,">$abname\_mir_overview_$level.out");
   print OT "";
   close OT;
   open (OT,">$abname\_mir_parameter_$level.out");
   print OT "";
   close OT;
   open (IN,"$infile");
   @data=<IN>;
   close IN;
   shift (@data);
   #generate the parameter
   foreach (@data){
	  $_=~s/\n//;
	  %para=();
	  @tmp=split("\t",$_);
	  for($m=1;$m<=@pmt;$m+=1){
		 $para{$pmt[$m-1]}=$tmp[$m-1];
	  }
	  print "Predicting to $para{'miid'}.";
	  $tmp=$para{'miid'};
	  $tmp=~s/-//g;
	  $m=length($para{'miid'})-length($tmp);
	  $tmp=$para{'miid'};
	  while ($m>1){
		 $tmp=~s/-//;
		 $m-=1;
	  }
	  $tmp2=$tmp;
	  $tmp2=~s/\_//g;
	  $m=length($tmp)-length($tmp2);
	  while ($m>1){
		 $tmp=~s/\_//;
		 $m-=1;
	  }
	  $serial=substr($tmp,index($tmp,"-")+1,index($tmp,"_")-index($tmp,"-")-1);
	  $rst=predict(%para);
	  if ($rst eq "1"){
		 print "...TRUE\n";
		 $rstcount{'candidate'}+=1;
		 if (($serial-$rstcount{'preserial'})>=2){
			$rstcount{'area'}+=1;
		 }
		 $rstcount{'preserial'}=$serial;
		 if ($para{'priid'} ne $rstcount{'prepriid'}){
			$rstcount{'pri'}+=1;
			$rstcount{'prepriid'}=$para{'priid'};
		 }
		 out(%para);
	  }
	  else{
		 print "...FALSE\n";
	  }
   }
}
