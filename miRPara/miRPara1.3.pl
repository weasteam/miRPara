#! /usr/bin/perl -w
#mirpara.pl
#################################USE############################################
use warnings;
use strict;
use Bio::Seq;#read and output sequences
use Bio::SeqIO;#read and output sequences
use Getopt::Long;
use Algorithm::SVM;
use Algorithm::SVM::DataSet;
use File::chdir;#to change the perl environment
##################################GOCOBLE#######################################
our $abname="noname";
our $species="overall";
our $level=7;
our $share="/home/weasteam/biology/WIV/Simon_Rayner/RNAi/datas/program/miRPara/test/";#share data,,add"/" at last
our $prilength=60;
our %rstcount;
our %mirbase;
my %options;
my $path;
my $infile="";
my $pri="";
my @splitname;
my @basedata;
my @base;
################################PROGROM START###################################
################################################################################
################################################################################
GetOptions \%options, 'version|v' => sub { version('miRPara.pl') },
						'help|h' => sub { usage() },
						'name|n=s' => \$abname,
						'species|s=s' => \$species,
						'level|l=i'=>\$level,
						'pri'=>\$pri,
						'prilengthlimit|p=i' => \$prilength or die("Type miRPara.pl -h for help!\n");
#################################data prepare###################################
our $svmmodel = new Algorithm::SVM(Model => $share."models/$species\_$level.model");
##infile
$infile=shift or die("Type miRPara.pl -h for help!\n");
if (length($infile) eq 0) {
   print STDERR "Error: data not specified\nRun 'miRPara.pl -h' for help\n";
   exit;
}
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
##blast data
open (IN,$share."mirbase.dat");#read the blast files
@basedata=<IN>;
close IN;
foreach (@basedata){
   $_=~ s/\n+//g;
   if ($_ ne ""){
	  @base=split(":",$_);
	  $mirbase{lc($base[0])}="$base[1]";
   }
}
@basedata=();
@base=();
##abname
if ($abname eq "noname"){#reset the abname
   if (index($infile,".") ne "-1"){
	  $abname=substr($infile,0,index($infile,"."));
   }
   else{
	  $abname=$infile;
   }
}
##extra parameter
if ($pri eq "1"){
   goto step4;
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
candidate();#followed by parameter and pmt2svmdata
##################step5##################
step5:
head();
if ($rstcount{'area'} ne 0){
   system "rm tmp*";#delete all tmp files
}
print "\nAll Done.\n";
##################STOP HERE##################
sub usage () {
    print <<EOF;
Usage: miRPara.pl [options] file [file]

Options:
-v, --version
-h, --help
-n, --name=<abbreviated name (3 characters) of the species>
-s, --species=<overall, animal, plant or virus> (defaults as overall)
-l, --Level=<1..10>(defaults to 7)
-p, --prilengthlimit=<limit to the pri-miRNA length> (defaults to 60)
-pri, --Consider all input sequences as pri-miRNAs

EOF
print 'Report bugs to weasteam@gmail.com or raynere@wh.iov.cn', "\n";
    exit;
}
sub version ($) {
	print "\n";
    print "$_[0] (miRPare) 1.3\n";
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
   my $infile;#the file name
   my $splitlength;#the output length
   my $overlap;#the overlap length
   my $seq_obj;#the sequence object
   my $seq;#the sequence of the file
   my $seqlength;#the length of the input sequences
   my $seqname;#the name of new sequences
   my $splitseq;#teh splited sequences
   my $split_obj;#the sequence of splited sequences
   my $outputseq;#the sequence of output sequences
   my $start;#start of a short sequence
   my $end;#end of a short sequence
   my $m=1;#the number to read different seqs
   my @splitname;
   ($infile,$splitlength,$overlap)=@_;#input the file
   open (OUT,">$abname\_splitseq.fasta");#CLEAR THE DATA
   print OUT "";
   close OUT;
   $seq_obj = Bio::SeqIO->new(-file => "$infile");
   while ($seq = $seq_obj->next_seq){#get the sequence
	  $seqlength = $seq->length();#get the length of the whole sequence
	  $start=1;#start of a short sequence
	  $end = $start -1 + $splitlength;#end of a short sequence
	  if ($seqlength<=$splitlength) {#if the sequence if smaller than the span-$splitlength, output the file directly
		 $seqname="$abname\_$m\_1-$seqlength";
		 push(@splitname,$seqname);
		 $splitseq=lc($seq->subseq(1,$seqlength));
		 $splitseq=~ s/t/u/g;
		 $split_obj = Bio::Seq->new(-seq => $splitseq,
				 -display_id => $seqname);
		 $outputseq = Bio::SeqIO->new(-file => ">>$abname\_splitseq.fasta",
				 -format => 'fasta' );
		 $outputseq->write_seq($split_obj);
	  }
	  else {#if the sequence is longer than the span-$splitlength, be cutted with $overlap overlapped
		 while ($start <= $seqlength) {
			if ($end >=$seqlength) {#in the case that the last loop
			   $end = $seqlength;
			}
			$seqname="$abname\_$m\_$start-$end";#define a fine name of new seq
			push(@splitname,$seqname);
			$splitseq=lc($seq->subseq($start,$end));
			$splitseq=~ s/t/u/g;
			$split_obj = Bio::Seq->new(-seq => $splitseq,
					-display_id => $seqname);#extract the seq
			$outputseq = Bio::SeqIO->new(-file => ">>$abname\_splitseq.fasta",
					-format => 'fasta' );#print the seq
			$outputseq->write_seq($split_obj);
			if ($end ne $seqlength) {#redefine the start and the end
			   $start = $end -$overlap +1;
			   $end = $start -1 + $splitlength;
			}
			else {#in the case that the last loop
			   $start = $end +1;
			}
		 }
	  }
	  $m+=1;
   }
   $seq_obj="";
   $seq="";
   $splitseq="";
   $split_obj="";
   $outputseq="";
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
	my $seq_obj;
	my $seq;
	my $id;
    ($parameter)=@_;#-l for the long file, and -s for the small file
    if ($parameter eq "long"){
	  $seq_obj = Bio::SeqIO->new(-file => "$abname\_splitseq.fasta");
	  while ($seq = $seq_obj->next_seq){
		 $id=$seq->display_id();
		 print "RNA folding with UNAFold to $id\n";
		 open (OUT,">tmp.fasta");
		 print OUT ">tmp\n";
		 print OUT $seq->seq;
		 close OUT;
		 system "UNAFold.pl tmp.fasta > error.log";
		 system "ct2out <tmp.fasta_1.ct> $id\.str";#run ct2out
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
   my $serial=1;#the number of hairpins
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
   (@splitname)=@_;
   open (NAME,">$abname\_hairpin.fasta");#clean the fasta file
   print NAME "";
   close NAME;
   foreach (@splitname){
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
			   $seq=hairpin2seq($hairpin{1},$hairpin{2},$hairpin{3},$hairpin{4});
			   foreach (@sequence){#whether it is repeat to provious seq
				  if ($seq eq $_){
					 $rp=1
				  }
			   }
			   if ($rp eq 0){
				  open (HAIRPIN,">>$abname\_hairpin.fasta");
				  print HAIRPIN ">$abname-mir-$serial","\n";
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
   @data=();
   %hairpin=();
   @sequence=();
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
   my $id="";#the id of the sequence
   my @strdata;#the data of *.out files
   my $seq_obj;#the object of sequence
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
   open (OUTPUT, ">$abname\.dat");#tab separate files
   print OUTPUT "";
   close OUTPUT;
   open (OT,">$abname\_mir_overview.out");
   print OT "";
   close OT;
   open (OT,">$abname\_mir_parameter.out");
   print OT "";
   close OT;
   $seq_obj = Bio::SeqIO->new(-file => "$abname\_hairpin.fasta");
   while ($seq = $seq_obj->next_seq){
	  $rstcount{'pri'}+=1;
	  $id=$seq->display_id();
	  open (OUT,">tmp.fasta");
	  print OUT ">$id\n";
	  print OUT $seq->seq;
	  close OUT;
	  callunafd("tmp");
	  open (DATA,"tmp.str");#read the name of the long sequences
	  @strdata=<DATA>;
	  close DATA;
	  ($m,$para{'primfe'},@strdata)=strtrim(1,@strdata);#read the minimal free energy
	  $para{'pristr_1'}=lc($strdata[0]);#read the second structure
	  $para{'pristr_2'}=lc($strdata[1]);
	  $para{'pristr_3'}=lc($strdata[2]);
	  $para{'pristr_4'}=lc($strdata[3]);
	  $para{'pristr_1'}=~ s/\n//g;
	  $para{'pristr_2'}=~ s/\n//g;
	  $para{'pristr_3'}=~ s/\n//g;
	  $para{'pristr_4'}=~ s/\n//g;
	  $para{'priseq'} = lc($seq->seq);#get the sequence
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
			while (($m<=24) and ($position+$m)<=length($seq5)){
			   $para{'id'}="$abname-MIR-$serial\_$m";
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
			   $para{'prilenth'}=length($para{'priseq'});#prilength
			   $para{'milength'}=$m;#milength
			   $s=join("",$para{'prestr_1'},$para{'prestr_2'},$para{'prestr_3'},$para{'prestr_4'});
			   $s=~ s/[\s\\-]+//g;
			   $para{'prelength'}=length($s);
			   #seq
			   $para{'preseq'}=substr($para{'priseq'},$position-1,$para{'prelength'});
			   $para{'miseq'}=substr($seq5,$position-1,$m);
			   #strand
			   $para{'strand'}='5';
			   print "Calculating miRNA candidate: $para{'id'} at $id\n";
			   open (OUT,">>mirpara.log");
			   print OUT "Calculating miRNA candidate: $para{'id'} at $id\n";
			   close OUT;
			   parameter(%para);#sent the information to calculate the parameters
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
			   $para{'id'}="$abname-MIR-$serial\_$m";
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
			   $para{'prilenth'}=length($para{'priseq'});#prilength
			   $para{'milength'}=$m;#milength
			   $s=join("",$para{'prestr_1'},$para{'prestr_2'},$para{'prestr_3'},$para{'prestr_4'});
			   $s=~ s/[\s\\-]+//g;
			   $para{'prelength'}=length($s);
			   #seq
			   $para{'preseq'}=substr($para{'priseq'},$m-$position-$para{'prelength'},$para{'prelength'});
			   $para{'miseq'}=substr($seq3,-$position,$m);
			   #strand
			   $para{'strand'}='3';
			   print "Calculating miRNA candidate: $para{'id'} at $id\n";
			   open (OUT,">>mirpara.log");
			   print OUT "Calculating miRNA candidate: $para{'id'} at $id\n";
			   close OUT;
			   parameter(%para);#sent the information to calculate the parameters
			   $m +=1;
		   }
		   $serial +=1;
		   $position-=1;
		 }
	  }#no budding stem
	  line1:
   }#each seq
   %para=();
   @strdata=();
   $s="";
   $seq5="";
   $seq3="";
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
   (%para)=@_;#receive the parameter
   #length
   $para{'prilength'}=length($para{'priseq'});
   $para{'prelength'}=length($para{'preseq'});
   $para{'milength'}=length($para{'miseq'});
   #premfe
   open (OUTPUT,">tmp.fasta");#output the pre seq
   print OUTPUT ">tmp\n";
   print OUTPUT $para{'preseq'};
   close OUTPUT;
   callunafd('tmp');#call unafold to generate the mfe
   open (DATA,"tmp.str");
   @strdata=<DATA>;
   close DATA;
   ($m,$para{'premfe'},@strdata)=strtrim(1,@strdata);
   #gc
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
   $s1=$para{'pristr_1'};
   chop($s1);#get rid of the black in the terminal
   my $terminalloopstem=0;
   until (substr($s1,-1,1) eq " "){#get the length of terminalloopstem
	  chop($s1);
	  $terminalloopstem++;
   }
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
   $rst=pmt2svmdat(%para);
   #output
   @pmt=pmt("all");
   foreach (@pmt){
	  $pmt.= "$para{$_}\t";
   }
   chop($pmt);
   open (DATAOUT, ">>$abname\_candidate.pmt");
   print DATAOUT "$pmt\n";
   close DATAOUT;
   if ($rst eq "1"){
	  $rstcount{'candidate'}+=1;
	  out(%para);
   }
   %para=();
   @strdata=();
   $s="";
   $s1="";
   $s2="";
   $pmt="";
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
	$s=join("",$seq1,$seq2);
	$s=~ s/\s+//g;
	$u=lc($s);
	$u=~ s/[ut]+//g;
	$a=lc($s);
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
		$ntcontent=sprintf '%4.4f', (length($seq)-length($s))/length($seq);
	}
	if ((lc($nt) eq "t") or (lc($nt) eq "u")){
		$s=~ s/[tu]+//g;
		$ntcontent=sprintf '%4.4f', (length($seq)-length($s))/length($seq);
	}
	if (lc($nt) eq "c"){
		$s=~ s/[c]+//g;
		$ntcontent=sprintf '%4.4f', (length($seq)-length($s))/length($seq);
	}
	if (lc($nt) eq "g"){
		$s=~ s/[g]+//g;
		$ntcontent=sprintf '%4.4f', (length($seq)-length($s))/length($seq);
	}
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
	$seq2=~ s/[-]+//g;
	$seq2=~ s/\s+/>/g;
	if (substr($seq1,-1,1) eq ">"){#get rid of the final > which might black splict.
		chop($seq1);
	}
	if (substr($seq2,-1,1) eq ">"){
		chop($seq2);
	}
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
	if (($unpairedbases+length($s)) ne 0) {
		$rate=sprintf '%4.4f',($unpairedbases)/($unpairedbases+length($s));
	}
	else{
		$rate=sprintf '%4.4f',1;
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
   my $testdata;
   my $rst=0;
   my @svmdata;
   my $label;
   my $testname;
   my @row;
   my ($data)=@_;
   $data =~ s/\n//g;
   @row=split(" ",$data);
   $label=substr($row[0],0,1);
   $testname=substr($row[0],2);
   @svmdata=split(",",$row[1]);
   $testdata = new Algorithm::SVM::DataSet(Label => 2,
									Data  => [@svmdata]);
   $rst = $svmmodel->predict($testdata);
   @svmdata=();
   @row=();
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
	  $value=0;#if nothing
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
   my $pmtdata;
   if (exists($mirbase{lc($para{'miseq'})})){
	  $blast=$mirbase{lc($para{'miseq'})};
   }
   $id=$para{'id'}.(" " x (20-length($para{'id'})));
   $miseq=uc($para{'miseq'}).(" " x (30-length($para{'miseq'})));
   open (OT,">>$abname\_mir_overview.out");
   print OT $id.$miseq.$blast."\n";
   close OT;
   $pmtdata =display(%para);
   open (OT,">>$abname\_mir_parameter.out");
   print OT "$pmtdata";
   close OT;
}
sub display{
	#>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	#===========The display was used to display the data files
	#===========Usage: display(%para)
	#===========The data will be returned
	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	my %para;
	my $data;
	my @p;
	my $l;
	(%para)=@_;
	@p=pmt('all');
	$l=length($para{'id'});
	$data=("=" x ((84-$l)/2)).">$para{'id'}<".("=" x ((84-$l)/2))."\n";
	foreach (@p){
		$data.=uc($_).":".(" " x (30-length($_))).$para{$_}."\n";
	}
	$data.=("-" x 86)."\n\n";
	%para=();
	@p=();
	return $data;
}
sub pmt2svmdat{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The pmt2svmdat.pl was used to generate the pmt data to svm data
   #===========Usage: perl pmt2xls.pl *.pmt
   #			Note: only *.pmt files were avaiable
   #===========The data were pre-filtered with range99
   #===========A each parameter value will be put in a *.dat
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my @svmparameter;#svm parameter data
   my @rangeparameter;#range parameter
   my $dat;#positive svm data
   my $value;#each single value
   my $rangeresult;
   my %para;
   my %rangeup;
   my %rangedown;
   my $rst=0;
   my $count;
   #print "Converting to SVM data";
   (%para)=@_;
   @svmparameter=pmt('svm');
   @rangeparameter=pmt("range");
   %rangeup=range("up");
   %rangedown=range("down");
   $rangeresult=0;
   $count=@rangeparameter;
   foreach (@rangeparameter){
	  $rangeresult +=compare($para{$_},$rangeup{$_},$rangedown{$_});
   }
   if (($rangeresult eq $count) and (uc($para{'buddingstem'}) eq "NO")){
	  $value=$para{"id"};
	  $dat="2_$value ";#clear the data
	  foreach (@svmparameter){
		 $value=nt2number($para{$_});
		 $dat .="$value,";
	  }
	  chop($dat);
	  $rst=predict($dat);
	  open (OUTPUT, ">>$abname\.dat");
	  print OUTPUT "$dat\n";
	  close OUTPUT;
   }
   @svmparameter=();
   @rangeparameter=();
   $dat="";
   $rangeresult="";
   %para=();
   return $rst;
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
	if ($n1 eq "NULL"){
		$n1="-";
	}
	if ($n2 eq "-"){
		$n2="";
	}
	if ($n3 eq "-"){
		$n3="";
	}
	if ($n2 ne ""){
		if (($n3 ne "") and ($n1 ne "-")){
			if (($n1 >=$n2) and ($n1<=$n3)){
			$result=1;
			}
			else{
				$result=0;
			}
		}
		else{
			if ($n1 ne "-"){
				if ($n1 >=$n2){
					$result=1;
				}
				else{
					$result=0;
				}
			}
			else{
				$result=1;
			}
		}
	}
	else{
		if (($n3 ne "")and ($n1 ne "-")){
			if (($n1<=$n3) and ($n1 ne "-")){
				$result=1;
			}
			else{
				$result=0;
			}
		}
		else{
			$result=1;
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
	  @parameter=("mistart","miend","milength","prelength","prilength",
				  "length_basalsegment","length_lowerstem","length_upperstem",
				  "length_topstem","length_terminalloop","mipairs","prepairs",
				  "pripairs","premfe","primfe","migc","pregc","prigc",
				  "mintcontent_a","mintcontent_c","mintcontent_g","mintcontent_u",
				  "prentcontent_a","prentcontent_c","prentcontent_g",
				  "prentcontent_u","printcontent_a","printcontent_c",
				  "printcontent_g","printcontent_u","miinternalloop",
				  "preinternalloop","priinternalloop","internalloop_lowerstem",
				  "internalloop_topstem","miinternalloopnumber",
				  "preinternalloopnumber","priinternalloopnumber",
				  "internalloopnumber_lowerstem","internalloopnumber_topstem",
				  "miunpairedbases","preunpairedbases","priunpairedbases",
				  "unpairedbases_lowerstem","unpairedbases_topstem",
				  "miunpairedrate","preunpairedrate","priunpairedrate",
				  "unpairedrate_lowerstem","unpairedrate_topstem","migu","pregu",
				  "prigu","strand","firstbase","penultimateposition",
				  "terminalnucleotide","loopstart","stability","upperstart",
				  "upperend","basalend");#62
   }
   elsif (lc($class) eq "svm"){
	  @parameter=("firstbase","prepairs","prelength","mipairs","mistart",
				   "miinternalloop","terminalnucleotide","prigc");#8
   }
   elsif (lc($class) eq "all"){
	  @parameter=("id",
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
	  "mistart","miend","stability");#74
   }
   else{
	  @parameter=();
   }
   return @parameter;
}
sub summary{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The summary was used to analysis the number of pri and mi RNAs
   #===========Usage: summary($abname)
   #===========the number of pri and miRNA will be return.
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my @data;
   my $m;#tmp number
   my $previous="";#previous candidate
   my @tmp;
   open (IN,"$abname\_mir_overview.out");
   @data=<IN>;
   close IN;
   $m=0;
   if ($rstcount{'candidate'} ne 0){
	  if (index($data[0],"-miR-") ne "-1"){
		 $rstcount{'area'}=@data;#number of miRNA
	  }
	  elsif (index($data[0],"-MIR-") ne "-1"){
		 $rstcount{'area'}=0;
		 foreach (@data){
			@tmp=split(" ",$_);
			$_=$tmp[0];
			if ($_ ne ""){
			   if ((lc(substr($_,0,3)) eq lc(substr($previous,0,3))) or ($previous eq "")){
				  if (($m eq "0") or (($m ne "0") and ((substr($_,index($_,'-MIR-')+5,index($_,'_')-index($_,'-MIR-')-5)-$m)>1))){
					 $rstcount{'area'}+=1;
				  }
				  $m=substr($_,index($_,'-MIR-')+5,index($_,'_')-index($_,'-MIR-')-5);
			   }
			   else{
				  $rstcount{'area'}+=1;
			   }
			   $previous=$_;
			}
		 }
	 }
   }
   else{
	 $rstcount{'area'}=0;
   }
   @data=();
}
sub range{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The range was used to provide range information
   #===========Usage: range(parameter)
   #			--up all parameters
   #			--down the range parameters
   #===========the corespond value will be return
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #<[1]mistart><[2]miend><[3]milength><[4]prelength><[5]prilength>
   #<[6]length_basalsegment><[7]length_lowerstem><[8]length_upperstem>
   #<[9]length_topstem><[10]length_terminalloop><[11]mipairs><[12]prepairs>
   #<[13]pripairs><[14]premfe><[15]primfe><[16]migc><[17]pregc><[18]prigc>
   #<[19]mintcontent_a><[20]mintcontent_c><[21]mintcontent_g><[22]mintcontent_u>
   #<[23]prentcontent_a><[24]prentcontent_c><[25]prentcontent_g><[26]prentcontent_u>
   #<[27]printcontent_a><[28]printcontent_c><[29]printcontent_g><[30]printcontent_u>
   #<[31]miinternalloop><[32]preinternalloop><[33]priinternalloop>
   #<[34]internalloop_lowerstem><[35]internalloop_topstem>
   #<[36]miinternalloopnumber><[37]preinternalloopnumber><[38]priinternalloopnumber>
   #<[39]internalloopnumber_lowerstem><[40]internalloopnumber_topstem>
   #<[41]miunpairedbases><[42]preunpairedbases><[43]priunpairedbases>
   #<[44]unpairedbases_lowerstem><[45]unpairedbases_topstem><[46]miunpairedrate>
   #<[47]preunpairedrate><[48]priunpairedrate><[49]unpairedrate_lowerstem>
   #<[50]unpairedrate_topstem><[51]migu><[52]pregu><[53]prigu><[54]strand>
   #<[55]firstbase><[56]penultimateposition><[57]terminalnucleotide>
   #<[58]loopstart><[59]stability><[60]upperstart><[61]upperend><[62]basalend>
   my ($id)=@_;
   my @range;
   my %range;
   my @rangeparameter=pmt("range");#the range parameter
   my $m=0;
   if (lc($id) eq "up"){
	  @range=("1","20","17","40","50","0","-1","17","-1","3","0","0","10",
				   "-254.83","-453.6","0.16","0.19","0.21","0","0","0","0",
				   "0.05","0.06","0.08","0.07","0.05","0.07","0.1","0.07","0",
				   "0","0","-1","-1","0","0","0","-1","-1","0","0","0","-1","-1",
				   "0","0","0","-1","-1","0","0","0","-","-","-","-","20.74",
				   "-1","1","19","0");
   }
   elsif (lc($id) eq "down"){
	  @range=("370.77","390.77","27.26","492.25","632","92.29","130.1","35",
			  "129.52","355.67","24","116.58","209.26","-5.92","-11.55","0.87",
			  "0.85","0.83","0.6","0.62","0.7","0.6","0.44","0.44","0.5","0.49",
			  "0.43","0.44","0.51","0.46","22.26","60.81","19.26","16","40.29",
			  "11.26","29","35.26","24","25.26","44","116.03","69.52","47.81",
			  "88.84","1","1","0.41","0.67","0.81","7.26","13","15","-","-","-",
			  "-","227.48","5","168.03","188.03","113.03");
   }
   while ($m<=61){
	  $range{$rangeparameter[$m]}=$range[$m];
	  $m+=1;
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
   summary();
   $head.= "miRNA predicting result by miRPara1.0\n";
   $head.= "By Yonggan Wu and Simon Rayner\n";
   $head.= 'Report bugs to weasteam@gmail.com or raynere@wh.iov.cn'."\n";
   $head.= "Wuhan Institute of Virology, Chinese Academic of Science.\n";
   $head.= "Wuhan, 430071, China\n\n";
   $head.= ("-" x 86)."\n";
   $head.= "Your data was predicted at level: $level (Total 10 levels)\n";
   $head.= "The number of pri-miRNAs, miRNA areas and miRNA candidates are:\n";
   $head.= "pri-miRNA           $rstcount{'pri'}\n";
   $head.= "miRNA candidates    $rstcount{'candidate'}\n";
   $head.= "miRNA region        $rstcount{'area'}\n";
   open (IN,"$abname\_mir_overview.out");
   @data=<IN>;
   close IN;
   unshift (@data,$head,("-" x 86)."\n",
			"Name                miRNA sequences               blast in miRBase12.0\n",
			("-" x 86)."\n");
   open (OUT,">$abname\_mir_overview.out");
   foreach (@data){
	  print OUT $_;
   }
   if ($rstcount{'area'} eq 0){
	  print OUT "No miRNA avaiable in your sequence.\n";
	  print OUT "You can try again with lower strict level.\n";
   }
   close OUT;
   open (IN,"$abname\_mir_parameter.out");
   @data=<IN>;
   close IN;
   unshift (@data,$head);
   open (OUT,">$abname\_mir_parameter.out");
   foreach (@data){
	  print OUT $_;
   }
   if ($rstcount{'area'} eq 0){
	  print OUT "\nNo miRNA avaiable in your sequence.\n";
	  print OUT "You can try again with lower strict level.\n";
   }
   close OUT;
}
