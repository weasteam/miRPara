#! /usr/bin/perl
#	by Yonggan Wu (weasteam@gmail.com)
#	by Simon Rayner ()
#	Version: 01
#	The script will prepare the data for SVM training from miRBase release
#	Input:
#		1) a miRBase folder (optional)
#	Requirements:
#		1) tar
#	Output:
#		
################################################################################
use strict;
use warnings;
use File::chdir;
use Cwd;
my $dir = getcwd;
$CWD = $dir;
my $usage="\n\nperl $0 <miRBase Version> <Full miRPara address>\n\n";
my $version=shift or die $usage;
if ($version<13){
	die "Error: sorry, we do not support miRBase before verison 13.0\n";
}
my $mirpara=shift or die $usage;
my $folder="./miRBase/$version";
my @files=("*THIS_IS_RELEASE_*","organisms.txt*","miRNA.dat*","hairpin.fa*","mature.fa*");
#check the given folder parameter
if (-e $folder){# if found the data
	$folder=~s/\/$//;
	my @tmp=glob("$folder/$files[0]");
	print "Checking required files...";
	for (my $i=1;$i<@files;$i++){
		my @tmp=glob("$folder/$files[$i]");
		if (@tmp eq 0){
			print "Warning: Not found $folder/$files[$i]!\nGoing to download from miRBase...\n";
			system "wget -P ./miRBase/$version/ ftp://mirbase.org/pub/mirbase/$version/$_";
		}
	}
	print "All Found\n";
}
else{#if not found
	system "mkdir -p ./miRBase/$version/";
	foreach (@files){
		system "wget -P ./miRBase/$version/ ftp://mirbase.org/pub/mirbase/$version/$_";
	}
}
#prepare the evidence infor
#use pre-miRNA accession number
my %evi;
if (not -e "./miRBase/$version/miRNA.dat"){
	system "gunzip ./miRBase/$version/miRNA.dat.gz";
}
print "Extracting the evidence from miRNA.dat...";
open (IN,"./miRBase/$version/miRNA.dat");
my $ac;
my $evi="";
my %evi_group;
while (<IN>){
	if ($_=~/^AC/){
		$_=~s/;\n//;
		$_=~s/AC\s+//;
		$ac=$_;
	}
	if ($_=~/^FT                   \/evidence/){
		$_=~s/FT                   \/evidence//;
		$_=~s/\n//;
		$_=~s/=//;
		if ($evi eq "" or $evi ne "experimental"){
			$evi=$_;
		}
	}
	if ($_=~/^SQ/){
		$evi{$ac}=$evi;
		$evi_group{$evi}=1;
		$evi="";
	}
}
close IN;
print "Done\n";
#prepare the species information;
print "Preparing the species information from organisms.txt...";
if (not -e "./miRBase/$version/organisms.txt"){
	system "gunzip ./miRBase/$version/organisms.txt.gz";
}
my %spe;
my %spe_group;
my %belong;
open (IN,"./miRBase/$version/organisms.txt");
while (<IN>){
	$_=~s/\n//;
	my @tmp=split("\t",$_);
	my @dat=split(";",$tmp[@tmp-1]);
	$spe{$tmp[0]}=$dat[0];
	$spe_group{$dat[0]}=1;
	$belong{$tmp[0]}=$dat[0];
	$belong{$dat[0]}=$dat[0];
}
close IN;
print "Done\n";
#finally, output the pre-miRNA sequences
print "Outputing the grouped pre-miRNA sequences...";
if (not -e "grouped_hairpin/spe" or not -e "grouped_hairpin"){
	system "mkdir -p grouped_hairpin/spe";
}
if (not -e "./miRBase/$version/hairpin.fa"){
	system "gunzip ./miRBase/$version/hairpin.fa.gz";
}
my %fh;
my @evi=keys(%evi_group);
my @spe_group=keys(%spe_group);
my @spe=keys (%spe);
foreach (@evi){
	$evi=$_;
	open ($fh{"$evi\_all"},">./grouped_hairpin/miRBase_$version\_$evi\_all.fasta");
	foreach (@spe_group){
		open ($fh{"$evi\_$_"},">./grouped_hairpin/miRBase_$version\_$evi\_$_.fasta");
	}
	foreach (@spe){
		open ($fh{"$evi\_$_"},">./grouped_hairpin/spe/miRBase_$version\_$evi\_$_.fasta");
	}
}

open (IN,"./miRBase/$version/hairpin.fa");
my $spe_group;
my $spe;
while (<IN>){
	if ($_=~/^>/){
		my @tmp=split(" ",$_);
		if (exists $evi{$tmp[1]}){
			$evi=$evi{$tmp[1]};
		}
		else{
			$evi="";
		}
		my @dat=split("-",$tmp[0]);
		$dat[0]=~s/>//;
		$spe=$dat[0];
		if (exists $spe{$dat[0]}){
			$spe_group=$spe{$dat[0]};
		}
		else{
			$spe_group="";
		}
	}
	if (exists $fh{"$evi\_$spe_group"}){
		print {$fh{"$evi\_$spe_group"}} $_;
	}
	if (exists $fh{"$evi\_all"}){
		print {$fh{"$evi\_all"}} $_;
	}
	if (exists $fh{"$evi\_$spe"}){
		print {$fh{"$evi\_$spe"}} $_;
	}
}
close IN;
foreach (@evi){
	$evi=$_;
	close $fh{"$evi\_all"};
	foreach (@spe_group){
		close $fh{"$evi\_$_"};
	}
	foreach (@spe){
		close $fh{"$evi\_$_"};
	}
}
print "Done\n";
################################################################################
#calculate the parameters
print "Calculating the parameters...\n";
my @expfiles=glob("./grouped_hairpin/miRBase_$version\_experimental_*.fasta");
my @spefiles=glob("./grouped_hairpin/spe/miRBase_$version\_experimental_*.fasta");
@files=();
push (@files,@expfiles);
push (@files,@spefiles);
foreach (@files){
	my $file=$_;
	system "perl $mirpara --pmt $file";
}
if (not -e "calculated_parameters"){
	mkdir "calculated_parameters";
}
system "mv ./grouped_hairpin/*.pmt ./calculated_parameters";
system "mv ./grouped_hairpin/spe/*.pmt ./calculated_parameters";
################################################################################
#generate the training parameters
print "generating the training parameters...\n";
my $out="training_parameters_$version";
if (not -e $out){
	system "mkdir $out";
}
my @pmt=glob("./calculated_parameters/miRBase_$version\_*.pmt");
#unzip the data
if (not -e "./miRBase/$version/mature.fa"){
	system "gunzip ./miRBase/$version/mature.fa.gz";
}
#get the mature sequences
print "reading the mature sequences...\n";
my %mature;
my $title;
open (IN,"./miRBase/$version/mature.fa");
while (<IN>){
	$_=~s/\n//;
	if ($_=~/^>/){
		my @tmp=split(" ",$_);
		$title=$tmp[0];
		$title=~s/>//;
		$title=lc($title);
	}
	else{
		$mature{$_}.="$title,";
	}
}
close IN;
#output the sequences
foreach (@pmt){
	my $pmt=$_;
	print "extract positive for $pmt...\n";
	my %start=();#upper start
	my %end=();#upper end
	my $outpmt=$pmt;
	$outpmt=~s/\.\/calculated_parameters\/miRBase_$version\_//;
	open (OUT,">./training_parameters_$version/$outpmt");
	open (IN,$pmt);
	while (<IN>){
		$_=~s/\n//;
		if ($_=~/^#/){
			print OUT "$_\n";
		}
		else{
			my @tmp=split("\t",$_);
			my @dat=split("_",$tmp[0]);
			my $id=lc($dat[0]);
			if (exists $mature{$tmp[2]}){
				if ($mature{$tmp[2]}=~/$id/){
					print OUT "$_\n";
					$start{$id}=$tmp[74];
					$end{$id}=$tmp[75];
				}
			}
		}
	}
	close IN;
	close OUT;
	#create the random sequences
	$outpmt=~s/\.pmt//;
	my $i;
	my @data=();
	my $mark="";
	my %fh;
	open (IN,$pmt);
	while (<IN>){
		$_=~s/\n//;
		if ($_=~/^#/){
			$title=$_;
			for ($i=1;$i<=20;$i++){
				open ($fh{"file_$i"},">./training_parameters/$outpmt\_random_$i.pmt");
				print {$fh{"file_$i"}} "$title\n";
			}
		}
		else{
			my @tmp=split("\t",$_);
			my @dat=split("_",$tmp[0]);
			my $id=$dat[0];
			print "calculating negative for $id...\n";
			#if ($id eq "aae-mir-929-2"){
			#	print "";
			#}
			if ($id ne $mark and $mark ne ""){
				for ($i=1;$i<=20;$i++){
					my $tmp=@data;
					#if the totally pmt is less than the levels
					if ($tmp<=$i){
						my @tmp=@data;
						foreach (@tmp){
							print {$fh{"file_$i"}} "$_\n";
						}
					}
					else{#the the pmt rows is more than levels, get one from random
						my %random=();
						my $j;
						my $t=0;
						for ($j=1;$j<=$i;$j++){
							getrand:
							$t++;
							my $rand=int(rand($tmp));#need to check here
							if (not exists $random{$rand}){
								$random{$rand}=1;
								print {$fh{"file_$i"}} "$data[$rand]\n";#need to check here
							}
							else{
								if ($t>500){
									die "Warning: too much runs for $outpmt at level $i, please check\n";
								}
								goto getrand;
							}
						}
					}
				}
				@data=();
			}
			if (exists $start{$id}){
				if (abs($tmp[74]-$start{$id})>5 and abs($tmp[75]-$end{$id})>5){
					push @data,$_;
				}
			}
			$mark=$id;
		}
	}
	close IN;
	for ($i=1;$i<=20;$i++){
		close $fh{"file_$i"};
	}
}
################################################################################
#model training
$out="models_$version";
if (-e $out){
	check:
	print "Warning: $out is exist, do you want to remove the folder? (yes/no)";
	my $check=<>;
	$check=~s/\n//;
	if (lc($check) eq "yes" or $check eq ""){
		print "removing the $out folder...\n";
		system "rm -rf $out";
	}
	elsif (lc($check) eq "no"){
		die "Please remove $out folder!\n";
	}
	else{
		goto check;
	}
}
system "mkdir $out";
if (not -e "training_parameters_$version"){
	die "Error: not found training_parameters_$version\n";
}
#collect the species that going to be created
print "collect the species that going to be created...\n";
my @group;
@files=glob("training_parameters_$version/experimental_*.pmt");
foreach (@files){
	next if /random/;
	my @tmp;
	open (IN,$_);
	@tmp=<IN>;
	close IN;
	if (@tmp>1){
		my $group=$_;
		$group=~s/training_parameters_$version\/experimental_//;
		$group=~s/\.pmt//;
		push @group,$group;
	}
}
#get the species belonging
print "get the species belonging...\n";
$belong{"all"}="Overall";
#start training
my $label;
#list the parameters
my %pmt;
my @tmp=pmt('mirbase');
my $tmp=0;
foreach (@tmp){
	$pmt{$_}=$tmp;
	$tmp+=1;
}
foreach (@group){
	my $group=$_;
	my @data=();
	my @svmpmt=pmt($belong{$group});
	if (@svmpmt eq 0){
		@svmpmt=pmt($belong{"all"});
	}
	#read the positive data
    open (IN,"training_parameters_$version/experimental_$group.pmt");
	while (<IN>){
		next if /^#/;
		$_=~s/\n//g;
		my @tmp=split("\t",$_);
		$label="1";
		my @svmone=();
		foreach (@svmpmt){
			my $value=nt2number($tmp[$pmt{$_}]);
			push(@svmone,$value);
		}
		my $dataset = new Algorithm::SVM::DataSet(Label => $label,Data  => [@svmone]);
		push (@data,$dataset);
	}
	close IN;
	my $i;
	for ($i=1;$i<=20;$i++){
		print "################################################################\n";
		print "training $group at level $i\n";
		print "################################################################\n";
		#read the negative data
		open (IN,"training_parameters_$version/experimental_$group\_random\_$i.pmt");
		while (<IN>){
			next if /^#/;
			$_=~s/\n//g;
			my @tmp=split("\t",$_);
			$label="0";
			my @svmone=();
			foreach (@svmpmt){
				my $value=nt2number($tmp[$pmt{$_}]);
				push(@svmone,$value);
			}
			my $dataset = new Algorithm::SVM::DataSet(Label => $label,Data  => [@svmone]);
			push (@data,$dataset);
		}
		close IN;
		my $svm = new Algorithm::SVM(Type => 'C-SVC',Kernel => 'radial');
		$svm->train(@data);
		#$accu = $svm->validate($cross);
       	#$c=$svm->coef0();
    	#$g=$svm->gamma();
        $svm->save("./models_$version/$group\_$i.model");
	}
}
sub pmt{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The pmt was used to provide pmt files
   #===========Usage: pmt(parameter)
   #			--all all parameters
   #			--range the range parameters
   #			--svmoverall the svm parameters of all species
   #			--svmoanimal the svm parameters of animal
   #			--svmplant the svm parameters of plants
   #			--svmvirus the svm parameters of all virus
   #			--display the display parameters
   #===========the parameters will be return
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my $class;
   my @parameter;
   ($class)=@_;
   if ($class eq "Viruses"){
   		$class="Metazoa";
   }
   if ($class eq "mirbase"){
	  @parameter=("miacc","miid","priacc","priid",
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
	  "mistart","miend","upperstart","upperend","stability");#79
   }
   elsif ($class eq "Overall"){
        @parameter=("unpairedrate_lowerstem","prelength","internalloopnumber_lowerstem","length_upperstem",
        "miinternalloop","firstbase","mintcontent_a","migc","pregc","prentcontent_u",
        "prentcontent_a","internalloop_topstem","preunpairedrate","mipairs","prepairs",
        "internalloopnumber_topstem","unpairedrate_topstem","mintcontent_c","mistart",
        "miunpairedrate","mintcontent_g","terminalnucleotide","prentcontent_c",
        "prentcontent_g","mintcontent_u");#25pmt
   }
#   elsif (lc($class) eq "svmanimal"){
#        @parameter=("internalloop_topstem","internalloopnumber_topstem","length_topstem","length_upperstem",
#        "migc","miinternalloop","miinternalloopnumber","mintcontent_a","mintcontent_c","mintcontent_g",
#        "mintcontent_u","mistart","miunpairedrate","penultimateposition","pregc","prelength",
#        "prentcontent_a","prentcontent_c","prentcontent_g","prentcontent_u","preunpairedrate",
#        "stability","unpairedrate_lowerstem","unpairedrate_topstem");#24pmt
#   }
   elsif ($class eq "Metazoa"){
        @parameter=("internalloop_topstem","internalloopnumber_topstem","length_topstem","length_upperstem",
        "migc","miinternalloop","miinternalloopnumber","mintcontent_a","mintcontent_c","mintcontent_g",
        "mintcontent_u","mistart","miunpairedrate","penultimateposition","pregc","prelength",
        "prentcontent_a","prentcontent_c","prentcontent_g","prentcontent_u","preunpairedrate",
        "stability","unpairedrate_topstem");#24pmt
   }
   elsif ($class eq "Viridiplantae"){
        @parameter=("firstbase","internalloop_topstem","internalloop_lowerstem","internalloopnumber_lowerstem",
        "length_upperstem","migc","migu","miinternalloop","miinternalloopnumber","mintcontent_a",
        "mintcontent_g","mintcontent_u","mipairs","penultimateposition","pregc","prentcontent_a",
        "prentcontent_c","prentcontent_g","prentcontent_u","preunpairedrate","stability",
        "unpairedrate_lowerstem","unpairedrate_topstem","upperstart");#24pmt
	}
   else{
	  @parameter=();
   }
   return @parameter;
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
