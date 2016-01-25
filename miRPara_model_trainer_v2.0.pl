#! /usr/bin/perl
#	by Yonggan Wu (weasteam@gmail.com)
#	by Simon Rayner ()
#	Version: 2.0 (2012-10-02 15:44:02 )
#	The script will prepare the data for SVM training from miRBase release
#	Changes since last version: calculate the parameters first and then group
################################################################################
use strict;
use warnings;
use File::chdir;
use Cwd;
use Algorithm::SVM;
use Algorithm::SVM::DataSet;
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
		if ($evi eq "" or $evi ne "experimental"){#if two miRNAs, will select if one one of them is from experimental
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
my %speg;
open (IN,"./miRBase/$version/organisms.txt");
while (<IN>){
	$_=~s/\n//;
	my @tmp=split("\t",$_);
	my @dat=split(";",$tmp[@tmp-1]);
	$spe{$tmp[0]}=$dat[0];
	$speg{$dat[0]}=$dat[0];
}
close IN;
print "Done\n";
#finally, output the pre-miRNA sequences
print "Outputing the experimentally varified sequences...";
if (not -e "experimental"){
	system "mkdir experimental";
}
if (not -e "./miRBase/$version/hairpin.fa"){
	system "gunzip ./miRBase/$version/hairpin.fa.gz";
}
my %fh;
my @evi=keys(%evi_group);
foreach (@evi){
	$evi=$_;
	open ($fh{$evi},">./experimental/miRBase_$version\_$evi\_overall.fasta");
}
open (IN,"./miRBase/$version/hairpin.fa");
my $spe;
while (<IN>){
	if ($_=~/^>/){
		my @tmp=split(" ",$_);
		if (exists $evi{$tmp[1]}){
			$evi=$evi{$tmp[1]};
		}
		else{
			$evi="11";
		}
	}
	if (exists $fh{$evi}){
		print {$fh{$evi}} $_;
	}
}
close IN;
foreach (@evi){
	$evi=$_;
	close $fh{$evi};
}
print "Done\n";
#goto step3;
################################################################################
#calculate the parameters
print "Calculating the parameters...";
system "perl $mirpara --pmt ./experimental/miRBase_$version\_experimental_overall.fasta";
print "Done\n";
step3:
################################################################################
#generate the training parameters
print "Generating the training parameters...\n";
my $out="training_parameters_$version";
if (not -e $out){
	system "mkdir $out";
}
#unzip the data
if (not -e "./miRBase/$version/mature.fa"){
	system "gunzip ./miRBase/$version/mature.fa.gz";
}
#get the mature sequences
print "Reading the mature sequences...\n";
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
#goto step4;
#output the sequences
print "Done\nExtracting positive parameters...\n";
my %start=();#upper start
my %end=();#upper end
my @spe=keys(%spe);
push @spe,keys(%speg);
%fh=();
foreach (@spe){
	$spe=$_;
	open ($fh{$spe},">./training_parameters_$version/$spe.pmt");
}
open ($fh{"overall"},">./training_parameters_$version/overall.pmt");
open (IN,"./experimental/miRBase_$version\_experimental_overall.pmt");
while (<IN>){
	$_=~s/\n//;
	if ($_=~/^#/){
		my $title=$_;
		foreach (@spe){
			print {$fh{$_}} "$title\n";
		}
		print {$fh{"overall"}} "$title\n";
	}
	else{
		my @tmp=split("\t",$_);
		my @dat=split("_",$tmp[0]);
		my $id=lc($dat[0]);
		@dat=split("-",$id);
		my $spe=$dat[0];
		if (exists $mature{$tmp[2]}){
			if ($mature{$tmp[2]}=~/$id/){
				print "Extracting positive parameters for $id\n";
				$start{$id}=$tmp[74];
				$end{$id}=$tmp[75];
				print {$fh{"overall"}} "$_\n";
				print {$fh{$spe}} "$_\n";
				print {$fh{$spe{$spe}}} "$_\n";
			}
		}
	}
}
close IN;
foreach (@spe){
	close $fh{$_};
}
close $fh{"overall"};
print "Done\n";
#create the random sequences
my $i;
my @data=();
my $mark="";
@spe=();
push @spe,keys(%spe);
%fh=();
foreach (@spe){
	my $spe=$_;
	open (IN,"./experimental/miRBase_$version\_experimental_overall.pmt");
	while (<IN>){
		$_=~s/\n//;
		if ($_=~/^#/){
			$title=$_;
			for ($i=1;$i<=20;$i++){
				open ($fh{"$spe\_$i"},">./training_parameters_$version/$spe\_random_$i.pmt");
				print {$fh{"$spe\_$i"}} "$title\n";
				if (not exists $fh{"overall_check_$i"}){
					open ($fh{"overall\_$i"},">./training_parameters_$version/overall\_random_$i.pmt");
					print {$fh{"overall\_$i"}} "$title\n";
					$fh{"overall_check_$i"}=1;
				}
				if (not exists $fh{"$spe{$spe}_check_$i"}){
					open ($fh{"$spe{$spe}\_$i"},">./training_parameters_$version/$spe{$spe}\_random_$i.pmt");
					print {$fh{"$spe{$spe}\_$i"}} "$title\n";
					$fh{"$spe{$spe}_check_$i"}=1;
				}
			}
		}
		else{
			my @tmp=split("\t",$_);
			my @dat=split("_",$tmp[0]);
			my $id=$dat[0];
			#print "Calculating negative for $id...\n";
			#if ($id eq "aae-mir-929-2"){
			#	print "";
			#}
			if ($id ne $mark and $mark ne ""){
				if ($mark=~/$spe/){
					print "Extracting negative parameters for $mark\n";
				}
				for ($i=1;$i<=20;$i++){
					my $tmp=@data;
					my @dat=split("-",$mark);
					my $s=$dat[0];
					#if the totally pmt is less than the levels
					if ($tmp<=$i){
						my @tmp=@data;
						foreach (@tmp){
							if ($s eq $spe){
								print {$fh{"$spe\_$i"}} "$_\n";
							}
							if ($fh{"overall_check_$i"} eq 1){
								print {$fh{"overall\_$i"}} "$_\n";
							}
							if (exists $fh{"$spe{$s}_check_$i"}){
								if ($fh{"$spe{$s}_check_$i"} eq 1){
									print {$fh{"$spe{$s}\_$i"}} "$_\n";
								}
							}
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
								if ($s eq $spe){
									print {$fh{"$spe\_$i"}} "$data[$rand]\n";#need to check here
								}
								if ($fh{"overall_check_$i"} eq 1){
									print {$fh{"overall\_$i"}} "$data[$rand]\n";
								}
								if (exists $fh{"$spe{$s}_check_$i"} ){
									if ($fh{"$spe{$s}_check_$i"} eq 1){
										print {$fh{"$spe{$s}\_$i"}} "$data[$rand]\n";
									}
								}
							}
							else{
								if ($t>500){
									die "Warning: too much runs for $s at level $i, please check\n";
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
		close $fh{"$spe\_$i"};
		if ($fh{"overall_check_$i"} eq 1){
			close $fh{"overall\_$i"};
			$fh{"overall_check_$i"}=0;
		}
		if ($fh{"$spe{$spe}_check_$i"} eq 1){
			close $fh{"$spe{$spe}\_$i"};
			$fh{"$spe{$spe}_check_$i"}=0;
		}
	}
}
#die "done for the third step";
step4:
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
@files=glob("./training_parameters_$version/*.pmt");
foreach (@files){
	next if /random/;
	my @tmp;
	open (IN,$_);
	@tmp=<IN>;
	close IN;
	if (@tmp>1){
		my $group=$_;
		$group=~s/\.\/training_parameters_$version\///;
		$group=~s/\.pmt//;
		push @group,$group;
	}
}
#get the species belonging
print "get the species belonging...\n";
$spe{"overall"}="Overall";
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
	my @svmpmt=pmt($spe{$group});
	if (@svmpmt eq 0){
		@svmpmt=pmt($spe{"Overall"});
	}
	#read the positive data
    open (IN,"training_parameters_$version/$group.pmt");
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
		open (IN,"training_parameters_$version/$group\_random\_$i.pmt");
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
