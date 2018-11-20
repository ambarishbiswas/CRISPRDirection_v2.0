#!/usr/bin/perl -W

#--------------------------------------------------------------------------------------------------------------------------------------------
#Program : 	CRISPRDirection_v1.0
#Author: 	Ambarish Biswas
#Contact:	Chris Brown [chris.brown@otago.ac.nz], Peter Fineran [peter.fineran@otago.ac.nz] or Ambarish Biswas [ambarish.biswas@otago.ac.nz]
#--------------------------------------------------------------------------------------------------------------------------------------------



use strict;

#-------------------------- change the parameters you like ----------------------------------------------------------------------------
my $check_motif_in_repeat=1; 
	my $Motif_match_score=4.50;			# Default: Score 4.5 (acts as a filter); If you choose to use this method as the other methods make the score 0.50 [i.e. PPV (1) - 0.50];
	my $motif="ATTGAAA.?";				# Default: "ATTGAAA.?" where '.?' is equal to N, another example is "A.?.?TGAAA[C|G]" which is same as ANNTGAAAC or ANNTGAAAG 

my $check_A_and_T_ratio_in_repeat=1;
	my $A_and_T_ratio_score=0.37;		# Default: PPV (0.87) - 0.50
	
my $check_similarity_with_reference_repeat=1; 
	my $Similarity_score=0.50;			# Default: Score 0.50 [i.e. PPV (1) - 0.50]; If you want to make it act as a filter use Score 4.5 or higher; 	
	my $allowed_no_of_mismatches=3;		# Default 3		#Note: We showed for up to allowed 6 bases mismatch there were no false predictions (Refer paper/supplement)
	
my $check_secondary_structure_of_repeat=1;		
	my $MFE_score=0.37;					# Default: PPV (0.87) - 0.50		
	my $MFE_cutoff=1;					# Default 1
	my $MFE_minimum_difference=0.25;	# Default 0.25
	my $MFE_exclude_bases=5;			# Default 5	

my $check_array_degeneracy=1;	 
	my $Array_degeneracy_score=0.41;	# Default: PPV (0.91) - 0.50
	my $permitted_mutation_per_array=0;	# Default 0

my $check_AT_distribution_in_flanks=1; 
	my $AT_distribution_score=0.27;		# Default: PPV (0.77) - 0.50
	my $AT_distribution_window=60;		# Default: 60
	my $AT_distribution_minimum_percentage_difference=10;	# Default: 10 


my $check_longer_leader=1;
	my $Longer_leader_score=0.18;		# Default: PPV (0.68) - 0.50



#--------------------- PILER-CR and CRT parameters-----------------------------------------------------------------------

my $piler_cr_args="";	# "-minarray 2 -minrepeat 18 -maxrepeat 67 -minspacer 18 -maxspacer 72";        		#------------ give suitable PILER-CR arguments or leave blank for defaults 
my $crt_args	="";	# "-minNR 2 -minRL 18 -maxRL 67 -minSL 18 -maxSL 72";			 		#------------ give suitable CRT arguments or leave blank for defaults	
	

#---------------------------------------------------------------------------------------


















############## do not change anyting below this, unless you know what you doing #######################################################



#---------------- check RNAfold installation ------------------------------------------------------------

my $RNAfold=`which RNAfold >&1 2>&1`; chomp $RNAfold; $RNAfold=~s/\r//g;

my $RNAfold_found="Found in $RNAfold";

if(not defined $RNAfold or $RNAfold eq "" or $RNAfold=~/:/)
{
print qq~   
   NOTE:
	RNAfold was not found in system path. Please download and install RNAfold (version 2.0.7 or higher)
	from http://www.tbi.univie.ac.at/\~ronny/RNA/index.html#download
	Latest version [V2.1.3]: [ http://www.tbi.univie.ac.at/\~ronny/RNA/packages/source/ViennaRNA-2.1.3.tar.gz ].
	
	The related method 'Analysis of RNA secondary structure' will not be available without RNAfold.
	
~;

	$check_secondary_structure_of_repeat=0; #-- skipping this method as RNAfold is not available
	$RNAfold_found="Not found.";
}
else{	
		#--- check version >= 2.0.7 ----------------------------------------
		my $RNAfold_version=`$RNAfold -V >&1 2>&1`;chomp $RNAfold_version; $RNAfold_version=~s/\r//g;
		
		$RNAfold_version=~s/RNAfold //; $RNAfold_version=~s/\.//g;
		
		if(not $RNAfold_version>=207)
			{
				$check_secondary_structure_of_repeat=0; #-- skipping this method as RNAfold is not available
				$RNAfold_found="Older version found.";
				
				print "\n    Note: It seems you have older version of RNAfold, please update your version from http://www.tbi.univie.ac.at/~ronny/RNA/index.html#download\n\n";
			}	
}


#------------ Check seqret installation----------------------------------------------------------------------------

my $seqret_path=`which seqret >&1 2>&1`; chomp $seqret_path; $seqret_path=~s/\r//g;

my $seqret_found="Found in $seqret_path";

if(not defined $seqret_path or $seqret_path eq "" or $seqret_path=~/:/)
{
print qq~   
   NOTE:
	seqret was not found in system path. Please download and install seqret [ http://emboss.sourceforge.net/download/ ]
	version EMBOSS:6.3.1 tools or higher required.	
	
	'seqret' is required for extracting genomic sequence from GBK files.
	
~;
	$seqret_found="Not found.";
}
#--------------------------------------------------------------------------------------------------------------------

#------------ Check water installation----------------------------------------------------------------------------

my $water_path=`which water >&1 2>&1`; chomp $water_path; $water_path=~s/\r//g;

my $water_found="Found in $water_path";

if(not defined $water_path or $water_path eq "" or $water_path=~/:/)
{
print qq~   
   NOTE:
	Water was not found in the system path. Please download and install water [ http://emboss.sourceforge.net/download/ ]
	version EMBOSS:6.3.1 tools or higher required.	
	
	'water' is required for comparing reference known repeats and CRISPR array repeats.
	
~;
	$water_found="Not found.";
}
#--------------------------------------------------------------------------------------------------------------------


if($RNAfold_found=~/Not/ and $water_found=~/Not/)
	{
		
print qq~		
   NOTE:
	You need to install the Emboss tools (version 6.3.1 or higher) and Vienna RNA package version 2.0.7 or higher. The other_executables folder
	has only the executables of RNAfold, water and seqret. Check each of them by the following commands from the current directory:
	
	./other_executables/RNAfold -h
	
	./other_executables/water -h
	
	./other_executables/seqret -h
	
	If you see all the programs are running successfully, then copying them to a location like '/usr/local/bin' [require root access] will solve the issue.
	Alternately export the 'other_executables' folder to your path will also work. To know how to export a directory refer this http://www.troubleshooters.com/linux/prepostpath.htm
~;		
	}


sub show_help()
	{
		my $piler_cr_found="Not found in current directory.";
		if(-e 'pilercr1.06/pilercr'){$piler_cr_found="Found in current directory.";}
		
		my $crt_found="Not found in current directory.";
		if(-e 'CRT1.2-CLI.jar'){$crt_found="Found in current directory.";}
		
		my $water_found="Not found in current directory.";
		if(-e 'water'){$water_found="Found in current directory.";}
		
		my $seqret_found="Not found in current directory.";
		if(-e 'seqret'){$seqret_found="Found in current directory.";}
		
		print qq~
CRISPRDirection Version:1.0 help:
Syntax:	
	perl CRISPRDirection.pl -f test.fa
	perl CRISPRDirection.pl -f test.fa -o output.txt

Options:
-------------------------------------------------------------------------------------------------------------------------- 
	-f		FASTA	A FASTA formatted file containing the sequence 		[e.g. -f test.fa]
	-r		FASTA	A FASTA formatted file containing only repeat(s) 	[e.g. -r test_repeats.fa]
	-g		Genbank	A .gbk file that contains sequence and annotations 	[e.g. -g NC_017625.gbk]
	-o		TEXT	A text file that will contain the output 		[e.g. -o output.txt] 
			(This option is optional, if not provided the output will be printed on standard input/output)
	-h/-help	HELP	shows this help texts	
--------------------------------------------------------------------------------------------------------------------------

CRISPRDirection dependencies:
	PILER-CR	$piler_cr_found.	
	CRT		$crt_found [Requires java runtime environment 1.5 or higher].	
	water		$water_found [Comes with EMBOSS:6.3.1 tools]
	seqret		$seqret_found [Comes with EMBOSS:6.3.1 tools]
	RNAfold		$RNAfold_found [Comes with Vienna RNA package]

NOTE: 
	Please make sure that the 'pilercr1.06/pilercr', 'CRT1.2-CLI.jar', 'RNAfold', 'water' and 'seqret' are in the current 
	directory and have execution permissions. The 'tmp' folder should have read and write permissions. An easy way to 
	do that is by issuing the command 'chmod -R 755 . && chmod 777 tmp' from the current directory.

	CRISPRDirection.pl should run under any unix based operating system that has a working 'perl' executable [comes with 
	default installations under all *nix based operating systems]. However, Mac OS users needs a working RNAfold (comes with 
	vienna RNA package).
	 
	To test different methods as specified in the literature, open the CRISPRDirection.pl program with any text editor 
	[e.g. gedit in RHEL/Fedora/CentOS, or vi in any *nix OS, or notepad in Windows OS] and change the parameters in the 
	top most section of the script. To toggle individual methods, locate the '\$check_' prefix and change the value to 
	1 (i.e. the method will be applied) or 0 (i.e. the method will not be applied). 
	
	Example:
	
		\$check_motif_in_repeat=0;	[ Default was 1, making it 0 will turn the method off.] 
		
	To change the parameter(s) of a particular method change the nested variables under that particular method.
	Example:
	
		\$check_array_degeneracy=1;	 
			\$Array_degeneracy_score=0.41;
			\$permitted_mutation_per_array=0;
	
	Making 	'\$permitted_mutation_per_array=2;' will instruct the program to allow maximum 2 bases as permitted mutations
	per CRISPR array.
	
	Finally, for version updates and bug fix refer to http://bioanalysis.otago.ac.nz/CRISPRDirection.
	 	
~;

return 1;
	}






my $idnumber=time()."_".int(rand(10000000));

my $all_gene_positions_file=$idnumber."_CDS_positions.txt";
open(WR,">tmp/$all_gene_positions_file");close(WR);

my $input_sequence_file=$idnumber."_fasta_sequence.txt";
my $output_file;
my @arr_sequence_files;
my %hash_id_lookup_table;


my $repeats_file=$idnumber."_repeats.txt";;
my $process_repeats_only=0;

my $array_seq_range="NA";

#------ step 0: get the input file ----------------------------------------
if(not defined $ARGV[0]){&show_help();exit;}

for(my $i=0;$i<=$#ARGV;$i++)
	{
		
		#--------------------------- input options ------------------------------------------------------------
		if($ARGV[$i]=~/-g/)
			{
				#----- input is a gbk file : first extract sequence from it------------------------------------
				my $filename=$ARGV[$i+1];			
				
				my $accession=`grep -w "ACCESSION  " $filename >&1`;  #---- do not change grep condition blindly, the No. of gapa are important
				$accession=~s/ACCESSION  //;chomp $accession; $accession=~s/\r//g;$accession=~s/^\s+//;
				
				
				my $defination=`grep -w "DEFINITION " $filename >&1`;
				$defination=~s/DEFINITION //; chomp $defination;$defination=~s/\r//g;$defination=~s/^\s+//;
				
				#print "\$accession=$accession\n\$defination=$defination\n";exit;
				#----------- call seqret to extract the fasta sequence-----------------------------------------
				
				unless(-e "tmp/$accession\.fna")
					{
						my $tmp_file=$idnumber."_tmp_seq.txt";
						my $ret_text=`seqret -sequence $filename -outseq tmp/$tmp_file >/dev/null 2>&1`;
						
						my $file_size=-s "tmp/$tmp_file";
						if($file_size<=0){print "Error: The GBK file does not contain any genomic sequence. Please download a full .gbk file\n"; exit;}
						#---- now open the $tmp_file and convert it to supported .fna style
						my $seq="";
						open(RD,"tmp/$tmp_file") or print "$!";
						while(my $line=<RD>)
							{
								if($line!~/>/)
									{
										chomp $line;
										$line=~s/\r//g;
										$line=uc($line);
										$seq=$seq.$line;
									}
							}
						close(RD);
						
						#----- now write the sequence to $accession.fna file
						open(WR,">tmp/$accession\.fna") or print "$!";
						print WR ">$accession\n$seq\n";
						close(WR);
						
						unlink("tmp/$tmp_file");	
					}
				push(@arr_sequence_files,$accession);
				$hash_id_lookup_table{$accession}="$accession-$defination";	
				##------------ get the accession ---------------------------------------------------------------
				#&process_fasta_file($tmp_file,\@arr_sequence_files,\%hash_id_lookup_table);		
				
						
				#----------- now get the gene positions -------------------------------------------------------				
				&process_gbk_file($filename,$all_gene_positions_file);
				#----------------------------------------------------------------------------------------------

			}
		elsif($ARGV[$i]=~/-f/)
			{
				#----- input is a fasta file ------------------------------------------------------------------
				my $tmp_file=$ARGV[$i+1];
				
				#--- now copy the sequence from $source_fasta_file to tmp/$input_sequence_file
				&process_fasta_file($tmp_file,\@arr_sequence_files,\%hash_id_lookup_table);
			}
		elsif($ARGV[$i]=~/-r/)
			{
				#----- input is a fasta repeat(s) only file ---------------------------------------------------
				my $tmp_file=$ARGV[$i+1];
				
				#--- now copy the sequence from $source_fasta_file to tmp/$input_sequence_file ----------------
				#&process_repeats_file('$tmp_file','tmp/$input_sequence_file',\%hash_of_original_identifiers);
				#----------------------------------------------------------------------------------------------
				
				system("cat $tmp_file >tmp/$repeats_file");
				
				$process_repeats_only=1;
				 
			}	
		
		
		#--------- output options -----------------------------------------------------------------------------
		elsif($ARGV[$i]=~/-o/)
			{
				#----- input is a fasta repeat(s) only file ---------------------------------------------------
				$output_file=$ARGV[$i+1];
				
				
				if(not defined $ARGV[$i+1])
					{
						print "\nWrong input. A filename is expected after -o [e.g. -o output.txt]\n\n";
						exit;
					}
				else{
						open(WR,">$output_file") or print "$!";close(WR);
					}	
			}
		
		#----- optional, will be removed later ----------------------------------------------------------------
		#elsif($ARGV[$i]=~/-array_seq_range/)
		#	{
		#		#----- input is a fasta repeat(s) only file ---------------------------------------------------
		#		$array_seq_range=$ARGV[$i+1];
		#		#open(WR,">$output_file") or print "$!";close(WR);
		#		
		#	}	
		elsif($ARGV[$i]=~/-h/ or $ARGV[$i]=~/-help/)
			{
				#----- CRISPRDirection help ---------------------------------------------------
				&show_help();
				
				exit;
			}	
	}




#--------- open library of repeats with confirmed direction -----------------------
my %lib_of_repeats_with_confirmed_direction;				
open(LIB_REP,"Ref_lib_files/verified_repeats_with_family.txt") or print "$!";
my @arr_ref_lib=<LIB_REP>;
close(LIB_REP);
				
foreach my $ref_line(@arr_ref_lib)
	{
		if($ref_line=~/^#/){next;}  # skip the comment lines if any
		
		chomp $ref_line;$ref_line=~s/\r//g;
		my($rep,$fam_type)=split('\t',$ref_line);
		my $translated_rep=$rep; $translated_rep=~tr/T/U/;
		$lib_of_repeats_with_confirmed_direction{$rep}=$fam_type;
		$lib_of_repeats_with_confirmed_direction{$translated_rep}=$fam_type;
	}








#------- Step 1: check if only repeats to process --------------------------

if($process_repeats_only==1)
	{
		
		open(RD,"tmp/$repeats_file") or print "$!";
		my @arr_rd=<RD>;
		close(RD);
		
		my $id; 
		my %hash_of_repeats;
		my $id_index=1;
		foreach my $line(@arr_rd)
			{
				chomp $line;$line=~s/\r//g;
				if($line=~/^>/){$line=~s/>//;$id=$line;}
				else{
						my $new_id="Repeat_".$id_index;
						$hash_of_repeats{$new_id}=$line;
						$hash_id_lookup_table{$new_id}=$id;
						$id_index++;
					}
			}
		
		#---------   now check the direction for each repeat and report ------------
		my $check_direction=1;
		my $no_of_repeats_to_process=keys %hash_of_repeats;
		
		foreach my $accession(sort keys %hash_of_repeats)
			{
				
				$no_of_repeats_to_process--;  
				print "Currently processing: $hash_id_lookup_table{$accession} \tRemaining: $no_of_repeats_to_process\n\n";	#next; 
		
		
		
				#---------------

				
				my $model_repeat=$hash_of_repeats{$accession};				
				
				if($check_direction==1)
					{
						
						#-------------------------------------------------------------------------------
							
						my ($matching_reference_repeat,$repeat_family,$array_direction_MEMO,$case_found,$array_direction);
						($matching_reference_repeat,$model_repeat,$array_direction,$repeat_family,$array_direction_MEMO,$case_found)=&check_repeats_direction($check_motif_in_repeat,$motif,$check_A_and_T_ratio_in_repeat,$check_similarity_with_reference_repeat,$allowed_no_of_mismatches,$check_secondary_structure_of_repeat,$MFE_cutoff,$MFE_minimum_difference,$MFE_exclude_bases,$check_array_degeneracy,$permitted_mutation_per_array,$check_AT_distribution_in_flanks,$AT_distribution_window,$AT_distribution_minimum_percentage_difference,$check_longer_leader,$Motif_match_score,$A_and_T_ratio_score,$Similarity_score,$MFE_score,$accession,$model_repeat,\%lib_of_repeats_with_confirmed_direction);
						
						#print "$matching_reference_repeat,$model_repeat,$array_direction,$repeat_family,$array_direction_MEMO,$case_found\n";

		
					
						if(defined $output_file and $output_file ne "")
							{				
								open(FPA,">>$output_file") or print "$!";
								flock(FPA,2);
								print FPA ">$hash_id_lookup_table{$accession}\n$model_repeat\n\n";								
								print FPA "# Directional analysis summary from each method: \n";
								my @arr_memo=split(';',$array_direction_MEMO);
								foreach my $memo(@arr_memo)
									{
										$memo=~s/\t+/, /g;$memo=~s/^\s+//;
										print FPA "# $memo\n";
									}				
								print FPA "\n";	
								close(FPA);
								#-----------------------------------------------------------------------------------

							}
						else{
								print ">$hash_id_lookup_table{$accession}\n$model_repeat\n\n";
								
								#print "$array_direction_MEMO\n";
								my @arr_memo2=split(';',$array_direction_MEMO);
								foreach my $memo(@arr_memo2)
									{
										$memo=~s/\t+/, /g;$memo=~s/^\s+//;
										print "# $memo\n";
									}
								print "\n\n";	
							}
					
					}
				
				
			}
		
		
		
		#------ exit ----------------------------------------------------------------
		#----- clean tmp/ folder --------------- 
		unlink<tmp/*\.mfe>;
		unlink<tmp/mr_*\.txt>;
		exit;
	}





#----- step 2: predict CRISPRs using PILER-CR and CRT ---------------------


my $all_array_file=$idnumber."_all_predicted_arrays.txt";
open(MAIN,">tmp/$all_array_file") or print "$!";

my $file_index=0;
foreach my $seq_file(@arr_sequence_files)
	{
			$file_index++;
			$seq_file=$seq_file."\.fna";
			
			#-------- run pilercr -------------------------------------------------------------------------------------------------
			my $pilercr_output_file=$idnumber."_pilercr_output_".$file_index.".txt";
			open(WR,">tmp/$pilercr_output_file") or print "$!";	close(WR); system("chmod 777 tmp/$pilercr_output_file");
			system("./pilercr1.06/pilercr -noinfo -quiet $piler_cr_args -in tmp/$seq_file -out tmp/$pilercr_output_file");
			
			#-------- run CRT -----------------------------------------------------------------------------------------------------
			my $crt_output_file=$idnumber."_crt_output_".$file_index.".txt";
			open(WR,">tmp/$crt_output_file") or print "$!";	close(WR);	system("chmod 777 tmp/$crt_output_file");	
			system("java -cp CRT1.2-CLI.jar crt $crt_args tmp/$seq_file tmp/$crt_output_file >/dev/null 2>&1");

			#-------- convert CRT to Piler-cr format and append to PILER-CR output file -------------------------------------------
			&convert_crt_arrays_to_pilercr_and_append($crt_output_file,$pilercr_output_file);
			
			
			#---- now write the output to the main compiled output file -----------------------------------------------------------	
			open(RD,"tmp/$pilercr_output_file") or print "$!";
			while( my $line=<RD> )
				{
					chomp $line;$line=~s/\r+//g;
					print MAIN "$line\n";
				}
			close(RD);	
			
			#----- delete individual output files ---------------------------------------------------------------------------------
			unlink("tmp/$pilercr_output_file");
			unlink("tmp/$crt_output_file");
	}
close(MAIN);











#----- step 3: Process the CRISPR arrays to get the longest one -----------
my %hash_of_pcr_arrays;
	
#open(RD,"compiled_pilercr_output.txt") or print "$!";
open(RD,"tmp/$all_array_file") or print "$!";
my @arr_compiled_pilercr_output=<RD>;
close(RD);


print "\nGoing to collect the arrays ...";
my @array_cpr_lines=@arr_compiled_pilercr_output;

#my %hash_of_nr_pilercr_repeats;
my %hash_of_all_pilercr_repeats;




#my $cnt_2=1;
my $cnt_2_2=1;
my $cnt_2_3=1;



my $array_index=0;	
my $total_records=0;

#open(WR,">pilercr_repeats.txt");

#while(my $line=<RD>)
for(my $i=0;$i<=$#arr_compiled_pilercr_output;$i++)
	{		
		my $line=$arr_compiled_pilercr_output[$i];
		
		if($line=~/^>/ and $arr_compiled_pilercr_output[$i+2]=~/Pos/)
			{
				
				my($accession,$repeat_sequence2,$no_of_spacers2);
				my($pcr_crispr_start,$pcr_crispr_stop,$repeat_start,$repeat_stop);
				my $array_line_start=$i;
				if($line=~/\|/)
					{
						my @tmp_arr_id=split('\|',$line);
						my @current_array;													
						$accession=$tmp_arr_id[3]; $accession=~ s/\.\d+$//;
					}
				else{
						$accession=$line; chomp $accession; $accession=~ s/\r+//g;$accession=~ s/>//g;
					}	
				
				
				$array_index++;

				
				if($arr_compiled_pilercr_output[$i+2]!~/\s+Pos\s+Repeat/){next;}
					
				
				
				my $coord_line=$arr_compiled_pilercr_output[$i+4]; 	#   1372383      26    96.2      43  ACTGGAGAGT    ......C...................    CAATTTAGAATGGCTACAAGCCGATGGTAATCAGCTAAGTCGG
				
				$coord_line=~s/^\s+//;
				
				my @arr_tmp=split('\s+',$coord_line);
				
				$pcr_crispr_start=$arr_tmp[0];
				$repeat_start=$arr_tmp[0];
				
			
				#------- get the last repeat -------------------------------------------------------------------------------------------
				my $no_of_perfect_repeats=0;
				my $j=$i+4;
				while($arr_compiled_pilercr_output[$j]!~/=======/)
					{
						#if($arr_compiled_pilercr_output[$j]=~/^\s+(\d+)\s+(\d+)\s+\S+          \S+/)#  10095660      36   100.0          TGGCAAACCA    ....................................    GTACTT
						#	{
						#		$pcr_crispr_stop=$1+$2-1;
						#		
						#	}
							
						if($arr_compiled_pilercr_output[$j+1]=~/^====/)#  10095660      36   100.0          TGGCAAACCA    ....................................    GTACTT
							{
								my $last_line=$arr_compiled_pilercr_output[$j];
								$last_line=~s/^\s+//;$last_line=~s/\s+/\t/g;
								if($last_line=~/^(\d+)\t(\d+)\t/)
									{
										$pcr_crispr_stop=$1+$2-1;
									}
									
								#--- get the repeat
								my @arr_tmp1=split('\t',$last_line);
								if($arr_tmp1[$#arr_tmp1-1]=~/\.\.\./){$arr_tmp1[$#arr_tmp1-1]=~s/\.+//g;my $c_repeat=$arr_tmp1[$#arr_tmp1-1];if($c_repeat eq ""){$no_of_perfect_repeats++;}}	
								elsif($arr_tmp1[$#arr_tmp1]=~/\.\.\./){$arr_tmp1[$#arr_tmp1]=~s/\.+//g;my $c_repeat=$arr_tmp1[$#arr_tmp1];if($c_repeat eq ""){$no_of_perfect_repeats++;}}
							}
						else{
								my $c_line=$arr_compiled_pilercr_output[$j];$c_line=~s/^\s+//;$c_line=~s/\s+/\t/g;
								my @arr_tmp1=split('\t',$c_line);
								my $c_repeat=$arr_tmp1[$#arr_tmp1-1];
								$c_repeat=~s/\.+//g;
								if($c_repeat eq ""){$no_of_perfect_repeats++;}
							}		
	
						$j++;	
					}
		
				#if($no_of_perfect_repeats<2){next;}
		
		#-------- now get the model repeat seq --------------------------
		my $pcr_crispr_range=$pcr_crispr_start."-".$pcr_crispr_stop;
		my $range="$pcr_crispr_start-$pcr_crispr_stop";
				
		my $array_line_stop=$j+1;		
				
		#----- now store the whole array in @original_array		
		
		
		my @original_pcr_array;
		
		for(my $k=$array_line_start;$k<=$array_line_stop;$k++)
			{
				#print "$arr_compiled_pilercr_output[$k]";
				chomp $array_cpr_lines[$k]; $array_cpr_lines[$k]=~ s/\r+//g;
				push(@original_pcr_array,$array_cpr_lines[$k]);
			}	
		
		my $orig_array=join("\n",@original_pcr_array);
		
		
		
		if($array_seq_range!~/NA/)
			{
				my ($filter_start,$filter_stop)=split('-',$array_seq_range);
				my $filter_mid_point=$filter_start+int(($filter_stop-$filter_start)/2);
				if($filter_mid_point>$pcr_crispr_start and $filter_mid_point<$pcr_crispr_stop)
					{
						#next;
					}
				else{
						next;
					}	
			}
		#----- remove any existing array that may be shorter than this array --------------------
		foreach my $existing_array_range(sort keys %{$hash_of_pcr_arrays{$accession}})
			{
				my($existing_array_start,$existing_array_stop)=split('-',$existing_array_range);
				if($existing_array_start>=$pcr_crispr_start and $existing_array_stop<=$pcr_crispr_stop)
					{
						delete $hash_of_pcr_arrays{$accession}{$existing_array_range};
						
						$total_records--;
					}
	
			}
			
		
		#--- now store the array --------------------------------------------------------------		
		$hash_of_pcr_arrays{$accession}{$range}=$orig_array;
		$total_records++;
	
		#print "$i - $array_line_stop: \$hash_of_pcr_arrays{$accession}{$range}\n$hash_of_pcr_arrays{$accession}{$range}\n\n";
		$i=$array_line_stop;
		#next;
		
		
		
	}
 else{next;}
	
 #	if($single_array_test>20){last;}	
}		
		
print ".. Done.\n\n";		
#exit;		
		
	

				






#----- step 4: pass the hash of CRISPR arrays to check_orientation ---------

print "Total Arrays to process: $total_records\n";

foreach my $accession(sort keys %hash_of_pcr_arrays)
 {	
	 

	#if($accession!~/NC_000913/){next;}
	
	my %hash_of_arrays_per_accession;
	my %hash_of_original_arrays_per_accession;
	#my %hash_of_questionable_arrays_per_accession;
	
	#-------------------------------------------------------------------
	foreach my $range(sort keys %{$hash_of_pcr_arrays{$accession}})
	  {
		
		my $key=$accession."-".$range;

		$total_records--;  
		print "Currently processing: $hash_id_lookup_table{$accession} Range:$range \tRemaining: $total_records\n\n";	#next; 
		
		
		#----get the reference repeat for this model_repeat
		#my $model_repeats_reference_string="NA";

		  
		my @current_array; 
		my @modified_array;
		my $model_repeat;
		#------------------------- store the whole array in @current_array after fixing the coords  ---------------------------------------------------------------------------
		my($range_start,$range_stop)=split("-",$range);
		my $array_start_position=$range_start;
		my $array_stop_position=$range_stop;
		
		$hash_of_original_arrays_per_accession{$range_start}{$range_stop}=$hash_of_pcr_arrays{$accession}{$range};
		
		
		my @original_array=split('\n',$hash_of_pcr_arrays{$accession}{$range});
		
		#-------------- get the model repeat ---------------------------------------------------------		
		my $last_line=$original_array[$#original_array];$last_line=~s/\s+/\t/g;
		my @tmp_arr_last_line=split('\t',$last_line);
		$model_repeat=$tmp_arr_last_line[$#tmp_arr_last_line];
		my $model_repeat_rev_comp=$model_repeat;$model_repeat_rev_comp=reverse $model_repeat_rev_comp;$model_repeat_rev_comp=~tr/ACGT/TGCA/;
		#---------------------------------------------------------------------------------------------
		

		 
		
		
		#$model_repeats_reference_string=~s/,/;/g;
		
		my $first_occurrence_of_gap=0;
		my $coord_diff=0;
		my $index=0;	
		my $skip_line=0;
		
		my $avg_spacer_length=0;
		my $total_spacer_length=0;
		my $no_of_spacers=0;
		
		
		
		for(my $k=0;$k<4;$k++)
			{
				chomp $original_array[$k]; $original_array[$k]=~ s/\r+//g;
				if($k<(2)){push(@current_array,$original_array[$k]);}				
				elsif($k==(2)){push(@current_array,"Position\tRepeat\tSpacer\tComment");}
				elsif($k==(3)){push(@current_array,"========\t======\t======\t=======");}
			}	
		
		for(my $k=4;$k<=$#original_array;$k++)
			{
				#print "$arr_compiled_pilercr_output[$k]";
				chomp $original_array[$k]; $original_array[$k]=~ s/\r+//g;
				#push(@original_array,$array_cpr_lines[$k]);
				
				
				#print "$array_cpr_lines[$k]\n";
				#-------------------------------------------------------------------

				#if($k<(4)){next;}
				
				#---- now make a simpler array with 4 columns ------------------------------------
				my $current_line=$original_array[$k]; 
				chomp $current_line; $current_line=~s/\r+//g;$current_line=~s/^\s+//;$current_line=~s/\s+$//;$current_line=~s/\s+/\t/g;
				
				if($current_line=~/^=/){push(@current_array,"========\t======\t======\t=======");$skip_line=1;}
				
				
				
				my @tmp_array=split('\t',$current_line);
				my($r_start,$r_length,$s_length,$r_seq,$s_seq);	
				my $insertion_base_and_position="";	
							

				
				
				
				if($skip_line==0)   #---if skip_line is not set to skip, continue
				{	#next;}	
				
						
				#print "$current_line\n";
							
				if($original_array[$k+1]!~/===/) #--- first repeat line, remember the case when there are no left flank (array starts at first base), check NC_014374
					{
						$r_start=$tmp_array[0];
						$r_seq=$tmp_array[$#tmp_array-1];
						$s_seq=$tmp_array[$#tmp_array];
						
						$total_spacer_length=$total_spacer_length+length($s_seq);
						$no_of_spacers++;
					}
				else{
						$r_start=$tmp_array[0];
						if($tmp_array[$#tmp_array-1]=~/\.\.\./)
							{
								$r_seq=$tmp_array[$#tmp_array-1];
							}
						else{
								$r_seq=$tmp_array[$#tmp_array];    # sometime pilercr doesn't give right flank 
							}	
								
						$s_seq=0;
					}	
	
				#print "\n\n\$r_start=$r_start\n";		
				#--------- check for no_of_gaps in the repeat ---------
				my $no_of_gaps=0;
					
				if($r_seq=~/-/)
					{
						$no_of_gaps= $r_seq=~ s/-/-/g;					
					}
				#---- address the coord difference
					
				if($#current_array>3)
					{
						#print "\$#current_array=$#current_array\n";	
						my $last_record=$current_array[$#current_array];		
						#print "\n\$last_record=$last_record\n";	
													
						my($last_rec_start,$last_repeat,$last_spacer,$comment)=split('\t',$last_record,);
								
						$last_repeat=~s/-//g;
						$last_spacer=~s/-//g;
								
						$r_start=$last_rec_start+length($last_repeat)+length($last_spacer);
						
	
					}
				
				#----- create a rec line and save it into current array --
				
				my $rec_line="$r_start\t$r_seq\t$s_seq\t";
				#print "$rec_line\n";
				push(@current_array,$rec_line);
				$index++;
				}
				
			} 
		
		if($no_of_spacers==0){print "@original_array\n\n";exit;}
		else{	
			  $avg_spacer_length=int($total_spacer_length/$no_of_spacers);	
		    }

				
		
		my $atleast_one_operation_performed=0;
		
		my $potential_alternate_repeat="NA";
		

		
		#---------- check direction of the array ---------------------------------------------

		my $direction_found=0;
		my $array_direction="NA";
		my $repeat_family="NA";
		my $matching_reference_repeat="NA";
		my $array_direction_MEMO="NA";
		
		my $all_gene_positions_folder="tmp";
		#my $all_gene_positions_file="all_gene_and_crispr_positions.txt";
		
		my $check_direction=1;
		
		if($check_direction==1)
			{
				
				#-------------------------------------------------------------------------------
				my $case_found=0;	
				my ($l_flank,$r_flank);
				($matching_reference_repeat,$model_repeat,$array_direction,$repeat_family,$array_direction_MEMO,$case_found)=&check_array_direction($check_motif_in_repeat,$motif,$check_A_and_T_ratio_in_repeat,$check_similarity_with_reference_repeat,$allowed_no_of_mismatches,$check_secondary_structure_of_repeat,$MFE_cutoff,$MFE_minimum_difference,$MFE_exclude_bases,$check_array_degeneracy,$permitted_mutation_per_array,$check_AT_distribution_in_flanks,$AT_distribution_window,$AT_distribution_minimum_percentage_difference,$check_longer_leader,$Motif_match_score,$A_and_T_ratio_score,$Similarity_score,$MFE_score,$Array_degeneracy_score,$AT_distribution_score,$Longer_leader_score,$array_start_position,$array_stop_position,$accession,$model_repeat,$all_gene_positions_folder,$all_gene_positions_file,\%lib_of_repeats_with_confirmed_direction,\@current_array,\@modified_array);
				
				#print "$matching_reference_repeat,$model_repeat,$array_direction,$repeat_family,$array_direction_MEMO,$case_found\n";
				undef @modified_array;	
			}
		
		
		#----- Step 5: report the orientation ------------------------------------- 

		if(defined $output_file and $output_file ne "")
			{		
				$hash_of_pcr_arrays{$accession}{$range}=~s/$accession/$hash_id_lookup_table{$accession}/gi;
						
				open(FPA,">>$output_file") or print "$!";
				flock(FPA,2);
				print FPA "$hash_of_pcr_arrays{$accession}{$range}\n";
						
				print FPA "\n";
					#print FPA "# Ref: $ref_array_direction [$model_repeats_reference_string] \t Predicted: $array_direction\n";
						
				print FPA "# Directional analysis summary from each method: \n";
				my @arr_memo=split(';',$array_direction_MEMO);
				foreach my $memo(@arr_memo)
					{
						$memo=~s/\t+/, /g;$memo=~s/^\s+//;
						print FPA "# $memo\n";
					}				
				print FPA "\n";	
				close(FPA);
				#-----------------------------------------------------------------------------------

			}
		else{
				$hash_of_pcr_arrays{$accession}{$range}=~s/$accession/$hash_id_lookup_table{$accession}/gi;
				
				print "$hash_of_pcr_arrays{$accession}{$range}\n\n";
				
				#print "$array_direction_MEMO\n";
				my @arr_memo2=split(';',$array_direction_MEMO);
				foreach my $memo(@arr_memo2)
					{
						$memo=~s/\t+/, /g;$memo=~s/^\s+//;
						print "# $memo\n";
					}
				print "\n\n";	
			}
		
		
		undef @current_array;	
		
				
	}

 }
 


#----- clean tmp/ folder --------------- 
#unlink<tmp/*\.mfe>;
unlink<tmp/mr_*\.txt>;

	







##################################################################### subs #############################################################



sub process_gbk_file()
	{
		my($input_gbk_file,$output_gene_position_file)=@_;
		
		my @arr_tmp_1=split("/",$input_gbk_file);
		
		my $accession=$arr_tmp_1[$#arr_tmp_1]; $accession=~s/\.\S+$//;
		
		#print "Accession : $accession [$input_gbk_file] \$input_gbk_file=$input_gbk_file \n";
		
		#-------------------------------------------------------------------------------
		open(RD,"$input_gbk_file") or print "$!";
		my @arr_gbk_file=<RD>;
		close(RD);
		
		#------  extract fasta seq ------------------------------------------------------
		#system("seqret -sequence $input_gbk_file -outseq $output_fasta_file");
		#------ extract the gene positions ----------------------------------------------
		open(TAB,">>tmp/$output_gene_position_file") or print "";
		flock(TAB,2);
		
		my @arr_gene_info_line=`grep -nH -w 'CDS' $input_gbk_file >&1`;
		#print "grep -nH -w '     CDS' tmp/$input_gbk_file\n @arr_gene_info_line";exit;
		
		foreach my $gene_info_line(@arr_gene_info_line)
			{
				if($gene_info_line=~/(\d+\.\.\d+)/)
					{
						my($start,$stop)=split('\.\.',$1);
						
						print TAB "$accession\tCDS\t$start\t$stop\tNA\tNA\n";
					}
			}
		
		
		
		#-------- extract Cas genes ------------------------------------------------------
		my @arr_lines=`grep -nH 'CRISPR' $input_gbk_file >&1`;
		
		my %hash_of_crispr_regions;
		foreach my $rec_line(@arr_lines)
			{
				#print "\n\n$rec_line";
				my($t_file,$line_number,$tmp_1)=split(':',$rec_line);
				
				my $cas_gene_start_stop;
				my $cas_gene_start;
				my $cas_gene_stop;
				my $cas_gene_found=0;
				my $i=$line_number-1;
				
				my $skip_this_rec=0;
				while($i>0)
					{
						my $current_line=$arr_gbk_file[$i];chomp $current_line; $current_line=~s/\r+$//;

						if($current_line=~/^     \S+\s+/)
							{
								#print "\t$i :$current_line\n";								
								#---------------------------------
								$current_line=~s/>//g;$current_line=~s/<//g;
								if($current_line=~/(\d+\.\.\d+)/)
									{
										$cas_gene_start_stop=$1;
										
										($cas_gene_start,$cas_gene_stop)=split('\.\.',$1);	
																			
										#print "\tCas gene start-stop=$cas_gene_start-$cas_gene_stop\n";
										$cas_gene_found=1;
									}
								else{
										$skip_this_rec=1;  #wrong rec, should not proceed further.
										last;
									}	
							}
							
							
						if($cas_gene_found!=0){last;}		
						$i--;
					}
				if($skip_this_rec==1){next;}
				if(not $cas_gene_start or not $cas_gene_stop){next;}
					
				#--------- now get the Cas family -----------------------------------	
				my $cas_gene_family_found=0;
				my $cas_det_line="";	
							
				my $j=$line_number-1; #---- as grep line number starts with 1 -------
				while($j<$#arr_gbk_file)
					{
						my $current_line=$arr_gbk_file[$j];chomp $current_line; $current_line=~s/^\s+//; $current_line=~s/\r+//g;
						my $next_line=$arr_gbk_file[$j+1];$next_line=~s/^\s+//;
						
						$cas_det_line=$cas_det_line.$current_line;
						
						if($next_line=~/^\//)
							{
								#$cas_det_line=$cas_det_line.$current_line;
								last;
							}
						elsif($next_line=~/^ORIGIN/)
							{
								last;
							}	
	
	
							
							
						if($cas_gene_family_found!=0){last;}		
						$j++;
					}
					
				$cas_det_line=~s/\/note=//g;$cas_det_line=~s/\"//g;$cas_det_line=~s/\'//g;		
				#print "\t$cas_det_line\n";
				
				#---- now extract the family type ----------
				my $cas_family="NA";
				my $cas_family_type="NA";
				my @arr_family_det=split(';',$cas_det_line);
				if($arr_family_det[0]=~/(C\D\D\d+)/)
					{
						$cas_family=$1; $cas_family=~s/\.\./_/g;  $cas_family=~s/^\s+//; $cas_family=~s/^,//;
						
						if($cas_family eq "")
							{
								$cas_family="NA";
							}
						#print "\tFamily:$cas_family\n";
					}	
				if(defined $arr_family_det[1] and $arr_family_det[1]=~/$cas_family(\S+)/)
					{
						$cas_family_type=$1; $cas_family_type=~s/^_//; $cas_family_type=~s/_/,/g;  $cas_family_type=~s/^\s+//; $cas_family_type=~s/^,//;
						
						if($cas_family_type eq "")
							{
								$cas_family_type="NA";
							}
						#print "\tType:$cas_family_type\n";
					}
				#print "\n";
				
				#------- write the CRISPR position, Family and type to the table
				
				#print TAB "$accession\tCRISPR\t$cas_gene_start\t$cas_gene_stop\t$cas_family\t$cas_family_type\n";
				
				#----------- now store the record in %hash_of_crispr_regions ------------------------------------
				if(defined $hash_of_crispr_regions{$cas_gene_start}{$cas_gene_stop})
					{
						my($fam,$type)=split('\t',$hash_of_crispr_regions{$cas_gene_start}{$cas_gene_stop});
						
						my $better_fam="NA"; my $better_type="NA";
						#------- check existing record ---------------------------------------------------------
						
						if($cas_family eq "NA" and $fam ne "NA"){$better_fam=$fam;}
						elsif($fam eq "NA" and $cas_family ne "NA"){$better_fam=$cas_family;}
						elsif($fam eq "NA" and $cas_family eq "NA"){$better_fam="NA";}
						
						if($cas_family_type eq "NA" and $type ne "NA"){$better_type=$type;}
						elsif($type eq "NA" and $cas_family_type ne "NA"){$better_type=$cas_family_type;}
						elsif($type eq "NA" and $cas_family_type eq "NA"){$better_type="NA";}
						
						$hash_of_crispr_regions{$cas_gene_start}{$cas_gene_stop}="$better_fam\t$better_type";
					}
				else{
						$hash_of_crispr_regions{$cas_gene_start}{$cas_gene_stop}="$cas_family\t$cas_family_type";
					}	
			}
			
		#---------------------- now print the CRISPR gene positions --------------------------------------------	
		foreach my $start_p(sort{$a<=>$b} keys %hash_of_crispr_regions)
			{
				foreach my $stop_p(sort{$a<=>$b} keys %{$hash_of_crispr_regions{$start_p}})
					{					
						print TAB "$accession\tCRISPR\t$start_p\t$stop_p\t$hash_of_crispr_regions{$start_p}{$stop_p}\n";
					}	
			}	
			
		close(TAB);
		return 1;
	}

sub process_fasta_file()
	{
		my($input_fasta_file,$arr_sequence_files,$hash_id_lookup_table)=@_;
		
		
		
		#----------- now open the user fasta file and extract the individual sequences, and pass them to first pilercr, then to CRISPRDetect -----
		open(RD,"$input_fasta_file") or print "$!: $input_fasta_file not found\n";


		my $seq_index=0;
		my $last_seq_id;
		my $seq="";
		while( my $line=<RD> )
			{				
				chomp $line; $line=~ s/\r//; $line=~ s/^\s+//; 
				#print "$line\n";
				
				if($line=~/^>/)
					{
						#print "matched\n";
						if($seq_index>0)
							{
								$seq=uc($seq);
								open(WR,">tmp/$last_seq_id\.fna") or print "$!";
								print WR ">$last_seq_id\n";
								print WR "$seq\n";
								close(WR);
								
								push(@$arr_sequence_files,$last_seq_id);
								$last_seq_id="";
								$seq="";
							}						
						
						my $o_id=$line;	chomp $o_id;$o_id=~ s/\r//; $o_id=~ s/^>//;
													#$o_id=~ s/[^a-zA-Z0-9_.-]/_/;
													
						my $time=time().int(rand(1000000000000));
					
						$last_seq_id="SEQUENCE_".$seq_index.$time;	
						$hash_id_lookup_table->{$last_seq_id}=$o_id;
						
						#print WR ">$new_id\n";
												
						$seq_index++;
					}
				else{
						#---------- get the sequence ---
						#print "matched\n";
						
						$line=uc($line);
						$seq=$seq.$line;
					}																	 								 						 
				#print WR "$line\n";
												
			 }
		close(RD); 
		#------- now write the last sequence
		$seq=uc($seq);
		
		#print "\$last_seq_id=$last_seq_id\n";
		open(WR,">tmp/$last_seq_id\.fna") or print "$!";
		print WR ">$last_seq_id\n";
		print WR "$seq\n";
		close(WR);
			 
			push(@$arr_sequence_files,$last_seq_id);
			 
			
		return 1;
	}
	
	
sub process_repeats_file()
	{
		my($input_repeat_file,$output_repeat_file,$hash_of_original_identifiers)=@_;
		
		return 1;
	}	
	
	
	
sub fill_string_with_gaps()
	{
		my($object,$length,$side)=@_;
		
		if($object eq "NA"){$object="";}
		
		my $gap_filled_string=$object;
		
		for(my $i=0;$i<=$length;$i++)
			{
				if(length($gap_filled_string)==$length){last;}
				
				if($side eq "RIGHT")
					{
						$gap_filled_string=$gap_filled_string." ";
					}
				else{
						$gap_filled_string=" ".$gap_filled_string;
					}		
			}
		
		return($gap_filled_string);
	}



sub convert_crt_arrays_to_pilercr_and_append()
	{
		my($crt_output_file,$pilercr_output_file)=@_;
		
		#print "$crt_output_file,$pilercr_output_file,$out_put_dir\n";
		my $highest_array_index=0;
		my %pilercr_array_ranges;	
		
		my $skip=1;
		
		if($skip==1)
		{	
		#----- open the PILER-CR output file and get the array index --------
		open(RD,"tmp/$pilercr_output_file") or print "$!";
		my @arr_pilercr=<RD>;
		close(RD);
		my $highest_array_index=0;
		my $acc;
		for(my $i=0;$i<$#arr_pilercr;$i++)
			{
				my $line=$arr_pilercr[$i];
				if($line=~/^Array/ and $arr_pilercr[$i+1]=~/^>/)
					{
						if($line=~/^Array (\d+)/)
							{
								if($1>$highest_array_index){$highest_array_index=$1;}
							}
						if($arr_pilercr[$i+1]=~/^>(\S+)/)
							{
								$acc=$1;
								chomp $acc; $acc=~s/\r+//g;
							}	
						#------ now get the start and stop ---------------------------
						my $j=$i+5;
						my $array_start;my $array_stop;
						while($arr_pilercr[$j]!~/====/)
							{
								if($j==$i+5)
									{
										my $first_line=$arr_pilercr[$j];$first_line=~s/^\s+//;
										if($first_line=~/^(\d+) /){$array_start=$1;}
									}
								if($arr_pilercr[$j+1]=~/====/)
									{
										#print "$array_start-$arr_pilercr[$j+1]<br>";	#exit;
										my $last_line=$arr_pilercr[$j];$last_line=~s/^\s+//;$last_line=~s/\s+/ /g;
										
										if($last_line=~/^(\d+) (\d+)/){$array_stop=$1+$2;}
										if($array_stop>$array_start)
										{
											my $range="$array_start-$array_stop";#print "$array_start-$array_stop<br>";
											$pilercr_array_ranges{$acc}{$range}=1;
										}
										last;
									}
								$j++;		
							}	
						
						#-----------------------------------------------------------------	
						$i=$j;	
					}
			}
			
		#print "\$highest_array_index=$highest_array_index<BR>";#exit;	
		}
		
		
		#---- now read the arrays predicted by CRT and convert them to piler-CR format --------------
		open(CRT,"tmp/$crt_output_file") or print "$!";
		my @arr_crt=<CRT>;
		close(CRT);
		
		my $accession="";
		for(my $i=0;$i<=$#arr_crt;$i++)
			{
				my $line=$arr_crt[$i]; chomp $line; $line=~s/\r+//g;
				
				if($line=~/ORGANISM:  (\S+)/)
					{
						$accession=$1;
						$highest_array_index=0;
						#print "\$accession=$accession\t\$i=$i\n";
					}
				
				elsif($arr_crt[$i]=~/^CRISPR/ and $arr_crt[$i]=~/Range:/)
					{
						#print "\$accession=$accession\t$line\t\$highest_array_index=$highest_array_index\n";#next;
						my $range;my $range_center; my $crt_array_length;
						
						if($arr_crt[$i]=~/Range: (\d+) - (\d+)/)
							{
								$range=$1."-".$2;
								$range_center=$1+int(($2-$1)/2);
								$crt_array_length=$2-$1+1;
							}	
						#--------- check if the array center belong to any piler-cr array already ----
						my $array_exist=0;
						if(defined $range)
							{
								foreach my $pilercr_range(keys %{$pilercr_array_ranges{$accession}})
									{
										#print "$pilercr_range\t$range_center<br>";
										my($r_start,$r_stop)=split('-',$pilercr_range);
										my $pilercr_array_length=$r_stop-$r_start+1;
										if($crt_array_length<$pilercr_array_length){$array_exist=1;last;}
										#if($range_center>=$r_start and $range_center<=$r_stop){$array_exist=1;last;}
									}
							}
						if($array_exist==1){next;}	#exit;	
						#-------- extract all the repeats, and get the model repeat ------------------
						my $j=$i+3;
						my $c_line="";
						my @arr_repeats;
						my %hash_of_repeats;
						while($arr_crt[$j]!~/^---/)
							{
								#print "\$arr_crt[$j]=$arr_crt[$j]";
								my @arr_line=split('\t+',$arr_crt[$j]);
								
								push(@arr_repeats,$arr_line[1]);
								
								if($hash_of_repeats{$arr_line[1]})
									{
										$hash_of_repeats{$arr_line[1]}=$hash_of_repeats{$arr_line[1]}+1;
									}
								else{
										$hash_of_repeats{$arr_line[1]}=1;
									}
								
								$j++;
							}
							


						
						
						
						#---------- now get the model repeat ----------------------------------------
						my $model_repeat="";
						my @arr_of_model_repeats;
						my $heighest_repeat_occurrence_score=0;
						foreach my $repeat(sort{$hash_of_repeats{$b}<=>$hash_of_repeats{$a}} keys %hash_of_repeats)
							{
								if(not $arr_of_model_repeats[0])
									{
											#push(@arr_of_model_repeats,$repeat);
										$arr_of_model_repeats[0]=$repeat;
										$heighest_repeat_occurrence_score=$hash_of_repeats{$repeat};
										$model_repeat=$repeat;
									}
								elsif($hash_of_repeats{$repeat}>=$heighest_repeat_occurrence_score)
									{
										push(@arr_of_model_repeats,$repeat);
										$heighest_repeat_occurrence_score=$hash_of_repeats{$repeat};
										$model_repeat=$repeat;
									}
							}
							
						#-------------- skip this array if no_of_repeats are less than 3-------------	
						#if($#arr_repeats<2){next;}
						
						#---- now skip the arrays below a certain degenaracy cutoff -----------------
						
							#--- check if there are at least 2 repeats with no degeneracy -----------
	
						#if($heighest_repeat_occurrence_score<2){next;}

						#--------- check if there is just one repeat in the array
						
						if($#arr_of_model_repeats>=0)
							{
								#----- now select the best model repeat by scoring the array
								my %tmp_hash1;
								foreach my $repeat(@arr_of_model_repeats)
									{

										my $array_degeneracy_score=&get_array_degeneracy_score($model_repeat,$repeat,\@arr_repeats);
										
										#print "\$array_degeneracy_score=$array_degeneracy_score\n";
										
										$tmp_hash1{$repeat}=$array_degeneracy_score;
									}
								
								foreach my $repeat(sort{$tmp_hash1{$b}<=>$tmp_hash1{$a}}keys %tmp_hash1)
									{
										#print "$repeat\t$tmp_hash1{$repeat}\n";
										
										$model_repeat=$repeat;
										last;
									}												
							}
						else{
								$model_repeat=$arr_of_model_repeats[0];
							}
						
						
						#----
						$model_repeat=~s/\s+$//;
						#print "\$model_repeat=$model_repeat\n";
					
						
						#----- now that the model repeat is found, convert the whole array on the go
						open(WR,">>tmp/$pilercr_output_file") or print "$!";
						$highest_array_index++;
						print WR "\n";
						print WR "Array $highest_array_index\n";
						print WR ">$accession\n";
						print WR "\n";
						print WR "       Pos  Repeat     %id  Spacer  Left flank    Repeat                               Spacer\n";
						
						my $tmp_equal_sign="";
						for(my $m=0;$m<length($model_repeat);$m++)
							{
								$tmp_equal_sign=$tmp_equal_sign."=";
							}
						print WR "==========  ======  ======  ======  ==========    $tmp_equal_sign    ======\n";
						
						
						my $k=$i+3;
						while($arr_crt[$k]!~/^---/)
							{
								my @arr_line=split('\t+',$arr_crt[$k]);
								
								my $start=$arr_line[0]; 										my $gapped_start=&fill_string_with_gaps($start,10,"LEFT");
								my $repeat_len=length($arr_line[1]);							my $gapped_repeat_len=&fill_string_with_gaps($repeat_len,6,"LEFT");
								my $percent_id="100.0";											my $gapped_percent_id=&fill_string_with_gaps($percent_id,6,"LEFT");
								my $spacer_len=length($arr_line[2]);							my $gapped_spacer_len=&fill_string_with_gaps($spacer_len,6,"LEFT");
								my $left_flank="NNNNNNNNNN";									my $gapped_left_flank=&fill_string_with_gaps($left_flank,10,"LEFT");
								my $repeat=&change_bases_to_dots($arr_line[1],$model_repeat);	my $gapped_repeat=&fill_string_with_gaps($repeat,length($model_repeat),"RIGHT");
								my $spacer="";
								if($arr_crt[$k+1]!~/^---/)
									{
										$spacer=$arr_line[2];
									}
								else{
										$spacer="NNNNNNNNNN";
										$spacer_len=0;
																								$gapped_spacer_len=&fill_string_with_gaps($spacer_len,6,"LEFT");
									}						
																								my $gapped_spacer=&fill_string_with_gaps($spacer,10,"RIGHT");			
								
								
								print WR "$gapped_start  $gapped_repeat_len  $gapped_percent_id  $gapped_spacer_len  $gapped_left_flank    $gapped_repeat    $gapped_spacer\n";
								
								$k++;
							}
						
						#---- now get the last line ----
						my $no_of_repeats;
						my $avg_repeat_len;
						my $avg_p_id=" ";
						my $avg_spacer_len;
						my $avg_left_flank=" ";
						if($arr_crt[$k+1]=~/Repeats:/)
							{
								my @arr_tmp_line=split('\t+',$arr_crt[$k+1]);
								if($arr_tmp_line[0]=~/Repeats: (\d+)/){$no_of_repeats=$1;}else{$no_of_repeats=0;}
								if($arr_tmp_line[1]=~/Length: (\d+)/){$avg_repeat_len=$1;}else{$avg_repeat_len=0;}
								if($arr_tmp_line[2]=~/Length: (\d+)/){$avg_spacer_len=$1;}else{$avg_spacer_len=0;}
							}
						
						my $gapped_no_of_repeats=&fill_string_with_gaps($no_of_repeats,10,"LEFT");	
						my $gapped_avg_repeat_len=&fill_string_with_gaps($avg_repeat_len,6,"LEFT");
						my $gapped_avg_p_id=&fill_string_with_gaps($avg_p_id,6,"LEFT");
						my $gapped_avg_spacer_len=&fill_string_with_gaps($avg_spacer_len,6,"LEFT");
						my $gapped_avg_left_flank=&fill_string_with_gaps($avg_left_flank,10,"LEFT");
						
						print WR "==========  ======  ======  ======  ==========    $tmp_equal_sign\n";						
						print WR "$gapped_no_of_repeats  $gapped_avg_repeat_len  $gapped_avg_p_id  $gapped_avg_spacer_len  $gapped_avg_left_flank    $model_repeat\n";
						print WR "\n\n";
						#print WR "SUMMARY BY SIMILARITY\n";
						close(WR);
						#--- end of current array ----
						$i=$k;
						#exit;
					}
			}
		
		
		return 1;
	}



	
	
sub get_array_degeneracy_score()
	{
		my($model_repeat,$current_repeat,$ref_arr_repeats)=@_;
		
		if(not $model_repeat)
			{
				#print "model_repeat missing\n";
				return;
			}
		#print "$model_repeat,$current_repeat,$ref_arr_repeats\n";
		
		my $arr_degen_score=0;
		
		foreach my $repeat(@{$ref_arr_repeats})
			{
				#print "\$repeat=$repeat\n";
				my $new_r_seq=&change_dots_to_bases($repeat,$model_repeat);
				
				my $sim_score=&get_similarity_score($new_r_seq,$current_repeat);
				
				my $neg_score=length($current_repeat)-$sim_score;
				
				$arr_degen_score=$arr_degen_score-$neg_score;
			}
			
		return($arr_degen_score);	
	}
	
sub change_dots_to_bases()
	{
		my($r_string,$mr_string)=@_;
		
		if(not $r_string){return("");}
		#print "\$bf_string,\$mr_string $bf_string,$mr_string\n";
		my $return_string;
		
		my @arr_1=split('',$r_string);
		my @arr_2=split('',$mr_string);
		
		for(my $i=0;$i<=$#arr_2;$i++)
			{
				if($arr_1[$i] eq ".")
					{
						$return_string=$return_string.$arr_2[$i];
					}
				else{
						$return_string=$return_string.$arr_1[$i];
					}	
			}
		
		return($return_string);
	}


sub get_similarity_score()
	{
		my($current_repeat,$model_repeat)=@_;
		
			
		my $similarity_score=0;
		
		
		
		my @arr_1=split('',$current_repeat);
		my @arr_2=split('',$model_repeat);
		
		
		
		for(my $i=0;$i<=$#arr_2;$i++)   
			{
				if($arr_1[$i] eq $arr_2[$i])
					{
						$similarity_score=$similarity_score+1;	
					}	
			}
		
		return $similarity_score;
	}

sub change_bases_to_dots()
	{
		my($bf_string,$mr_string)=@_;
		
		if(not $bf_string){return("");}
		#print "\$bf_string,\$mr_string $bf_string,$mr_string\n";
		my $return_string;
		
		my @arr_1=split('',$bf_string);
		my @arr_2=split('',$mr_string);
		
		for(my $i=0;$i<=$#arr_2;$i++)
			{
				if($arr_1[$i] eq $arr_2[$i])
					{
						$return_string=$return_string.".";
					}
				else{
						$return_string=$return_string.$arr_1[$i];
					}	
			}
		
		return($return_string);
	}



sub get_matching_reference_repeat_and_direction()
	{
		
		my($accession,$model_repeat,$allowed_no_of_mismatches,$lib_of_repeats_with_confirmed_direction)=@_;
		
		#print "\nInside sub get-matching_reference_repeat_and_direction with $accession,$model_repeat :\n";
		
		my $matching_reference_repeat="NA";
		my $array_direction="";
		my $observed_percent_similarity=0;
		
		my $case_found=0;
		my $match_found=0;
		my $ref_repeat_family="NA";
	
		my $allowed_percent_similarity=100-int(($allowed_no_of_mismatches*100)/length($model_repeat)); 				#----- here it has nothing to do with 95% similarity
		my $minimum_length_distribution=100-int(($allowed_no_of_mismatches*100)/length($model_repeat));

		#------------ first check in the lib of repeats with the model_repeat 
		foreach my $ref_repeat(keys %{$lib_of_repeats_with_confirmed_direction})
			{
				$ref_repeat=~tr/U/T/;

				
						my $model_repeat_rc=$model_repeat; $model_repeat_rc=reverse $model_repeat_rc; $model_repeat_rc=~tr/ACGT/TGCA/;
						
						
						
						if($ref_repeat=~/$model_repeat/ or $model_repeat=~/$ref_repeat/)
							{
								if($ref_repeat=~/$model_repeat/)
									{
										$observed_percent_similarity=int((length($model_repeat)/int(length($ref_repeat)))*100);
									}
								elsif($model_repeat=~/$ref_repeat/)
									{
										$observed_percent_similarity=int((length($ref_repeat)/int(length($model_repeat)))*100);
									}	
								
								if($observed_percent_similarity<80){next;}
									
								$array_direction="F";
								$ref_repeat_family=$lib_of_repeats_with_confirmed_direction->{$ref_repeat};
								
								$match_found=1;
								$case_found=1;
								$matching_reference_repeat=$ref_repeat;
								#print "\n\nForward: [$accession] The $model_repeat matched with $ref_repeat and belongs to group: $lib_of_repeats_with_confirmed_direction->{$ref_repeat} .\n";
								
								#$observed_percent_similarity=100;
								last;
							}
						elsif($ref_repeat=~/$model_repeat_rc/ or $model_repeat_rc=~/$ref_repeat/)
							{
								
								if($ref_repeat=~/$model_repeat_rc/)
									{
										$observed_percent_similarity=int((length($model_repeat_rc)/int(length($ref_repeat)))*100);
									}
								elsif($model_repeat_rc=~/$ref_repeat/)
									{
										$observed_percent_similarity=int((length($ref_repeat)/int(length($model_repeat_rc)))*100);
									}	
								
								if($observed_percent_similarity<80){next;}
								
								$array_direction="R";
								$ref_repeat_family=$lib_of_repeats_with_confirmed_direction->{$ref_repeat};
								
								$match_found=1;
								$case_found=1;
								$matching_reference_repeat=$ref_repeat;
								
								#print "\n\nReverse: [$accession] The $model_repeat_rc matched with $ref_repeat and belongs to group: $lib_of_repeats_with_confirmed_direction->{$ref_repeat} .\n";	
								
								#$observed_percent_similarity=100;					
								last;
							}	
							
			}

		if($match_found==0)
			{
				foreach my $ref_repeat(sort{length($lib_of_repeats_with_confirmed_direction->{$b})<=>length($lib_of_repeats_with_confirmed_direction->{$a})}keys %{$lib_of_repeats_with_confirmed_direction})
					{
						$ref_repeat=~tr/U/T/;
						
						my $top_line="";
						my $bottom_line="";	
				
				
						my $model_repeat_rc=$model_repeat; $model_repeat_rc=reverse $model_repeat_rc; $model_repeat_rc=~tr/ACGT/TGCA/;	
									
						
									
						#while(not $best_bottom_line)
						#{
						
						my $time=time();	
						$time=$time.int(rand(1000000000));	
						my $file_seq_1=$time.$accession."_".$allowed_percent_similarity."model_repeat_bases.txt";				
						my $file_seq_2=$time.$accession."_".$allowed_percent_similarity."repeat_bases.txt";
						my $file_seq_3=$time.$accession."_".$allowed_percent_similarity."rev_model_repeat_bases.txt";
						
						my $outfile1=$time.$accession."_".$allowed_percent_similarity."output1.txt";
						my $outfile2=$time.$accession."_".$allowed_percent_similarity."output2.txt";
										
						system("echo '$model_repeat' >tmp/$file_seq_1");
						system("echo '$ref_repeat' >tmp/$file_seq_2");					
						system("echo '$model_repeat_rc' >tmp/$file_seq_3");
										
										
						my $ret_msg1=`water -asequence tmp/$file_seq_1 -bsequence tmp/$file_seq_2 -gapopen 10 -gapextend 9 -outfile tmp/$outfile1 >/dev/null 2>&1`;
						my $ret_msg2=`water -asequence tmp/$file_seq_3 -bsequence tmp/$file_seq_2 -gapopen 10 -gapextend 9 -outfile tmp/$outfile2 >/dev/null 2>&1`;				
										#---- now open the output file and get the alignment, store the alignment in a compiled file
						
						
						
						#------------------- check forward orientation ---------------------------------------------------------------				
						open(RD1,"tmp/$outfile1") or print "$!";
						my @arr_rd1=<RD1>;
						close(RD1);
						
						my $tl_count1=0;
						my $bl_count1=0;
						
						foreach my $line(@arr_rd1)
							{
												#print $line,"\n";
								chomp $line;$line=~s/\r//g;
								if(not $line or $line=~/#/){next;}
								elsif($line=~/\d+/)
									{
										#print "$line\n";
										$line=~s/^\s+//;$line=~s/\s+/\t/g;
										my($start,$seq,$stop)=split('\t',$line);
										if($tl_count1==0){$top_line=$line;$tl_count1++;}
										elsif($bl_count1==0){$bottom_line=$line;$bl_count1++;}
									}						
							}
										
										#print "\n";	
						unlink("tmp/$file_seq_1");
						unlink("tmp/$file_seq_2");
						unlink("tmp/$file_seq_3");															
						
						unlink("tmp/$outfile1");							
										
						#print "$top_line\n$bottom_line\n\n";	
						
						my($top_start1,$top_seq1,$top_stop1)=split('\t',$top_line);
						#	$top_start1=$top_start1-1;
						#	$top_stop1=$top_stop1-1;
							
						my($bottom_start1,$bottom_seq1,$bottom_stop1)=split('\t',$bottom_line);
						#	$bottom_start1=$bottom_start1-1;
						#	$bottom_stop1=$bottom_stop1-1;
						
						if($top_seq1=~/-/ or $bottom_seq1=~/-/){next;}
						
						my $similarity_score1=&get_similarity_score($top_seq1,$bottom_seq1);
						
						#$observed_percent_similarity=$similarity_score1;
						if(length($ref_repeat)>length($model_repeat))
							{
								$observed_percent_similarity=int(($similarity_score1/int(length($ref_repeat)))*100);
							}
						else{
								$observed_percent_similarity=int(($similarity_score1/int(length($model_repeat)))*100);
							}
						#$observed_percent_similarity=int(($similarity_score1/int(length($model_repeat)))*100);
						
						#print "$ref_repeat ($observed_percent_similarity): $similarity_score1>length($model_repeat)*.95",$similarity_score1>length($model_repeat)*.95,"\n";#exit;
						#if((length($top_seq1)>=length($model_repeat)*0.95 or length($model_repeat)>=length($top_seq1)) and (length($ref_repeat)-$similarity_score1)<=$allowed_no_of_mismatches and $observed_percent_similarity>=80)#>=int(length($ref_repeat)*.95))
						if((length($top_seq1)>=length($model_repeat)*($minimum_length_distribution/100) or length($model_repeat)>=length($top_seq1)) and (length($ref_repeat)-$similarity_score1)<=$allowed_no_of_mismatches and $observed_percent_similarity>=$allowed_percent_similarity)#>=int(length($ref_repeat)*.95))
							{
								$array_direction="F";
								$match_found=1;	
								$case_found=1;
								$matching_reference_repeat=$ref_repeat;
								$ref_repeat_family=$lib_of_repeats_with_confirmed_direction->{$ref_repeat};
								
								
								
								#$observed_percent_similarity=$similarity_score1;
								#$observed_percent_similarity=int(($similarity_score1/int(length($model_repeat)))*100);
								
								#print "\n\nQ: Forward: [$accession] The $model_repeat matched with $ref_repeat and belongs to group: $lib_of_repeats_with_confirmed_direction->{$ref_repeat} .\n";
								#print "\$observed_percent_similarity=$observed_percent_similarity\n";
								#print "\t$top_line\n\t$bottom_line\n\n";
								
								
							}
						else{
								#unlink("tmp/$outfile1");
							}	
						
						if($match_found==1){last;}
						
						#-------------- now check reverse orientation ---------------------------------------------
						
						my $top_line2="";
						my $bottom_line2="";
						
						open(RD2,"tmp/$outfile2") or print "$!";
						my @arr_rd2=<RD2>;
						close(RD2);
						
						my $tl_count2=0;
						my $bl_count2=0;
						
						foreach my $line(@arr_rd2)
							{
												#print $line,"\n";
								chomp $line;$line=~s/\r//g;
								if(not $line or $line=~/#/){next;}
								elsif($line=~/\d+/)
									{
										#print "$line\n";
										$line=~s/^\s+//;$line=~s/\s+/\t/g;
										my($start,$seq,$stop)=split('\t',$line);
										if($tl_count2==0){$top_line2=$line;$tl_count2++;}
										elsif($bl_count2==0){$bottom_line2=$line;$bl_count2++;}
									}						
							}
										
										#print "\n";						
						unlink("tmp/$outfile2");							
										
						#print "$top_line\n$bottom_line\n\n";	
						
						my($top_start2,$top_seq2,$top_stop2)=split('\t',$top_line2);
						#	$top_start1=$top_start1-1;
						#	$top_stop1=$top_stop1-1;
							
						my($bottom_start2,$bottom_seq2,$bottom_stop2)=split('\t',$bottom_line2);
						#	$bottom_start1=$bottom_start1-1;
						#	$bottom_stop1=$bottom_stop1-1;
						
						
						if($top_seq2=~/-/ or $bottom_seq2=~/-/){next;}
						
						my $similarity_score2=&get_similarity_score($top_seq2,$bottom_seq2);
						#$observed_percent_similarity=$similarity_score2;
						if(length($ref_repeat)>length($model_repeat))
							{
								$observed_percent_similarity=int(($similarity_score2/int(length($ref_repeat)))*100);
							}
						else{
								$observed_percent_similarity=int(($similarity_score2/int(length($model_repeat)))*100);
							}		
						
						#if((length($top_seq2)>length($model_repeat)*0.95 or length($model_repeat)>=length($top_seq2))and (length($ref_repeat)-$similarity_score2)<=$allowed_no_of_mismatches and $observed_percent_similarity>=80)#$similarity_score2>=int(length($ref_repeat)*.95))
						if((length($top_seq2)>length($model_repeat)*($minimum_length_distribution/100) or length($model_repeat)>=length($top_seq2))and (length($ref_repeat)-$similarity_score2)<=$allowed_no_of_mismatches and $observed_percent_similarity>=$allowed_percent_similarity)#$similarity_score2>=int(length($ref_repeat)*.95))
							{					
								
								$array_direction="R";
								$match_found=1;	
								$case_found=1;
								$matching_reference_repeat=$ref_repeat;		
								$ref_repeat_family=$lib_of_repeats_with_confirmed_direction->{$ref_repeat};	
										
								
								
								
								#print "\n\nQ: Reverse: [$accession] The $model_repeat matched with $model_repeat_rc and belongs to group: $lib_of_repeats_with_confirmed_direction->{$ref_repeat} .\n";								
								#print "\$observed_percent_similarity=$observed_percent_similarity\n";
								#print "\t$top_line2\n\t$bottom_line2\n\n";					
								
							}
						else{
								#unlink("tmp/$outfile2");
							}
						if($match_found==1){last;}
					}
		
			}
			
			
		return($matching_reference_repeat,$ref_repeat_family,$array_direction,$observed_percent_similarity);	
	}

	
sub get_confidence_score()
	{
		my($total_score,$achieved_score)=@_;
		
		my $confidence="NA";
		
		#-- check if the obtained score is >60% of the total_score -----
		if($achieved_score>($total_score*0.66) and $achieved_score>0.5)   #--0.5 is minimum cutoff to assign any prediction HIGH confidence
			{
				$confidence="HIGH";
			}
		elsif($achieved_score>($total_score*0.33))
			{
				$confidence="MEDIUM";
			}
		else{
				$confidence="LOW";
			}		
		
		return($confidence);
	}	
	
	

sub check_array_direction()
	{
		my($check_motif_in_repeat,$motif,$check_A_and_T_ratio_in_repeat,$check_similarity_with_reference_repeat,$allowed_no_of_mismatches,$check_secondary_structure_of_repeat,$MFE_cutoff,$MFE_minimum_difference,$MFE_exclude_bases,$check_array_degeneracy,$permitted_mutation_per_array,$check_AT_distribution_in_flanks,$AT_distribution_window,$AT_distribution_minimum_percentage_difference,$check_longer_leader,$Motif_match_score,$A_and_T_ratio_score,$Similarity_score,$MFE_score,$Array_degeneracy_score,$AT_distribution_score,$Longer_leader_score,$array_start_position,$array_stop_position,$accession,$model_repeat,$all_gene_positions_folder,$all_gene_positions_file,$lib_of_repeats_with_confirmed_direction,$current_array,$modified_array)=@_;
		#print "\tGoing to check for consensus sequence with \$model_repeat=$model_repeat and ref:$model_repeats_reference_string for $accession:\n\n";
		
		
		
		my $array_direction="NA";
		my $repeat_family="NA";
		my $matching_reference_repeat="NA";
		my $match_found=0;
		my $observed_percent_similarity=0;
			
		#($matching_reference_repeat,$repeat_family,$array_direction,$observed_percent_similarity)=split(';',$model_repeats_reference_string);		
		
		
		

		
		
		
		
		#-------------------------- some common parameters -----------------------------------------------------			
		my $case_found=0;
		my $new_model_repeat="";
		my $potential_alternate_repeat="NA";
		my $array_direction_MEMO="";
		

		
		
		#---------- now open the sequence file and get the sequence string -------------------------------------------------
		open(SEQ,"tmp/$accession\.fna") or print "$!";
		my @arr_seq=<SEQ>;
		close(SEQ);
		my $species_seq=$arr_seq[1]; chomp $species_seq;$species_seq=~s/\r//g;
		
		
		#--------------------------------------------- first check the model_repeats for matching to ref repeat, motif and A_and_T count ----------------------------------------------------------
		
		
		#---- check motif -------------------------------------------------
		if($check_motif_in_repeat==1)
			{					
				
				my $model_repeat_rc=$model_repeat; $model_repeat_rc= reverse $model_repeat_rc; $model_repeat_rc=~tr/ACGT/TGCA/;
				
				my $suggested_direction="NA";
										
				if($model_repeat=~/$motif/ and $model_repeat_rc=~/$motif/){$suggested_direction="NA";}
				elsif($model_repeat=~/$motif/ ){$suggested_direction="F";}
				elsif($model_repeat_rc=~/$motif/){$suggested_direction="R";}

				#-----  change all .? to N
				$motif=~s/\.\?/\(N\)/g;
				if($suggested_direction!~/NA/)
					{
						$array_direction_MEMO=$array_direction_MEMO." Motif $motif match prediction:         $suggested_direction Score: $Motif_match_score/$Motif_match_score;";
					}
				else{
						$array_direction_MEMO=$array_direction_MEMO." Motif $motif match prediction:         $suggested_direction Score: 0/$Motif_match_score;";
					}		
				
			}
		else{
				#-----  change all .? to N
				$motif=~s/\.\?/\(N\)/g;
				$array_direction_MEMO=$array_direction_MEMO." Motif $motif matched prediction:       NA ;";
			}	


		#---- check A&T -------------------------------------------------
		if($check_A_and_T_ratio_in_repeat==1)
			{				
				#--- replace all U with T if present -----------
				$model_repeat=~tr/U/T/;

				##################### check direction using no. of As and Ts: mainly applicable to Archaeal genomes #######################		
												
				my $no_of_As_in_repeat=$model_repeat=~s/A/A/gi;
				my $no_of_Ts_in_repeat=$model_repeat=~s/T/T/gi;	
				my $A_and_T_suggestion="NA";
											
				my $at_p_d=(($no_of_As_in_repeat+$no_of_Ts_in_repeat)/length($model_repeat))*100; $at_p_d=sprintf("%.2f",$at_p_d);
				my $gc_p_d=((length($model_repeat) -$no_of_As_in_repeat-$no_of_Ts_in_repeat)/length($model_repeat))*100; $gc_p_d=sprintf("%.2f",$gc_p_d);
				
				#---- skip if AT% is higher than GC% in the repeat -----------------
				if($at_p_d>$gc_p_d)
					{
						$array_direction_MEMO=$array_direction_MEMO." A,T distribution in repeat prediction:     NA [Repeat is AT rich:$at_p_d%AT];";
					}
				else{	
						#------------------------------------------------------------------
											
						if($no_of_As_in_repeat>$no_of_Ts_in_repeat)
							{
								$A_and_T_suggestion="F";
							}
						elsif($no_of_Ts_in_repeat>$no_of_As_in_repeat)
							{
								$A_and_T_suggestion="R";
							}							

				
						$array_direction_MEMO=$array_direction_MEMO." A,T distribution in repeat prediction:     $A_and_T_suggestion [$no_of_As_in_repeat,$no_of_Ts_in_repeat] Score: $A_and_T_ratio_score/$A_and_T_ratio_score;";
					}
			}
		else{
						$array_direction_MEMO=$array_direction_MEMO." A,T distribution in repeat prediction:     NA ;";
			}			
		
		
		#------- check reference-repeats ---------------------------------------------
		if($check_similarity_with_reference_repeat==1)
			{
				
				($matching_reference_repeat,$repeat_family,$array_direction,$observed_percent_similarity)=&get_matching_reference_repeat_and_direction($accession,$model_repeat,$allowed_no_of_mismatches,$lib_of_repeats_with_confirmed_direction);
				
				#print "$matching_reference_repeat,$repeat_family,$array_direction,$observed_percent_similarity\n";
				#----------------------------append the suggestion in MEMO string -------------------------------------------	
				if($matching_reference_repeat!~/NA/)
					{
						$array_direction_MEMO=$array_direction_MEMO." Reference repeat match prediction:         $array_direction [matched $matching_reference_repeat with $observed_percent_similarity\% identity] Score: $Similarity_score/$Similarity_score;";
					}
				else{
						$array_direction_MEMO=$array_direction_MEMO." Reference repeat match prediction:         NA ;";
					}		
				#-------------------------------------------now check for the AT richness/longer leader/degeneracy -------------------------------
		
			}
		else{
						$array_direction_MEMO=$array_direction_MEMO." Reference repeat match prediction:         NA ;";
			}
		
						
		
	
		#------- check MFE ---------------------------------------------
		if($check_secondary_structure_of_repeat==1)
			{		

				my $forward_strand=$model_repeat;
				my $MFE_suggestion_memo="";
					
				if($MFE_exclude_bases>0)
					{
						my $string_N="";
						for(1..$MFE_exclude_bases){$string_N=$string_N."N";}
						
						$forward_strand=~s/^\S{$MFE_exclude_bases}/$string_N/;
						$forward_strand=~s/\S{$MFE_exclude_bases}$/$string_N/;						
					}	
				my $other_strand=$forward_strand; $other_strand=reverse $other_strand; $other_strand=~tr/ACGT/TGCA/;	

				
				
				my $mfe_1=&run_rnafold($forward_strand);				
				my $mfe_2=&run_rnafold($other_strand);
				

				#print "\n";
				
				my $MFE_suggestion="NA";	
				
				#---- check if meets minimum valid threshhold ---
				my $valid_case=0;
				if( (abs( abs($mfe_1) - abs($mfe_2) )>=$MFE_minimum_difference) and ( abs($mfe_1)>=$MFE_cutoff or abs($mfe_2)>=$MFE_cutoff))
					{
						$valid_case=1;
					}				
				
				if($mfe_1<$mfe_2 and $valid_case==1){$MFE_suggestion="F";}
				elsif($mfe_2<$mfe_1 and $valid_case==1){$MFE_suggestion="R";}
				
								
				
				if($MFE_suggestion!~/NA/)
					{
						$MFE_suggestion_memo=$MFE_suggestion_memo."$MFE_suggestion [$mfe_1,$mfe_2] Score: $MFE_score/$MFE_score";
					}
				else{
						$MFE_suggestion_memo=$MFE_suggestion_memo."$MFE_suggestion [$mfe_1,$mfe_2] Score: 0/$MFE_score";
					}		
				
				
				#}
				
				$MFE_suggestion_memo=~s/\t$//;
				
				$array_direction_MEMO=$array_direction_MEMO." Secondary structural analysis prediction:  $MFE_suggestion_memo;";
				#------------------------------------------------------------------------------------------
				#system("echo '$accession\t$array_start_position\t$array_stop_position\t$array_direction [$observed_percent_similarity]\t$MFE_suggestion [$mfe_1,$mfe_2]\t$model_repeat' >>all_repeats_MFE_suggestions.txt");
				#system("echo '$accession\t$array_start_position\t$array_stop_position\t$array_direction [$observed_percent_similarity]\t$MFE_suggestion_memo\t$model_repeat' >>all_repeats_MFE_suggestions.txt");
			}
		else{
				$array_direction_MEMO=$array_direction_MEMO." Secondary structural analysis prediction:  NA ;";
			}		
	
	
	
	
	
	
	
		#------------------- next analyze the other predictions using array information -------------------------------------------------------------------
		

		#-------- check degeneracy -------------------------------------			
		my $degeneracy_suggestion="NA";				
		if($check_array_degeneracy==1 and $#{$current_array}>4)
			{		
						
						#------------------------- check degeneracy -------------------------------------------------------------------
						my $old_top_degeneracy=0;   #keep a backup for printing -----------
						my $old_bottom_degeneracy=0;
						my $top_degeneracy=0;
						my $middle_degeneracy=0;
						my $bottom_degeneracy=0;
						my $degeneracy_in_first_repeat=0;
						my $degeneracy_in_last_repeat=0;
						
						my @arr_repeats;
						my @arr_comments;
							
						for(my $k1=4;$k1<=$#{$current_array}-1;$k1++)
							{
								#print "@{$current_array-[0]}\n";
								my $current_line=$$current_array[$k1]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
											
								#print "C: $current_line\n";
											
								my @tmp_array=split('\t',$current_line);
								my($r_start,$r_length,$s_length,$r_seq,$s_seq);
								my $comment="";		
								$r_start=$tmp_array[0];
								$r_seq=$tmp_array[1]; my $r_seq_1=$r_seq;$r_seq_1=~s/-//g;
								$s_seq=$tmp_array[2]; my $s_seq_1=$s_seq;$s_seq_1=~s/-//g;
								
								my $no_of_insertions_and_deletions=0;
								my $dotless_r_seq=$r_seq;$dotless_r_seq=~s/\.+//g;
								
								if($tmp_array[3])
									{
										$comment=$tmp_array[3];$comment=~s/^\s+//;
										if($comment!~/^Del/)
											{
												my @tmp_array=split(' ',$comment);
										
												#($cur_insertion_bases,$cur_insertion_positions)=split(' ',$cur_comment);
												my $cur_insertion_bases=$tmp_array[0];
												
												$cur_insertion_bases=~s/,//g;
												$no_of_insertions_and_deletions=length($cur_insertion_bases);
												
												#if(defined $tmp_array[2] and $tmp_array[2]=~/Del/){$no_of_insertions_and_deletions++;}
											}
										else{
												#$no_of_insertions_and_deletions=1;
											}	
									}
								#--------------------------------------------------------							
								push(@arr_repeats,$r_seq);	
								push(@arr_comments,$comment);	
								
								#------ new scoring system ------------------------------
										#my $current_repeat=$r_seq; 
										#my $dotless_current_repeat=$current_repeat; $dotless_current_repeat=~s/\.//g; $$dotless_current_repeat=~s/-+/-/g;
								if($k1==4)
									{
										my $next_line=$$current_array[$k1+1]; chomp $next_line; $next_line=~s/\r+//g; $next_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
										my @tmp_array2=split('\t',$next_line);
										my $next_r_seq=$tmp_array2[1];
										#if($r_seq eq $next_r_seq){$degeneracy_in_first_repeat=$no_of_insertions_and_deletions;}
										#else{
										#		$degeneracy_in_first_repeat=length($current_repeat)+$no_of_insertions_and_deletions;
										#	}
										
										#if($r_seq ne $next_r_seq or $no_of_insertions_and_deletions>0)
										if($dotless_r_seq ne "" or $no_of_insertions_and_deletions>0)
											{
												#print "\tCR:$r_seq\n\tNR:$next_r_seq\n";
												$degeneracy_in_first_repeat=1;
											}
										
									}
								elsif($k1==$#{$current_array}-1)
									{
										my $prev_line=$$current_array[$k1-1]; chomp $prev_line; $prev_line=~s/\r+//g; $prev_line=~s/^\s+//;	#$current_line=~s/\s+/\t/g;
										my @tmp_array2=split('\t',$prev_line);
										my $prev_r_seq=$tmp_array2[1];
										#if($r_seq eq $prev_r_seq){$degeneracy_in_last_repeat=$no_of_insertions_and_deletions;}
										#else{
										#		$degeneracy_in_last_repeat=length($current_repeat)+$no_of_insertions_and_deletions;
										#	}
										
										
										#if($r_seq ne $prev_r_seq or $no_of_insertions_and_deletions>0){$degeneracy_in_last_repeat=1;}
										if($dotless_r_seq ne "" or $no_of_insertions_and_deletions>0){$degeneracy_in_last_repeat=1;}
									}								
							}
							
							
						my $skip_old_scoring=1;
						if($skip_old_scoring==1)
						{	
							
						my $devider;
						
						if($#arr_repeats<=4){$devider=0.40;}else{$devider=0.33;}   # for shorter arrays which may have only 5 repeats this works well
							
						for(my $i=0;$i<=$#arr_repeats;$i++)
							{
								my $current_repeat=$arr_repeats[$i]; 
								my $dotless_current_repeat=$current_repeat;  $dotless_current_repeat=~s/\.//g; $dotless_current_repeat=~s/-+/-/g;
								my $current_comment;
								my $no_of_insertions_and_deletions=0;
								
								if($arr_comments[$i])
									{
										$current_comment=$arr_comments[$i]; $current_comment=~s/^\s+//g;
										
										if($current_comment!~/^Del/)
											{
												my @tmp_array=split(' ',$current_comment);
										
												#($cur_insertion_bases,$cur_insertion_positions)=split(' ',$cur_comment);
												my $cur_insertion_bases=$tmp_array[0];
												
												$cur_insertion_bases=~s/,//g;
												$no_of_insertions_and_deletions=length($cur_insertion_bases);
												
												#if(defined $tmp_array[2] and $tmp_array[2]=~/Del/){$no_of_insertions_and_deletions++;}
											}
										else{
												#$no_of_insertions_and_deletions=1; # do not score deletion from comment,as this will be scored in repeat anyway
											}
									}
								
								if($i<= int($#arr_repeats*$devider))
									{
										#---- check if the current repeat is same as the next repeat or not
										my $prev_r_seq=$arr_repeats[$i-1];
										my $next_r_seq=$arr_repeats[$i+1];
										
										if($i==0) 
											{
												if(length($dotless_current_repeat)==0 or $current_repeat eq $next_r_seq){$top_degeneracy=$top_degeneracy+$no_of_insertions_and_deletions;}
												else{
														$top_degeneracy=$top_degeneracy+length($dotless_current_repeat)+$no_of_insertions_and_deletions;
													}
											}
										else{
												if(length($dotless_current_repeat)==0 or $current_repeat eq $next_r_seq  or $current_repeat eq $prev_r_seq){$top_degeneracy=$top_degeneracy+$no_of_insertions_and_deletions;}
												else{
														$top_degeneracy=$top_degeneracy+length($dotless_current_repeat)+$no_of_insertions_and_deletions;
													}
											}
										
										#$top_degeneracy=$top_degeneracy+length($current_repeat);
										
											#----- check if the first repeat has any degeneracy---
											#if($i==0 and length($current_repeat)>0){$degeneracy_in_first_repeat=1;}
									}
								elsif($i>=int(($#arr_repeats*$devider)*2))
									{
										my $prev_r_seq=$arr_repeats[$i-1];
										my $next_r_seq=$arr_repeats[$i+1];
															#-- then compare with the previous repeat, else compare with the next repeat,
															# advantage: when, the first repeat of this region comes for checking, it should be checked with next one, instead of the previous one, as previous one can be a perfect repeat, See the example:
															# >NC_012724

															#	   Pos  Repeat     %id  Spacer  Left flank    Repeat                               Spacer
															#==========  ======  ======  ======  ==========    ============================    ======
															#   1995500      28   100.0      32  NNNNNNNNNN    ............................    CTACCGGCGCGAACTGGAGCGCCTGCTGCTCT           
															#   1995560      28   100.0      33  NNNNNNNNNN    ............................    CACGAGCGCACACGGCACCGGCCAGTGCGTGCA           
															#   1995621      28   100.0      33  NNNNNNNNNN    ............................    GCGCTTCAACAGGATTCTCTGACGTCGCCCTGA           
															#   1995682      28   100.0      33  NNNNNNNNNN    ............................    CTCGAACGGGGACGCGGCATGACTTACCGCGAG           
															#   1995743      28   100.0      32  NNNNNNNNNN    ................C......AT...    CGCGTCATGCGAAACGGATACGTCGAGATCCC           
															#   1995803      28   100.0       1  NNNNNNNNNN    ................C......AT...    NNNNNNNNNN
															#==========  ======  ======  ======  ==========    ============================
															#		 6      28              32                TTTCTAAGCTGCCTACACGGCAGCGAAC

										if($i==$#arr_repeats) 
											{
												if(length($dotless_current_repeat)==0 or $current_repeat eq $prev_r_seq){$bottom_degeneracy=$bottom_degeneracy+$no_of_insertions_and_deletions;}
												else{
														#print "CR:$current_repeat\nPR:$prev_r_seq\n\n";
														$bottom_degeneracy=$bottom_degeneracy+length($dotless_current_repeat)+$no_of_insertions_and_deletions;
													}
											}
										else{
												if(length($dotless_current_repeat)==0 or $current_repeat eq $next_r_seq or $current_repeat eq $prev_r_seq){$bottom_degeneracy=$bottom_degeneracy+$no_of_insertions_and_deletions;}
												else{
														#print "CR:$current_repeat\nPR:$prev_r_seq\n\n";
														$bottom_degeneracy=$bottom_degeneracy+length($dotless_current_repeat)+$no_of_insertions_and_deletions;
													}
											}
										#$bottom_degeneracy=$bottom_degeneracy+length($current_repeat);
											#if($i==$#arr_repeats and length($current_repeat)>0){$degeneracy_in_last_repeat=1;}
									}
								else{
										my $next_r_seq=$arr_repeats[$i+1];
										if(length($dotless_current_repeat)==0 or $current_repeat eq $next_r_seq){$middle_degeneracy=$middle_degeneracy+$no_of_insertions_and_deletions;}
										else{
												$middle_degeneracy=$middle_degeneracy+length($dotless_current_repeat)+$no_of_insertions_and_deletions;
											}
										#$middle_degeneracy=$middle_degeneracy+length($dotless_current_repeat)+$no_of_insertions_and_deletions;
										
									}	
										
								#print "$arr_repeats[$i]\n";
							}

						my $total_array_degeneracy=$top_degeneracy+$middle_degeneracy+$bottom_degeneracy;
						
						#---- add one extra point if the top/bottom first repeat shows degeneracy 
						$old_top_degeneracy=$top_degeneracy;   #keep a backup for printing -----------
						$old_bottom_degeneracy=$bottom_degeneracy;
						$top_degeneracy=$top_degeneracy+$degeneracy_in_first_repeat;			
						$bottom_degeneracy=$bottom_degeneracy+$degeneracy_in_last_repeat;
						#------------------------------------------------------------------------
						
						if($total_array_degeneracy<=$permitted_mutation_per_array){$degeneracy_suggestion="NA";}			
						elsif($bottom_degeneracy>$top_degeneracy){$degeneracy_suggestion="F";}
						elsif($top_degeneracy>$bottom_degeneracy){$degeneracy_suggestion="R";}
						elsif($bottom_degeneracy==$top_degeneracy)
							{
								#if($degeneracy_in_first_repeat!=1 and $degeneracy_in_last_repeat==1){$degeneracy_suggestion="F";}
								#elsif($degeneracy_in_first_repeat==1 and $degeneracy_in_last_repeat!=1){$degeneracy_suggestion="R";}
								#else{
										$degeneracy_suggestion="NA";
								#	}	
							}
						else{$degeneracy_suggestion="NA";}
						}
						
						#----- new scoring -------------------
						#if($degeneracy_in_last_repeat>$degeneracy_in_first_repeat){$degeneracy_suggestion="F";}
						#elsif($degeneracy_in_first_repeat>$degeneracy_in_last_repeat){$degeneracy_suggestion="R";}
						#else{$degeneracy_suggestion="NA";}
						#----------------------------append the suggestion in MEMO string -------------------------------------------	
						if($degeneracy_suggestion!~/NA/)
							{
								$array_direction_MEMO=$array_direction_MEMO." Array degeneracy analysis prediction:      $degeneracy_suggestion [$old_top_degeneracy-$old_bottom_degeneracy] Score: $Array_degeneracy_score/$Array_degeneracy_score;";
							}
						else{
								$array_direction_MEMO=$array_direction_MEMO." Array degeneracy analysis prediction:      $degeneracy_suggestion [$old_top_degeneracy-$old_bottom_degeneracy] Score: 0/$Array_degeneracy_score;";
							}		
						#$array_direction_MEMO=$array_direction_MEMO." $degeneracy_suggestion [$degeneracy_in_first_repeat-$degeneracy_in_last_repeat],";						
						#--------------- record it to a file fore directional analysis ----------------------------------------------
						#system("echo '$accession\t$array_start_position\t$array_stop_position\t$array_direction [$observed_percent_similarity]\t$degeneracy_suggestion [$old_top_degeneracy-$old_bottom_degeneracy]\t$degeneracy_in_first_repeat\t$degeneracy_in_last_repeat' >>array_degeneracy_predictions.txt");
			}
		else{
				$array_direction_MEMO=$array_direction_MEMO." Array degeneracy analysis prediction:      NA ;";
			}					

				
					
		#-------- check AT richness ------------------------------------			
		my $at_distribution="";		
			
		if($check_AT_distribution_in_flanks==1)
			{
				my $at_richness_suggestion="NA";					
				
				
				
						
				#------------------------------- Window for AT richness and at_richness_suggestion -----------------------------------						
						
				my $total_percent_at_in_ls=0;
				my $total_percent_at_in_ts=0;
						
				#if($AT_distribution_window>=75)
				#	{							
				#	for(my $i=15;$i<=175;$i=$i+15)
				#		{
				#	$AT_distribution_window=$i;
				#-------- now get the AT richness in 150base window before and after the array -------------------------------
				my $leading_seq=substr($species_seq,($array_start_position-1-$AT_distribution_window),$AT_distribution_window);
				my $trailing_seq=substr($species_seq,$array_stop_position-1,$AT_distribution_window);
								
				my $total_a_in_ls=$leading_seq=~s/A/A/gi;
				my $total_t_in_ls=$leading_seq=~s/T/T/gi;
				my $total_at_in_ls=$total_a_in_ls+$total_t_in_ls;
				$total_percent_at_in_ls=sprintf("%.1f",($total_at_in_ls/$AT_distribution_window)*100);
								
								
				my $total_a_in_ts=$trailing_seq=~s/A/A/gi;
				my $total_t_in_ts=$trailing_seq=~s/T/T/gi;
				my $total_at_in_ts=$total_a_in_ts+$total_t_in_ts;
				$total_percent_at_in_ts=sprintf("%.1f",($total_at_in_ts/$AT_distribution_window)*100);
								
				my $higher_richness;
				if($total_percent_at_in_ls>$total_percent_at_in_ts*(1+$AT_distribution_minimum_percentage_difference/100)) # minimum 10% higher
					{
						$higher_richness="F";
					}
				elsif($total_percent_at_in_ts>$total_percent_at_in_ls*(1+$AT_distribution_minimum_percentage_difference/100))
					{
						$higher_richness="R";
					}	
				else{
						$higher_richness="NA";
					}	
				
				if($higher_richness !~/NA/)
					{				
						$at_distribution=$at_distribution."$higher_richness [$total_percent_at_in_ls-$total_percent_at_in_ts]\%AT Score: $AT_distribution_score/$AT_distribution_score";
					}
				else{
						$at_distribution=$at_distribution."$higher_richness [$total_percent_at_in_ls-$total_percent_at_in_ts]\%AT Score: 0/$AT_distribution_score";
					}				
				#	}# end of for
								
						
				if($total_percent_at_in_ls>$total_percent_at_in_ts)
					{
						$at_richness_suggestion="F";
					}
				elsif($total_percent_at_in_ls<$total_percent_at_in_ts)
					{
						$at_richness_suggestion="R";
					}
				else{
						$at_richness_suggestion="NA";
					}								
						#	}
						
				#----------------------------append the suggestion in MEMO string ------------------------	
				#$array_direction_MEMO=$array_direction_MEMO." $at_richness_suggestion [$total_percent_at_in_ls-$total_percent_at_in_ts]\%AT,";
				$array_direction_MEMO=$array_direction_MEMO." AT richness analysis in flanks prediction: $at_distribution;";
				$at_distribution=~s/\t$//;
						#system("echo '$accession\t$array_start_position\t$array_stop_position\t$array_direction [$observed_percent_similarity]\t$at_distribution' >>at_richness_distribution.txt");		
			}#--- end of AT_richness block						
		else{
				$array_direction_MEMO=$array_direction_MEMO." AT richness analysis in flanks prediction: NA ;";
			}					
		
		
		
		#------ check longer leader ------------------------------------				
		if($check_longer_leader==1 and $all_gene_positions_file ne "NA")
			{	
						
				my $longer_leader_suggestion="NA";
				my $longer_l_suggestions="";	
						
				#---------------------------------- get all the gene positions from the all_gene_and_crispr_positions.txt file --------------------
				my @arr_gene_positions=`grep -w '$accession' $all_gene_positions_folder/$all_gene_positions_file >&1`;
				#----------------------------------------------------------------------------------------------------------------------------------
						
						
				if($#arr_gene_positions<=0)
					{
						$array_direction_MEMO=$array_direction_MEMO." Longer leader analysis prediction:         NA ;";
						#next;
					}
				else{

						#------------ now get the previous gene stop and next gene start position -----------------------------------------
						my $previous_gene_stop=0;
						my $next_gene_start=length($species_seq)+1;
						foreach my $gene_det_line(@arr_gene_positions)
								{
									my @arr_tmp1=split('\t',$gene_det_line);
											
									if($arr_tmp1[3]>$previous_gene_stop and $arr_tmp1[3]<$array_start_position)
										{
											$previous_gene_stop=$arr_tmp1[3];
										}
											
									if($arr_tmp1[2]<$next_gene_start and $arr_tmp1[2]>$array_stop_position){$next_gene_start=$arr_tmp1[2];}
								}
						my $left_len=  $array_start_position - $previous_gene_stop + 1;
						my $right_len= $next_gene_start - $array_stop_position + 1;

						
						
						
						
							
						#------------------------------- Longer leader suggestion -----------------------------------------------------------						
						
						my $i=1;
						#for(my $i=0;$i<=1;$i=$i+0.2)  # check and record  the % differences starting with 0% (thats minimum 1 base), 20%, 40%, 60%, 80% and 100%
						#{
							$longer_leader_suggestion="NA";
							if($left_len<$right_len)
								{								
									if($right_len>int($left_len*(1+$i)) and ($left_len>75 or $right_len>75)) #---- there should be a length difference of at least 1.5 times or 75 bases [discuss with chris/peter ]
										{
											$longer_leader_suggestion="R";
										}
									else{
											$longer_leader_suggestion="NA";
										}		
								}
							elsif($left_len>$right_len)
								{
									if($left_len>int($right_len*(1+$i)) and ($left_len>75 or $right_len>75)) #---- there should be a length difference of at least 1.5 times or 75 bases [discuss with chris/peter ]
										{
											$longer_leader_suggestion="F";
										}
									else{
											$longer_leader_suggestion="NA";
										}									
								}
							my $current_diff=(1+$i)*100; $current_diff=$current_diff."%";
							#print "$accession $array_direction [$observed_percent_similarity] |--- $previous_gene_stop ---- $array_start_position ---- $array_stop_position ---- $next_gene_start ----| [$left_len - $right_len] \n";	
							
							if($longer_leader_suggestion!~/NA/)
								{
									$longer_l_suggestions=$longer_l_suggestions."$longer_leader_suggestion [$left_len,$right_len] Score: $Longer_leader_score/$Longer_leader_score";
								}
							else{
									$longer_l_suggestions=$longer_l_suggestions."$longer_leader_suggestion [$left_len,$right_len] Score: 0/$Longer_leader_score";
								}		
							#------------------------------
							
							
						#}
						
						$array_direction_MEMO=$array_direction_MEMO." Longer leader analysis prediction:         $longer_l_suggestions;";
						#system("echo '$accession\t$array_start_position\t$array_stop_position\t$array_direction [$observed_percent_similarity]\tDist: [$left_len-$right_len]\t$longer_l_suggestions' >>longer_leader_prediction.txt");
						
					}
			}
		else{
				$array_direction_MEMO=$array_direction_MEMO." Longer leader analysis prediction:         NA ;";
			}		
						
						
			
		#----------------- now process the MEMO field and get the Final direction ----------------------------------------------------------
		# $A_and_T_ratio_score,$Similarity_score,$MFE_score,$Array_degeneracy_score,$AT_distribution_score,$Longer_leader_score,
		my @arr_memo=split(';',$array_direction_MEMO);
		my $score_F=0;
		my $score_R=0;
		my $total_confidence=0;
		foreach my $analysis(@arr_memo)
			{
				if($analysis=~/Motif/i)
					{
						my $c_score=$Motif_match_score;
											
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1]; $direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
							#$direction=~s/\s+$//;#if($direction=~/ \[\S+\]?$/){$direction=~s/ \[\S+\]?$//;}
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}	
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}					
					}
				elsif($analysis=~/A,T distribution/i)
					{
						my $c_score=$A_and_T_ratio_score;
						
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1];$direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}	
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}					
					}
				elsif($analysis=~/Reference repeat/i)
					{
						my $c_score=$Similarity_score;
						
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1];$direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}	
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}					
					}
				elsif($analysis=~/Secondary/i)
					{
						my $c_score=$MFE_score;
						
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1];$direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}						
					}
				
				
				elsif($analysis=~/Array degeneracy/i)
					{
						my $c_score=$Array_degeneracy_score;
						
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1];$direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}		
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}				
					}
				elsif($analysis=~/AT richness/i)
					{
						my $c_score=$AT_distribution_score;
						
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1];$direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}						
					}
				elsif($analysis=~/Longer leader/i)
					{
						my $c_score=$Longer_leader_score;					
						
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1];$direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}	
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}					
					}						
			}
		
		my $number_of_F=$array_direction_MEMO=~s/F/F/g;
		my $number_of_R=$array_direction_MEMO=~s/R/R/g;
		my $final_array_direction="NA";
		my $final_array_direction_MEMO;	
		
		my $confidence="NA";
		if($score_F>$score_R)
			{
				$final_array_direction="F";
				my $score_diff=$score_F-$score_R;
				$confidence=&get_confidence_score($total_confidence,$score_diff);
			}	
		if($score_F<$score_R)
			{
				$final_array_direction="R";
				my $score_diff=$score_R-$score_F;
				$confidence=&get_confidence_score($total_confidence,$score_diff);
			}
		
		#my $difference_in_prediction=abs(abs($score_F)-abs($score_R));		
		#if($difference_in_prediction>0.66){$confidence="HIGH";}
		#elsif($difference_in_prediction>0.33){$confidence="NORMAL";}
		#else{$confidence="POOR";}
		#$array_direction_MEMO=~s/,$//;
		
		$array_direction_MEMO=$array_direction_MEMO."  ; Final direction:         $final_array_direction [$score_F,$score_R   Confidence: $confidence] ;";
		#print "Final_array_direction\t$final_array_direction\n\n";					
			
		
		
		
		#my $l_flank=substr($species_seq,$array_start_position-1-100,100);
		#my $r_flank=substr($species_seq,$array_stop_position-1,100);
		#-----------------------------------------------------------------------------------------------------------------------------------	
		return($matching_reference_repeat,$model_repeat,$array_direction,$repeat_family,$array_direction_MEMO,$case_found);	
	
	}



sub check_repeats_direction()
	{
		my($check_motif_in_repeat,$motif,$check_A_and_T_ratio_in_repeat,$check_similarity_with_reference_repeat,$allowed_no_of_mismatches,$check_secondary_structure_of_repeat,$MFE_cutoff,$MFE_minimum_difference,$MFE_exclude_bases,$check_array_degeneracy,$permitted_mutation_per_array,$check_AT_distribution_in_flanks,$AT_distribution_window,$AT_distribution_minimum_percentage_difference,$check_longer_leader,$Motif_match_score,$A_and_T_ratio_score,$Similarity_score,$MFE_score,$accession,$model_repeat,$lib_of_repeats_with_confirmed_direction)=@_;
		#print "\tGoing to check for consensus sequence with \$model_repeat=$model_repeat and ref:$model_repeats_reference_string for $accession:\n\n";
		


		my $array_direction="NA";
		my $repeat_family="NA";
		my $matching_reference_repeat="NA";
		my $match_found=0;
		my $observed_percent_similarity=0;
			
		#($matching_reference_repeat,$repeat_family,$array_direction,$observed_percent_similarity)=split(';',$model_repeats_reference_string);		
		
		
		
		#my $check_motif_in_repeat=1;my $motif="ATTGAAA.?";
		#my $check_A_and_T_ratio_in_repeat=1;
		#my $check_similarity_with_reference_repeat=1; my $allowed_no_of_mismatches=3;
		#my $check_secondary_structure_of_repeat=1;		
		#	my $MFE_cutoff=1;		# Default 1
		#	my $MFE_minimum_difference=0.25;	# Default 0.25
		#	my $MFE_exclude_bases=5;				# Default 5
	

		
		
		
		
		#-------------------------- some common parameters -----------------------------------------------------			
		my $case_found=0;
		my $new_model_repeat="";
		my $potential_alternate_repeat="NA";
		my $array_direction_MEMO="";
		

		
		
		#---------- now open the sequence file and get the sequence string -------------------------------------------------
		#open(SEQ,"tmp/$accession\.fna") or print "$!";
		#my @arr_seq=<SEQ>;
		#close(SEQ);
		#my $species_seq=$arr_seq[1]; chomp $species_seq;$species_seq=~s/\r//g;
		
		
		#--------------------------------------------- first check the model_repeats for matching to ref repeat, motif and A_and_T count ----------------------------------------------------------
		
		
		#---- check motif -------------------------------------------------
		if($check_motif_in_repeat==1)
			{					
				
				my $model_repeat_rc=$model_repeat; $model_repeat_rc= reverse $model_repeat_rc; $model_repeat_rc=~tr/ACGT/TGCA/;				
				my $suggested_direction="NA";
										
				if($model_repeat=~/$motif/ and $model_repeat_rc=~/$motif/){$suggested_direction="NA";}
				elsif($model_repeat=~/$motif/ ){$suggested_direction="F";}
				elsif($model_repeat_rc=~/$motif/){$suggested_direction="R";}

				#-----  change all .? to N
				$motif=~s/\.\?/\(N\)/g;
				if($suggested_direction!~/NA/)
					{
						$array_direction_MEMO=$array_direction_MEMO." Motif $motif match prediction:         $suggested_direction Score: $Motif_match_score/$Motif_match_score;";
					}
				else{
						$array_direction_MEMO=$array_direction_MEMO." Motif $motif match prediction:         $suggested_direction Score: 0/$Motif_match_score;";
					}		
				
			}
		else{
				$motif=~s/\.\?/\(N\)/g;
				$array_direction_MEMO=$array_direction_MEMO." Motif $motif matched prediction:       NA ;";
			}	


		#---- check A&T ratio -------------------------------------------------
		if($check_A_and_T_ratio_in_repeat==1)
			{				
				#--- replace all U with T if present -----------
				$model_repeat=~tr/U/T/;

				##################### check direction using no. of As and Ts: mainly applicable to Archaeal genomes #######################		
												
				my $no_of_As_in_repeat=$model_repeat=~s/A/A/gi;
				my $no_of_Ts_in_repeat=$model_repeat=~s/T/T/gi;	
				my $A_and_T_suggestion="NA";
											
				my $at_p_d=(($no_of_As_in_repeat+$no_of_Ts_in_repeat)/length($model_repeat))*100; $at_p_d=sprintf("%.2f",$at_p_d);
				my $gc_p_d=((length($model_repeat) -$no_of_As_in_repeat-$no_of_Ts_in_repeat)/length($model_repeat))*100; $gc_p_d=sprintf("%.2f",$gc_p_d);
											
				if($no_of_As_in_repeat>$no_of_Ts_in_repeat)
					{
						$A_and_T_suggestion="F";
					}
				elsif($no_of_Ts_in_repeat>$no_of_As_in_repeat)
					{
						$A_and_T_suggestion="R";
					}							

				
				$array_direction_MEMO=$array_direction_MEMO." A,T distribution in repeat prediction:     $A_and_T_suggestion [$no_of_As_in_repeat,$no_of_Ts_in_repeat] Score: $A_and_T_ratio_score/$A_and_T_ratio_score;";
				
			}
		else{
				 $array_direction_MEMO=$array_direction_MEMO." A,T distribution in repeat prediction:     NA ;";
			}			
		
		
		#------- check reference-repeats ---------------------------------------------
		if($check_similarity_with_reference_repeat==1)
			{
				
				($matching_reference_repeat,$repeat_family,$array_direction,$observed_percent_similarity)=&get_matching_reference_repeat_and_direction($accession,$model_repeat,$allowed_no_of_mismatches,$lib_of_repeats_with_confirmed_direction);
				
				#print "$matching_reference_repeat,$repeat_family,$array_direction,$observed_percent_similarity\n";
				#----------------------------append the suggestion in MEMO string -------------------------------------------	
				if($matching_reference_repeat!~/NA/)
					{
						$array_direction_MEMO=$array_direction_MEMO." Reference repeat match prediction:         $array_direction [matched $matching_reference_repeat with $observed_percent_similarity\% identity] Score: $Similarity_score/$Similarity_score;";
					}
				else{
						$array_direction_MEMO=$array_direction_MEMO." Reference repeat match prediction:         NA ;";
					}		
				#-------------------------------------------now check for the AT richness/longer leader/degeneracy -------------------------------
		
			}
		else{
						$array_direction_MEMO=$array_direction_MEMO." Reference repeat match prediction:         NA ;";
			}
		
						
		
	
		#------- check MFE ---------------------------------------------
		if($check_secondary_structure_of_repeat==1)
			{				
						

				my $forward_strand=$model_repeat;
				my $MFE_suggestion_memo="";
					
				if($MFE_exclude_bases>0)
					{
						my $string_N="";
						for(1..$MFE_exclude_bases){$string_N=$string_N."N";}
						
						$forward_strand=~s/^\S{$MFE_exclude_bases}/$string_N/;
						$forward_strand=~s/\S{$MFE_exclude_bases}$/$string_N/;						
					}
				
				
				
				my $other_strand=$forward_strand; $other_strand=reverse $other_strand; $other_strand=~tr/ACGT/TGCA/;
				
				my $mfe_1=&run_rnafold($forward_strand);				
				my $mfe_2=&run_rnafold($other_strand);
				
				#print "\n";
				
				my $MFE_suggestion="NA";	
				my $valid_case=0;
				if( (abs( abs($mfe_1) - abs($mfe_2) )>=$MFE_minimum_difference) and ( abs($mfe_1)>=$MFE_cutoff or abs($mfe_2)>=$MFE_cutoff))
					{
						$valid_case=1;
					}				
				
				if($mfe_1<$mfe_2 and $valid_case==1){$MFE_suggestion="F";}
				elsif($mfe_2<$mfe_1 and $valid_case==1){$MFE_suggestion="R";}
				

								
				
				if($MFE_suggestion!~/NA/)
					{
						$MFE_suggestion_memo=$MFE_suggestion_memo."$MFE_suggestion [$mfe_1,$mfe_2] Score: $MFE_score/$MFE_score";
					}
				else{
						$MFE_suggestion_memo=$MFE_suggestion_memo."$MFE_suggestion [$mfe_1,$mfe_2] Score: 0/$MFE_score";
					}		
				
				
				#}
				
				$MFE_suggestion_memo=~s/\t$//;
				
				$array_direction_MEMO=$array_direction_MEMO." Secondary structural analysis prediction:  $MFE_suggestion_memo;";
				#------------------------------------------------------------------------------------------
				
			}
		else{
				$array_direction_MEMO=$array_direction_MEMO." Secondary structural analysis prediction:  NA ;";
			}		
	
	
	
	
	
	
	
	
				
						
			
		#----------------- now process the MEMO field and get the Final direction ----------------------------------------------------------
		# $A_and_T_ratio_score,$Similarity_score,$MFE_score,$Array_degeneracy_score,$AT_distribution_score,$Longer_leader_score,
		my @arr_memo=split(';',$array_direction_MEMO);
		my $score_F=0;
		my $score_R=0;
		my $total_confidence=0;
		foreach my $analysis(@arr_memo)
			{
				if($analysis=~/Motif/i)
					{
						my $c_score=$Motif_match_score;
											
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1]; $direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
							#$direction=~s/\s+$//;#if($direction=~/ \[\S+\]?$/){$direction=~s/ \[\S+\]?$//;}
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}	
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}					
					}
				elsif($analysis=~/A,T distribution/i)
					{
						my $c_score=$A_and_T_ratio_score;
						
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1];$direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}	
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}					
					}
				elsif($analysis=~/Reference repeat/i)
					{
						my $c_score=$Similarity_score;
						
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1];$direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}	
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}					
					}
				elsif($analysis=~/Secondary/i)
					{
						my $c_score=$MFE_score;
						
						#----------------------------------------------------
						my @arr_tmp=split(':',$analysis);
						my $direction=$arr_tmp[1];$direction=~s/^\s+//; 
						my @arr_tmp2=split(' ',$direction);
						$direction=$arr_tmp2[0];
						
						if($direction =~ /F/){$score_F=$score_F+$c_score;}
						elsif($direction =~ /R/){$score_R=$score_R+$c_score;}
						
						#-----------------------------------------------------
						if($direction !~ /NA/)
							{
								$total_confidence=$total_confidence+$c_score;
							}						
					}
				
				
				
			}
		
		#my $number_of_F=$array_direction_MEMO=~s/F/F/g;
		#my $number_of_R=$array_direction_MEMO=~s/R/R/g;
		my $final_array_direction="NA";
		my $final_array_direction_MEMO;	
		
		my $confidence="NA";
		if($score_F>$score_R)
			{
				$final_array_direction="F";
				my $score_diff=$score_F-$score_R;
				$confidence=&get_confidence_score($total_confidence,$score_diff);
			}	
		if($score_F<$score_R)
			{
				$final_array_direction="R";
				my $score_diff=$score_R-$score_F;
				$confidence=&get_confidence_score($total_confidence,$score_diff);
			}
		#$array_direction_MEMO=~s/,$//;
		
		#print "\n\nFinal direction:         $final_array_direction [$score_F,$score_R   Confidence: $prediction_sensitivity] ;\n\n";
		$array_direction_MEMO=$array_direction_MEMO."  ; Final direction: $final_array_direction [$score_F,$score_R   Confidence: $confidence] ;";
		#print "Final_array_direction\t$final_array_direction\n\n";					
			
		
		#-----------------------------------------------------------------------------------------------------------------------------------	
		return($matching_reference_repeat,$model_repeat,$final_array_direction,$repeat_family,$array_direction_MEMO,$case_found);	
	}	
	
	

sub run_rnafold()
	{
		my($sequence)=shift(@_);
		my $mfe=0;
		
	
		my $time = time()."_".rand(0000);
		my $rep_seq_file="mr_".$time.".txt";				
		my $out_file="MFE_".$time.".mfe";
		
		open(WR,">tmp/$rep_seq_file");close(WR);	system("chmod 777 tmp/$rep_seq_file");
		open(WR,">tmp/$out_file");close(WR);		system("chmod 777 tmp/$out_file");
				
		my $id=">Seq_$time";				
		my $string=$id."\n".$sequence;
				
		system("echo '$string' >tmp/$rep_seq_file");		
		
		
		##---- check RNAfold installation ----------------------------------------------------
		#my $RNAfold=`which RNAfold >&1 2>&1`; chomp $RNAfold; $RNAfold=~s/\r//g;
		#----- run RNAfold -------------------------------------------------------------------		
		system("RNAfold --noLP --noPS < tmp/$rep_seq_file >tmp/$out_file");					#print "\n$RNAfold --noPS < tmp/$rep_seq_file >tmp/$out_file\n";
		#-------------------------------------------------------------------------------------		
				
				
		open(MFE,"tmp/$out_file") or print "$!: $out_file not found\n";
		flock(MFE,2);
		my @arr_mfe=<MFE>;
		close(MFE); 
				

		foreach my $l(@arr_mfe)
			{						
				chomp $l;$l=~s/\r+//g;
				if($l=~/\)$/)
					{
								#print ": $l\n";
						my @arr_l=split(' ',$l);
						$mfe=$arr_l[$#arr_l];
						$mfe=~s/\(//g;
						$mfe=~s/\)//g;
						$mfe=~s/\s+//g;
								#last;
					}
			}
			
		unlink("tmp/$rep_seq_file");
		unlink("tmp/$out_file");	
					
		return($mfe);
	}	
	
	
