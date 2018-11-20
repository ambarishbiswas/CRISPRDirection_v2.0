#!/usr/bin/perl -W

#--- Note: This program reads the related properties from different methods (e.g: for AT richness in flanks methods, the %AT distribution in different windows)
#---       and calculated the PPV for each of them. The Best PPVs without compromising the specificity is selected from each method, and applied to the final method
#---



use strict;

my $q=10;		# number of tests
my $r=100;      # number of arrays per set


#------ change the following parameters to weight higher or lower a specific method
my $at_window_length=60;  			# Default 60 
my $minimum_percent_at_diff=10;		# Default 10 

my $permitted_no_of_bases=0;		# Default 0

my $minimum_energy_threshold=1;		# Default 1
my $minimum_energy_difference=0.25;	# Default 0.25
my $mfe_target_base=5;				# Default 5




my %hash_of_records_used_for_training;






#-------------------- dont modify anything below this unless you know what you are doing :) --------------------------------------------------------------






#---------------------- read files -----------------------------------------------------------------------------------------------------------------------
my %hash_of_analysis_results;

#---- make a fasta file of all unique repeats ---------------------
my %hash_tmp;
open(RD,"all_repeats_from_sample_set.txt") or print "$!";
my @arr_tmp=<RD>;
close(RD);
for (my $i=0;$i<=$#arr_tmp;$i++)
	{
		if($arr_tmp[$i]=~/>/)
			{
				my $header_l=$arr_tmp[$i]; chomp $header_l; $header_l=~s/\r//g;
				my $seq=$arr_tmp[$i+1];	chomp $seq; $seq=~s/\r//g;
				my $seq_rc=$seq;$seq_rc= reverse $seq_rc;$seq_rc=~tr/ACGT/TGCA/;
				if(not defined $hash_tmp{$seq} and not defined $hash_tmp{$seq_rc})
					{
						#print "$header_l\n$seq\n";
						$hash_tmp{$seq}=1;
					}
			}
	}
print "\n";	


#----------- get all bacterial accessions -----------------------------	
	
my %hash_of_bacterial_genomes;		
my %hash_of_bacterial_acc;
open(RD,"all_bacterial_sequence_list.txt") or print "$!";
my @arr_tmp1=<RD>;
close(RD);
foreach my $line(@arr_tmp1)
	{
		if($line=~/\/(\S+)\.gbk/)
			{
				$hash_of_bacterial_acc{$1}=1;
				#print "$1\n";
				
				my @arr_t1=split('\/',$line);
				$hash_of_bacterial_genomes{$arr_t1[0]}=1;
			}
	}
my $total_bct_genomes=keys %hash_of_bacterial_genomes;	
#print "Total bacterial genomes=$total_bct_genomes\n";	


#----------- get all archaeal accessions -----------------------------	

my %hash_of_archaeal_genomes;	
my %hash_of_archaeal_acc;
open(RD,"all_archaeal_sequence_list.txt") or print "$!";
my @arr_tmp2=<RD>;
close(RD);
foreach my $line(@arr_tmp2)
	{
		if($line=~/\/(\S+)\.gbk/)
			{
				$hash_of_archaeal_acc{$1}=1;
				#print "$1\n";
				my @arr_t1=split('\/',$line);
				$hash_of_archaeal_genomes{$arr_t1[0]}=1;
			}
	}
my $total_archaeal_genomes=keys %hash_of_archaeal_genomes;	
#print "Total archaeal genomes=$total_archaeal_genomes\n";				
#exit;	


#------------------------------------------------------------------

my %hash_of_ref_repeat_direction;
######------------- open the nr file for the accessions -------------------------------------------------------------
open(RD,"sample_space_acc_repeat_and_their_reference.txt") or print "$!";
#open(RD,"array_degeneracy_predictions.txt") or print "$!";
my @arr_file=<RD>;
close(RD);

my %hash_of_all_records;
my %hash_of_ref_repeats;
my %hash_of_all_repeats_and_ref_direction;

my $index=1;
foreach my $line(@arr_file)
	{
			chomp $line; $line=~s/\r+//g;
			
			my @arr_tmp=split('\t',$line);
			
			#my($accession,$start,$stop)=split('-',$arr_tmp[0]);
			my $key=$arr_tmp[0];
										#my $key="$arr_tmp[0]-$arr_tmp[1]-$arr_tmp[2]";
			$hash_of_all_records{$index}=$key;
			
			$hash_of_ref_repeats{$key}=	$arr_tmp[2];
			
						
			my @arr_ref_det=split(',',$arr_tmp[2]);		
			
			
			$hash_of_ref_repeat_direction{$key}=$arr_ref_det[2];
			
			$hash_of_all_repeats_and_ref_direction{$key}{$arr_tmp[1]}=$arr_ref_det[2];
			
			#print "\$hash_of_all_repeats_and_ref_direction{$key}{$arr_tmp[1]}=$hash_of_all_repeats_and_ref_direction{$key}{$arr_tmp[1]}\n";
			$index++;
	}		
my $total_records=keys %hash_of_all_records;
#exit;


#----------- now take 50% of the total arrays for training ------------------------------------------------------- 
my $total_training_sets=0;

while($total_training_sets<int($total_records/2))
	{
		#--- now randomly pick a record from %hash_of_all_records
		my $random_no=int(rand($total_records-1));
							
		$random_no=$random_no+1;
								
		my $current_record=$hash_of_all_records{$random_no};
								
		#----- skip if not part of training set -------------------------------------
		if(not defined $hash_of_records_used_for_training{$current_record})
			{
				$hash_of_records_used_for_training{$current_record}=1;
				$total_training_sets++;
			}
		#----------------------------------------------------------------------------		
		
	}

#--- print the training set ---------------------------------------
print "Records used in training:\n";
my $rec_no=1;
foreach my $training_rec(keys %hash_of_records_used_for_training)
	{
		print "$rec_no:$training_rec\n";
		$rec_no++;
	}
#exit;
#-------------------------------------------------------------------------------------------------------------------







#---------------------- checks and PPVs ------------------------------------------------------------------------------------------------------------------




for(my $s=9;$s<=9;$s++)   #  $s=0 will apply all methods, $s=1 to $s=6 will show summary of all methods excluding one at a time, $s=7 and $s=8 will show  repeat specific and array specific methods summary, $s=9 to $s=14 will show individual analysis summary
{

my $check_motif=1;
my $check_A_and_T_in_repeats=1;

my $check_at_richness=1;
my $check_degeneracy=1;
my $check_mfe=1;
my $check_longer_leader=1;
my $check_final_analysis=1;
my $check_final_analysis_of_all_arrays=1;


if($s>=9)
	{
		$check_motif=0;
		$check_A_and_T_in_repeats=0;
		$check_mfe=0;
		
		$check_at_richness=0;
		$check_degeneracy=0;		
		$check_longer_leader=0;
		
		#my $check_final_analysis=1;
	}


my $mean_motif_PPV=5;  #   setting it to 5 will make it a filter 
my $mean_A_T_distribution_in_repeats_PPV=0;
my $mean_at_richness_ppv=0;
my $mean_degeneracy_ppv=0;
my $mean_mfe_ppv=0;
my $mean_longer_leader_ppv=0;


if($s==1)
	{
		$check_motif=0;		
		print "\n********************************************************   Analysis without Motif in repeats *********************************************\n\n";
	}
elsif($s==2)
	{
		$check_A_and_T_in_repeats=0;
		print "\n********************************************************   Analysis without A and T distribution in arrays *********************************************\n\n";
	}	
elsif($s==3)
	{
		$check_degeneracy=0;
		print "\n********************************************************   Analysis without degenerracy in repeats *********************************************\n\n";
	}	
		
elsif($s==4)
	{
		$check_at_richness=0;
		print "\n********************************************************   Analysis without AT richness in arrays *********************************************\n\n";
	}	

elsif($s==5)
	{
		$check_mfe=0;
		print "\n********************************************************   Analysis without MFE in repeats *********************************************\n\n";
	}
elsif($s==6)
	{
		$check_longer_leader=0;
		print "\n********************************************************   Analysis without longer leader analysis in arrays *********************************************\n\n";
	}	
	
elsif($s==7)
	{
		$check_at_richness=0;
		$check_degeneracy=0;
		$check_longer_leader=0;		
		print "\n********************************************************   Analysis without AT richness, Degeneracy, and longer leader methods *********************************************\n\n";
	}
elsif($s==8)
	{
		$check_motif=0;
		$check_A_and_T_in_repeats=0;
		$check_mfe=0;
		print "\n********************************************************   Analysis without Motif,A&T ratio, and MFE analysis in arrays *********************************************\n\n";
	}	
	
		
if($s==9)
	{
		$check_motif=1;		
		print "\n********************************************************   Analysis with only Motif method*********************************************\n\n";
	}		

if($s==10)
	{
		$check_A_and_T_in_repeats=1;		
		print "\n********************************************************   Analysis with only check_A_and_T_in_repeats method*********************************************\n\n";
	}
if($s==11)
	{
		$check_mfe=1;		
		print "\n********************************************************   Analysis with only MFE method *********************************************\n\n";
	}	
	
if($s==12)
	{
		$check_at_richness=1;		
		print "\n********************************************************   Analysis with only check_at_richness method*********************************************\n\n";
	}		

if($s==13)
	{
		$check_degeneracy=1;		
		print "\n********************************************************   Analysis with only check_degeneracy method*********************************************\n\n";
	}
if($s==14)
	{
		$check_longer_leader=1;		
		print "\n********************************************************   Analysis with only check_longer_leader method *********************************************\n\n";
	}							
############## no need to change anything below ######################################################################





#---- run whole process 100 times and store their corresponding Standard Deviation and Standard error ---------------



#---- initialize the PPVs ---------------------------------------------------------------------------

my $at_richness_ppv=0;
my $degeneracy_ppv=0;
my $mfe_ppv=0;
my $longer_leader_ppv=0;



###################################### training methods ######################################################
my %hash_of_repeat_keys_and_motif_found; 
if($check_motif==1)
	{
		print "\nMotif analysis for training set:\n";
		
		my %hash_of_motif_varients;
		
		#my $motif="ATTGAAA[C|G]";
		
		
		my @arr_motif_bases=split('',"ATTGAAA");
		$hash_of_motif_varients{"ATTGAAA.?"}=0;
		for(my $i=0;$i<=$#arr_motif_bases;$i++)
			{
				my $string="";
				#---add the bases before $i
				for(my $j=0;$j<$i;$j++)
					{
						$string=$string.$arr_motif_bases[$j];
					}
				#--- for $i, add \s
				$string=$string.".";
				#--- now add the rest ---------
				for(my $k=$i+1;$k<=$#arr_motif_bases;$k++)
					{
						$string=$string.$arr_motif_bases[$k];
					}	
				
				$string=$string.".?";	
				#print "$string\n";	
				
				$hash_of_motif_varients{$string}=$i+1;
			}
		$hash_of_motif_varients{"ATTGAAA[C|G]"}=10;	
		$hash_of_motif_varients{"[A|T][A|T][A|T]GAAA.?"}=11;
		#$hash_of_motif_varients{"GAAA.?"}=12;	
		print "\n\n";
							#exit;
							#$hash_of_motif_varients{"ATTGAAA[C|G]"}=1;
							#$hash_of_motif_varients{"ATTGAAA[A|C|G|T]"}=2;
							#$hash_of_motif_varients{"A[A|C|G|T]TGAAA[A|C|G|T]"}=3;
							#$hash_of_motif_varients{"AT[A|C|G|T]GAAA[A|C|G|T]"}=4;
							#$hash_of_motif_varients{"ATTGAAA[A|C|G|T][A|C|G|T]"}=5;
		
		my %hash_of_motif_and_matches;
		my %hash_of_motif_and_false_matches;
		my %hash_of_repeats_and_referrence_direction;

		
		
		#---- get 100 random keys from the current 230 arrays-----------
		my $total_training_sets=0;
		my %training_sets;
		while($total_training_sets<100)
			{
				#--- now randomly pick a record from %hash_of_all_records
				my $random_no=int(rand(229));
									
				$random_no=$random_no+1;
										
				my $current_record=$hash_of_all_records{$random_no};
										
				#----- skip if not part of training set -------------------------------------
				if(defined $hash_of_records_used_for_training{$current_record} and not defined $training_sets{$current_record})
					{
						$training_sets{$current_record}=1;
						$total_training_sets++;
					}
				#----------------------------------------------------------------------------		
				
			}
		
		my $rec_count=1;
		foreach my $key(sort keys %hash_of_all_repeats_and_ref_direction)
			{
				#----- skip if not part of training set -------------------------------------
				if(not defined $hash_of_records_used_for_training{$key}){next;}
				elsif(not defined $training_sets{$key}){next;}
				else{
						#print "$rec_count:$key\n";
						$rec_count++;
					}
				#----------------------------------------------------------------------------
								
				foreach my $repeat(sort keys %{$hash_of_all_repeats_and_ref_direction{$key}})
					{
						#--- replace all U with T if present -----------
						$repeat=~tr/U/T/;
						my $repeat_rev_comp=$repeat; $repeat_rev_comp= reverse $repeat_rev_comp; $repeat_rev_comp=~tr/ACGT/TGCA/;
						
						#---- extract last 8 bases from the repeat and its other strand -----------------------
						my $sub_repeat=substr($repeat,-9);
						my $sub_repeat_rev_comp=substr($repeat_rev_comp,-9);
						#--- now check if the repeat or the other strand contains the motif in the last 10 bases 
						my $suggested_direction="NA";
						
						#---- now get the reference of this repeat and check if it's correct ----------
						my $ref_direction=$hash_of_all_repeats_and_ref_direction{$key}{$repeat};
						#----- search the full motif -----------------------------------------------------------
						
						foreach my $motif(sort{$hash_of_motif_varients{$a}<=>$hash_of_motif_varients{$b}} keys %hash_of_motif_varients)
							{
								$suggested_direction="NA";
								
								if($repeat=~/$motif/ and $repeat_rev_comp=~/$motif/){$suggested_direction="NA";}
								elsif($repeat=~/$motif/ ){$suggested_direction="F";}
								elsif($repeat_rev_comp=~/$motif/){$suggested_direction="R";}
								
								if($suggested_direction ne "NA")
									{
										if($suggested_direction eq $ref_direction)
											{
												$hash_of_motif_and_matches{$motif}{$key}=$suggested_direction;
												
												if($motif=~/^ATTGAAA\.\?/)
													{
														#print "$motif\n";
														$hash_of_repeat_keys_and_motif_found{$key}=$suggested_direction;
													}	
											}
										else{
												$hash_of_motif_and_false_matches{$motif}{$key}=$suggested_direction;
											}	
										
										#last;
										
										#print "True prediction: $key       \t$motif\t$ref_direction\t$suggested_direction\t$repeat\n";
									}
		
							}

					
						#print "$suggested_direction\t$ref_direction\n";	
							
						if($suggested_direction!~/NA/ and $suggested_direction ne $ref_direction)
							{
								#print "False prediction: $key       \t$motif\t$ref_direction\t$suggested_direction\t$repeat\n";
							}
						elsif($suggested_direction!~/NA/)
							{
								#print "True prediction: $key       \t$motif\t$ref_direction\t$suggested_direction\t$repeat\n";
							}					
					}
			}

		#------- now print the ---------------------------------------------------
		print "Program identified the following CRISPR arrays to contain ATTGAAAN at the 3' end of the array:\n";
		foreach my $motif(sort{$hash_of_motif_varients{$a}<=>$hash_of_motif_varients{$b}} keys %hash_of_motif_and_matches)
			{					
				if($motif!~/^ATTGAAA\.\?/){next;}
				
				#print "\n$motif : TP in the following arrays:\n";
				my $count=0;
				foreach my $key(keys %{$hash_of_motif_and_matches{$motif}})
					{
						
						#print "\t$count:$key      \t$motif\t$hash_of_ref_repeat_direction{$key}\t$hash_of_motif_and_matches{$motif}{$key}\n";
						
						$count++;	
					}
				print "\t$motif : Total true prediction = $count\n";	
				
				#print "$motif : FP in the following arrays:\n";
				my $count2=0;
				foreach my $key(keys %{$hash_of_motif_and_false_matches{$motif}})
					{

						#print "\t$count2:$key      \t$motif\t$hash_of_ref_repeat_direction{$key}\t$hash_of_motif_and_false_matches{$motif}{$key}\n";

						$count2++;	
					}
				print "\t$motif : Total false predictions = $count2\n";	
				
				my $count3=100-$count-$count2;
				print "\t$motif : Total no predictions = $count3\n";	
			}	
		my $tmp_tot1=keys %hash_of_repeat_keys_and_motif_found;
		#print "\n\tTotal: $tmp_tot1 arrays\n\n";	
		
	}

if($check_motif==1)
	{
		
		print "\nMotif analysis for test set:\n";
		
		my %hash_of_motif_varients;
		
		#my $motif="ATTGAAA[C|G]";
		
		
		my @arr_motif_bases=split('',"ATTGAAA");
		$hash_of_motif_varients{"ATTGAAA.?"}=0;
		for(my $i=0;$i<=$#arr_motif_bases;$i++)
			{
				my $string="";
				#---add the bases before $i
				for(my $j=0;$j<$i;$j++)
					{
						$string=$string.$arr_motif_bases[$j];
					}
				#--- for $i, add \s
				$string=$string.".";
				#--- now add the rest ---------
				for(my $k=$i+1;$k<=$#arr_motif_bases;$k++)
					{
						$string=$string.$arr_motif_bases[$k];
					}	
				
				$string=$string.".?";	
				#print "$string\n";	
				
				$hash_of_motif_varients{$string}=$i+1;
			}
		$hash_of_motif_varients{"ATTGAAA[C|G]"}=10;	
		$hash_of_motif_varients{"[A|T][A|T][A|T]GAAA.?"}=11;
		#$hash_of_motif_varients{"GAAA.?"}=12;	
		print "\n\n";
							#exit;
							#$hash_of_motif_varients{"ATTGAAA[C|G]"}=1;
							#$hash_of_motif_varients{"ATTGAAA[A|C|G|T]"}=2;
							#$hash_of_motif_varients{"A[A|C|G|T]TGAAA[A|C|G|T]"}=3;
							#$hash_of_motif_varients{"AT[A|C|G|T]GAAA[A|C|G|T]"}=4;
							#$hash_of_motif_varients{"ATTGAAA[A|C|G|T][A|C|G|T]"}=5;
		
		my %hash_of_motif_and_matches;
		my %hash_of_motif_and_false_matches;
		my %hash_of_repeats_and_referrence_direction;

		
		
		#---- get 100 random keys from the current 230 arrays-----------
		my $total_training_sets=0;
		my %training_sets;
		while($total_training_sets<100)
			{
				#--- now randomly pick a record from %hash_of_all_records
				my $random_no=int(rand(229));
									
				$random_no=$random_no+1;
										
				my $current_record=$hash_of_all_records{$random_no};
										
				#----- skip if not part of training set -------------------------------------
				if(defined $hash_of_records_used_for_training{$current_record} and not defined $training_sets{$current_record})
					{
						$training_sets{$current_record}=1;
						$total_training_sets++;
					}
				#----------------------------------------------------------------------------		
				
			}
		
		my $rec_count=1;
		foreach my $key(sort keys %hash_of_all_repeats_and_ref_direction)
			{
				#----- skip if not part of training set -------------------------------------
				if(not defined $hash_of_records_used_for_training{$key}){next;}
				elsif(not defined $training_sets{$key}){next;}
				else{
						#print "$rec_count:$key\n";
						$rec_count++;
					}
				#----------------------------------------------------------------------------
								
				foreach my $repeat(sort keys %{$hash_of_all_repeats_and_ref_direction{$key}})
					{
						#--- replace all U with T if present -----------
						$repeat=~tr/U/T/;
						my $repeat_rev_comp=$repeat; $repeat_rev_comp= reverse $repeat_rev_comp; $repeat_rev_comp=~tr/ACGT/TGCA/;
						
						#---- extract last 8 bases from the repeat and its other strand -----------------------
						my $sub_repeat=substr($repeat,-9);
						my $sub_repeat_rev_comp=substr($repeat_rev_comp,-9);
						#--- now check if the repeat or the other strand contains the motif in the last 10 bases 
						my $suggested_direction="NA";
						
						#---- now get the reference of this repeat and check if it's correct ----------
						my $ref_direction=$hash_of_all_repeats_and_ref_direction{$key}{$repeat};
						#----- search the full motif -----------------------------------------------------------
						
						foreach my $motif(sort{$hash_of_motif_varients{$a}<=>$hash_of_motif_varients{$b}} keys %hash_of_motif_varients)
							{
								$suggested_direction="NA";
								
								if($repeat=~/$motif/ and $repeat_rev_comp=~/$motif/){$suggested_direction="NA";}
								elsif($repeat=~/$motif/ ){$suggested_direction="F";}
								elsif($repeat_rev_comp=~/$motif/){$suggested_direction="R";}
								
								if($suggested_direction ne "NA")
									{
										if($suggested_direction eq $ref_direction)
											{
												$hash_of_motif_and_matches{$motif}{$key}=$suggested_direction;
												
												if($motif=~/^ATTGAAA\.\?/)
													{
														#print "$motif\n";
														$hash_of_repeat_keys_and_motif_found{$key}=$suggested_direction;
													}	
											}
										else{
												$hash_of_motif_and_false_matches{$motif}{$key}=$suggested_direction;
											}	
										
										#last;
										
										#print "True prediction: $key       \t$motif\t$ref_direction\t$suggested_direction\t$repeat\n";
									}
		
							}

					
						#print "$suggested_direction\t$ref_direction\n";	
							
						if($suggested_direction!~/NA/ and $suggested_direction ne $ref_direction)
							{
								#print "False prediction: $key       \t$motif\t$ref_direction\t$suggested_direction\t$repeat\n";
							}
						elsif($suggested_direction!~/NA/)
							{
								#print "True prediction: $key       \t$motif\t$ref_direction\t$suggested_direction\t$repeat\n";
							}					
					}
			}

		#------- now print the ---------------------------------------------------
		print "Program identified the following CRISPR arrays to contain ATTGAAAN at the 3' end of the array:\n";
		foreach my $motif(sort{$hash_of_motif_varients{$a}<=>$hash_of_motif_varients{$b}} keys %hash_of_motif_and_matches)
			{					
				if($motif!~/^ATTGAAA\.\?/){next;}
				
				#print "\n$motif : TP in the following arrays:\n";
				my $count=0;
				foreach my $key(keys %{$hash_of_motif_and_matches{$motif}})
					{
						
						#print "\t$count:$key      \t$motif\t$hash_of_ref_repeat_direction{$key}\t$hash_of_motif_and_matches{$motif}{$key}\n";
						
						$count++;	
					}
				print "\t$motif : Total true prediction = $count\n";	
				
				#print "$motif : FP in the following arrays:\n";
				my $count2=0;
				foreach my $key(keys %{$hash_of_motif_and_false_matches{$motif}})
					{

						#print "\t$count2:$key      \t$motif\t$hash_of_ref_repeat_direction{$key}\t$hash_of_motif_and_false_matches{$motif}{$key}\n";

						$count2++;	
					}
				print "\t$motif : Total false predictions = $count2\n";	
				
				my $count3=100-$count-$count2;
				print "\t$motif : Total no predictions = $count3\n";
			}	
		my $tmp_tot1=keys %hash_of_repeat_keys_and_motif_found;
		#print "\n\tTotal: $tmp_tot1 arrays\n\n";
		
		
		
	}



#exit;


#---------- check A_and_T distribution in repeats ----------------------------------------------
my %hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats;
if($check_A_and_T_in_repeats==1)
	{
		
		print "\nProgram is going to check A and T ratio in training set (repeats):\n";
	
		my %hash_of_A_and_T_results;
		
		
		
		for(my $check=0;$check<=2;$check++)
			{		
				#print "\n";
		
		
				for(my $p=1;$p<=$q;$p++)
					{	
						
						#print "The following records are used in A&T ratio training:$p:\n";
						my $rec_count=1;			
						#-------now pick 100 keys from the sample --------------------	

					my %hash_of_current_sample;
					my %hash_of_ref_repeat_direction_current_sample;
					my %hash_of_unique_repeats;
					for(1..$r*90)
						{
							my $found=0;
							while($found!=1)
							{
								my $random_no=int(rand($total_records-1));
							
								$random_no=$random_no+1;
								
								my $current_record=$hash_of_all_records{$random_no};
								
								#----- skip if not part of training set -------------------------------------
								if(not defined $hash_of_records_used_for_training{$current_record}){next;}
								#----------------------------------------------------------------------------
								
								
								if(not defined $hash_of_current_sample{$current_record})
									{
										my $key=$current_record; 
										my $rec_looks_good=0;
										
										
										
										#------ skip the ones with motif ---------------------
										#if($hash_of_repeat_keys_and_motif_found{$key}){next;}
										#------------------------------------------------------
										
										#$key=~s/-/\t/g;
											#print "$key\n"; 
										#----- check if the prediction is not NA -------------------------------------------------------------------
										my $rec_line=`grep -w '$key' sample_space_acc_repeat_and_their_reference.txt >&1`; chomp $rec_line;$rec_line=~s/\r//g;
											#print "grep '$key' sample_space_acc_repeat_and_their_reference.txt\n";exit;
										my @arr_tmp1=split('\t',$rec_line);
										
										my $repeat=$arr_tmp1[1]; chomp $repeat;$repeat=~s/\s+//g;
										
										
										my $no_of_As_in_repeat=$repeat=~s/A/A/gi;
										my $no_of_Ts_in_repeat=$repeat=~s/T/T/gi;	
										#my $A_and_T_suggestion="NA";
										
										my $at_p_d=(($no_of_As_in_repeat+$no_of_Ts_in_repeat)/length($repeat))*100; $at_p_d=sprintf("%.2f",$at_p_d);
										my $gc_p_d=((length($repeat) -$no_of_As_in_repeat-$no_of_Ts_in_repeat)/length($repeat))*100; $gc_p_d=sprintf("%.2f",$gc_p_d);
										#------ get 100 valid cases -----------------------
										
										if($check==0)
											{
												$rec_looks_good=1;
											}
										elsif($check==1 and $at_p_d>$gc_p_d)
											{
												
												$rec_looks_good=1;
												
												
												#--- check if the AT% is > GC%
											}
										elsif($check==2 and $at_p_d<=$gc_p_d)
											{												
												$rec_looks_good=1;					
												
											}	
										
										#---------------------------------------------------------------------	
										#if($hash_of_unique_repeats{$repeat}){next;}
										#----------------------------------------------------------------------
										if($rec_looks_good==1)
											{
												#print "$rec_count:$current_record\n";
												$rec_count++;
												
												$hash_of_current_sample{$current_record}=1;
												$hash_of_unique_repeats{$repeat}=1;
												##print "$current_record\n";
												$found=1;
											}
									}
							}
							
							my $current_no_of_entries=keys %hash_of_current_sample;
							if($current_no_of_entries>=$r)
								{
									my $no_of_unique_repeats_in_sample=keys %hash_of_unique_repeats;
									#print "$p:\$current_no_of_entries=$current_no_of_entries with \$no_of_unique_repeats_in_sample=$no_of_unique_repeats_in_sample\n";
									last;
								}
						}
							
					#print "OK... sample selected\n";
							
					my $A_and_T_TP_count=0;
					my $A_and_T_FP_count=0;
					my $A_and_T_NP_count=0;
						
					foreach my $key(sort keys %hash_of_all_repeats_and_ref_direction)
						{
								
									
									#if($hash_of_repeat_keys_and_motif_found{$key}){next;}
									#--------- redundant/non-redundant toggle ------------------
									if(not $hash_of_current_sample{$key}){next;}#print "OK\n";
									#-----------------------------------------------------------
									foreach my $repeat(sort keys %{$hash_of_all_repeats_and_ref_direction{$key}})
										{
											
											#---- now get the reference of this repeat and check if it's correct ----------
											my $ref_direction=$hash_of_all_repeats_and_ref_direction{$key}{$repeat};
										
											#--- replace all U with T if present -----------
											$repeat=~tr/U/T/;

											##################### check direction using no. of As and Ts: mainly applicable to Archaeal genomes #######################		
												
											my $no_of_As_in_repeat=$repeat=~s/A/A/gi;
											my $no_of_Ts_in_repeat=$repeat=~s/T/T/gi;	
											my $A_and_T_suggestion="NA";
											
											my $at_p_d=(($no_of_As_in_repeat+$no_of_Ts_in_repeat)/length($repeat))*100; $at_p_d=sprintf("%.2f",$at_p_d);
											my $gc_p_d=((length($repeat) -$no_of_As_in_repeat-$no_of_Ts_in_repeat)/length($repeat))*100; $gc_p_d=sprintf("%.2f",$gc_p_d);
											
											
											
											
											
											if($check==0)#----- first check A and T count without checking the AT distribution ----
												{
													if($no_of_As_in_repeat>$no_of_Ts_in_repeat)
														{
															$A_and_T_suggestion="F";
														}
													elsif($no_of_Ts_in_repeat>$no_of_As_in_repeat)
														{
															$A_and_T_suggestion="R";
														}	
													
													#print "$A_and_T_suggestion\n";
													#------------ save the prediction in a hash -------------------------------------------------
													#$hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats{$key}=$A_and_T_suggestion;	
													#--------------------------------------------------------------------------------------------
													
													if($A_and_T_suggestion ne "NA")
														{								
															if($A_and_T_suggestion eq $ref_direction)
																{
																	$A_and_T_TP_count++;
																}
															else{
																	$A_and_T_FP_count++;
																}	
														}
													else{
															$A_and_T_NP_count++;
														}
												}
											elsif($check==1 and $at_p_d>$gc_p_d)
												{								
													#print "$repeat    \t is AT rich\t $ref_direction [$at_p_d,$gc_p_d] \t";
														
													if($no_of_As_in_repeat>$no_of_Ts_in_repeat)
														{
															$A_and_T_suggestion="F";
														}
													elsif($no_of_Ts_in_repeat>$no_of_As_in_repeat)
														{
															$A_and_T_suggestion="R";
														}	
													
													#print "$A_and_T_suggestion\n";
													#------------ save the prediction in a hash -------------------------------------------------
													#$hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats{$key}=$A_and_T_suggestion;	
													#--------------------------------------------------------------------------------------------
													
													if($A_and_T_suggestion ne "NA")
														{								
															if($A_and_T_suggestion eq $ref_direction)
																{
																	$A_and_T_TP_count++;
																}
															else{
																	$A_and_T_FP_count++;
																}	
														}
													else{
															$A_and_T_NP_count++;
														}
												}
											elsif($check==2 and $at_p_d<=$gc_p_d)
												{
													#print "$repeat    \t is GC rich\t $ref_direction [$at_p_d,$gc_p_d] \t"; 
													
													if($no_of_As_in_repeat>$no_of_Ts_in_repeat)
														{
															$A_and_T_suggestion="F";
														}
													elsif($no_of_Ts_in_repeat>$no_of_As_in_repeat)
														{
															$A_and_T_suggestion="R";
														}	
													
													#print "$A_and_T_suggestion\n";
													#------------ save the prediction in a hash -------------------------------------------------
													$hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats{$key}=$A_and_T_suggestion;	
													#--------------------------------------------------------------------------------------------
													
													if($A_and_T_suggestion ne "NA")
														{								
															if($A_and_T_suggestion eq $ref_direction)
																{
																	$A_and_T_TP_count++;
																}
															else{
																	$A_and_T_FP_count++;
																}	
														}
													else{
															$A_and_T_NP_count++;
														}
												}
											else{
													#print "\t$key\t with $repeat\t [$at_p_d,$gc_p_d]\n";
												}			
										}
						}

							
						
						#-------------- print Summary of direcional result based on the distribution of As and Ts in repeats ------------------------------------------
						my $A_and_T_PPV=$A_and_T_TP_count/($A_and_T_TP_count+$A_and_T_FP_count); $A_and_T_PPV=sprintf("%.2f",$A_and_T_PPV);
						my $A_and_T_total=$A_and_T_TP_count+$A_and_T_FP_count+$A_and_T_NP_count;
							#print "\n\n Summary of direcional result based on the distribution of As and Ts in repeats:\n";
						#print "CHECK=$check\tSet=$p\tTotal:$A_and_T_total\tTP:$A_and_T_TP_count\tFP:$A_and_T_FP_count\tNP:$A_and_T_NP_count\tPPV:$A_and_T_PPV\n";
						
						if($check==2){$mean_A_T_distribution_in_repeats_PPV=$A_and_T_PPV;}
						
						$hash_of_A_and_T_results{$check}{'TOTAL'}{$p}=$A_and_T_total;
						$hash_of_A_and_T_results{$check}{'PPV'}{$p}=$A_and_T_PPV;
						$hash_of_A_and_T_results{$check}{'TP'}{$p}=$A_and_T_TP_count;
						$hash_of_A_and_T_results{$check}{'FP'}{$p}=$A_and_T_FP_count;
						$hash_of_A_and_T_results{$check}{'NP'}{$p}=$A_and_T_NP_count;
						
				}
			#print "Finishes 100 tests for $check: \n";	
			
		}	#--- end of $check ------
		
	
	#--------- now calculate the SD and SE ------------------------------------------
	#print "Summary of Potential for RNA secondary structure analysis:\n";
	print "Test-type\tAvg_total_arrays_processed\tTP\ttp_sd\ttp_se\tFP\tfp_sd\tfp_se\tNP\tnp_sd\tnp_se\tPPV\tppv_sd\tppv_se\n";
	foreach my $check(sort{$a <=> $b} keys %hash_of_A_and_T_results)
		{
			my $avg_total_arrays_processed=0;		
			my $ppv_avg=0;	my $ppv_sd=0;	my $ppv_se=0;		
			my $tp_avg=0;	my $tp_sd=0;	my $tp_se=0;		
			my $fp_avg=0;	my $fp_sd=0;	my $fp_se=0;		
			my $np_avg=0;	my $np_sd=0;	my $np_se=0;
			
			
			foreach my $type(sort{$a cmp $b} keys %{$hash_of_A_and_T_results{$check}})
				{				
					#print "$type=$type\n";
					my @arr_data;
					
					foreach my $pass(sort{$a <=> $b} keys %{$hash_of_A_and_T_results{$check}{$type}})
						{
							push(@arr_data,$hash_of_A_and_T_results{$check}{$type}{$pass});
							#print "@arr_data\n";
							#print "$analysis\t$type\t$pass\t$hash_of_analysis_results{$analysis}{$type}{$pass}\n";
						}
						
					#--- now calculate the mean, SD and SE ------------------
					my $ave;my $std; my $std_err=0;
					
					$ave = &average(\@arr_data);
					
					if($type!~/PPV/){$ave=sprintf("%.0f",$ave);}
						
					if($type!~/TOTAL/)
						{
							$std = &stdev(\@arr_data);	
							
							if($std !=0 and $ave!=0)
								{				
									$std_err=($std/($ave ** 0.5)); $std_err=sprintf("%.2f",$std_err);
								}
							else{
									$std_err=0;
								}	
						}
					if($type=~/TOTAL/)
						{
							$avg_total_arrays_processed=int($ave);
							#print "Total=$avg_total_arrays_processed\n";
							next;
						}
					elsif($type=~/PPV/){$ppv_avg=$ave;$ppv_sd=$std;$ppv_se=$std_err;}
					elsif($type=~/TP/){$tp_avg=$ave;$tp_sd=$std;$tp_se=$std_err;}
					elsif($type=~/FP/){$fp_avg=$ave;$fp_sd=$std;$fp_se=$std_err;}
					elsif($type=~/NP/){$np_avg=$ave;$np_sd=$std;$np_se=$std_err;}
					#print "$check\t$type\t$ave\t$std\n";
					
					if($type=~/PPV/ and $check==$mfe_target_base){$mean_mfe_ppv=$ppv_avg;}
				}
			
			my $test_type="";
			if($check==0){$test_type="All repeats";}
			elsif($check==1){$test_type="Repeats with AT%>GC%";}
			elsif($check==2){$test_type="Repeats with AT%<=GC%";}
			
			print "$test_type\t$avg_total_arrays_processed\t$tp_avg\t$tp_sd\t$tp_se\t$fp_avg\t$fp_sd\t$fp_se\t$np_avg\t$np_sd\t$np_se\t$ppv_avg\t$ppv_sd\t$ppv_se\n";	
		}
		
	
	} # end of $check_A_and_T_in_repeats==1


if($check_A_and_T_in_repeats==1)
	{
		
		print "\nProgram is going to check A and T ratio in test set (repeats):\n";
	
		my %hash_of_A_and_T_results;
		
		
		
		for(my $check=0;$check<=2;$check++)
			{		
				#print "\n";
		
		
				for(my $p=1;$p<=$q;$p++)
					{	
						
						#print "The following records are used in A&T ratio training:$p:\n";
						my $rec_count=1;			
						#-------now pick 100 keys from the sample --------------------	

					my %hash_of_current_sample;
					my %hash_of_ref_repeat_direction_current_sample;
					my %hash_of_unique_repeats;
					for(1..$r*90)
						{
							my $found=0;
							while($found!=1)
							{
								my $random_no=int(rand($total_records-1));
							
								$random_no=$random_no+1;
								
								my $current_record=$hash_of_all_records{$random_no};
								
								#----- skip if not part of training set -------------------------------------
								if(not defined $hash_of_records_used_for_training{$current_record}){next;}
								#----------------------------------------------------------------------------
								
								
								if(not defined $hash_of_current_sample{$current_record})
									{
										my $key=$current_record; 
										my $rec_looks_good=0;
										
										
										
										#------ skip the ones with motif ---------------------
										#if($hash_of_repeat_keys_and_motif_found{$key}){next;}
										#------------------------------------------------------
										
										#$key=~s/-/\t/g;
											#print "$key\n"; 
										#----- check if the prediction is not NA -------------------------------------------------------------------
										my $rec_line=`grep -w '$key' sample_space_acc_repeat_and_their_reference.txt >&1`; chomp $rec_line;$rec_line=~s/\r//g;
											#print "grep '$key' sample_space_acc_repeat_and_their_reference.txt\n";exit;
										my @arr_tmp1=split('\t',$rec_line);
										
										my $repeat=$arr_tmp1[1]; chomp $repeat;$repeat=~s/\s+//g;
										
										
										my $no_of_As_in_repeat=$repeat=~s/A/A/gi;
										my $no_of_Ts_in_repeat=$repeat=~s/T/T/gi;	
										#my $A_and_T_suggestion="NA";
										
										my $at_p_d=(($no_of_As_in_repeat+$no_of_Ts_in_repeat)/length($repeat))*100; $at_p_d=sprintf("%.2f",$at_p_d);
										my $gc_p_d=((length($repeat) -$no_of_As_in_repeat-$no_of_Ts_in_repeat)/length($repeat))*100; $gc_p_d=sprintf("%.2f",$gc_p_d);
										#------ get 100 valid cases -----------------------
										
										if($check==0)
											{
												$rec_looks_good=1;
											}
										elsif($check==1 and $at_p_d>$gc_p_d)
											{
												
												$rec_looks_good=1;
												
												
												#--- check if the AT% is > GC%
											}
										elsif($check==2 and $at_p_d<=$gc_p_d)
											{												
												$rec_looks_good=1;					
												
											}	
										
										#---------------------------------------------------------------------	
										#if($hash_of_unique_repeats{$repeat}){next;}
										#----------------------------------------------------------------------
										if($rec_looks_good==1)
											{
												#print "$rec_count:$current_record\n";
												$rec_count++;
												
												$hash_of_current_sample{$current_record}=1;
												$hash_of_unique_repeats{$repeat}=1;
												##print "$current_record\n";
												$found=1;
											}
									}
							}
							
							my $current_no_of_entries=keys %hash_of_current_sample;
							if($current_no_of_entries>=$r)
								{
									my $no_of_unique_repeats_in_sample=keys %hash_of_unique_repeats;
									#print "$p:\$current_no_of_entries=$current_no_of_entries with \$no_of_unique_repeats_in_sample=$no_of_unique_repeats_in_sample\n";
									last;
								}
						}
							
					#print "OK... sample selected\n";
							
					my $A_and_T_TP_count=0;
					my $A_and_T_FP_count=0;
					my $A_and_T_NP_count=0;
						
					foreach my $key(sort keys %hash_of_all_repeats_and_ref_direction)
						{
								
									
									#if($hash_of_repeat_keys_and_motif_found{$key}){next;}
									#--------- redundant/non-redundant toggle ------------------
									if(not $hash_of_current_sample{$key}){next;}#print "OK\n";
									#-----------------------------------------------------------
									foreach my $repeat(sort keys %{$hash_of_all_repeats_and_ref_direction{$key}})
										{
											
											#---- now get the reference of this repeat and check if it's correct ----------
											my $ref_direction=$hash_of_all_repeats_and_ref_direction{$key}{$repeat};
										
											#--- replace all U with T if present -----------
											$repeat=~tr/U/T/;

											##################### check direction using no. of As and Ts: mainly applicable to Archaeal genomes #######################		
												
											my $no_of_As_in_repeat=$repeat=~s/A/A/gi;
											my $no_of_Ts_in_repeat=$repeat=~s/T/T/gi;	
											my $A_and_T_suggestion="NA";
											
											my $at_p_d=(($no_of_As_in_repeat+$no_of_Ts_in_repeat)/length($repeat))*100; $at_p_d=sprintf("%.2f",$at_p_d);
											my $gc_p_d=((length($repeat) -$no_of_As_in_repeat-$no_of_Ts_in_repeat)/length($repeat))*100; $gc_p_d=sprintf("%.2f",$gc_p_d);
											
											
											
											
											
											if($check==0)#----- first check A and T count without checking the AT distribution ----
												{
													if($no_of_As_in_repeat>$no_of_Ts_in_repeat)
														{
															$A_and_T_suggestion="F";
														}
													elsif($no_of_Ts_in_repeat>$no_of_As_in_repeat)
														{
															$A_and_T_suggestion="R";
														}	
													
													#print "$A_and_T_suggestion\n";
													#------------ save the prediction in a hash -------------------------------------------------
													#$hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats{$key}=$A_and_T_suggestion;	
													#--------------------------------------------------------------------------------------------
													
													if($A_and_T_suggestion ne "NA")
														{								
															if($A_and_T_suggestion eq $ref_direction)
																{
																	$A_and_T_TP_count++;
																}
															else{
																	$A_and_T_FP_count++;
																}	
														}
													else{
															$A_and_T_NP_count++;
														}
												}
											elsif($check==1 and $at_p_d>$gc_p_d)
												{								
													#print "$repeat    \t is AT rich\t $ref_direction [$at_p_d,$gc_p_d] \t";
														
													if($no_of_As_in_repeat>$no_of_Ts_in_repeat)
														{
															$A_and_T_suggestion="F";
														}
													elsif($no_of_Ts_in_repeat>$no_of_As_in_repeat)
														{
															$A_and_T_suggestion="R";
														}	
													
													#print "$A_and_T_suggestion\n";
													#------------ save the prediction in a hash -------------------------------------------------
													#$hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats{$key}=$A_and_T_suggestion;	
													#--------------------------------------------------------------------------------------------
													
													if($A_and_T_suggestion ne "NA")
														{								
															if($A_and_T_suggestion eq $ref_direction)
																{
																	$A_and_T_TP_count++;
																}
															else{
																	$A_and_T_FP_count++;
																}	
														}
													else{
															$A_and_T_NP_count++;
														}
												}
											elsif($check==2 and $at_p_d<=$gc_p_d)
												{
													#print "$repeat    \t is GC rich\t $ref_direction [$at_p_d,$gc_p_d] \t"; 
													
													if($no_of_As_in_repeat>$no_of_Ts_in_repeat)
														{
															$A_and_T_suggestion="F";
														}
													elsif($no_of_Ts_in_repeat>$no_of_As_in_repeat)
														{
															$A_and_T_suggestion="R";
														}	
													
													#print "$A_and_T_suggestion\n";
													#------------ save the prediction in a hash -------------------------------------------------
													$hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats{$key}=$A_and_T_suggestion;	
													#--------------------------------------------------------------------------------------------
													
													if($A_and_T_suggestion ne "NA")
														{								
															if($A_and_T_suggestion eq $ref_direction)
																{
																	$A_and_T_TP_count++;
																}
															else{
																	$A_and_T_FP_count++;
																}	
														}
													else{
															$A_and_T_NP_count++;
														}
												}
											else{
													#print "\t$key\t with $repeat\t [$at_p_d,$gc_p_d]\n";
												}			
										}
						}

							
						
						#-------------- print Summary of direcional result based on the distribution of As and Ts in repeats ------------------------------------------
						my $A_and_T_PPV=$A_and_T_TP_count/($A_and_T_TP_count+$A_and_T_FP_count); $A_and_T_PPV=sprintf("%.2f",$A_and_T_PPV);
						my $A_and_T_total=$A_and_T_TP_count+$A_and_T_FP_count+$A_and_T_NP_count;
							#print "\n\n Summary of direcional result based on the distribution of As and Ts in repeats:\n";
						#print "CHECK=$check\tSet=$p\tTotal:$A_and_T_total\tTP:$A_and_T_TP_count\tFP:$A_and_T_FP_count\tNP:$A_and_T_NP_count\tPPV:$A_and_T_PPV\n";
						
						if($check==2){$mean_A_T_distribution_in_repeats_PPV=$A_and_T_PPV;}
						
						$hash_of_A_and_T_results{$check}{'TOTAL'}{$p}=$A_and_T_total;
						$hash_of_A_and_T_results{$check}{'PPV'}{$p}=$A_and_T_PPV;
						$hash_of_A_and_T_results{$check}{'TP'}{$p}=$A_and_T_TP_count;
						$hash_of_A_and_T_results{$check}{'FP'}{$p}=$A_and_T_FP_count;
						$hash_of_A_and_T_results{$check}{'NP'}{$p}=$A_and_T_NP_count;
						
				}
			#print "Finishes 100 tests for $check: \n";	
			
		}	#--- end of $check ------
		
	
	#--------- now calculate the SD and SE ------------------------------------------
	#print "Summary of Potential for RNA secondary structure analysis:\n";
	print "Test-type\tAvg_total_arrays_processed\tTP\ttp_sd\ttp_se\tFP\tfp_sd\tfp_se\tNP\tnp_sd\tnp_se\tPPV\tppv_sd\tppv_se\n";
	foreach my $check(sort{$a <=> $b} keys %hash_of_A_and_T_results)
		{
			my $avg_total_arrays_processed=0;		
			my $ppv_avg=0;	my $ppv_sd=0;	my $ppv_se=0;		
			my $tp_avg=0;	my $tp_sd=0;	my $tp_se=0;		
			my $fp_avg=0;	my $fp_sd=0;	my $fp_se=0;		
			my $np_avg=0;	my $np_sd=0;	my $np_se=0;
			
			
			foreach my $type(sort{$a cmp $b} keys %{$hash_of_A_and_T_results{$check}})
				{				
					#print "$type=$type\n";
					my @arr_data;
					
					foreach my $pass(sort{$a <=> $b} keys %{$hash_of_A_and_T_results{$check}{$type}})
						{
							push(@arr_data,$hash_of_A_and_T_results{$check}{$type}{$pass});
							#print "@arr_data\n";
							#print "$analysis\t$type\t$pass\t$hash_of_analysis_results{$analysis}{$type}{$pass}\n";
						}
						
					#--- now calculate the mean, SD and SE ------------------
					my $ave;my $std; my $std_err=0;
					
					$ave = &average(\@arr_data);
					
					if($type!~/PPV/){$ave=sprintf("%.0f",$ave);}
						
					if($type!~/TOTAL/)
						{
							$std = &stdev(\@arr_data);	
							
							if($std !=0 and $ave!=0)
								{				
									$std_err=($std/($ave ** 0.5)); $std_err=sprintf("%.2f",$std_err);
								}
							else{
									$std_err=0;
								}	
						}
					if($type=~/TOTAL/)
						{
							$avg_total_arrays_processed=int($ave);
							#print "Total=$avg_total_arrays_processed\n";
							next;
						}
					elsif($type=~/PPV/){$ppv_avg=$ave;$ppv_sd=$std;$ppv_se=$std_err;}
					elsif($type=~/TP/){$tp_avg=$ave;$tp_sd=$std;$tp_se=$std_err;}
					elsif($type=~/FP/){$fp_avg=$ave;$fp_sd=$std;$fp_se=$std_err;}
					elsif($type=~/NP/){$np_avg=$ave;$np_sd=$std;$np_se=$std_err;}
					#print "$check\t$type\t$ave\t$std\n";
					
					if($type=~/PPV/ and $check==$mfe_target_base){$mean_mfe_ppv=$ppv_avg;}
				}
			
			my $test_type="";
			if($check==0){$test_type="All repeats";}
			elsif($check==1){$test_type="Repeats with AT%>GC%";}
			elsif($check==2){$test_type="Repeats with AT%<=GC%";}
			
			print "$test_type\t$avg_total_arrays_processed\t$tp_avg\t$tp_sd\t$tp_se\t$fp_avg\t$fp_sd\t$fp_se\t$np_avg\t$np_sd\t$np_se\t$ppv_avg\t$ppv_sd\t$ppv_se\n";	
		}
		
	
	} # end of $check_A_and_T_in_repeats==1



#----- then check AT-richness ---------------------------------------------------------------------
my %hash_of_at_richness_records;
if($check_at_richness==1)
{
	print "\n\nProgram is going to check the distribution of AT in different windows in flanking regions of training set:\n";
	my %hash_of_atr_results;

	for(my $p=1;$p<=$q;$p++)
		{

	
		#print "The following records are used in training:$p:\n";
		my $rec_count=1;
		#-------now pick 100 keys from the sample --------------------	

		my %hash_of_current_sample;
		my %hash_of_ref_repeat_direction_current_sample;
		for(1..$r*90)
			{
				my $found=0;
				while($found!=1)
				{
					my $random_no=int(rand($total_records-1));
				
					$random_no=$random_no+1;
					
					my $current_record=$hash_of_all_records{$random_no};
					
					#----- skip if not part of training set -------------------------------------
					if(not defined $hash_of_records_used_for_training{$current_record}){next;}
		
					#----------------------------------------------------------------------------
													
					if(not defined $hash_of_current_sample{$current_record})
						{
													
							my $key=$current_record; 
							
							#print "$rec_count:$current_record\n";
							#$rec_count++;
									
							#------ skip the ones with motif ---------------------
							#if($hash_of_repeat_keys_and_motif_found{$key}){ next;}
							#------------------------------------------------------
							
							
							$key=~s/-/\t/g;
							
							#print "$key\n"; 
							
							#----- check if the prediction is not NA  for window 135 -------------------------------------------------------------------
							my $rec_line=`grep '$key' at_richness_distribution.txt >&1`;
							#print "grep '$key' at_richness_distribution.txt\n$rec_line\n";
							my @arr_tmp1=split('\t',$rec_line);
							
							if($arr_tmp1[9]=~/NA/){next;}
							
							$hash_of_ref_repeat_direction_current_sample{$current_record}=$hash_of_ref_repeat_direction{$current_record};
							$hash_of_current_sample{$current_record}=1;
							#print "$current_record\n";
							$found=1;
						}
				}
				my $current_no_of_entries=keys %hash_of_current_sample;
				if($current_no_of_entries>=$r)
					{
						#print "\$current_no_of_entries=$current_no_of_entries;";
						last;
					}
			}
		my $t_sample_size=keys %hash_of_current_sample;
		#print "\$t_sample_size=$t_sample_size\n";
		#exit;
		######### check AT_richness#############

			
		for (my $i=1;$i<15;$i++)
			{
				
					my ($total1,$window_1,$ppv1,$tp1,$fp1,$np1,$total2,$window_2,$ppv2,$tp2,$fp2,$np2,$total3,$window_3,$ppv3,$tp3,$fp3,$np3,$total4,$window_4,$ppv4,$tp4,$fp4,$np4,$total5,$window_5,$ppv5,$tp5,$fp5,$np5,$total6,$window_6,$ppv6,$tp6,$fp6,$np6,$total7,$window_7,$ppv7,$tp7,$fp7,$np7,$total8,$window_8,$ppv8,$tp8,$fp8,$np8,$total9,$window_9,$ppv9,$tp9,$fp9,$np9,$total10,$window_10,$ppv10,$tp10,$fp10,$np10,$total11,$window_11,$ppv11,$tp11,$fp11,$np11);

					($total1,$window_1,$ppv1,$tp1,$fp1,$np1,$total2,$window_2,$ppv2,$tp2,$fp2,$np2,$total3,$window_3,$ppv3,$tp3,$fp3,$np3,$total4,$window_4,$ppv4,$tp4,$fp4,$np4,$total5,$window_5,$ppv5,$tp5,$fp5,$np5,$total6,$window_6,$ppv6,$tp6,$fp6,$np6,$total7,$window_7,$ppv7,$tp7,$fp7,$np7,$total8,$window_8,$ppv8,$tp8,$fp8,$np8,$total9,$window_9,$ppv9,$tp9,$fp9,$np9,$total10,$window_10,$ppv10,$tp10,$fp10,$np10,$total11,$window_11,$ppv11,$tp11,$fp11,$np11)=&at_richness_analysis($i,\%hash_of_current_sample,\%hash_of_at_richness_records);
						

					$hash_of_atr_results{$i}{$window_1}{'TOTAL'}{$p}=$total1;
					$hash_of_atr_results{$i}{$window_1}{'PPV'}{$p}=$ppv1;
					$hash_of_atr_results{$i}{$window_1}{'TP'}{$p}=$tp1;
					$hash_of_atr_results{$i}{$window_1}{'FP'}{$p}=$fp1;
					$hash_of_atr_results{$i}{$window_1}{'NP'}{$p}=$np1;
					
					$hash_of_atr_results{$i}{$window_2}{'TOTAL'}{$p}=$total2;
					$hash_of_atr_results{$i}{$window_2}{'PPV'}{$p}=$ppv2;
					$hash_of_atr_results{$i}{$window_2}{'TP'}{$p}=$tp2;
					$hash_of_atr_results{$i}{$window_2}{'FP'}{$p}=$fp2;
					$hash_of_atr_results{$i}{$window_2}{'NP'}{$p}=$np2;
					
					$hash_of_atr_results{$i}{$window_3}{'TOTAL'}{$p}=$total3;
					$hash_of_atr_results{$i}{$window_3}{'PPV'}{$p}=$ppv3;
					$hash_of_atr_results{$i}{$window_3}{'TP'}{$p}=$tp3;
					$hash_of_atr_results{$i}{$window_3}{'FP'}{$p}=$fp3;
					$hash_of_atr_results{$i}{$window_3}{'NP'}{$p}=$np3;
					
					$hash_of_atr_results{$i}{$window_4}{'TOTAL'}{$p}=$total4;
					$hash_of_atr_results{$i}{$window_4}{'PPV'}{$p}=$ppv4;
					$hash_of_atr_results{$i}{$window_4}{'TP'}{$p}=$tp4;
					$hash_of_atr_results{$i}{$window_4}{'FP'}{$p}=$fp4;
					$hash_of_atr_results{$i}{$window_4}{'NP'}{$p}=$np4;
					
					$hash_of_atr_results{$i}{$window_5}{'TOTAL'}{$p}=$total5;
					$hash_of_atr_results{$i}{$window_5}{'PPV'}{$p}=$ppv5;
					$hash_of_atr_results{$i}{$window_5}{'TP'}{$p}=$tp5;
					$hash_of_atr_results{$i}{$window_5}{'FP'}{$p}=$fp5;
					$hash_of_atr_results{$i}{$window_5}{'NP'}{$p}=$np5;
					
					$hash_of_atr_results{$i}{$window_6}{'TOTAL'}{$p}=$total6;
					$hash_of_atr_results{$i}{$window_6}{'PPV'}{$p}=$ppv6;
					$hash_of_atr_results{$i}{$window_6}{'TP'}{$p}=$tp6;
					$hash_of_atr_results{$i}{$window_6}{'FP'}{$p}=$fp6;
					$hash_of_atr_results{$i}{$window_6}{'NP'}{$p}=$np6;
					
					$hash_of_atr_results{$i}{$window_7}{'TOTAL'}{$p}=$total7;
					$hash_of_atr_results{$i}{$window_7}{'PPV'}{$p}=$ppv7;
					$hash_of_atr_results{$i}{$window_7}{'TP'}{$p}=$tp7;
					$hash_of_atr_results{$i}{$window_7}{'FP'}{$p}=$fp7;
					$hash_of_atr_results{$i}{$window_7}{'NP'}{$p}=$np7;
					
					$hash_of_atr_results{$i}{$window_8}{'TOTAL'}{$p}=$total8;
					$hash_of_atr_results{$i}{$window_8}{'PPV'}{$p}=$ppv8;
					$hash_of_atr_results{$i}{$window_8}{'TP'}{$p}=$tp8;
					$hash_of_atr_results{$i}{$window_8}{'FP'}{$p}=$fp8;
					$hash_of_atr_results{$i}{$window_8}{'NP'}{$p}=$np8;
					
					$hash_of_atr_results{$i}{$window_9}{'TOTAL'}{$p}=$total9;
					$hash_of_atr_results{$i}{$window_9}{'PPV'}{$p}=$ppv9;
					$hash_of_atr_results{$i}{$window_9}{'TP'}{$p}=$tp9;
					$hash_of_atr_results{$i}{$window_9}{'FP'}{$p}=$fp9;
					$hash_of_atr_results{$i}{$window_9}{'NP'}{$p}=$np9;				
					
					$hash_of_atr_results{$i}{$window_10}{'TOTAL'}{$p}=$total10;
					$hash_of_atr_results{$i}{$window_10}{'PPV'}{$p}=$ppv10;
					$hash_of_atr_results{$i}{$window_10}{'TP'}{$p}=$tp10;
					$hash_of_atr_results{$i}{$window_10}{'FP'}{$p}=$fp10;
					$hash_of_atr_results{$i}{$window_10}{'NP'}{$p}=$np10;		
				
					$hash_of_atr_results{$i}{$window_11}{'TOTAL'}{$p}=$total11;
					$hash_of_atr_results{$i}{$window_11}{'PPV'}{$p}=$ppv11;
					$hash_of_atr_results{$i}{$window_11}{'TP'}{$p}=$tp11;
					$hash_of_atr_results{$i}{$window_11}{'FP'}{$p}=$fp11;
					$hash_of_atr_results{$i}{$window_11}{'NP'}{$p}=$np11;								
			}

	} # end of $p


	print "\nSummary of AT richness analysis:\n\n";

#--------- now calculate the SD and SE ------------------------------------------
	#print "Window\tMinimum % difference in AT\tavg_total_arrays_processed\tppv_avg\tppv_sd\tppv_se\ttp_avg\ttp_sd\ttp_se\tfp_avg\tfp_sd\tfp_se\tnp_avg\tnp_sd\tnp_se\n";
	
	foreach my $min_percent_diff(sort{$a <=> $b} keys %hash_of_atr_results)
		{
			print "For min_percent_difference=$min_percent_diff\n";
			print "Window size (nt)\tAvg_total_arrays_processed\tTP\ttp_sd\ttp_se\tFP\tfp_sd\tfp_se\tNP\tnp_sd\tnp_se\tPPV\tppv_sd\tppv_se\n";
			foreach my $window(sort{$a <=> $b} keys %{$hash_of_atr_results{$min_percent_diff}})
				{
					my $avg_total_arrays_processed=0;		
					my $ppv_avg=0;	my $ppv_sd=0;	my $ppv_se=0;		
					my $tp_avg=0;	my $tp_sd=0;	my $tp_se=0;		
					my $fp_avg=0;	my $fp_sd=0;	my $fp_se=0;		
					my $np_avg=0;	my $np_sd=0;	my $np_se=0;
					
					
					foreach my $type(sort{$a cmp $b} keys %{$hash_of_atr_results{$min_percent_diff}{$window}})
						{				
							#print "$type=$type\n";
							my @arr_data;
							
							foreach my $pass(sort{$a <=> $b} keys %{$hash_of_atr_results{$min_percent_diff}{$window}{$type}})
								{
									push(@arr_data,$hash_of_atr_results{$min_percent_diff}{$window}{$type}{$pass});
									#print "@arr_data\n";
									#print "$analysis\t$type\t$pass\t$hash_of_analysis_results{$analysis}{$type}{$pass}\n";
								}
								
							#--- now calculate the mean, SD and SE ------------------
							my $ave;my $std; my $std_err=0;
							
							$ave = &average(\@arr_data);
							
							if($type!~/PPV/){$ave=sprintf("%.0f",$ave);}
								
							if($type!~/TOTAL/)
								{
									$std = &stdev(\@arr_data);	
									
									if($std !=0 and $ave!=0)
										{				
											$std_err=($std/($ave ** 0.5)); $std_err=sprintf("%.2f",$std_err);
										}
								}
							if($type=~/TOTAL/)
								{
									$avg_total_arrays_processed=int($ave);
									#print "Total=$avg_total_arrays_processed\n";
									next;
								}
							elsif($type=~/PPV/){$ppv_avg=$ave;$ppv_sd=$std;$ppv_se=$std_err;}
							elsif($type=~/TP/){$tp_avg=$ave;$tp_sd=$std;$tp_se=$std_err;}
							elsif($type=~/FP/){$fp_avg=$ave;$fp_sd=$std;$fp_se=$std_err;}
							elsif($type=~/NP/){$np_avg=$ave;$np_sd=$std;$np_se=$std_err;}
							#print "$min_percent_diff\t$type\t$ave\t$std\n";
							
							if($type=~/PPV/ and $min_percent_diff==$minimum_percent_at_diff and $window==$at_window_length){$mean_at_richness_ppv=$ppv_avg;}
						}
						
					print "$window\t$avg_total_arrays_processed\t$tp_avg\t$tp_sd\t$tp_se\t$fp_avg\t$fp_sd\t$fp_se\t$np_avg\t$np_sd\t$np_se\t$ppv_avg\t$ppv_sd\t$ppv_se\n";	
				}	
			print "\n";
		}



}# end of at richness

if($check_at_richness==1)
{
	print "\n\nProgram is going to check the distribution of AT in different windows in flanking regions of test set:\n";
	my %hash_of_atr_results;

	for(my $p=1;$p<=$q;$p++)
		{

	
		#print "The following records are used in training:$p:\n";
		my $rec_count=1;
		#-------now pick 100 keys from the sample --------------------	

		my %hash_of_current_sample;
		my %hash_of_ref_repeat_direction_current_sample;
		for(1..$r*90)
			{
				my $found=0;
				while($found!=1)
				{
					my $random_no=int(rand($total_records-1));
				
					$random_no=$random_no+1;
					
					my $current_record=$hash_of_all_records{$random_no};
					
					#----- skip if not part of training set -------------------------------------
					if(not defined $hash_of_records_used_for_training{$current_record}){next;}
		
					#----------------------------------------------------------------------------
													
					if(not defined $hash_of_current_sample{$current_record})
						{
													
							my $key=$current_record; 
							
							#print "$rec_count:$current_record\n";
							#$rec_count++;
									
							#------ skip the ones with motif ---------------------
							#if($hash_of_repeat_keys_and_motif_found{$key}){ next;}
							#------------------------------------------------------
							
							
							$key=~s/-/\t/g;
							
							#print "$key\n"; 
							
							#----- check if the prediction is not NA  for window 135 -------------------------------------------------------------------
							my $rec_line=`grep '$key' at_richness_distribution.txt >&1`;
							#print "grep '$key' at_richness_distribution.txt\n$rec_line\n";
							my @arr_tmp1=split('\t',$rec_line);
							
							if($arr_tmp1[9]=~/NA/){next;}
							
							$hash_of_ref_repeat_direction_current_sample{$current_record}=$hash_of_ref_repeat_direction{$current_record};
							$hash_of_current_sample{$current_record}=1;
							#print "$current_record\n";
							$found=1;
						}
				}
				my $current_no_of_entries=keys %hash_of_current_sample;
				if($current_no_of_entries>=$r)
					{
						#print "\$current_no_of_entries=$current_no_of_entries;";
						last;
					}
			}
		my $t_sample_size=keys %hash_of_current_sample;
		#print "\$t_sample_size=$t_sample_size\n";
		#exit;
		######### check AT_richness#############

			
		for (my $i=1;$i<15;$i++)
			{
				
					my ($total1,$window_1,$ppv1,$tp1,$fp1,$np1,$total2,$window_2,$ppv2,$tp2,$fp2,$np2,$total3,$window_3,$ppv3,$tp3,$fp3,$np3,$total4,$window_4,$ppv4,$tp4,$fp4,$np4,$total5,$window_5,$ppv5,$tp5,$fp5,$np5,$total6,$window_6,$ppv6,$tp6,$fp6,$np6,$total7,$window_7,$ppv7,$tp7,$fp7,$np7,$total8,$window_8,$ppv8,$tp8,$fp8,$np8,$total9,$window_9,$ppv9,$tp9,$fp9,$np9,$total10,$window_10,$ppv10,$tp10,$fp10,$np10,$total11,$window_11,$ppv11,$tp11,$fp11,$np11);

					($total1,$window_1,$ppv1,$tp1,$fp1,$np1,$total2,$window_2,$ppv2,$tp2,$fp2,$np2,$total3,$window_3,$ppv3,$tp3,$fp3,$np3,$total4,$window_4,$ppv4,$tp4,$fp4,$np4,$total5,$window_5,$ppv5,$tp5,$fp5,$np5,$total6,$window_6,$ppv6,$tp6,$fp6,$np6,$total7,$window_7,$ppv7,$tp7,$fp7,$np7,$total8,$window_8,$ppv8,$tp8,$fp8,$np8,$total9,$window_9,$ppv9,$tp9,$fp9,$np9,$total10,$window_10,$ppv10,$tp10,$fp10,$np10,$total11,$window_11,$ppv11,$tp11,$fp11,$np11)=&at_richness_analysis($i,\%hash_of_current_sample,\%hash_of_at_richness_records);
						

					$hash_of_atr_results{$i}{$window_1}{'TOTAL'}{$p}=$total1;
					$hash_of_atr_results{$i}{$window_1}{'PPV'}{$p}=$ppv1;
					$hash_of_atr_results{$i}{$window_1}{'TP'}{$p}=$tp1;
					$hash_of_atr_results{$i}{$window_1}{'FP'}{$p}=$fp1;
					$hash_of_atr_results{$i}{$window_1}{'NP'}{$p}=$np1;
					
					$hash_of_atr_results{$i}{$window_2}{'TOTAL'}{$p}=$total2;
					$hash_of_atr_results{$i}{$window_2}{'PPV'}{$p}=$ppv2;
					$hash_of_atr_results{$i}{$window_2}{'TP'}{$p}=$tp2;
					$hash_of_atr_results{$i}{$window_2}{'FP'}{$p}=$fp2;
					$hash_of_atr_results{$i}{$window_2}{'NP'}{$p}=$np2;
					
					$hash_of_atr_results{$i}{$window_3}{'TOTAL'}{$p}=$total3;
					$hash_of_atr_results{$i}{$window_3}{'PPV'}{$p}=$ppv3;
					$hash_of_atr_results{$i}{$window_3}{'TP'}{$p}=$tp3;
					$hash_of_atr_results{$i}{$window_3}{'FP'}{$p}=$fp3;
					$hash_of_atr_results{$i}{$window_3}{'NP'}{$p}=$np3;
					
					$hash_of_atr_results{$i}{$window_4}{'TOTAL'}{$p}=$total4;
					$hash_of_atr_results{$i}{$window_4}{'PPV'}{$p}=$ppv4;
					$hash_of_atr_results{$i}{$window_4}{'TP'}{$p}=$tp4;
					$hash_of_atr_results{$i}{$window_4}{'FP'}{$p}=$fp4;
					$hash_of_atr_results{$i}{$window_4}{'NP'}{$p}=$np4;
					
					$hash_of_atr_results{$i}{$window_5}{'TOTAL'}{$p}=$total5;
					$hash_of_atr_results{$i}{$window_5}{'PPV'}{$p}=$ppv5;
					$hash_of_atr_results{$i}{$window_5}{'TP'}{$p}=$tp5;
					$hash_of_atr_results{$i}{$window_5}{'FP'}{$p}=$fp5;
					$hash_of_atr_results{$i}{$window_5}{'NP'}{$p}=$np5;
					
					$hash_of_atr_results{$i}{$window_6}{'TOTAL'}{$p}=$total6;
					$hash_of_atr_results{$i}{$window_6}{'PPV'}{$p}=$ppv6;
					$hash_of_atr_results{$i}{$window_6}{'TP'}{$p}=$tp6;
					$hash_of_atr_results{$i}{$window_6}{'FP'}{$p}=$fp6;
					$hash_of_atr_results{$i}{$window_6}{'NP'}{$p}=$np6;
					
					$hash_of_atr_results{$i}{$window_7}{'TOTAL'}{$p}=$total7;
					$hash_of_atr_results{$i}{$window_7}{'PPV'}{$p}=$ppv7;
					$hash_of_atr_results{$i}{$window_7}{'TP'}{$p}=$tp7;
					$hash_of_atr_results{$i}{$window_7}{'FP'}{$p}=$fp7;
					$hash_of_atr_results{$i}{$window_7}{'NP'}{$p}=$np7;
					
					$hash_of_atr_results{$i}{$window_8}{'TOTAL'}{$p}=$total8;
					$hash_of_atr_results{$i}{$window_8}{'PPV'}{$p}=$ppv8;
					$hash_of_atr_results{$i}{$window_8}{'TP'}{$p}=$tp8;
					$hash_of_atr_results{$i}{$window_8}{'FP'}{$p}=$fp8;
					$hash_of_atr_results{$i}{$window_8}{'NP'}{$p}=$np8;
					
					$hash_of_atr_results{$i}{$window_9}{'TOTAL'}{$p}=$total9;
					$hash_of_atr_results{$i}{$window_9}{'PPV'}{$p}=$ppv9;
					$hash_of_atr_results{$i}{$window_9}{'TP'}{$p}=$tp9;
					$hash_of_atr_results{$i}{$window_9}{'FP'}{$p}=$fp9;
					$hash_of_atr_results{$i}{$window_9}{'NP'}{$p}=$np9;				
					
					$hash_of_atr_results{$i}{$window_10}{'TOTAL'}{$p}=$total10;
					$hash_of_atr_results{$i}{$window_10}{'PPV'}{$p}=$ppv10;
					$hash_of_atr_results{$i}{$window_10}{'TP'}{$p}=$tp10;
					$hash_of_atr_results{$i}{$window_10}{'FP'}{$p}=$fp10;
					$hash_of_atr_results{$i}{$window_10}{'NP'}{$p}=$np10;		
				
					$hash_of_atr_results{$i}{$window_11}{'TOTAL'}{$p}=$total11;
					$hash_of_atr_results{$i}{$window_11}{'PPV'}{$p}=$ppv11;
					$hash_of_atr_results{$i}{$window_11}{'TP'}{$p}=$tp11;
					$hash_of_atr_results{$i}{$window_11}{'FP'}{$p}=$fp11;
					$hash_of_atr_results{$i}{$window_11}{'NP'}{$p}=$np11;								
			}

	} # end of $p


	print "\nSummary of AT richness analysis:\n\n";

#--------- now calculate the SD and SE ------------------------------------------
	#print "Window\tMinimum % difference in AT\tavg_total_arrays_processed\tppv_avg\tppv_sd\tppv_se\ttp_avg\ttp_sd\ttp_se\tfp_avg\tfp_sd\tfp_se\tnp_avg\tnp_sd\tnp_se\n";
	
	foreach my $min_percent_diff(sort{$a <=> $b} keys %hash_of_atr_results)
		{
			print "For min_percent_difference=$min_percent_diff\n";
			print "Window size (nt)\tAvg_total_arrays_processed\tTP\ttp_sd\ttp_se\tFP\tfp_sd\tfp_se\tNP\tnp_sd\tnp_se\tPPV\tppv_sd\tppv_se\n";
			foreach my $window(sort{$a <=> $b} keys %{$hash_of_atr_results{$min_percent_diff}})
				{
					my $avg_total_arrays_processed=0;		
					my $ppv_avg=0;	my $ppv_sd=0;	my $ppv_se=0;		
					my $tp_avg=0;	my $tp_sd=0;	my $tp_se=0;		
					my $fp_avg=0;	my $fp_sd=0;	my $fp_se=0;		
					my $np_avg=0;	my $np_sd=0;	my $np_se=0;
					
					
					foreach my $type(sort{$a cmp $b} keys %{$hash_of_atr_results{$min_percent_diff}{$window}})
						{				
							#print "$type=$type\n";
							my @arr_data;
							
							foreach my $pass(sort{$a <=> $b} keys %{$hash_of_atr_results{$min_percent_diff}{$window}{$type}})
								{
									push(@arr_data,$hash_of_atr_results{$min_percent_diff}{$window}{$type}{$pass});
									#print "@arr_data\n";
									#print "$analysis\t$type\t$pass\t$hash_of_analysis_results{$analysis}{$type}{$pass}\n";
								}
								
							#--- now calculate the mean, SD and SE ------------------
							my $ave;my $std; my $std_err=0;
							
							$ave = &average(\@arr_data);
							
							if($type!~/PPV/){$ave=sprintf("%.0f",$ave);}
								
							if($type!~/TOTAL/)
								{
									$std = &stdev(\@arr_data);	
									
									if($std !=0 and $ave!=0)
										{				
											$std_err=($std/($ave ** 0.5)); $std_err=sprintf("%.2f",$std_err);
										}
								}
							if($type=~/TOTAL/)
								{
									$avg_total_arrays_processed=int($ave);
									#print "Total=$avg_total_arrays_processed\n";
									next;
								}
							elsif($type=~/PPV/){$ppv_avg=$ave;$ppv_sd=$std;$ppv_se=$std_err;}
							elsif($type=~/TP/){$tp_avg=$ave;$tp_sd=$std;$tp_se=$std_err;}
							elsif($type=~/FP/){$fp_avg=$ave;$fp_sd=$std;$fp_se=$std_err;}
							elsif($type=~/NP/){$np_avg=$ave;$np_sd=$std;$np_se=$std_err;}
							#print "$min_percent_diff\t$type\t$ave\t$std\n";
							
							if($type=~/PPV/ and $min_percent_diff==$minimum_percent_at_diff and $window==$at_window_length){$mean_at_richness_ppv=$ppv_avg;}
						}
						
					print "$window\t$avg_total_arrays_processed\t$tp_avg\t$tp_sd\t$tp_se\t$fp_avg\t$fp_sd\t$fp_se\t$np_avg\t$np_sd\t$np_se\t$ppv_avg\t$ppv_sd\t$ppv_se\n";	
				}	
			print "\n";
		}



}# end of at richness





#----- then check degeneracy  ---------------------------------------------------------------------
my %hash_of_degeneracy_records;
if($check_degeneracy==1)
	{
		print "\nProgram is going to process CRISPR array degeneracy records from the training sets:\n";
		my %hash_of_degen_res;
		
		for(my $p=1;$p<=$q;$p++)
			{
				#print "The following records are used in training:$p:\n";
				my $rec_count=1;
				#-------now pick 100 keys from the sample --------------------	

				my %hash_of_current_sample;
				my %hash_of_ref_repeat_direction_current_sample;
				for(1..$r*90)
					{
						my $found=0;
						while($found!=1)
						{
							my $random_no=int(rand($total_records-1));
						
							$random_no=$random_no+1;
							
							my $current_record=$hash_of_all_records{$random_no};
							
							#----- skip if not part of training set -------------------------------------
							if(not defined $hash_of_records_used_for_training{$current_record}){next;}

							#----------------------------------------------------------------------------
							
							
							if(not defined $hash_of_current_sample{$current_record})
								{
									my $key=$current_record; 
									
									#print "$rec_count:$current_record\n";
									#$rec_count++;
									
									
									#------ skip the ones with motif ---------------------
									#if($hash_of_repeat_keys_and_motif_found{$key}){next;}
									#------------------------------------------------------
							
									$key=~s/-/\t/g;
									#print "$key\n"; 
									#----- check if the prediction is not NA -------------------------------------------------------------------
									my $rec_line=`grep '$key' array_degeneracy_predictions.txt >&1`;
									my @arr_tmp1=split('\t',$rec_line);
									if($arr_tmp1[4]=~/NA \[0\-0\]/)
										{
											#print "$arr_tmp1[4]\n";
											next;
										}
									
									$hash_of_ref_repeat_direction_current_sample{$current_record}=$hash_of_ref_repeat_direction{$current_record};
									$hash_of_current_sample{$current_record}=1;
									#print "$current_record\n";
									$found=1;
								}
						}
						my $current_no_of_entries=keys %hash_of_current_sample;
						if($current_no_of_entries>=$r)
							{
								#print "\$current_no_of_entries=$current_no_of_entries;";
								last;
							}
					}

				###### check degeneracy ###############


		
				for(my $k=0;$k<=2;$k++)
					{
						my $degeneracy_total=0;
						my $degeneracy_tp=0;
						my $degeneracy_fp=0;
						my $degeneracy_np=0;
						
						($degeneracy_total,$degeneracy_ppv,$degeneracy_tp,$degeneracy_fp,$degeneracy_np)=&array_degeneracy_analysis($k,$permitted_no_of_bases,\%hash_of_current_sample,\%hash_of_degeneracy_records);
						
						$hash_of_degen_res{$k}{'TOTAL'}{$p}=$degeneracy_total;
						$hash_of_degen_res{$k}{'PPV'}{$p}=$degeneracy_ppv;
						$hash_of_degen_res{$k}{'TP'}{$p}=$degeneracy_tp;
						$hash_of_degen_res{$k}{'FP'}{$p}=$degeneracy_fp;
						$hash_of_degen_res{$k}{'NP'}{$p}=$degeneracy_np;
						#print "\t$degeneracy_total,$degeneracy_ppv,$degeneracy_tp,$degeneracy_fp,$degeneracy_np\n";
					}
			

			}

	#--------- now calculate the SD and SE ------------------------------------------
	print "Summary of array degeneracy analysis:\n";
	print "Allowed_mutation\tAvg_total_arrays_processed\tTP\ttp_sd\ttp_se\tFP\tfp_sd\tfp_se\tNP\tnp_sd\tnp_se\tPPV\tppv_sd\tppv_se\n";
	foreach my $allowed_mutation(sort{$a <=> $b} keys %hash_of_degen_res)
		{
			my $avg_total_arrays_processed=0;		
			my $ppv_avg=0;	my $ppv_sd=0;	my $ppv_se=0;		
			my $tp_avg=0;	my $tp_sd=0;	my $tp_se=0;		
			my $fp_avg=0;	my $fp_sd=0;	my $fp_se=0;		
			my $np_avg=0;	my $np_sd=0;	my $np_se=0;
			
			
			foreach my $type(sort{$a cmp $b} keys %{$hash_of_degen_res{$allowed_mutation}})
				{				
					#print "$type=$type\n";
					my @arr_data;
					
					foreach my $pass(sort{$a <=> $b} keys %{$hash_of_degen_res{$allowed_mutation}{$type}})
						{
							push(@arr_data,$hash_of_degen_res{$allowed_mutation}{$type}{$pass});
							#print "@arr_data\n";
							#print "$analysis\t$type\t$pass\t$hash_of_analysis_results{$analysis}{$type}{$pass}\n";
						}
						
					#--- now calculate the mean, SD and SE ------------------
					my $ave;my $std; my $std_err=0;
					
					$ave = &average(\@arr_data);
					
					if($type!~/PPV/){$ave=sprintf("%.0f",$ave);}
						
					if($type!~/TOTAL/)
						{
							$std = &stdev(\@arr_data);	
							
							if($std !=0)
								{				
									$std_err=($std/($ave ** 0.5)); $std_err=sprintf("%.2f",$std_err);
								}
						}
					if($type=~/TOTAL/)
						{
							$avg_total_arrays_processed=int($ave);
							#print "Total=$avg_total_arrays_processed\n";
							next;
						}
					elsif($type=~/PPV/){$ppv_avg=$ave;$ppv_sd=$std;$ppv_se=$std_err;}
					elsif($type=~/TP/){$tp_avg=$ave;$tp_sd=$std;$tp_se=$std_err;}
					elsif($type=~/FP/){$fp_avg=$ave;$fp_sd=$std;$fp_se=$std_err;}
					elsif($type=~/NP/){$np_avg=$ave;$np_sd=$std;$np_se=$std_err;}
					#print "$allowed_mutation\t$type\t$ave\t$std\n";
					
					if($type=~/PPV/ and $allowed_mutation==$permitted_no_of_bases){$mean_degeneracy_ppv=$ppv_avg;}
				}
				
			print "$allowed_mutation\t$avg_total_arrays_processed\t$tp_avg\t$tp_sd\t$tp_se\t$fp_avg\t$fp_sd\t$fp_se\t$np_avg\t$np_sd\t$np_se\t$ppv_avg\t$ppv_sd\t$ppv_se\n";	
		}	


}

if($check_degeneracy==1)
	{
		print "\nProgram is going to process CRISPR array degeneracy records from the test sets:\n";
		my %hash_of_degen_res;
		
		for(my $p=1;$p<=$q;$p++)
			{
				#print "The following records are used in training:$p:\n";
				my $rec_count=1;
				#-------now pick 100 keys from the sample --------------------	

				my %hash_of_current_sample;
				my %hash_of_ref_repeat_direction_current_sample;
				for(1..$r*90)
					{
						my $found=0;
						while($found!=1)
						{
							my $random_no=int(rand($total_records-1));
						
							$random_no=$random_no+1;
							
							my $current_record=$hash_of_all_records{$random_no};
							
							#----- skip if not part of training set -------------------------------------
							if(not defined $hash_of_records_used_for_training{$current_record}){next;}

							#----------------------------------------------------------------------------
							
							
							if(not defined $hash_of_current_sample{$current_record})
								{
									my $key=$current_record; 
									
									#print "$rec_count:$current_record\n";
									#$rec_count++;
									
									
									#------ skip the ones with motif ---------------------
									#if($hash_of_repeat_keys_and_motif_found{$key}){next;}
									#------------------------------------------------------
							
									$key=~s/-/\t/g;
									#print "$key\n"; 
									#----- check if the prediction is not NA -------------------------------------------------------------------
									my $rec_line=`grep '$key' array_degeneracy_predictions.txt >&1`;
									my @arr_tmp1=split('\t',$rec_line);
									if($arr_tmp1[4]=~/NA \[0\-0\]/)
										{
											#print "$arr_tmp1[4]\n";
											next;
										}
									
									$hash_of_ref_repeat_direction_current_sample{$current_record}=$hash_of_ref_repeat_direction{$current_record};
									$hash_of_current_sample{$current_record}=1;
									#print "$current_record\n";
									$found=1;
								}
						}
						my $current_no_of_entries=keys %hash_of_current_sample;
						if($current_no_of_entries>=$r)
							{
								#print "\$current_no_of_entries=$current_no_of_entries;";
								last;
							}
					}

				###### check degeneracy ###############


		
				for(my $k=0;$k<=2;$k++)
					{
						my $degeneracy_total=0;
						my $degeneracy_tp=0;
						my $degeneracy_fp=0;
						my $degeneracy_np=0;
						
						($degeneracy_total,$degeneracy_ppv,$degeneracy_tp,$degeneracy_fp,$degeneracy_np)=&array_degeneracy_analysis($k,$permitted_no_of_bases,\%hash_of_current_sample,\%hash_of_degeneracy_records);
						
						$hash_of_degen_res{$k}{'TOTAL'}{$p}=$degeneracy_total;
						$hash_of_degen_res{$k}{'PPV'}{$p}=$degeneracy_ppv;
						$hash_of_degen_res{$k}{'TP'}{$p}=$degeneracy_tp;
						$hash_of_degen_res{$k}{'FP'}{$p}=$degeneracy_fp;
						$hash_of_degen_res{$k}{'NP'}{$p}=$degeneracy_np;
						#print "\t$degeneracy_total,$degeneracy_ppv,$degeneracy_tp,$degeneracy_fp,$degeneracy_np\n";
					}
			

			}

	#--------- now calculate the SD and SE ------------------------------------------
	print "Summary of array degeneracy analysis:\n";
	print "Allowed_mutation\tAvg_total_arrays_processed\tTP\ttp_sd\ttp_se\tFP\tfp_sd\tfp_se\tNP\tnp_sd\tnp_se\tPPV\tppv_sd\tppv_se\n";
	foreach my $allowed_mutation(sort{$a <=> $b} keys %hash_of_degen_res)
		{
			my $avg_total_arrays_processed=0;		
			my $ppv_avg=0;	my $ppv_sd=0;	my $ppv_se=0;		
			my $tp_avg=0;	my $tp_sd=0;	my $tp_se=0;		
			my $fp_avg=0;	my $fp_sd=0;	my $fp_se=0;		
			my $np_avg=0;	my $np_sd=0;	my $np_se=0;
			
			
			foreach my $type(sort{$a cmp $b} keys %{$hash_of_degen_res{$allowed_mutation}})
				{				
					#print "$type=$type\n";
					my @arr_data;
					
					foreach my $pass(sort{$a <=> $b} keys %{$hash_of_degen_res{$allowed_mutation}{$type}})
						{
							push(@arr_data,$hash_of_degen_res{$allowed_mutation}{$type}{$pass});
							#print "@arr_data\n";
							#print "$analysis\t$type\t$pass\t$hash_of_analysis_results{$analysis}{$type}{$pass}\n";
						}
						
					#--- now calculate the mean, SD and SE ------------------
					my $ave;my $std; my $std_err=0;
					
					$ave = &average(\@arr_data);
					
					if($type!~/PPV/){$ave=sprintf("%.0f",$ave);}
						
					if($type!~/TOTAL/)
						{
							$std = &stdev(\@arr_data);	
							
							if($std !=0)
								{				
									$std_err=($std/($ave ** 0.5)); $std_err=sprintf("%.2f",$std_err);
								}
						}
					if($type=~/TOTAL/)
						{
							$avg_total_arrays_processed=int($ave);
							#print "Total=$avg_total_arrays_processed\n";
							next;
						}
					elsif($type=~/PPV/){$ppv_avg=$ave;$ppv_sd=$std;$ppv_se=$std_err;}
					elsif($type=~/TP/){$tp_avg=$ave;$tp_sd=$std;$tp_se=$std_err;}
					elsif($type=~/FP/){$fp_avg=$ave;$fp_sd=$std;$fp_se=$std_err;}
					elsif($type=~/NP/){$np_avg=$ave;$np_sd=$std;$np_se=$std_err;}
					#print "$allowed_mutation\t$type\t$ave\t$std\n";
					
					if($type=~/PPV/ and $allowed_mutation==$permitted_no_of_bases){$mean_degeneracy_ppv=$ppv_avg;}
				}
				
			print "$allowed_mutation\t$avg_total_arrays_processed\t$tp_avg\t$tp_sd\t$tp_se\t$fp_avg\t$fp_sd\t$fp_se\t$np_avg\t$np_sd\t$np_se\t$ppv_avg\t$ppv_sd\t$ppv_se\n";	
		}	


}


#----- next check MFE ---------------------------------------------------------------------
my %hash_of_mfe_records;
if($check_mfe==1)
{
	print "\n\nProgram is going to process the records of RNA secondary structure analysis from the training set:\n";
	my %hash_of_mfe_results;
	
	for(my $p=1;$p<=$q;$p++)
		{
		#print "\$p=$p\n";
		#print "The following records are used in training:$p:\n";
		my $rec_count=1;
		#-------now pick 100 keys from the sample --------------------	

		my %hash_of_current_sample;
		my %hash_of_ref_repeat_direction_current_sample;
		my %hash_of_unique_repeats;
		
		for(1..$r*40)
			{
				my $found=0;
				
				my $no_of_tries=0;
				while($found!=1)
				{
					my $random_no=int(rand($total_records-1));
				
					$random_no=$random_no+1;
					
					my $current_record=$hash_of_all_records{$random_no};
					
					#----- skip if not part of training set -------------------------------------
					if(not defined $hash_of_records_used_for_training{$current_record}){next;}
					#----------------------------------------------------------------------------
					
					if(not defined $hash_of_current_sample{$current_record})
						{
							my $key=$current_record; 
							
							#print "$rec_count:$current_record\n";
							#$rec_count++;
									
							#------ skip the ones with motif ---------------------
							#if($hash_of_repeat_keys_and_motif_found{$key}){next;}
							#------------------------------------------------------
							
							$key=~s/-/\t/g;
							#print "$key\n"; 
							#----- check if the prediction is not NA -------------------------------------------------------------------
							my $rec_line=`grep '$key' all_repeats_MFE_suggestions.txt >&1`; chomp $rec_line;$rec_line=~s/\r//g;
							my @arr_tmp1=split('\t',$rec_line);
							
							my $current_repeat=$arr_tmp1[$#arr_tmp1];
							my $current_repeat_rev_comp=$current_repeat; $current_repeat_rev_comp=reverse $current_repeat_rev_comp;$current_repeat_rev_comp=~tr/ACGT/TGCA/;
							
															#if($arr_tmp1[10]=~/NA/ or $hash_of_unique_repeats{$current_repeat} or $hash_of_unique_repeats{$current_repeat_rev_comp}){next;}
							#if($arr_tmp1[10]=~/NA/ or $hash_of_unique_repeats{$current_repeat}){next;}
							
							$hash_of_ref_repeat_direction_current_sample{$current_record}=$hash_of_ref_repeat_direction{$current_record};
							$hash_of_current_sample{$current_record}=1;
							$hash_of_unique_repeats{$current_repeat}=1;
							#print "$current_record\n";
							$found=1;
						}
						
					$no_of_tries++;
					
					if($no_of_tries>1000){last;}	
				}
				
				my $current_no_of_entries=keys %hash_of_current_sample;
				if($current_no_of_entries>=$r)
					{
						my $no_of_unique_repeats_in_sample=keys %hash_of_unique_repeats;
						#print "$p:\$current_no_of_entries=$current_no_of_entries with \$no_of_unique_repeats_in_sample=$no_of_unique_repeats_in_sample\n";
						last;
					}
			}

		###### check MFE ####################

			for(my $j=0;$j<=11;$j++)
				{
					my $mfe_total=0;
					my $mfe_tp=0;
					my $mfe_fp=0;
					my $mfe_np=0;
					
					($mfe_total,$mfe_ppv,$mfe_tp,$mfe_fp,$mfe_np)=&repeats_structural_analysis($j,$minimum_energy_threshold,\%hash_of_current_sample,\%hash_of_mfe_records);
					
					$hash_of_mfe_results{$j}{'TOTAL'}{$p}=$mfe_total;
					$hash_of_mfe_results{$j}{'PPV'}{$p}=$mfe_ppv;
					$hash_of_mfe_results{$j}{'TP'}{$p}=$mfe_tp;
					$hash_of_mfe_results{$j}{'FP'}{$p}=$mfe_fp;
					$hash_of_mfe_results{$j}{'NP'}{$p}=$mfe_np;
					
					#print "$mfe_total,$mfe_ppv,$mfe_tp,$mfe_fp,$mfe_np\n";
				}


		} # end of 100 passes


	#--------- now calculate the SD and SE ------------------------------------------
	print "Summary of Potential for RNA secondary structure analysis:\n";
	print "Bases_altered\tavg_total_arrays_processed\tTP\ttp_sd\ttp_se\tFP\tfp_sd\tfp_se\tNP\tnp_sd\tnp_se\tPPV\tppv_sd\tppv_se\n";
	foreach my $bases_altered(sort{$a <=> $b} keys %hash_of_mfe_results)
		{
			my $avg_total_arrays_processed=0;		
			my $ppv_avg=0;	my $ppv_sd=0;	my $ppv_se=0;		
			my $tp_avg=0;	my $tp_sd=0;	my $tp_se=0;		
			my $fp_avg=0;	my $fp_sd=0;	my $fp_se=0;		
			my $np_avg=0;	my $np_sd=0;	my $np_se=0;
			
			
			foreach my $type(sort{$a cmp $b} keys %{$hash_of_mfe_results{$bases_altered}})
				{				
					#print "$type=$type\n";
					my @arr_data;
					
					foreach my $pass(sort{$a <=> $b} keys %{$hash_of_mfe_results{$bases_altered}{$type}})
						{
							push(@arr_data,$hash_of_mfe_results{$bases_altered}{$type}{$pass});
							#print "@arr_data\n";
							#print "$analysis\t$type\t$pass\t$hash_of_analysis_results{$analysis}{$type}{$pass}\n";
						}
						
					#--- now calculate the mean, SD and SE ------------------
					my $ave;my $std; my $std_err=0;
					
					$ave = &average(\@arr_data);
					
					if($type!~/PPV/){$ave=sprintf("%.0f",$ave);}
						
					if($type!~/TOTAL/)
						{
							$std = &stdev(\@arr_data);	
							
							if($std !=0 and $ave!=0)
								{				
									$std_err=($std/($ave ** 0.5)); $std_err=sprintf("%.2f",$std_err);
								}
							else{
									$std_err=0;
								}	
						}
					if($type=~/TOTAL/)
						{
							$avg_total_arrays_processed=int($ave);
							#print "Total=$avg_total_arrays_processed\n";
							next;
						}
					elsif($type=~/PPV/){$ppv_avg=$ave;$ppv_sd=$std;$ppv_se=$std_err;}
					elsif($type=~/TP/){$tp_avg=$ave;$tp_sd=$std;$tp_se=$std_err;}
					elsif($type=~/FP/){$fp_avg=$ave;$fp_sd=$std;$fp_se=$std_err;}
					elsif($type=~/NP/){$np_avg=$ave;$np_sd=$std;$np_se=$std_err;}
					#print "$bases_altered\t$type\t$ave\t$std\n";
					
					if($type=~/PPV/ and $bases_altered==$mfe_target_base){$mean_mfe_ppv=$ppv_avg;}
				}
				
			print "$bases_altered\t$avg_total_arrays_processed\t$tp_avg\t$tp_sd\t$tp_se\t$fp_avg\t$fp_sd\t$fp_se\t$np_avg\t$np_sd\t$np_se\t$ppv_avg\t$ppv_sd\t$ppv_se\n";	
		}	




}

if($check_mfe==1)
{
	print "\n\nProgram is going to process the records of RNA secondary structure analysis from the test set:\n";
	my %hash_of_mfe_results;
	
	for(my $p=1;$p<=$q;$p++)
		{
		#print "\$p=$p\n";
		#print "The following records are used in training:$p:\n";
		my $rec_count=1;
		#-------now pick 100 keys from the sample --------------------	

		my %hash_of_current_sample;
		my %hash_of_ref_repeat_direction_current_sample;
		my %hash_of_unique_repeats;
		
		for(1..$r*40)
			{
				my $found=0;
				
				my $no_of_tries=0;
				while($found!=1)
				{
					my $random_no=int(rand($total_records-1));
				
					$random_no=$random_no+1;
					
					my $current_record=$hash_of_all_records{$random_no};
					
					#----- skip if not part of training set -------------------------------------
					if(not defined $hash_of_records_used_for_training{$current_record}){next;}
					#----------------------------------------------------------------------------
					
					if(not defined $hash_of_current_sample{$current_record})
						{
							my $key=$current_record; 
							
							#print "$rec_count:$current_record\n";
							#$rec_count++;
									
							#------ skip the ones with motif ---------------------
							#if($hash_of_repeat_keys_and_motif_found{$key}){next;}
							#------------------------------------------------------
							
							$key=~s/-/\t/g;
							#print "$key\n"; 
							#----- check if the prediction is not NA -------------------------------------------------------------------
							my $rec_line=`grep '$key' all_repeats_MFE_suggestions.txt >&1`; chomp $rec_line;$rec_line=~s/\r//g;
							my @arr_tmp1=split('\t',$rec_line);
							
							my $current_repeat=$arr_tmp1[$#arr_tmp1];
							my $current_repeat_rev_comp=$current_repeat; $current_repeat_rev_comp=reverse $current_repeat_rev_comp;$current_repeat_rev_comp=~tr/ACGT/TGCA/;
							
															#if($arr_tmp1[10]=~/NA/ or $hash_of_unique_repeats{$current_repeat} or $hash_of_unique_repeats{$current_repeat_rev_comp}){next;}
							#if($arr_tmp1[10]=~/NA/ or $hash_of_unique_repeats{$current_repeat}){next;}
							
							$hash_of_ref_repeat_direction_current_sample{$current_record}=$hash_of_ref_repeat_direction{$current_record};
							$hash_of_current_sample{$current_record}=1;
							$hash_of_unique_repeats{$current_repeat}=1;
							#print "$current_record\n";
							$found=1;
						}
						
					$no_of_tries++;
					
					if($no_of_tries>1000){last;}	
				}
				
				my $current_no_of_entries=keys %hash_of_current_sample;
				if($current_no_of_entries>=$r)
					{
						my $no_of_unique_repeats_in_sample=keys %hash_of_unique_repeats;
						#print "$p:\$current_no_of_entries=$current_no_of_entries with \$no_of_unique_repeats_in_sample=$no_of_unique_repeats_in_sample\n";
						last;
					}
			}

		###### check MFE ####################

			for(my $j=0;$j<=11;$j++)
				{
					my $mfe_total=0;
					my $mfe_tp=0;
					my $mfe_fp=0;
					my $mfe_np=0;
					
					($mfe_total,$mfe_ppv,$mfe_tp,$mfe_fp,$mfe_np)=&repeats_structural_analysis($j,$minimum_energy_threshold,\%hash_of_current_sample,\%hash_of_mfe_records);
					
					$hash_of_mfe_results{$j}{'TOTAL'}{$p}=$mfe_total;
					$hash_of_mfe_results{$j}{'PPV'}{$p}=$mfe_ppv;
					$hash_of_mfe_results{$j}{'TP'}{$p}=$mfe_tp;
					$hash_of_mfe_results{$j}{'FP'}{$p}=$mfe_fp;
					$hash_of_mfe_results{$j}{'NP'}{$p}=$mfe_np;
					
					#print "$mfe_total,$mfe_ppv,$mfe_tp,$mfe_fp,$mfe_np\n";
				}


		} # end of 100 passes


	#--------- now calculate the SD and SE ------------------------------------------
	print "Summary of Potential for RNA secondary structure analysis:\n";
	print "Bases_altered\tavg_total_arrays_processed\tTP\ttp_sd\ttp_se\tFP\tfp_sd\tfp_se\tNP\tnp_sd\tnp_se\tPPV\tppv_sd\tppv_se\n";
	foreach my $bases_altered(sort{$a <=> $b} keys %hash_of_mfe_results)
		{
			my $avg_total_arrays_processed=0;		
			my $ppv_avg=0;	my $ppv_sd=0;	my $ppv_se=0;		
			my $tp_avg=0;	my $tp_sd=0;	my $tp_se=0;		
			my $fp_avg=0;	my $fp_sd=0;	my $fp_se=0;		
			my $np_avg=0;	my $np_sd=0;	my $np_se=0;
			
			
			foreach my $type(sort{$a cmp $b} keys %{$hash_of_mfe_results{$bases_altered}})
				{				
					#print "$type=$type\n";
					my @arr_data;
					
					foreach my $pass(sort{$a <=> $b} keys %{$hash_of_mfe_results{$bases_altered}{$type}})
						{
							push(@arr_data,$hash_of_mfe_results{$bases_altered}{$type}{$pass});
							#print "@arr_data\n";
							#print "$analysis\t$type\t$pass\t$hash_of_analysis_results{$analysis}{$type}{$pass}\n";
						}
						
					#--- now calculate the mean, SD and SE ------------------
					my $ave;my $std; my $std_err=0;
					
					$ave = &average(\@arr_data);
					
					if($type!~/PPV/){$ave=sprintf("%.0f",$ave);}
						
					if($type!~/TOTAL/)
						{
							$std = &stdev(\@arr_data);	
							
							if($std !=0 and $ave!=0)
								{				
									$std_err=($std/($ave ** 0.5)); $std_err=sprintf("%.2f",$std_err);
								}
							else{
									$std_err=0;
								}	
						}
					if($type=~/TOTAL/)
						{
							$avg_total_arrays_processed=int($ave);
							#print "Total=$avg_total_arrays_processed\n";
							next;
						}
					elsif($type=~/PPV/){$ppv_avg=$ave;$ppv_sd=$std;$ppv_se=$std_err;}
					elsif($type=~/TP/){$tp_avg=$ave;$tp_sd=$std;$tp_se=$std_err;}
					elsif($type=~/FP/){$fp_avg=$ave;$fp_sd=$std;$fp_se=$std_err;}
					elsif($type=~/NP/){$np_avg=$ave;$np_sd=$std;$np_se=$std_err;}
					#print "$bases_altered\t$type\t$ave\t$std\n";
					
					if($type=~/PPV/ and $bases_altered==$mfe_target_base){$mean_mfe_ppv=$ppv_avg;}
				}
				
			print "$bases_altered\t$avg_total_arrays_processed\t$tp_avg\t$tp_sd\t$tp_se\t$fp_avg\t$fp_sd\t$fp_se\t$np_avg\t$np_sd\t$np_se\t$ppv_avg\t$ppv_sd\t$ppv_se\n";	
		}	




}



#----- next check longer leader ---------------------------------------------------------------------
my %hash_of_longer_leader_records;
if($check_longer_leader==1)
{
	print "\n\nProgram is now processing the records of longer leader analysis from training set:\n";
	my %hash_of_ll_results;
	
	for(my $p=1;$p<=$q;$p++)
		{
			
		#print "The following records are used in training:$p:\n";
		my $rec_count=1;
		#-------now pick 100 keys from the sample --------------------	

		my %hash_of_current_sample;
		my %hash_of_ref_repeat_direction_current_sample;
		
		
		for(1..$r*90)
			{
				my $found=0;
				while($found!=1)
					{
						my $random_no=int(rand($total_records-1));
					
						$random_no=$random_no+1;
						
						my $current_record=$hash_of_all_records{$random_no};
								
						#----- skip if not part of training set -------------------------------------
						if(not defined $hash_of_records_used_for_training{$current_record}){next;}
						#----------------------------------------------------------------------------
								
						if(not defined $hash_of_current_sample{$current_record})
							{
								
								#print "$rec_count:$current_record\n";
								$rec_count++;
								
								my $key=$current_record; 
								
								#print "$rec_count:$current_record\n";
								#$rec_count++;		
								
								#------ skip the ones with motif ---------------------
								#if($hash_of_repeat_keys_and_motif_found{$key}){ next;}
								#------------------------------------------------------
							
								$key=~s/-/\t/g;
								#print "$key\n"; 
								#----- check if the prediction is not NA -------------------------------------------------------------------
								my $rec_line=`grep '$key' longer_leader_prediction.txt >&1`;
								my @arr_tmp1=split('\t',$rec_line);
								#if($arr_tmp1[10]=~/NA/){next;}
								
								$hash_of_ref_repeat_direction_current_sample{$current_record}=$hash_of_ref_repeat_direction{$current_record};
								$hash_of_current_sample{$current_record}=1;
								#print "$current_record\n";
								$found=1;
							}
					}
		
				my $current_no_of_entries=keys %hash_of_current_sample;
				if($current_no_of_entries>=$r)
					{
						#print "\$current_no_of_entries=$current_no_of_entries;";
						last;
					}
			}







	
		#####  check longer leader ############

		my $longer_leader_total=0;
		my $longer_leader_tp=0;
		my $longer_leader_fp=0;
		my $longer_leader_np=0;
		
		#my $lb=10000000000000000; my $ub=0;
		
		my($ll_total1,$ll_ppv1,$ll_tp1,$ll_fp1,$ll_np1,$ll_total2,$ll_ppv2,$ll_tp2,$ll_fp2,$ll_np2,$ll_total3,$ll_ppv3,$ll_tp3,$ll_fp3,$ll_np3,$ll_total4,$ll_ppv4,$ll_tp4,$ll_fp4,$ll_np4,$ll_total5,$ll_ppv5,$ll_tp5,$ll_fp5,$ll_np5,$ll_total6,$ll_ppv6,$ll_tp6,$ll_fp6,$ll_np6)=&analysis_of_longer_leader($longer_leader_total,$longer_leader_ppv,$longer_leader_tp,$longer_leader_fp,$longer_leader_np,\%hash_of_current_sample,\%hash_of_longer_leader_records);

		#print "\n\$lb=$lb\t\$ub=$ub\n";

		$hash_of_ll_results{'1'}{'TOTAL'}{$p}=$ll_total1;
		$hash_of_ll_results{'1'}{'PPV'}{$p}=$ll_ppv1;
		$hash_of_ll_results{'1'}{'TP'}{$p}=$ll_tp1;
		$hash_of_ll_results{'1'}{'FP'}{$p}=$ll_fp1;
		$hash_of_ll_results{'1'}{'NP'}{$p}=$ll_np1;
		

		$hash_of_ll_results{'1.2'}{'TOTAL'}{$p}=$ll_total2;
		$hash_of_ll_results{'1.2'}{'PPV'}{$p}=$ll_ppv2;
		$hash_of_ll_results{'1.2'}{'TP'}{$p}=$ll_tp2;
		$hash_of_ll_results{'1.2'}{'FP'}{$p}=$ll_fp2;
		$hash_of_ll_results{'1.2'}{'NP'}{$p}=$ll_np2;
		
		
		$hash_of_ll_results{'1.4'}{'TOTAL'}{$p}=$ll_total3;
		$hash_of_ll_results{'1.4'}{'PPV'}{$p}=$ll_ppv3;
		$hash_of_ll_results{'1.4'}{'TP'}{$p}=$ll_tp3;
		$hash_of_ll_results{'1.4'}{'FP'}{$p}=$ll_fp3;
		$hash_of_ll_results{'1.4'}{'NP'}{$p}=$ll_np3;
		
		
		$hash_of_ll_results{'1.6'}{'TOTAL'}{$p}=$ll_total4;
		$hash_of_ll_results{'1.6'}{'PPV'}{$p}=$ll_ppv4;
		$hash_of_ll_results{'1.6'}{'TP'}{$p}=$ll_tp4;
		$hash_of_ll_results{'1.6'}{'FP'}{$p}=$ll_fp4;
		$hash_of_ll_results{'1.6'}{'NP'}{$p}=$ll_np4;
		
		
		$hash_of_ll_results{'1.8'}{'TOTAL'}{$p}=$ll_total5;
		$hash_of_ll_results{'1.8'}{'PPV'}{$p}=$ll_ppv5;
		$hash_of_ll_results{'1.8'}{'TP'}{$p}=$ll_tp5;
		$hash_of_ll_results{'1.8'}{'FP'}{$p}=$ll_fp5;
		$hash_of_ll_results{'1.8'}{'NP'}{$p}=$ll_np5;
		
		
		$hash_of_ll_results{'2'}{'TOTAL'}{$p}=$ll_total6;
		$hash_of_ll_results{'2'}{'PPV'}{$p}=$ll_ppv6;
		$hash_of_ll_results{'2'}{'TP'}{$p}=$ll_tp6;
		$hash_of_ll_results{'2'}{'FP'}{$p}=$ll_fp6;
		$hash_of_ll_results{'2'}{'NP'}{$p}=$ll_np6;
	}

	print "\nSummary of longer leader analysis:\n";

#--------- now calculate the SD and SE ------------------------------------------
	print "Minimum length difference\tAvg_total_arrays_processed\tTP\ttp_sd\ttp_se\tFP\tfp_sd\tfp_se\tNP\tnp_sd\tnp_se\tPPV\tppv_sd\tppv_se\n";
	foreach my $min_length_diff(sort{$a <=> $b} keys %hash_of_ll_results)
		{
			my $avg_total_arrays_processed=0;		
			my $ppv_avg=0;	my $ppv_sd=0;	my $ppv_se=0;		
			my $tp_avg=0;	my $tp_sd=0;	my $tp_se=0;		
			my $fp_avg=0;	my $fp_sd=0;	my $fp_se=0;		
			my $np_avg=0;	my $np_sd=0;	my $np_se=0;
			
			
			foreach my $type(sort{$a cmp $b} keys %{$hash_of_ll_results{$min_length_diff}})
				{				
					#print "$type=$type\n";
					my @arr_data;
					
					foreach my $pass(sort{$a <=> $b} keys %{$hash_of_ll_results{$min_length_diff}{$type}})
						{
							push(@arr_data,$hash_of_ll_results{$min_length_diff}{$type}{$pass});
							#print "@arr_data\n";
							#print "$analysis\t$type\t$pass\t$hash_of_analysis_results{$analysis}{$type}{$pass}\n";
						}
						
					#--- now calculate the mean, SD and SE ------------------
					my $ave;my $std; my $std_err=0;
					
					$ave = &average(\@arr_data);
					
					if($type!~/PPV/){$ave=sprintf("%.0f",$ave);}	
					
					if($type!~/TOTAL/)
						{
							$std = &stdev(\@arr_data);	
							
							if($std !=0 and $ave!=0)
								{				
									$std_err=($std/($ave ** 0.5)); $std_err=sprintf("%.2f",$std_err);
								}
						}
					if($type=~/TOTAL/)
						{
							$avg_total_arrays_processed=int($ave);
							#print "Total=$avg_total_arrays_processed\n";
							next;
						}
					elsif($type=~/PPV/){$ppv_avg=$ave;$ppv_sd=$std;$ppv_se=$std_err;}
					elsif($type=~/TP/){$tp_avg=$ave;$tp_sd=$std;$tp_se=$std_err;}
					elsif($type=~/FP/){$fp_avg=$ave;$fp_sd=$std;$fp_se=$std_err;}
					elsif($type=~/NP/){$np_avg=$ave;$np_sd=$std;$np_se=$std_err;}
					#print "$min_length_diff\t$type\t$ave\t$std\n";
					
					if($type=~/PPV/ and $min_length_diff==2){$mean_longer_leader_ppv=$ppv_avg;}
				}
				
			print "$min_length_diff\t$avg_total_arrays_processed\t$tp_avg\t$tp_sd\t$tp_se\t$fp_avg\t$fp_sd\t$fp_se\t$np_avg\t$np_sd\t$np_se\t$ppv_avg\t$ppv_sd\t$ppv_se\n";	
		}	




}

if($check_longer_leader==1)
{
	print "\n\nProgram is now processing the records of longer leader analysis from test set:\n";
	my %hash_of_ll_results;
	
	for(my $p=1;$p<=$q;$p++)
		{
			
		#print "The following records are used in training:$p:\n";
		my $rec_count=1;
		#-------now pick 100 keys from the sample --------------------	

		my %hash_of_current_sample;
		my %hash_of_ref_repeat_direction_current_sample;
		
		
		for(1..$r*90)
			{
				my $found=0;
				while($found!=1)
					{
						my $random_no=int(rand($total_records-1));
					
						$random_no=$random_no+1;
						
						my $current_record=$hash_of_all_records{$random_no};
								
						#----- skip if not part of training set -------------------------------------
						if(not defined $hash_of_records_used_for_training{$current_record}){next;}
						#----------------------------------------------------------------------------
								
						if(not defined $hash_of_current_sample{$current_record})
							{
								
								#print "$rec_count:$current_record\n";
								$rec_count++;
								
								my $key=$current_record; 
								
								#print "$rec_count:$current_record\n";
								#$rec_count++;		
								
								#------ skip the ones with motif ---------------------
								#if($hash_of_repeat_keys_and_motif_found{$key}){ next;}
								#------------------------------------------------------
							
								$key=~s/-/\t/g;
								#print "$key\n"; 
								#----- check if the prediction is not NA -------------------------------------------------------------------
								my $rec_line=`grep '$key' longer_leader_prediction.txt >&1`;
								my @arr_tmp1=split('\t',$rec_line);
								#if($arr_tmp1[10]=~/NA/){next;}
								
								$hash_of_ref_repeat_direction_current_sample{$current_record}=$hash_of_ref_repeat_direction{$current_record};
								$hash_of_current_sample{$current_record}=1;
								#print "$current_record\n";
								$found=1;
							}
					}
		
				my $current_no_of_entries=keys %hash_of_current_sample;
				if($current_no_of_entries>=$r)
					{
						#print "\$current_no_of_entries=$current_no_of_entries;";
						last;
					}
			}







	
		#####  check longer leader ############

		my $longer_leader_total=0;
		my $longer_leader_tp=0;
		my $longer_leader_fp=0;
		my $longer_leader_np=0;
		
		#my $lb=10000000000000000; my $ub=0;
		
		my($ll_total1,$ll_ppv1,$ll_tp1,$ll_fp1,$ll_np1,$ll_total2,$ll_ppv2,$ll_tp2,$ll_fp2,$ll_np2,$ll_total3,$ll_ppv3,$ll_tp3,$ll_fp3,$ll_np3,$ll_total4,$ll_ppv4,$ll_tp4,$ll_fp4,$ll_np4,$ll_total5,$ll_ppv5,$ll_tp5,$ll_fp5,$ll_np5,$ll_total6,$ll_ppv6,$ll_tp6,$ll_fp6,$ll_np6)=&analysis_of_longer_leader($longer_leader_total,$longer_leader_ppv,$longer_leader_tp,$longer_leader_fp,$longer_leader_np,\%hash_of_current_sample,\%hash_of_longer_leader_records);

		#print "\n\$lb=$lb\t\$ub=$ub\n";

		$hash_of_ll_results{'1'}{'TOTAL'}{$p}=$ll_total1;
		$hash_of_ll_results{'1'}{'PPV'}{$p}=$ll_ppv1;
		$hash_of_ll_results{'1'}{'TP'}{$p}=$ll_tp1;
		$hash_of_ll_results{'1'}{'FP'}{$p}=$ll_fp1;
		$hash_of_ll_results{'1'}{'NP'}{$p}=$ll_np1;
		

		$hash_of_ll_results{'1.2'}{'TOTAL'}{$p}=$ll_total2;
		$hash_of_ll_results{'1.2'}{'PPV'}{$p}=$ll_ppv2;
		$hash_of_ll_results{'1.2'}{'TP'}{$p}=$ll_tp2;
		$hash_of_ll_results{'1.2'}{'FP'}{$p}=$ll_fp2;
		$hash_of_ll_results{'1.2'}{'NP'}{$p}=$ll_np2;
		
		
		$hash_of_ll_results{'1.4'}{'TOTAL'}{$p}=$ll_total3;
		$hash_of_ll_results{'1.4'}{'PPV'}{$p}=$ll_ppv3;
		$hash_of_ll_results{'1.4'}{'TP'}{$p}=$ll_tp3;
		$hash_of_ll_results{'1.4'}{'FP'}{$p}=$ll_fp3;
		$hash_of_ll_results{'1.4'}{'NP'}{$p}=$ll_np3;
		
		
		$hash_of_ll_results{'1.6'}{'TOTAL'}{$p}=$ll_total4;
		$hash_of_ll_results{'1.6'}{'PPV'}{$p}=$ll_ppv4;
		$hash_of_ll_results{'1.6'}{'TP'}{$p}=$ll_tp4;
		$hash_of_ll_results{'1.6'}{'FP'}{$p}=$ll_fp4;
		$hash_of_ll_results{'1.6'}{'NP'}{$p}=$ll_np4;
		
		
		$hash_of_ll_results{'1.8'}{'TOTAL'}{$p}=$ll_total5;
		$hash_of_ll_results{'1.8'}{'PPV'}{$p}=$ll_ppv5;
		$hash_of_ll_results{'1.8'}{'TP'}{$p}=$ll_tp5;
		$hash_of_ll_results{'1.8'}{'FP'}{$p}=$ll_fp5;
		$hash_of_ll_results{'1.8'}{'NP'}{$p}=$ll_np5;
		
		
		$hash_of_ll_results{'2'}{'TOTAL'}{$p}=$ll_total6;
		$hash_of_ll_results{'2'}{'PPV'}{$p}=$ll_ppv6;
		$hash_of_ll_results{'2'}{'TP'}{$p}=$ll_tp6;
		$hash_of_ll_results{'2'}{'FP'}{$p}=$ll_fp6;
		$hash_of_ll_results{'2'}{'NP'}{$p}=$ll_np6;
	}

	print "\nSummary of longer leader analysis:\n";

#--------- now calculate the SD and SE ------------------------------------------
	print "Minimum length difference\tAvg_total_arrays_processed\tTP\ttp_sd\ttp_se\tFP\tfp_sd\tfp_se\tNP\tnp_sd\tnp_se\tPPV\tppv_sd\tppv_se\n";
	foreach my $min_length_diff(sort{$a <=> $b} keys %hash_of_ll_results)
		{
			my $avg_total_arrays_processed=0;		
			my $ppv_avg=0;	my $ppv_sd=0;	my $ppv_se=0;		
			my $tp_avg=0;	my $tp_sd=0;	my $tp_se=0;		
			my $fp_avg=0;	my $fp_sd=0;	my $fp_se=0;		
			my $np_avg=0;	my $np_sd=0;	my $np_se=0;
			
			
			foreach my $type(sort{$a cmp $b} keys %{$hash_of_ll_results{$min_length_diff}})
				{				
					#print "$type=$type\n";
					my @arr_data;
					
					foreach my $pass(sort{$a <=> $b} keys %{$hash_of_ll_results{$min_length_diff}{$type}})
						{
							push(@arr_data,$hash_of_ll_results{$min_length_diff}{$type}{$pass});
							#print "@arr_data\n";
							#print "$analysis\t$type\t$pass\t$hash_of_analysis_results{$analysis}{$type}{$pass}\n";
						}
						
					#--- now calculate the mean, SD and SE ------------------
					my $ave;my $std; my $std_err=0;
					
					$ave = &average(\@arr_data);
					
					if($type!~/PPV/){$ave=sprintf("%.0f",$ave);}	
					
					if($type!~/TOTAL/)
						{
							$std = &stdev(\@arr_data);	
							
							if($std !=0 and $ave!=0)
								{				
									$std_err=($std/($ave ** 0.5)); $std_err=sprintf("%.2f",$std_err);
								}
						}
					if($type=~/TOTAL/)
						{
							$avg_total_arrays_processed=int($ave);
							#print "Total=$avg_total_arrays_processed\n";
							next;
						}
					elsif($type=~/PPV/){$ppv_avg=$ave;$ppv_sd=$std;$ppv_se=$std_err;}
					elsif($type=~/TP/){$tp_avg=$ave;$tp_sd=$std;$tp_se=$std_err;}
					elsif($type=~/FP/){$fp_avg=$ave;$fp_sd=$std;$fp_se=$std_err;}
					elsif($type=~/NP/){$np_avg=$ave;$np_sd=$std;$np_se=$std_err;}
					#print "$min_length_diff\t$type\t$ave\t$std\n";
					
					if($type=~/PPV/ and $min_length_diff==2){$mean_longer_leader_ppv=$ppv_avg;}
				}
				
			print "$min_length_diff\t$avg_total_arrays_processed\t$tp_avg\t$tp_sd\t$tp_se\t$fp_avg\t$fp_sd\t$fp_se\t$np_avg\t$np_sd\t$np_se\t$ppv_avg\t$ppv_sd\t$ppv_se\n";	
		}	




}

#exit;
#exit;


print "Summary of all individual analysis from training set:\n";
print qq~
	Motif PPV= $mean_motif_PPV 
	A and T ratio PPV=	$mean_A_T_distribution_in_repeats_PPV
	Structural analysis (MFE) PPV with $mfe_target_base bases neutralized from both end= $mean_mfe_ppv

	AT richness PPV in window $at_window_length with min. $minimum_percent_at_diff\% diff = $mean_at_richness_ppv
	Array degeneracy PPV with $permitted_no_of_bases permitted mutation= $mean_degeneracy_ppv
	Longer leader PPV with minimum double length difference=  $mean_longer_leader_ppv\n\n~;



#------------------------------ print individual PPV, TP and FPs with their SD and SE -----------------------------------------------------------------------


print "\n\nFinal Analysis using randonmy selected $r arrays and repeated $q times:\n";

#----- open each file and get the relative columns ----------------------------------
my %hash_of_longer_leader_records_final;
my %hash_of_at_richness_records_final;
my %hash_of_degeneracy_records_final;
my %hash_of_mfe_records_final;




#---- AT_richness file -------------------------------------
my $file2="at_richness_distribution.txt";

open(RD,"$file2") or print "$!";
my @arr_file2=<RD>;
close(RD);

foreach my $line(@arr_file2)
	{
		chomp $line; $line=~s/\r+//g;$line=~s/>//g;#$line=~s/ \[\S+\]//g;		
		my @arr_tmp=split('\t',$line);			
		my $key="$arr_tmp[0]-$arr_tmp[1]-$arr_tmp[2]";
		my $at_richness_suggestion;
		
		if($at_window_length==15){$at_richness_suggestion=$arr_tmp[4];}
		elsif($at_window_length==30){$at_richness_suggestion=$arr_tmp[5];}
		if($at_window_length==45){$at_richness_suggestion=$arr_tmp[6];}
		elsif($at_window_length==60){$at_richness_suggestion=$arr_tmp[7];}
		elsif($at_window_length==75){$at_richness_suggestion=$arr_tmp[8];}
		elsif($at_window_length==90){$at_richness_suggestion=$arr_tmp[9];}
		elsif($at_window_length==105){$at_richness_suggestion=$arr_tmp[10];}
		elsif($at_window_length==120){$at_richness_suggestion=$arr_tmp[11];}
		elsif($at_window_length==135){$at_richness_suggestion=$arr_tmp[12];}
		elsif($at_window_length==150){$at_richness_suggestion=$arr_tmp[13];}
		elsif($at_window_length==165){$at_richness_suggestion=$arr_tmp[14];}

		
		$at_richness_suggestion=~s/ \[(\S+)\]\%AT//g;	my($l_seq,$t_seq)=split('-',$1);		$at_richness_suggestion=~s/^\s+//;	
		if(abs($l_seq-$t_seq)<=$minimum_percent_at_diff){$at_richness_suggestion="NA";}	#if($at_richness_suggestion5=~/NA/){next;}
				#if($ref_repeat_suggestion ne $at_richness_suggestion6 and  $at_richness_suggestion6!~/NA/){$fp6++;}
			
		$hash_of_at_richness_records_final{$key}="$at_richness_suggestion [$l_seq-$t_seq]";
	}
#---- Degeneracy file -------------------------------------
my $file3="array_degeneracy_predictions.txt";

open(RD,"$file3") or print "$!";
my @arr_file3=<RD>;
close(RD);

foreach my $line(@arr_file3)
	{
		chomp $line; $line=~s/\r+//g;$line=~s/>//g;#$line=~s/ \[\S+\]//g;		
		my @arr_tmp=split('\t',$line);			
		my $key="$arr_tmp[0]-$arr_tmp[1]-$arr_tmp[2]";	
		
		my $degeneracy_in_first_repeat=$arr_tmp[5];
		my $degeneracy_in_last_repeat=$arr_tmp[6];
		
		my $top_degen; my $bottom_degen;my $total_arry_degeneracy=0;
					
		my $suggested_direction="NA";
					
					
		if($arr_tmp[4]=~/\[(\d+)-(\d+)\]/)
			{
				if($1+$2<=$permitted_no_of_bases)
					{
						$hash_of_degeneracy_records_final{$key}="$suggested_direction [$1-$2]";
						next;
					}
				
				$top_degen=$1 + $degeneracy_in_first_repeat;
				$bottom_degen=$2 + $degeneracy_in_last_repeat;
			}
					

		if($top_degen>$bottom_degen){$suggested_direction="R";}
		if($top_degen<$bottom_degen){$suggested_direction="F";}
							
		$hash_of_degeneracy_records_final{$key}="$suggested_direction [$top_degen-$bottom_degen]";
	}
#---- MFE file -------------------------------------
my $file4="all_repeats_MFE_suggestions.txt";

open(RD,"$file4") or print "$!";
my @arr_file4=<RD>;
close(RD);
foreach my $line(@arr_file4)
	{
		chomp $line; $line=~s/\r+//g;$line=~s/>//g;#$line=~s/ \[\S+\]//g;		
		my @arr_tmp=split('\t',$line);			
		my $key="$arr_tmp[0]-$arr_tmp[1]-$arr_tmp[2]";	
			
		if($arr_tmp[4+$mfe_target_base]=~/ \[(\S+)\]/)
			{
				my $mfe_scores=$1;
				my($forward_mfe,$reverse_mfe)=split(',',$mfe_scores);
						
				if(abs($forward_mfe)<$minimum_energy_threshold and abs($reverse_mfe)<$minimum_energy_threshold)
					{
						$hash_of_mfe_records_final{$key}="NA [$forward_mfe,$reverse_mfe]";
						next;
					}
				elsif(abs(abs($forward_mfe)-abs($reverse_mfe))<$minimum_energy_difference)
					{
						$hash_of_mfe_records_final{$key}="NA [$forward_mfe,$reverse_mfe]";
						next;
					}	
			}	
				
		$hash_of_mfe_records_final{$key}=$arr_tmp[4+$mfe_target_base];
	}



#---- Longer_leader file ---------------------------------		
my $file1="longer_leader_prediction.txt";

open(RD,"$file1") or print "$!";
my @arr_file1=<RD>;
close(RD);
foreach my $line(@arr_file1)
	{
		chomp $line; $line=~s/\r+//g;$line=~s/>//g;#$line=~s/ \[\S+\]//g;		
		my @arr_tmp=split('\t',$line);			
		my $key="$arr_tmp[0]-$arr_tmp[1]-$arr_tmp[2]";			
		$hash_of_longer_leader_records_final{$key}=$arr_tmp[10];
	}






#----------------------  first check for 100 random samples-----------------------------------

if($check_final_analysis==1)
{
	#print "Program is now applying the individual predictions to make a decision about the array direction for 100 random set for 100 times:\n";
	
	print "\n\nAnalysis of 100 random arrays from the training set:\n";
	
	my %hash_of_com_results;
	
	for(my $p=1;$p<=$q;$p++)
		{

			
		#print "Analysis of 100 random arrays from the remaining 230 arrays that were not part of the training set:\n";
		my $rec_count=1;
		#-------now pick 100 keys from the sample --------------------	

		my %hash_of_current_sample;
		my %hash_of_ref_repeat_direction_current_sample;
		
		
		for(1..$r*40)
			{
				my $found=0;
				while($found!=1)
					{
						my $random_no=int(rand($total_records-1));
					
						$random_no=$random_no+1;
						
						my $current_record=$hash_of_all_records{$random_no};

						#----- skip if not part of training set -------------------------------------
						if(not defined $hash_of_records_used_for_training{$current_record}){next;}

						#----------------------------------------------------------------------------
								
														
						if(not defined $hash_of_current_sample{$current_record})
							{
								my $key=$current_record;
								
								#print "$rec_count:$current_record\n";
								#$rec_count++;
								
								#------ skip the ones with motif ---------------------
								#if($hash_of_repeat_keys_and_motif_found{$key}){ next;}
								#------------------------------------------------------
								
																
								$key=~s/-/\t/g;
								#print "$key\n"; 
								#----- check if the prediction is not NA -------------------------------------------------------------------
								my $rec_line=`grep '$key' longer_leader_prediction.txt >&1`;
								my @arr_tmp1=split('\t',$rec_line);
								#if($arr_tmp1[10]=~/NA/){next;}
								
								$hash_of_ref_repeat_direction_current_sample{$current_record}=$hash_of_ref_repeat_direction{$current_record};
								$hash_of_current_sample{$current_record}=1;
								#print "$current_record\n";
								$found=1;
							}
					}
		
				my $current_no_of_entries=keys %hash_of_current_sample;
				if($current_no_of_entries>=$r)
					{
						#print "\$current_no_of_entries=$current_no_of_entries;";
						last;
					}
			}




			
			my $count_T=0;
			my $count_F=0;
			my $count_NA=0;
			my $total_c=0;			
						
			foreach my $key(sort keys %hash_of_ref_repeat_direction)		
				{	
					
					
					#--------- redundant/non-redundant toggle ------------------
					if(not $hash_of_current_sample{$key}){next;}
					#-----------------------------------------------------------
					
					my $score_F=0;
					my $score_R=0;
					my $score_NA=0;		
					
					
					my $ref_repeat=$hash_of_ref_repeat_direction{$key};
					
					if(not defined $hash_of_ref_repeat_direction{$key})
						{
							print "ref_repeat not found for $key\n";
						}	
					my $predicted_array_direction="NA";
					
					my $ref_direction=$hash_of_ref_repeat_direction{$key};




					#my $longer_leader_prediction="NA [1]";
					
					#my $longer_leader_prediction=$hash_of_longer_leader_records_final{$key};		
					#my $at_richness_prediction=$hash_of_at_richness_records_final{$key};
					#my $degeneracy_prediction=$hash_of_degeneracy_records_final{$key};
					#my $mfe_prediction=$hash_of_mfe_records_final{$key};
					


					

					#if($ref_repeat=~/(\S) \[/){$ref_direction=$1;}else{print "\n\$hash_of_ref_repeat_direction{$key}=$hash_of_ref_repeat_direction{$key}\n";next;}
					
					#print "\$ref_direction=$ref_direction\n";
					#------- apply individual predictions and their corresponding scoring parameter ---------
					
					#--------- motif toggle ------------------ -------------------------------------------------------------------
					if($check_motif==1 and defined $hash_of_repeat_keys_and_motif_found{$key})
						{
							my $motif_prediction=$hash_of_repeat_keys_and_motif_found{$key};	
							my $motif_prediction_score=$mean_motif_PPV-0.50;
							#------ MOTIF score ---------------------------------	
							if($motif_prediction!~/^NA/)
								{
									if($motif_prediction=~/(\S+)/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$motif_prediction_score;
												}
											else{
													$score_R=$score_R+$motif_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$motif_prediction_score;
								}
							
							#if($hash_of_repeat_keys_and_motif_found{$key}){$count_T++;$total_c++; next;}		 					
						}	
					
					#-------- A_and_T_toggle ----------------	
					if($check_A_and_T_in_repeats==1)
						{
							if($hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats{$key})
								{
									#------- analysis result of A and T distribution in repeats -----------------------------------------					
									my $A_T_in_repeat_prediction=$hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats{$key};
									my $A_T_distribution_in_repeats_score=$mean_A_T_distribution_in_repeats_PPV-0.50;
									if($A_T_in_repeat_prediction!~/NA/)
										{
											if($A_T_in_repeat_prediction eq "F")
												{
													$score_F=$score_F+$A_T_distribution_in_repeats_score;
												}
											else{
													$score_R=$score_R+$A_T_distribution_in_repeats_score;
												}
										}
								}
						}	
	
					if($check_mfe==1)
						{	
							my $mfe_prediction=$hash_of_mfe_records_final{$key};	
							my $mfe_prediction_score=$mean_mfe_ppv-0.50;
							
							#------ MFE score ---------------------------------	
							if($mfe_prediction!~/^NA/)
								{
									if($mfe_prediction=~/(\S+)/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$mfe_prediction_score;
												}
											else{
													$score_R=$score_R+$mfe_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$mfe_prediction_score;
								}
						}
					
											
	



					if($check_at_richness==1)
						{
							
							my $at_richness_prediction=$hash_of_at_richness_records_final{$key};
							my $at_richness_prediction_score=$mean_at_richness_ppv-0.50;
							
							#------ at_richness score ------------------------	
							if($at_richness_prediction!~/^NA/)
								{
									if($at_richness_prediction=~/(\S+) \[/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$at_richness_prediction_score;
												}
											elsif($1 eq "R")
												{
													$score_R=$score_R+$at_richness_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$at_richness_prediction_score;
								}
						}	
					
					if($check_degeneracy==1)
						{
							my $degeneracy_prediction=$hash_of_degeneracy_records_final{$key};
							my $degeneracy_prediction_score=$mean_degeneracy_ppv-0.50;
							
							#------- degeneracy score ------------------	
							if($degeneracy_prediction!~/^NA/)
								{
									if($degeneracy_prediction=~/(\S+) \[/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$degeneracy_prediction_score;
												}
											else{
													$score_R=$score_R+$degeneracy_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$degeneracy_prediction_score;
								}
							
						}
					
					
				
					if($check_longer_leader==1)
						{	
							my $longer_leader_prediction=$hash_of_longer_leader_records_final{$key};	
							my $longer_leader_prediction_score=$mean_longer_leader_ppv-0.50;
							
							#------- longer leader score -----------------------------------
							if($longer_leader_prediction!~/^NA/)
								{
									
									if($longer_leader_prediction=~/(\S+) \[/)
										{
											#print "$longer_leader_prediction\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$longer_leader_prediction_score;
												}
											else{
													$score_R=$score_R+$longer_leader_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$longer_leader_prediction_score;
								}							
						}	
	
										
					
					
					#----- now calculate direction ---------------------
					
					#if($score_NA>$score_F and $score_NA>$score_R){$predicted_array_direction="NA";}
					#els
					if($score_F>$score_R){$predicted_array_direction="F";}
					elsif($score_R>$score_F){$predicted_array_direction="R";}
					
					
					if($predicted_array_direction eq $ref_direction)
						{
							#print "$key \t$ref_repeat \t$longer_leader_prediction \t$at_richness_prediction \t$degeneracy_prediction \t $ref_direction-$predicted_array_direction\n";
							$count_T++;
						}
					elsif($predicted_array_direction ne "NA")
						{
							#print "$key    \t$ref_repeat \t$longer_leader_prediction \t$at_richness_prediction \t$degeneracy_prediction \t MFE:$mfe_prediction \t $ref_direction-$predicted_array_direction\n";
							$count_F++;
						}
					else{
							#print "$key    \t$ref_repeat \t$longer_leader_prediction \t$at_richness_prediction \t$degeneracy_prediction \t MFE:$mfe_prediction \t $ref_direction-$predicted_array_direction\n";
							$count_NA++;
							#
						}			
					$total_c++;				
				}
		
			my $final_PPV=0;
			if($count_T!=0)
				{
					$final_PPV=($count_T/($count_T+$count_F))*1; $final_PPV=sprintf("%.2f",$final_PPV);
				}	
				
			#print "Total: $total_c \tFinal PPV: $final_PPV [TP:$count_T FP:$count_F NP:$count_NA]\n\n";

			$hash_of_com_results{'COM'}{'TOTAL'}{$p}=$total_c;
			$hash_of_com_results{'COM'}{'PPV'}{$p}=$final_PPV;
			$hash_of_com_results{'COM'}{'TP'}{$p}=$count_T;
			$hash_of_com_results{'COM'}{'FP'}{$p}=$count_F;
			$hash_of_com_results{'COM'}{'NP'}{$p}=$count_NA;


		} # end of $p loop

	#--------- now calculate the SD and SE ------------------------------------------
	print "Analysis\tAvg_total_arrays_processed\tTP\ttp_sd\ttp_se\tFP\tfp_sd\tfp_se\tNP\tnp_sd\tnp_se\tPPV\tppv_sd\tppv_se\n";
	foreach my $analysis(sort keys %hash_of_com_results)
		{
			my $avg_total_arrays_processed=0;		
			my $ppv_avg=0;	my $ppv_sd=0;	my $ppv_se=0;		
			my $tp_avg=0;	my $tp_sd=0;	my $tp_se=0;		
			my $fp_avg=0;	my $fp_sd=0;	my $fp_se=0;		
			my $np_avg=0;	my $np_sd=0;	my $np_se=0;
			
			
			foreach my $type(sort{$a cmp $b} keys %{$hash_of_com_results{$analysis}})
				{				
					#print "$type=$type\n";
					my @arr_data;
					
					foreach my $pass(sort{$a <=> $b} keys %{$hash_of_com_results{$analysis}{$type}})
						{
							push(@arr_data,$hash_of_com_results{$analysis}{$type}{$pass});
							#print "@arr_data\n";
							#print "$analysis\t$type\t$pass\t$hash_of_analysis_results{$analysis}{$type}{$pass}\n";
						}
						
					#--- now calculate the mean, SD and SE ------------------
					my $ave;my $std; my $std_err=0;
					
					$ave = &average(\@arr_data);
					
					if($type!~/PPV/){$ave=sprintf("%.0f",$ave);}		
						
					if($type!~/TOTAL/)
						{
							$std = &stdev(\@arr_data);	
							
							if($std !=0 and $ave!=0)
								{				
									$std_err=($std/($ave ** 0.5)); $std_err=sprintf("%.2f",$std_err);
								}
							else{
									$std_err=0;
								}	
						}
					if($type=~/TOTAL/)
						{
							$avg_total_arrays_processed=int($ave);
							#print "Total=$avg_total_arrays_processed\n";
							next;
						}
					elsif($type=~/PPV/){$ppv_avg=$ave;$ppv_sd=$std;$ppv_se=$std_err;}
					elsif($type=~/TP/){$tp_avg=$ave;$tp_sd=$std;$tp_se=$std_err;}
					elsif($type=~/FP/){$fp_avg=$ave;$fp_sd=$std;$fp_se=$std_err;}
					elsif($type=~/NP/){$np_avg=$ave;$np_sd=$std;$np_se=$std_err;}
					#print "$analysis\t$type\t$ave\t$std\n";
					
					#if($type=~/PPV/){$mean_mfe_ppv=$ppv_avg;}
				}
				
			print "Combined\t$avg_total_arrays_processed\t$tp_avg\t$tp_sd\t$tp_se\t$fp_avg\t$fp_sd\t$fp_se\t$np_avg\t$np_sd\t$np_se\t$ppv_avg\t$ppv_sd\t$ppv_se\n";	
		}	



}


if($check_final_analysis==1)
{
	#print "Program is now applying the individual predictions to make a decision about the array direction for 100 random set for 100 times:\n";
	
	print "\n\nAnalysis of 100 random arrays from the test set:\n";
	
	my %hash_of_com_results;
	
	for(my $p=1;$p<=$q;$p++)
		{

			
		#print "Analysis of 100 random arrays from the remaining 230 arrays that were not part of the training set:\n";
		my $rec_count=1;
		#-------now pick 100 keys from the sample --------------------	

		my %hash_of_current_sample;
		my %hash_of_ref_repeat_direction_current_sample;
		
		
		for(1..$r*40)
			{
				my $found=0;
				while($found!=1)
					{
						my $random_no=int(rand($total_records-1));
					
						$random_no=$random_no+1;
						
						my $current_record=$hash_of_all_records{$random_no};

						#----- skip if not part of training set -------------------------------------
						if(not defined $hash_of_records_used_for_training{$current_record}){next;}

						#----------------------------------------------------------------------------
								
														
						if(not defined $hash_of_current_sample{$current_record})
							{
								my $key=$current_record;
								
								#print "$rec_count:$current_record\n";
								#$rec_count++;
								
								#------ skip the ones with motif ---------------------
								#if($hash_of_repeat_keys_and_motif_found{$key}){ next;}
								#------------------------------------------------------
								
																
								$key=~s/-/\t/g;
								#print "$key\n"; 
								#----- check if the prediction is not NA -------------------------------------------------------------------
								my $rec_line=`grep '$key' longer_leader_prediction.txt >&1`;
								my @arr_tmp1=split('\t',$rec_line);
								#if($arr_tmp1[10]=~/NA/){next;}
								
								$hash_of_ref_repeat_direction_current_sample{$current_record}=$hash_of_ref_repeat_direction{$current_record};
								$hash_of_current_sample{$current_record}=1;
								#print "$current_record\n";
								$found=1;
							}
					}
		
				my $current_no_of_entries=keys %hash_of_current_sample;
				if($current_no_of_entries>=$r)
					{
						#print "\$current_no_of_entries=$current_no_of_entries;";
						last;
					}
			}




			
			my $count_T=0;
			my $count_F=0;
			my $count_NA=0;
			my $total_c=0;			
						
			foreach my $key(sort keys %hash_of_ref_repeat_direction)		
				{	
					
					
					#--------- redundant/non-redundant toggle ------------------
					if(not $hash_of_current_sample{$key}){next;}
					#-----------------------------------------------------------
					
					my $score_F=0;
					my $score_R=0;
					my $score_NA=0;		
					
					
					my $ref_repeat=$hash_of_ref_repeat_direction{$key};
					
					if(not defined $hash_of_ref_repeat_direction{$key})
						{
							print "ref_repeat not found for $key\n";
						}	
					my $predicted_array_direction="NA";
					
					my $ref_direction=$hash_of_ref_repeat_direction{$key};




					#my $longer_leader_prediction="NA [1]";
					
					#my $longer_leader_prediction=$hash_of_longer_leader_records_final{$key};		
					#my $at_richness_prediction=$hash_of_at_richness_records_final{$key};
					#my $degeneracy_prediction=$hash_of_degeneracy_records_final{$key};
					#my $mfe_prediction=$hash_of_mfe_records_final{$key};
					


					

					#if($ref_repeat=~/(\S) \[/){$ref_direction=$1;}else{print "\n\$hash_of_ref_repeat_direction{$key}=$hash_of_ref_repeat_direction{$key}\n";next;}
					
					#print "\$ref_direction=$ref_direction\n";
					#------- apply individual predictions and their corresponding scoring parameter ---------
					
					#--------- motif toggle ------------------ -------------------------------------------------------------------
					if($check_motif==1 and defined $hash_of_repeat_keys_and_motif_found{$key})
						{
							my $motif_prediction=$hash_of_repeat_keys_and_motif_found{$key};	
							my $motif_prediction_score=$mean_motif_PPV-0.50;
							#------ MOTIF score ---------------------------------	
							if($motif_prediction!~/^NA/)
								{
									if($motif_prediction=~/(\S+)/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$motif_prediction_score;
												}
											else{
													$score_R=$score_R+$motif_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$motif_prediction_score;
								}
							
							#if($hash_of_repeat_keys_and_motif_found{$key}){$count_T++;$total_c++; next;}		 					
						}	
					
					#-------- A_and_T_toggle ----------------	
					if($check_A_and_T_in_repeats==1)
						{
							if($hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats{$key})
								{
									#------- analysis result of A and T distribution in repeats -----------------------------------------					
									my $A_T_in_repeat_prediction=$hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats{$key};
									my $A_T_distribution_in_repeats_score=$mean_A_T_distribution_in_repeats_PPV-0.50;
									if($A_T_in_repeat_prediction!~/NA/)
										{
											if($A_T_in_repeat_prediction eq "F")
												{
													$score_F=$score_F+$A_T_distribution_in_repeats_score;
												}
											else{
													$score_R=$score_R+$A_T_distribution_in_repeats_score;
												}
										}
								}
						}	
	
					if($check_mfe==1)
						{	
							my $mfe_prediction=$hash_of_mfe_records_final{$key};	
							my $mfe_prediction_score=$mean_mfe_ppv-0.50;
							
							#------ MFE score ---------------------------------	
							if($mfe_prediction!~/^NA/)
								{
									if($mfe_prediction=~/(\S+)/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$mfe_prediction_score;
												}
											else{
													$score_R=$score_R+$mfe_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$mfe_prediction_score;
								}
						}
					
											
	



					if($check_at_richness==1)
						{
							
							my $at_richness_prediction=$hash_of_at_richness_records_final{$key};
							my $at_richness_prediction_score=$mean_at_richness_ppv-0.50;
							
							#------ at_richness score ------------------------	
							if($at_richness_prediction!~/^NA/)
								{
									if($at_richness_prediction=~/(\S+) \[/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$at_richness_prediction_score;
												}
											elsif($1 eq "R")
												{
													$score_R=$score_R+$at_richness_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$at_richness_prediction_score;
								}
						}	
					
					if($check_degeneracy==1)
						{
							my $degeneracy_prediction=$hash_of_degeneracy_records_final{$key};
							my $degeneracy_prediction_score=$mean_degeneracy_ppv-0.50;
							
							#------- degeneracy score ------------------	
							if($degeneracy_prediction!~/^NA/)
								{
									if($degeneracy_prediction=~/(\S+) \[/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$degeneracy_prediction_score;
												}
											else{
													$score_R=$score_R+$degeneracy_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$degeneracy_prediction_score;
								}
							
						}
					
					
				
					if($check_longer_leader==1)
						{	
							my $longer_leader_prediction=$hash_of_longer_leader_records_final{$key};	
							my $longer_leader_prediction_score=$mean_longer_leader_ppv-0.50;
							
							#------- longer leader score -----------------------------------
							if($longer_leader_prediction!~/^NA/)
								{
									
									if($longer_leader_prediction=~/(\S+) \[/)
										{
											#print "$longer_leader_prediction\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$longer_leader_prediction_score;
												}
											else{
													$score_R=$score_R+$longer_leader_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$longer_leader_prediction_score;
								}							
						}	
	
										
					
					
					#----- now calculate direction ---------------------
					
					#if($score_NA>$score_F and $score_NA>$score_R){$predicted_array_direction="NA";}
					#els
					if($score_F>$score_R){$predicted_array_direction="F";}
					elsif($score_R>$score_F){$predicted_array_direction="R";}
					
					
					if($predicted_array_direction eq $ref_direction)
						{
							#print "$key \t$ref_repeat \t$longer_leader_prediction \t$at_richness_prediction \t$degeneracy_prediction \t $ref_direction-$predicted_array_direction\n";
							$count_T++;
						}
					elsif($predicted_array_direction ne "NA")
						{
							#print "$key    \t$ref_repeat \t$longer_leader_prediction \t$at_richness_prediction \t$degeneracy_prediction \t MFE:$mfe_prediction \t $ref_direction-$predicted_array_direction\n";
							$count_F++;
						}
					else{
							#print "$key    \t$ref_repeat \t$longer_leader_prediction \t$at_richness_prediction \t$degeneracy_prediction \t MFE:$mfe_prediction \t $ref_direction-$predicted_array_direction\n";
							$count_NA++;
							#
						}			
					$total_c++;				
				}
		
			my $final_PPV=0;
			if($count_T!=0)
				{
					$final_PPV=($count_T/($count_T+$count_F))*1; $final_PPV=sprintf("%.2f",$final_PPV);
				}	
				
			#print "Total: $total_c \tFinal PPV: $final_PPV [TP:$count_T FP:$count_F NP:$count_NA]\n\n";

			$hash_of_com_results{'COM'}{'TOTAL'}{$p}=$total_c;
			$hash_of_com_results{'COM'}{'PPV'}{$p}=$final_PPV;
			$hash_of_com_results{'COM'}{'TP'}{$p}=$count_T;
			$hash_of_com_results{'COM'}{'FP'}{$p}=$count_F;
			$hash_of_com_results{'COM'}{'NP'}{$p}=$count_NA;


		} # end of $p loop

	#--------- now calculate the SD and SE ------------------------------------------
	print "Analysis\tAvg_total_arrays_processed\tTP\ttp_sd\ttp_se\tFP\tfp_sd\tfp_se\tNP\tnp_sd\tnp_se\tPPV\tppv_sd\tppv_se\n";
	foreach my $analysis(sort keys %hash_of_com_results)
		{
			my $avg_total_arrays_processed=0;		
			my $ppv_avg=0;	my $ppv_sd=0;	my $ppv_se=0;		
			my $tp_avg=0;	my $tp_sd=0;	my $tp_se=0;		
			my $fp_avg=0;	my $fp_sd=0;	my $fp_se=0;		
			my $np_avg=0;	my $np_sd=0;	my $np_se=0;
			
			
			foreach my $type(sort{$a cmp $b} keys %{$hash_of_com_results{$analysis}})
				{				
					#print "$type=$type\n";
					my @arr_data;
					
					foreach my $pass(sort{$a <=> $b} keys %{$hash_of_com_results{$analysis}{$type}})
						{
							push(@arr_data,$hash_of_com_results{$analysis}{$type}{$pass});
							#print "@arr_data\n";
							#print "$analysis\t$type\t$pass\t$hash_of_analysis_results{$analysis}{$type}{$pass}\n";
						}
						
					#--- now calculate the mean, SD and SE ------------------
					my $ave;my $std; my $std_err=0;
					
					$ave = &average(\@arr_data);
					
					if($type!~/PPV/){$ave=sprintf("%.0f",$ave);}		
						
					if($type!~/TOTAL/)
						{
							$std = &stdev(\@arr_data);	
							
							if($std !=0 and $ave!=0)
								{				
									$std_err=($std/($ave ** 0.5)); $std_err=sprintf("%.2f",$std_err);
								}
							else{
									$std_err=0;
								}	
						}
					if($type=~/TOTAL/)
						{
							$avg_total_arrays_processed=int($ave);
							#print "Total=$avg_total_arrays_processed\n";
							next;
						}
					elsif($type=~/PPV/){$ppv_avg=$ave;$ppv_sd=$std;$ppv_se=$std_err;}
					elsif($type=~/TP/){$tp_avg=$ave;$tp_sd=$std;$tp_se=$std_err;}
					elsif($type=~/FP/){$fp_avg=$ave;$fp_sd=$std;$fp_se=$std_err;}
					elsif($type=~/NP/){$np_avg=$ave;$np_sd=$std;$np_se=$std_err;}
					#print "$analysis\t$type\t$ave\t$std\n";
					
					#if($type=~/PPV/){$mean_mfe_ppv=$ppv_avg;}
				}
				
			print "Combined\t$avg_total_arrays_processed\t$tp_avg\t$tp_sd\t$tp_se\t$fp_avg\t$fp_sd\t$fp_se\t$np_avg\t$np_sd\t$np_se\t$ppv_avg\t$ppv_sd\t$ppv_se\n";	
		}	



}



#------------ finally process all the 230 arrays used in  ------------------------------------------



if($check_final_analysis_of_all_arrays==1)
	{
	
	print "\n\nFinal prediction using all individual methods against first set of 230 arrays:\n";


	my $count_T=0;
	my $count_F=0;
	my $count_NA=0;
	my $total_c=0;	
	
	print "The following records are used:\n";		
	my $rec_index=1;					
	foreach my $key(sort keys %hash_of_ref_repeat_direction)		
		{			

			#----- skip the training set -------------------------------------
			if(not defined $hash_of_records_used_for_training{$key}){next;}
			else{
					print "$rec_index:$key\n";
					$rec_index++;
				}	
			#----------------------------------------------------------------------------

						
					my $score_F=0;
					my $score_R=0;
					my $score_NA=0;				
			
			
					my($accession,$start,$stop)=split('-',$key);
					my $species="Unknown";
					if($hash_of_archaeal_acc{$accession})
						{
							$species="Archaea";#next;
						}
					elsif($hash_of_bacterial_acc{$accession})
						{
							$species="Bacteria";#next;
						}	

					#---------------- get the reference direction for this array --------------------------------------------------	
					
					my $ref_repeat=$hash_of_ref_repeat_direction{$key};
					
					if(not defined $hash_of_ref_repeat_direction{$key})
						{
							print "ref_repeat not found for $key\n";
						}
						
					my $predicted_array_direction="NA";					
					my $ref_direction=$hash_of_ref_repeat_direction{$key};	
					
					
					#---- get the % AT of the genome --------------------
					my $skip_checking_genome_at_p=0;
					
					if($skip_checking_genome_at_p==1)
						{
							my $genome_at_percentage=0;
							my $percent_at_in_left_window=0;
							my $percent_at_in_right_window=0;
							
							my @arr_genome_lines=`grep -A 1 '>' tmp/$accession\.fna >&1`;
							foreach my $g_line(@arr_genome_lines)
								{
									if($g_line!~/>/)
										{
											chomp $g_line;$g_line=~s/\r//g;
											my $genome_length=length($g_line);
											
											my $no_of_at=$g_line=~s/[A|T]/[A|T]/gi;
											
											$genome_at_percentage=($no_of_at/$genome_length)*100;
											$genome_at_percentage=sprintf("%.2f",$genome_at_percentage);
											
											#---- get left window and right window sequences for AT% calculation --------------------
											my $left_window=substr($g_line,$start-$at_window_length,$at_window_length);
											my $right_window=substr($g_line,$stop,$at_window_length);
											
											my $no_of_at_in_left=$left_window=~s/[A|T]/[A|T]/gi;
											my $no_of_at_in_right=$right_window=~s/[A|T]/[A|T]/gi;
											
											$percent_at_in_left_window=sprintf("%.2f",(($no_of_at_in_left/length($left_window))*100));
											$percent_at_in_right_window=sprintf("%.2f",(($no_of_at_in_right/length($left_window))*100));
										}
								}
						
						}		

					
					
					
					
											

					
					#------- apply individual predictions and their corresponding scoring parameter ---------
					
					
					
					
					
					#--------- repeat specific analysis ------ -------------------------------------------------------------------
					if($check_motif==1 and defined $hash_of_repeat_keys_and_motif_found{$key})
						{
							my $motif_prediction=$hash_of_repeat_keys_and_motif_found{$key};	
							my $motif_prediction_score=$mean_motif_PPV-0.50;
							#------ MOTIF score ---------------------------------	
							if($motif_prediction!~/^NA/)
								{
									if($motif_prediction=~/(\S+)/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$motif_prediction_score;
												}
											else{
													$score_R=$score_R+$motif_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$motif_prediction_score;
								}
							
							#if($hash_of_repeat_keys_and_motif_found{$key}){$count_T++;$total_c++; next;}		 					
						}	
					
						
					if($check_A_and_T_in_repeats==1)
						{
							if($hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats{$key})
								{
									#------- analysis result of A and T distribution in repeats -----------------------------------------					
									my $A_T_in_repeat_prediction=$hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats{$key};
									my $A_T_distribution_in_repeats_score=$mean_A_T_distribution_in_repeats_PPV-0.50;
									if($A_T_in_repeat_prediction!~/NA/)
										{
											if($A_T_in_repeat_prediction eq "F")
												{
													$score_F=$score_F+$A_T_distribution_in_repeats_score;
												}
											else{
													$score_R=$score_R+$A_T_distribution_in_repeats_score;
												}
										}
								}
						}	
	
										
					if($check_mfe==1)
						{	
							my $mfe_prediction=$hash_of_mfe_records_final{$key};	
							my $mfe_prediction_score=$mean_mfe_ppv-0.50;
							
							#------ MFE score ---------------------------------	
							if($mfe_prediction!~/^NA/)
								{
									if($mfe_prediction=~/(\S+)/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$mfe_prediction_score;
												}
											else{
													$score_R=$score_R+$mfe_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$mfe_prediction_score;
								}
						}
					
		


					#--- array specific analysis -----------------------------------------------------------------------------------

					if($check_at_richness==1)
						{
							
							my $at_richness_prediction=$hash_of_at_richness_records_final{$key};
							my $at_richness_prediction_score=$mean_at_richness_ppv-0.50;
							
							#------ at_richness score ------------------------	
							if($at_richness_prediction!~/^NA/)
								{
									if($at_richness_prediction=~/(\S+) \[/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$at_richness_prediction_score;
												}
											elsif($1 eq "R")
												{
													$score_R=$score_R+$at_richness_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$at_richness_prediction_score;
								}
						}	
					
					if($check_degeneracy==1)
						{
							my $degeneracy_prediction=$hash_of_degeneracy_records_final{$key};
							my $degeneracy_prediction_score=$mean_degeneracy_ppv-0.50;
							
							#------- degeneracy score ------------------	
							if($degeneracy_prediction!~/^NA/)
								{
									if($degeneracy_prediction=~/(\S+) \[/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$degeneracy_prediction_score;
												}
											else{
													$score_R=$score_R+$degeneracy_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$degeneracy_prediction_score;
								}
							
						}
					
					
				
					if($check_longer_leader==1)
						{	
							my $longer_leader_prediction=$hash_of_longer_leader_records_final{$key};	
							my $longer_leader_prediction_score=$mean_longer_leader_ppv-0.50;
							
							#------- longer leader score -----------------------------------
							if($longer_leader_prediction!~/^NA/)
								{
									
									if($longer_leader_prediction=~/(\S+) \[/)
										{
											#print "$longer_leader_prediction\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$longer_leader_prediction_score;
												}
											else{
													$score_R=$score_R+$longer_leader_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$longer_leader_prediction_score;
								}							
						}	
	
					
					



						
	
					
					
					
					#----- now calculate direction ---------------------
					
					#if($score_NA>$score_F and $score_NA>$score_R){$predicted_array_direction="NA";}
					#els
					if($score_F>$score_R){$predicted_array_direction="F";}
					elsif($score_R>$score_F){$predicted_array_direction="R";}
					
					
					if($ref_direction eq $predicted_array_direction)
						{
							#print "$key \t$ref_repeat \t$longer_leader_prediction \t$at_richness_prediction \t$degeneracy_prediction \t $ref_direction-$predicted_array_direction\n";
							$count_T++;
						}
					elsif($predicted_array_direction ne "NA")
						{
							#print "$key [$species]   \t$ref_repeat \t$longer_leader_prediction \t$at_richness_prediction \t$degeneracy_prediction \t MFE:$mfe_prediction \t $ref_direction-$predicted_array_direction\t$hash_of_ref_repeats{$key}\n";
							
							#print "$key [$species]\n";
							
							$count_F++;
						}
					else{
							#print "$key    \t$ref_repeat \t$longer_leader_prediction \t$at_richness_prediction \t$degeneracy_prediction \t MFE:$mfe_prediction \t $ref_direction-$predicted_array_direction\n";
							$count_NA++;
							#
						}			
					$total_c++;				
		}
	
	my $final_PPV=0;
	if($count_T !=0)
		{
			$final_PPV=($count_T/($count_T+$count_F))*1; $final_PPV=sprintf("%.2f",$final_PPV);
		}		
	#print "Total: $total_c \tFinal PPV: $final_PPV [TP:$count_T FP:$count_F NP:$count_NA]\n\n";
	
	print "Analysis\tTotal\tTP\tFP\tNP\tPPV\n";
	print "Final\t$total_c\t$count_T\t$count_F\t$count_NA\t$final_PPV\n\n";

	}










####################################    ananlysis of remaining sets #####################################################################
print "\n\n############## Program will now apply the models on the second set of 230 arrays \n\n ";

print "The second set contains the following arrays:\n";
my $rec_index=1;
foreach my $key(sort keys %hash_of_ref_repeat_direction)		
	{			

		#----- skip the training set -------------------------------------
			if(defined $hash_of_records_used_for_training{$key}){next;}
			else{
					print "$rec_index:$key\n";
					$rec_index++;
				}
	}		

#----------------------  first check for 100 random samples from remaining set-----------------------------------


if($check_final_analysis==1)
{
	#print "Program is now applying the individual predictions to make a decision about the array direction for 100 random set for 100 times:\n";
	
	print "\n\nAnalysis of 100 randomly selected arrays from the remaining 230 arrays that were not part of the first set used in building models:\n";
	
	my %hash_of_com_results;
	
	for(my $p=1;$p<=$q;$p++)
		{

			
		#print "Analysis of 100 random arrays from the remaining 230 arrays that were not part of the training set:\n";
		my $rec_count=1;
		#-------now pick 100 keys from the sample --------------------	

		my %hash_of_current_sample;
		my %hash_of_ref_repeat_direction_current_sample;
		
		
		for(1..$r*40)
			{
				my $found=0;
				while($found!=1)
					{
						my $random_no=int(rand($total_records-1));
					
						$random_no=$random_no+1;
						
						my $current_record=$hash_of_all_records{$random_no};

						#----- skip if not part of training set -------------------------------------
						if(defined $hash_of_records_used_for_training{$current_record}){next;}

						#----------------------------------------------------------------------------
								
														
						if(not defined $hash_of_current_sample{$current_record})
							{
								my $key=$current_record;
								
								#print "$rec_count:$current_record\n";
								#$rec_count++;
								
								#------ skip the ones with motif ---------------------
								#if($hash_of_repeat_keys_and_motif_found{$key}){ next;}
								#------------------------------------------------------
								
																
								$key=~s/-/\t/g;
								#print "$key\n"; 
								#----- check if the prediction is not NA -------------------------------------------------------------------
								my $rec_line=`grep '$key' longer_leader_prediction.txt >&1`;
								my @arr_tmp1=split('\t',$rec_line);
								#if($arr_tmp1[10]=~/NA/){next;}
								
								$hash_of_ref_repeat_direction_current_sample{$current_record}=$hash_of_ref_repeat_direction{$current_record};
								$hash_of_current_sample{$current_record}=1;
								#print "$current_record\n";
								$found=1;
							}
					}
		
				my $current_no_of_entries=keys %hash_of_current_sample;
				if($current_no_of_entries>=$r)
					{
						#print "\$current_no_of_entries=$current_no_of_entries;";
						last;
					}
			}




			
			my $count_T=0;
			my $count_F=0;
			my $count_NA=0;
			my $total_c=0;			
						
			foreach my $key(sort keys %hash_of_ref_repeat_direction)		
				{	
					
					
					#--------- redundant/non-redundant toggle ------------------
					if(not $hash_of_current_sample{$key}){next;}
					#-----------------------------------------------------------
					
					my $score_F=0;
					my $score_R=0;
					my $score_NA=0;		
					
					
					my $ref_repeat=$hash_of_ref_repeat_direction{$key};
					
					if(not defined $hash_of_ref_repeat_direction{$key})
						{
							print "ref_repeat not found for $key\n";
						}	
					my $predicted_array_direction="NA";
					
					my $ref_direction=$hash_of_ref_repeat_direction{$key};




					#my $longer_leader_prediction="NA [1]";
					
					#my $longer_leader_prediction=$hash_of_longer_leader_records_final{$key};		
					#my $at_richness_prediction=$hash_of_at_richness_records_final{$key};
					#my $degeneracy_prediction=$hash_of_degeneracy_records_final{$key};
					#my $mfe_prediction=$hash_of_mfe_records_final{$key};
					


					

					#if($ref_repeat=~/(\S) \[/){$ref_direction=$1;}else{print "\n\$hash_of_ref_repeat_direction{$key}=$hash_of_ref_repeat_direction{$key}\n";next;}
					
					#print "\$ref_direction=$ref_direction\n";
					#------- apply individual predictions and their corresponding scoring parameter ---------
					
					#--------- motif toggle ------------------ -------------------------------------------------------------------
					if($check_motif==1 and defined $hash_of_repeat_keys_and_motif_found{$key})
						{
							my $motif_prediction=$hash_of_repeat_keys_and_motif_found{$key};	
							my $motif_prediction_score=$mean_motif_PPV-0.50;
							#------ MOTIF score ---------------------------------	
							if($motif_prediction!~/^NA/)
								{
									if($motif_prediction=~/(\S+)/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$motif_prediction_score;
												}
											else{
													$score_R=$score_R+$motif_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$motif_prediction_score;
								}
							
							#if($hash_of_repeat_keys_and_motif_found{$key}){$count_T++;$total_c++; next;}		 					
						}	
					
					#-------- A_and_T_toggle ----------------	
					if($check_A_and_T_in_repeats==1)
						{
							if($hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats{$key})
								{
									#------- analysis result of A and T distribution in repeats -----------------------------------------					
									my $A_T_in_repeat_prediction=$hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats{$key};
									my $A_T_distribution_in_repeats_score=$mean_A_T_distribution_in_repeats_PPV-0.50;
									if($A_T_in_repeat_prediction!~/NA/)
										{
											if($A_T_in_repeat_prediction eq "F")
												{
													$score_F=$score_F+$A_T_distribution_in_repeats_score;
												}
											else{
													$score_R=$score_R+$A_T_distribution_in_repeats_score;
												}
										}
								}
						}	
	
					if($check_mfe==1)
						{	
							my $mfe_prediction=$hash_of_mfe_records_final{$key};	
							my $mfe_prediction_score=$mean_mfe_ppv-0.50;
							
							#------ MFE score ---------------------------------	
							if($mfe_prediction!~/^NA/)
								{
									if($mfe_prediction=~/(\S+)/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$mfe_prediction_score;
												}
											else{
													$score_R=$score_R+$mfe_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$mfe_prediction_score;
								}
						}
					
											
	



					if($check_at_richness==1)
						{
							
							my $at_richness_prediction=$hash_of_at_richness_records_final{$key};
							my $at_richness_prediction_score=$mean_at_richness_ppv-0.50;
							
							#------ at_richness score ------------------------	
							if($at_richness_prediction!~/^NA/)
								{
									if($at_richness_prediction=~/(\S+) \[/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$at_richness_prediction_score;
												}
											elsif($1 eq "R")
												{
													$score_R=$score_R+$at_richness_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$at_richness_prediction_score;
								}
						}	
					
					if($check_degeneracy==1)
						{
							my $degeneracy_prediction=$hash_of_degeneracy_records_final{$key};
							my $degeneracy_prediction_score=$mean_degeneracy_ppv-0.50;
							
							#------- degeneracy score ------------------	
							if($degeneracy_prediction!~/^NA/)
								{
									if($degeneracy_prediction=~/(\S+) \[/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$degeneracy_prediction_score;
												}
											else{
													$score_R=$score_R+$degeneracy_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$degeneracy_prediction_score;
								}
							
						}
					
					
				
					if($check_longer_leader==1)
						{	
							my $longer_leader_prediction=$hash_of_longer_leader_records_final{$key};	
							my $longer_leader_prediction_score=$mean_longer_leader_ppv-0.50;
							
							#------- longer leader score -----------------------------------
							if($longer_leader_prediction!~/^NA/)
								{
									
									if($longer_leader_prediction=~/(\S+) \[/)
										{
											#print "$longer_leader_prediction\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$longer_leader_prediction_score;
												}
											else{
													$score_R=$score_R+$longer_leader_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$longer_leader_prediction_score;
								}							
						}	
	
										
					
					
					#----- now calculate direction ---------------------
					
					#if($score_NA>$score_F and $score_NA>$score_R){$predicted_array_direction="NA";}
					#els
					if($score_F>$score_R){$predicted_array_direction="F";}
					elsif($score_R>$score_F){$predicted_array_direction="R";}
					
					
					if($predicted_array_direction eq $ref_direction)
						{
							#print "$key \t$ref_repeat \t$longer_leader_prediction \t$at_richness_prediction \t$degeneracy_prediction \t $ref_direction-$predicted_array_direction\n";
							$count_T++;
						}
					elsif($predicted_array_direction ne "NA")
						{
							#print "$key    \t$ref_repeat \t$longer_leader_prediction \t$at_richness_prediction \t$degeneracy_prediction \t MFE:$mfe_prediction \t $ref_direction-$predicted_array_direction\n";
							$count_F++;
						}
					else{
							#print "$key    \t$ref_repeat \t$longer_leader_prediction \t$at_richness_prediction \t$degeneracy_prediction \t MFE:$mfe_prediction \t $ref_direction-$predicted_array_direction\n";
							$count_NA++;
							#
						}			
					$total_c++;				
				}
		
			my $final_PPV=0;
			if($count_T!=0)
				{
					$final_PPV=($count_T/($count_T+$count_F))*1; $final_PPV=sprintf("%.2f",$final_PPV);
				}	
				
			#print "Total: $total_c \tFinal PPV: $final_PPV [TP:$count_T FP:$count_F NP:$count_NA]\n\n";

			$hash_of_com_results{'COM'}{'TOTAL'}{$p}=$total_c;
			$hash_of_com_results{'COM'}{'PPV'}{$p}=$final_PPV;
			$hash_of_com_results{'COM'}{'TP'}{$p}=$count_T;
			$hash_of_com_results{'COM'}{'FP'}{$p}=$count_F;
			$hash_of_com_results{'COM'}{'NP'}{$p}=$count_NA;


		} # end of $p loop

	#--------- now calculate the SD and SE ------------------------------------------
	print "Analysis\tAvg_total_arrays_processed\tTP\ttp_sd\ttp_se\tFP\tfp_sd\tfp_se\tNP\tnp_sd\tnp_se\tPPV\tppv_sd\tppv_se\n";
	foreach my $analysis(sort keys %hash_of_com_results)
		{
			my $avg_total_arrays_processed=0;		
			my $ppv_avg=0;	my $ppv_sd=0;	my $ppv_se=0;		
			my $tp_avg=0;	my $tp_sd=0;	my $tp_se=0;		
			my $fp_avg=0;	my $fp_sd=0;	my $fp_se=0;		
			my $np_avg=0;	my $np_sd=0;	my $np_se=0;
			
			
			foreach my $type(sort{$a cmp $b} keys %{$hash_of_com_results{$analysis}})
				{				
					#print "$type=$type\n";
					my @arr_data;
					
					foreach my $pass(sort{$a <=> $b} keys %{$hash_of_com_results{$analysis}{$type}})
						{
							push(@arr_data,$hash_of_com_results{$analysis}{$type}{$pass});
							#print "@arr_data\n";
							#print "$analysis\t$type\t$pass\t$hash_of_analysis_results{$analysis}{$type}{$pass}\n";
						}
						
					#--- now calculate the mean, SD and SE ------------------
					my $ave;my $std; my $std_err=0;
					
					$ave = &average(\@arr_data);
					
					if($type!~/PPV/){$ave=sprintf("%.0f",$ave);}		
						
					if($type!~/TOTAL/)
						{
							$std = &stdev(\@arr_data);	
							
							if($std !=0 and $ave!=0)
								{				
									$std_err=($std/($ave ** 0.5)); $std_err=sprintf("%.2f",$std_err);
								}
							else{
									$std_err=0;
								}	
						}
					if($type=~/TOTAL/)
						{
							$avg_total_arrays_processed=int($ave);
							#print "Total=$avg_total_arrays_processed\n";
							next;
						}
					elsif($type=~/PPV/){$ppv_avg=$ave;$ppv_sd=$std;$ppv_se=$std_err;}
					elsif($type=~/TP/){$tp_avg=$ave;$tp_sd=$std;$tp_se=$std_err;}
					elsif($type=~/FP/){$fp_avg=$ave;$fp_sd=$std;$fp_se=$std_err;}
					elsif($type=~/NP/){$np_avg=$ave;$np_sd=$std;$np_se=$std_err;}
					#print "$analysis\t$type\t$ave\t$std\n";
					
					#if($type=~/PPV/){$mean_mfe_ppv=$ppv_avg;}
				}
				
			print "Combined\t$avg_total_arrays_processed\t$tp_avg\t$tp_sd\t$tp_se\t$fp_avg\t$fp_sd\t$fp_se\t$np_avg\t$np_sd\t$np_se\t$ppv_avg\t$ppv_sd\t$ppv_se\n";	
		}	



}



#------------ finally process all the remaining 230 arrays ------------------------------------------



if($check_final_analysis_of_all_arrays==1)
	{
	
	print "\n\nFinal prediction against all 230 arrays from the second set:\n";


	my $count_T=0;
	my $count_F=0;
	my $count_NA=0;
	my $total_c=0;	
	
	#print "The following records are used:\n";		
	my $rec_index=1;					
	foreach my $key(sort keys %hash_of_ref_repeat_direction)		
		{			

			#----- skip the training set -------------------------------------
			if(defined $hash_of_records_used_for_training{$key}){next;}
			else{
					#print "$rec_index:$key\n";
					$rec_index++;
				}	
			#----------------------------------------------------------------------------

						
					my $score_F=0;
					my $score_R=0;
					my $score_NA=0;				
			
			
					my($accession,$start,$stop)=split('-',$key);
					my $species="Unknown";
					if($hash_of_archaeal_acc{$accession})
						{
							$species="Archaea";#next;
						}
					elsif($hash_of_bacterial_acc{$accession})
						{
							$species="Bacteria";#next;
						}	

					#---------------- get the reference direction for this array --------------------------------------------------	
					
					my $ref_repeat=$hash_of_ref_repeat_direction{$key};
					
					if(not defined $hash_of_ref_repeat_direction{$key})
						{
							print "ref_repeat not found for $key\n";
						}
						
					my $predicted_array_direction="NA";					
					my $ref_direction=$hash_of_ref_repeat_direction{$key};	
					
					
					#---- get the % AT of the genome --------------------
					my $skip_checking_genome_at_p=0;
					
					if($skip_checking_genome_at_p==1)
						{
							my $genome_at_percentage=0;
							my $percent_at_in_left_window=0;
							my $percent_at_in_right_window=0;
							
							my @arr_genome_lines=`grep -A 1 '>' tmp/$accession\.fna >&1`;
							foreach my $g_line(@arr_genome_lines)
								{
									if($g_line!~/>/)
										{
											chomp $g_line;$g_line=~s/\r//g;
											my $genome_length=length($g_line);
											
											my $no_of_at=$g_line=~s/[A|T]/[A|T]/gi;
											
											$genome_at_percentage=($no_of_at/$genome_length)*100;
											$genome_at_percentage=sprintf("%.2f",$genome_at_percentage);
											
											#---- get left window and right window sequences for AT% calculation --------------------
											my $left_window=substr($g_line,$start-$at_window_length,$at_window_length);
											my $right_window=substr($g_line,$stop,$at_window_length);
											
											my $no_of_at_in_left=$left_window=~s/[A|T]/[A|T]/gi;
											my $no_of_at_in_right=$right_window=~s/[A|T]/[A|T]/gi;
											
											$percent_at_in_left_window=sprintf("%.2f",(($no_of_at_in_left/length($left_window))*100));
											$percent_at_in_right_window=sprintf("%.2f",(($no_of_at_in_right/length($left_window))*100));
										}
								}
						
						}		

					
					
					
					
											

					
					#------- apply individual predictions and their corresponding scoring parameter ---------
					
					
					
					
					
					#--------- repeat specific analysis ------ -------------------------------------------------------------------
					if($check_motif==1 and defined $hash_of_repeat_keys_and_motif_found{$key})
						{
							my $motif_prediction=$hash_of_repeat_keys_and_motif_found{$key};	
							my $motif_prediction_score=$mean_motif_PPV-0.50;
							#------ MOTIF score ---------------------------------	
							if($motif_prediction!~/^NA/)
								{
									if($motif_prediction=~/(\S+)/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$motif_prediction_score;
												}
											else{
													$score_R=$score_R+$motif_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$motif_prediction_score;
								}
							
							#if($hash_of_repeat_keys_and_motif_found{$key}){$count_T++;$total_c++; next;}		 					
						}	
					
						
					if($check_A_and_T_in_repeats==1)
						{
							if($hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats{$key})
								{
									#------- analysis result of A and T distribution in repeats -----------------------------------------					
									my $A_T_in_repeat_prediction=$hash_of_repeat_keys_and_A_and_T_distribution_analysis_in_repeats{$key};
									my $A_T_distribution_in_repeats_score=$mean_A_T_distribution_in_repeats_PPV-0.50;
									if($A_T_in_repeat_prediction!~/NA/)
										{
											if($A_T_in_repeat_prediction eq "F")
												{
													$score_F=$score_F+$A_T_distribution_in_repeats_score;
												}
											else{
													$score_R=$score_R+$A_T_distribution_in_repeats_score;
												}
										}
								}
						}	
	
										
					if($check_mfe==1)
						{	
							my $mfe_prediction=$hash_of_mfe_records_final{$key};	
							my $mfe_prediction_score=$mean_mfe_ppv-0.50;
							
							#------ MFE score ---------------------------------	
							if($mfe_prediction!~/^NA/)
								{
									if($mfe_prediction=~/(\S+)/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$mfe_prediction_score;
												}
											else{
													$score_R=$score_R+$mfe_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$mfe_prediction_score;
								}
						}
					
		


					#--- array specific analysis -----------------------------------------------------------------------------------

					if($check_at_richness==1)
						{
							
							my $at_richness_prediction=$hash_of_at_richness_records_final{$key};
							my $at_richness_prediction_score=$mean_at_richness_ppv-0.50;
							
							#------ at_richness score ------------------------	
							if($at_richness_prediction!~/^NA/)
								{
									if($at_richness_prediction=~/(\S+) \[/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$at_richness_prediction_score;
												}
											elsif($1 eq "R")
												{
													$score_R=$score_R+$at_richness_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$at_richness_prediction_score;
								}
						}	
					
					if($check_degeneracy==1)
						{
							my $degeneracy_prediction=$hash_of_degeneracy_records_final{$key};
							my $degeneracy_prediction_score=$mean_degeneracy_ppv-0.50;
							
							#------- degeneracy score ------------------	
							if($degeneracy_prediction!~/^NA/)
								{
									if($degeneracy_prediction=~/(\S+) \[/)
										{
											#print "$1\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$degeneracy_prediction_score;
												}
											else{
													$score_R=$score_R+$degeneracy_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$degeneracy_prediction_score;
								}
							
						}
					
					
				
					if($check_longer_leader==1)
						{	
							my $longer_leader_prediction=$hash_of_longer_leader_records_final{$key};	
							my $longer_leader_prediction_score=$mean_longer_leader_ppv-0.50;
							
							#------- longer leader score -----------------------------------
							if($longer_leader_prediction!~/^NA/)
								{
									
									if($longer_leader_prediction=~/(\S+) \[/)
										{
											#print "$longer_leader_prediction\n";
											if($1 eq "F")
												{
													$score_F=$score_F+$longer_leader_prediction_score;
												}
											else{
													$score_R=$score_R+$longer_leader_prediction_score;
												}
										}
								}
							else{
									$score_NA=$score_NA+$longer_leader_prediction_score;
								}							
						}	
	
					
					



						
	
					
					
					
					#----- now calculate direction ---------------------
					
					#if($score_NA>$score_F and $score_NA>$score_R){$predicted_array_direction="NA";}
					#els
					if($score_F>$score_R){$predicted_array_direction="F";}
					elsif($score_R>$score_F){$predicted_array_direction="R";}
					
					
					if($ref_direction eq $predicted_array_direction)
						{
							#print "$key \t$ref_repeat \t$longer_leader_prediction \t$at_richness_prediction \t$degeneracy_prediction \t $ref_direction-$predicted_array_direction\n";
							$count_T++;
						}
					elsif($predicted_array_direction ne "NA")
						{
							#print "$key [$species]   \t$ref_repeat \t$longer_leader_prediction \t$at_richness_prediction \t$degeneracy_prediction \t MFE:$mfe_prediction \t $ref_direction-$predicted_array_direction\t$hash_of_ref_repeats{$key}\n";
							
							#print "$key [$species]\n";
							
							$count_F++;
						}
					else{
							#print "$key    \t$ref_repeat \t$longer_leader_prediction \t$at_richness_prediction \t$degeneracy_prediction \t MFE:$mfe_prediction \t $ref_direction-$predicted_array_direction\n";
							$count_NA++;
							#
						}			
					$total_c++;				
		}
	
	my $final_PPV=0;
	if($count_T !=0)
		{
			$final_PPV=($count_T/($count_T+$count_F))*1; $final_PPV=sprintf("%.2f",$final_PPV);
		}		
	#print "Total: $total_c \tFinal PPV: $final_PPV [TP:$count_T FP:$count_F NP:$count_NA]\n\n";
	
	print "Analysis\tTotal\tTP\tFP\tNP\tPPV\n";
	print "Final\t$total_c\t$count_T\t$count_F\t$count_NA\t$final_PPV\n\n";

	}






}

























#######-------- subs -----------------------------------------------------------------------

sub at_richness_analysis()
	{
		my($i,$hash_of_current_sample,$hash_of_at_richness_records)=@_;
		
		#print "\n\nAnalysing AT richness in the flanks:\n";
	
		my $file="at_richness_distribution.txt";

		open(RD,"$file") or print "$!";
		my @arr_file=<RD>;
		close(RD);

		my $at_richness_total=0;
		my $at_richness_ppv=0;
		my $at_richness_tp=0;
		my $at_richness_fp=0;
		my $at_richness_np=0;
	#--------------------------------------------
	my %hash_of_at_richness_results;   #------ not needed for main analysis
#	for (my $i=.5;$i<15;$i=$i+.5)
#	{
		
		my $ref_rep_and_leader_false_prediction=0;

		my $ref_rep_and_leader_true_prediction=0;
		my $ref_rep_and_degeneracy_prediction=0;

		my $tp1=0;
		my $tp2=0;
		my $tp3=0;
		my $tp4=0;
		my $tp5=0;
		my $tp6=0;
		my $tp7=0;
		my $tp8=0;
		my $tp9=0;
		my $tp10=0;
		my $tp11=0;

		my $fp1=0;
		my $fp2=0;
		my $fp3=0;
		my $fp4=0;
		my $fp5=0;
		my $fp6=0;
		my $fp7=0;
		my $fp8=0;
		my $fp9=0;	
		my $fp10=0;	
		my $fp11=0;	
		
		
		
		
		
		my $total=0;	
		my $minimum_difference=$i;
		my $skip_line=0;
		foreach my $line(@arr_file)
			{
				#print "$line\n";next;
				#if($line!~/100\%/){next;}
				
				$skip_line=0;
				
				chomp $line; $line=~s/\r+//g;
				my @arr_tmp=split('\t+',$line);
				
				my $accession=$arr_tmp[0];
				my $key="$arr_tmp[0]-$arr_tmp[1]-$arr_tmp[2]";
				
				if(not $hash_of_current_sample->{$key}){next;}
				
				#----- store this record in hash, will be used to optimize scoring parameters --------------------------------
				if(not defined $hash_of_at_richness_records->{$key} and $i==1)
					{
						$hash_of_at_richness_records->{$key}=$arr_tmp[6];
					}
				else{
						#$skip_line=1;
						#next;
					}	
				#-------------------------------------------------------------------------------------------------------------	
				
				my $ref_repeat_suggestion=$arr_tmp[3]; 		$ref_repeat_suggestion=~s/ \[(\S+)\]//g; 	my $p_sim=$1;	$ref_repeat_suggestion=~s/^\s+//; #if($arr_tmp[3]!~/NA/){next;}
				#my $longer_leader_suggestion=$arr_tmp[4]; 	$longer_leader_suggestion=~s/ \[\S+\]//g;	$longer_leader_suggestion=~s/^\s+//;
				
				
				
				
				my $at_richness_suggestion1=$arr_tmp[4]; 	$at_richness_suggestion1=~s/ \[(\S+)\]\%AT//g;	my($l_seq1,$t_seq1)=split('-',$1);		$at_richness_suggestion1=~s/^\s+//;
						if(abs($l_seq1-$t_seq1)<=$minimum_difference){$at_richness_suggestion1="NA";}	#if($at_richness_suggestion5=~/NA/){next;}
						if($ref_repeat_suggestion ne $at_richness_suggestion1 and  $at_richness_suggestion1!~/NA/){$fp1++;}
						
				
				my $at_richness_suggestion2=$arr_tmp[5]; 	$at_richness_suggestion2=~s/ \[(\S+)\]\%AT//g;	my($l_seq2,$t_seq2)=split('-',$1);		$at_richness_suggestion2=~s/^\s+//;	
						if(abs($l_seq2-$t_seq2)<=$minimum_difference){$at_richness_suggestion2="NA";}	#if($at_richness_suggestion5=~/NA/){next;}
						if($ref_repeat_suggestion ne $at_richness_suggestion2 and  $at_richness_suggestion2!~/NA/){$fp2++;}
						
						
				my $at_richness_suggestion3=$arr_tmp[6]; 	$at_richness_suggestion3=~s/ \[(\S+)\]\%AT//g;	my($l_seq3,$t_seq3)=split('-',$1);		$at_richness_suggestion3=~s/^\s+//;	
						if(abs($l_seq3-$t_seq3)<=$minimum_difference){$at_richness_suggestion3="NA";}	#if($at_richness_suggestion5=~/NA/){next;}
						if($ref_repeat_suggestion ne $at_richness_suggestion3 and  $at_richness_suggestion3!~/NA/){$fp3++;}
						
						
				my $at_richness_suggestion4=$arr_tmp[7]; 	$at_richness_suggestion4=~s/ \[(\S+)\]\%AT//g;	my($l_seq4,$t_seq4)=split('-',$1);		$at_richness_suggestion4=~s/^\s+//;	
						if(abs($l_seq4-$t_seq4)<=$minimum_difference){$at_richness_suggestion4="NA";}	#if($at_richness_suggestion5=~/NA/){next;}
						if($ref_repeat_suggestion ne $at_richness_suggestion4 and  $at_richness_suggestion4!~/NA/){$fp4++;}
				
				#-----------------------------------------------
					
				my $at_richness_suggestion5=$arr_tmp[8]; 	$at_richness_suggestion5=~s/ \[(\S+)\]\%AT//g;	my($l_seq5,$t_seq5)=split('-',$1);	$at_richness_suggestion5=~s/^\s+//;		
						if(abs($l_seq5-$t_seq5)<=$minimum_difference){$at_richness_suggestion5="NA";}	#if($at_richness_suggestion5=~/NA/){next;}
						if($ref_repeat_suggestion ne $at_richness_suggestion5 and  $at_richness_suggestion5!~/NA/){$fp5++;}
				
				#----------------------------------------------
				
						
				my $at_richness_suggestion6=$arr_tmp[9]; 	$at_richness_suggestion6=~s/ \[(\S+)\]\%AT//g;	my($l_seq6,$t_seq6)=split('-',$1);		$at_richness_suggestion6=~s/^\s+//;	
						if(abs($l_seq6-$t_seq6)<=$minimum_difference){$at_richness_suggestion6="NA";}	#if($at_richness_suggestion5=~/NA/){next;}
						if($ref_repeat_suggestion ne $at_richness_suggestion6 and  $at_richness_suggestion6!~/NA/){$fp6++;}
						
						
				my $at_richness_suggestion7=$arr_tmp[10]; 	$at_richness_suggestion7=~s/ \[(\S+)\]\%AT//g;	my($l_seq7,$t_seq7)=split('-',$1);		$at_richness_suggestion7=~s/^\s+//;	
						if(abs($l_seq7-$t_seq7)<=$minimum_difference){$at_richness_suggestion7="NA";}	#if($at_richness_suggestion5=~/NA/){next;}
						if($ref_repeat_suggestion ne $at_richness_suggestion7 and  $at_richness_suggestion7!~/NA/){$fp7++;}
						
						
				my $at_richness_suggestion8=$arr_tmp[11]; 	$at_richness_suggestion8=~s/ \[(\S+)\]\%AT//g;	my($l_seq8,$t_seq8)=split('-',$1);		$at_richness_suggestion8=~s/^\s+//;	
						if(abs($l_seq8-$t_seq8)<=$minimum_difference){$at_richness_suggestion8="NA";}	#if($at_richness_suggestion5=~/NA/){next;}
						if($ref_repeat_suggestion ne $at_richness_suggestion8 and  $at_richness_suggestion8!~/NA/){$fp8++;}
						
						
				my $at_richness_suggestion9=$arr_tmp[12]; 	$at_richness_suggestion9=~s/ \[(\S+)\]\%AT//g;	my($l_seq9,$t_seq9)=split('-',$1);		$at_richness_suggestion9=~s/^\s+//;		
						if(abs($l_seq9-$t_seq9)<=$minimum_difference){$at_richness_suggestion9="NA";}	#if($at_richness_suggestion5=~/NA/){next;}
						if($ref_repeat_suggestion ne $at_richness_suggestion9 and  $at_richness_suggestion9!~/NA/){$fp9++;}
						
				
				my $at_richness_suggestion10=$arr_tmp[13]; 	$at_richness_suggestion10=~s/ \[(\S+)\]\%AT//g;	my($l_seq10,$t_seq10)=split('-',$1);		$at_richness_suggestion10=~s/^\s+//;		
						if(abs($l_seq10-$t_seq10)<=$minimum_difference){$at_richness_suggestion10="NA";}	#if($at_richness_suggestion5=~/NA/){next;}
						if($ref_repeat_suggestion ne $at_richness_suggestion10 and  $at_richness_suggestion10!~/NA/){$fp10++;}
						
				my $at_richness_suggestion11=$arr_tmp[14]; 	$at_richness_suggestion11=~s/ \[(\S+)\]\%AT//g;	my($l_seq11,$t_seq11)=split('-',$1);		$at_richness_suggestion11=~s/^\s+//;		
						if(abs($l_seq11-$t_seq11)<=$minimum_difference){$at_richness_suggestion11="NA";}	#if($at_richness_suggestion5=~/NA/){next;}
						if($ref_repeat_suggestion ne $at_richness_suggestion11 and  $at_richness_suggestion11!~/NA/){$fp11++;}	
				#my $degeneracy_suggestion=$arr_tmp[14]; 		$degeneracy_suggestion=~s/ \[\S+\]//g;		$degeneracy_suggestion=~s/^\s+//;
				
				
				
				#if($p_sim<100){next;}
				#if($line=~/NA/){next;}
				#--------------------------------------------------------------------------------------------------------------------------------
				#print "$arr_tmp[0]\t$ref_repeat_suggestion\t$at_richness_suggestion1\t$at_richness_suggestion2\t$at_richness_suggestion3\t$at_richness_suggestion4\t$at_richness_suggestion5\t$at_richness_suggestion6\t$at_richness_suggestion7\t$at_richness_suggestion8\t$at_richness_suggestion9\n";
				
				#if($ref_repeat_suggestion eq $longer_leader_suggestion){$ref_rep_and_leader_true_prediction++;};	#else{$ref_rep_and_leader_false_prediction++;}
					#else{$ref_rep_and_leader_false_prediction++;}
				#if($ref_repeat_suggestion eq $degeneracy_suggestion){$ref_rep_and_degeneracy_prediction++;};	#else{$ref_rep_and_leader_false_prediction++;}
				
				if($ref_repeat_suggestion eq $at_richness_suggestion1){$tp1++;};
				if($ref_repeat_suggestion eq $at_richness_suggestion2){$tp2++;};
				if($ref_repeat_suggestion eq $at_richness_suggestion3){$tp3++;};
				if($ref_repeat_suggestion eq $at_richness_suggestion4){$tp4++;};
				
				#---- for unfixed version, this one is the best--------------------------------------------------
				if($ref_repeat_suggestion eq $at_richness_suggestion5){$tp5++;};
				#------------------------------------------------------------------------------------------------
				
				if($ref_repeat_suggestion eq $at_richness_suggestion6){$tp6++;};
				if($ref_repeat_suggestion eq $at_richness_suggestion7){$tp7++;};
				if($ref_repeat_suggestion eq $at_richness_suggestion8){$tp8++;};
				if($ref_repeat_suggestion eq $at_richness_suggestion9){$tp9++;};
				if($ref_repeat_suggestion eq $at_richness_suggestion10){$tp10++;};
				if($ref_repeat_suggestion eq $at_richness_suggestion11){$tp11++;};
				$total++;
			}
			#print "-------------------------------------------------------------------------------------------------------\n";

			#my $p1=sprintf("%.2f",($ref_rep_and_leader_true_prediction/$tp)*100);
			#my $p3=sprintf("%.2f",($ref_rep_and_degeneracy_prediction/$tp)*100);
			
			if($skip_line==1){next;}
			
			my $ppv1=sprintf("%.2f",($tp1/($tp1+$fp1))*1);
			
			#---- for unfixed version, this one is the best-------------------------------------------------------------------------------------------------------
			my $ppv2=sprintf("%.2f",($tp2/($tp2+$fp2))*1);
			#-----------------------------------------------------------------------------------------------------------------------------------------------------
			
			my $ppv3=sprintf("%.2f",($tp3/($tp3+$fp3))*1);
			my $ppv4=sprintf("%.2f",($tp4/($tp4+$fp4))*1);
			my $ppv5=sprintf("%.2f",($tp5/($tp5+$fp5))*1);
			my $ppv6=sprintf("%.2f",($tp6/($tp6+$fp6))*1);
			my $ppv7=sprintf("%.2f",($tp7/($tp7+$fp7))*1);
			my $ppv8=sprintf("%.2f",($tp8/($tp8+$fp8))*1);
			my $ppv9=sprintf("%.2f",($tp9/($tp9+$fp9))*1);
			my $ppv10=sprintf("%.2f",($tp10/($tp10+$fp10))*1);
			my $ppv11=sprintf("%.2f",($tp11/($tp11+$fp11))*1);

			#print qq~Min.Diff: $minimum_difference:\tTotal:$tp\tPPV: $ppv1\% [TP:$tp1 FP:$fp1]\tPPV: $ppv2\% [TP:$tp2 FP:$fp2]\tPPV: $ppv3\% [TP:$tp3 FP:$fp3]\tPPV: $ppv4\% [TP:$tp4 FP:$fp4]\tPPV: $ppv5\% [TP:$tp5 FP:$fp5]\tPPV: $ppv6\% [TP:$tp6 FP:$fp6]\tPPV: $ppv7\% [TP:$tp7 FP:$fp7]\tPPV: $ppv8\% [TP:$tp8 FP:$fp8]\tPPV: $ppv9\% [TP:$tp9 FP:$fp9]\n~;
			my $np1= $total-$tp1-$fp1;
			my $np2= $total-$tp2-$fp2;
			my $np3= $total-$tp3-$fp3;
			my $np4= $total-$tp4-$fp4;
			my $np5= $total-$tp5-$fp5;
			my $np6= $total-$tp6-$fp6;
			my $np7= $total-$tp7-$fp7;
			my $np8= $total-$tp8-$fp8;
			my $np9= $total-$tp9-$fp9;
			my $np10= $total-$tp10-$fp10;
			my $np11= $total-$tp11-$fp11;
			
			#print qq~			


#			print qq~\n\n\n\nTotal:$total and Min.Diff $minimum_difference\%\tPPV\tTP\tFP\tNP
#Window 45nt\t$ppv1\t$tp1\t$fp1\t$np1
#Window 60nt\t$ppv2\t$tp2\t$fp2\t$np2
#Window 75nt\t$ppv3\t$tp3\t$fp3\t$np3
#Window 90nt\t$ppv4\t$tp4\t$fp4\t$np4
#Window 105nt\t$ppv5\t$tp5\t$fp5\t$np5
#Window 120nt\t$ppv6\t$tp6\t$fp6\t$np6
#Window 135nt\t$ppv7\t$tp7\t$fp7\t$np7
#Window 150nt\t$ppv8\t$tp8\t$fp8\t$np8
#Window 165nt\t$ppv9\t$tp9\t$fp9\t$np9\n~;
			
			my $window_1=15;
			my $window_2=30;
			my $window_3=45;
			my $window_4=60;
			my $window_5=75;
			my $window_6=90;
			my $window_7=105;
			my $window_8=120;
			my $window_9=135;
			my $window_10=150;
			my $window_11=165;
			
			
			$hash_of_at_richness_results{$window_1}{$minimum_difference}="$ppv1\t$tp1\t$fp1\t$np1";
			$hash_of_at_richness_results{$window_2}{$minimum_difference}="$ppv2\t$tp2\t$fp2\t$np2";
			$hash_of_at_richness_results{$window_3}{$minimum_difference}="$ppv3\t$tp3\t$fp3\t$np3";
			$hash_of_at_richness_results{$window_4}{$minimum_difference}="$ppv4\t$tp4\t$fp4\t$np4";
			$hash_of_at_richness_results{$window_5}{$minimum_difference}="$ppv5\t$tp5\t$fp5\t$np5";
			$hash_of_at_richness_results{$window_6}{$minimum_difference}="$ppv6\t$tp6\t$fp6\t$np6";
			$hash_of_at_richness_results{$window_7}{$minimum_difference}="$ppv7\t$tp7\t$fp7\t$np7";
			$hash_of_at_richness_results{$window_8}{$minimum_difference}="$ppv8\t$tp8\t$fp8\t$np8";
			$hash_of_at_richness_results{$window_9}{$minimum_difference}="$ppv9\t$tp9\t$fp9\t$np9";
			$hash_of_at_richness_results{$window_10}{$minimum_difference}="$ppv10\t$tp10\t$fp10\t$np10";
			$hash_of_at_richness_results{$window_11}{$minimum_difference}="$ppv11\t$tp11\t$fp11\t$np11";
			
			#if($i==1)  #---means 10% longer than other
			#	{
					$at_richness_total=$total;
					$at_richness_ppv=$ppv7;					
					$at_richness_tp=$tp7;
					$at_richness_fp=$fp7;
					$at_richness_np=$np7;
			#	}
#		}
		
		
		##print qq~\nWindow and Min. % Diff.\tPPV\tTP\tFP\tNP\n~;		
		#foreach my $window(sort{$a<=>$b}keys %hash_of_at_richness_results)
		#	{
		#		my $tp=0;
		#		foreach my $min_diff(sort{$a<=>$b}keys %{$hash_of_at_richness_results{$window}})
		#			{
		#				my $x_coord=$window+$tp;
		#				
		#				print "$x_coord\t$hash_of_at_richness_results{$window}{$min_diff}\n";
		#				$tp=$tp+0.25;
		#			}
		#	}
		#------------------------------------------------------------------------------------------------
		
		
		#return($at_richness_total,$at_richness_ppv,$at_richness_tp,$at_richness_fp,$at_richness_np);
		return($total,$window_1,$ppv1,$tp1,$fp1,$np1,$total,$window_2,$ppv2,$tp2,$fp2,$np2,$total,$window_3,$ppv3,$tp3,$fp3,$np3,$total,$window_4,$ppv4,$tp4,$fp4,$np4,$total,$window_5,$ppv5,$tp5,$fp5,$np5,$total,$window_6,$ppv6,$tp6,$fp6,$np6,$total,$window_7,$ppv7,$tp7,$fp7,$np7,$total,$window_8,$ppv8,$tp8,$fp8,$np8,$total,$window_9,$ppv9,$tp9,$fp9,$np9,$total,$window_10,$ppv10,$tp10,$fp10,$np10,$total,$window_11,$ppv11,$tp11,$fp11,$np11);
	}



sub array_degeneracy_analysis()
	{
		my($k,$permitted_no_of_bases,$hash_of_current_sample,$hash_of_degeneracy_records)=@_;
		
		#print "\n\nAnalysing array degeneracy:\n";
	
		my $file="array_degeneracy_predictions.txt";

		open(RD,"$file") or print "$!";
		my @arr_file=<RD>;
		close(RD);
		
		my $count1;
		
		#print "Total\tPPV\tTP\tFP\tNP\n";
		
		
		#for(my $k=0;$k<=2;$k++)
		#	{
				
			$count1=0;
			my $count2=0;
			my $count3=0;
			my $count4=0;
			foreach my $line(@arr_file)
				{
					chomp $line; $line=~s/\r+//g;
					
					my @arr_tmp=split('\t',$line);
					
					my $accession=$arr_tmp[0];
					my $key="$arr_tmp[0]-$arr_tmp[1]-$arr_tmp[2]";
					
					#--------- redundant/non-redundant toggle ------------------
					if(not $hash_of_current_sample->{$key}){next;}
					#-----------------------------------------------------------
					my $degeneracy_in_first_repeat=0;
					my $degeneracy_in_last_repeat=0;
					
					$degeneracy_in_first_repeat=$arr_tmp[5];
					$degeneracy_in_last_repeat=$arr_tmp[6];
					
					my $top_degen; my $bottom_degen;
					my $total_arry_degeneracy=0;
					
					my $suggested_direction="NA";
					
					
					if($arr_tmp[4]=~/\[(\d+)-(\d+)\]/)
						{
							$top_degen=$1 +$degeneracy_in_first_repeat;
							$bottom_degen=$2+ $degeneracy_in_last_repeat;
						}
					
					if(abs($top_degen+$bottom_degen)>0)
							{
								$total_arry_degeneracy=abs($top_degen+$bottom_degen);
								#($top_degen>0 and $bottom_degen>0)
								#$arr_tmp[4]="NA";
								$count1++;
							}
					if($top_degen>$bottom_degen){$suggested_direction="R";}
					if($top_degen<$bottom_degen){$suggested_direction="F";}
					
					#----- store this record in hash, will be used to optimize scoring parameters --------------------------------
					if(not defined $hash_of_degeneracy_records->{$key})
						{
							$hash_of_degeneracy_records->{$key}=$suggested_direction;  #$arr_tmp[4];
						}	
					else{
							#next;
						}	
					#-------------------------------------------------------------------------------------------------------------			
					if($arr_tmp[3]=~/ \[\S+\]/ ){$arr_tmp[3]=~s/ \[(\S+)\]//g; }
					if($arr_tmp[4]=~/ \[\S+\]/ ){$arr_tmp[4]=~s/ \[(\S+)\]//g; }
					

					if($arr_tmp[3] eq $suggested_direction and ($top_degen>$k or $bottom_degen>$k) and $total_arry_degeneracy>$permitted_no_of_bases)
						{
							$count2++;
						}
					elsif($suggested_direction ne "NA" and ($top_degen>$k or $bottom_degen>$k) and $total_arry_degeneracy>$permitted_no_of_bases)
						{

							$count3++;
						}
					#elsif($top_degen>$k or $bottom_degen>$k)
					else{
							#print "$line\n";
							$count4++;
							#next;
						}
					#else{next;}	

				}
				
			my $degen_ppv=($count2/($count2+$count3));
			   $degen_ppv=sprintf("%.2f",$degen_ppv);
			#my $total=$count2+$count3+$count4;
			#print "$k allowed mutations\t$degen_ppv\t$count2\t$count3\t$count4\n";	
			
			#if($k==0)
			#	{
					
					my $degeneracy_total=$count2+$count3+$count4;
					my $degeneracy_ppv=$degen_ppv;
					my $degeneracy_tp=$count2;
					my $degeneracy_fp=$count3;
					my $degeneracy_np=$count4;
			#	}
		#}

		#print "Total arrays :$count1\n\n\n";
		
		#--------------------------------------------------------------------------------------------------
		return($count1,$degeneracy_ppv,$degeneracy_tp,$degeneracy_fp,$degeneracy_np);
	}


sub repeats_structural_analysis()
	{
		my($j,$minimum_energy_threshold,$hash_of_current_sample,$hash_of_mfe_records)=@_;
		
		#print "\n\nAnalysing MFE:\n";
	
		my $file="all_repeats_MFE_suggestions.txt";

		open(RD,"$file") or print "$!";
		my @arr_file=<RD>;
		close(RD);
		
		my %hash_of_nr_repeats;	
		

		my $count1;
		#print "Total: \tPPV\tTP\tFP\tNP\n";
		
		#for(my $j=0;$j<=11;$j++)
		#{
			$count1=0;
			my $count2=0;
			my $count3=0;
			my $count4=0;
			my $mfe_ppv_tmp=0;
			my %hash_of_unique_repeats;
			
			my $i=4+$j; #--- Because the first 4 columns contains Accession, start, stop and reference direction
			
			foreach my $line(@arr_file)
			{
				chomp $line; $line=~s/\r$//;
				#print "\#:$line\n";
				
				my @arr_tmp=split('\t',$line);
				

				
				my $accession=$arr_tmp[0];
				my $key="$arr_tmp[0]-$arr_tmp[1]-$arr_tmp[2]";
				
				my $repeat=$arr_tmp[$#arr_tmp];
				my $reverse_comp_repeat=$repeat;	$reverse_comp_repeat=reverse $reverse_comp_repeat;$reverse_comp_repeat=~tr/ACGT/TGCA/;
				#--------- redundant/non-redundant toggle ------------------
				if(not $hash_of_current_sample->{$key}){next;}#print "OK\n";
				#-----------------------------------------------------------
		
				my $mfe_scores;my $forward_mfe; my $reverse_mfe;
				
				if($arr_tmp[$i]=~/ \[(\S+)\]/)
					{
						$mfe_scores=$1;
						($forward_mfe,$reverse_mfe)=split(',',$mfe_scores);
						
						if(abs($forward_mfe)<$minimum_energy_threshold and abs($reverse_mfe)<$minimum_energy_threshold)
							{
								$count1++;
								$count4++;
								next;
							}
					}
					
				#print "$arr_tmp[3]\n";
						
				if($arr_tmp[3]=~/ \[\S+\]/ ){$arr_tmp[3]=~s/ \[(\S+)\]//g; }
				if($arr_tmp[$i]=~/ \[\S+\]/ ){$arr_tmp[$i]=~s/ \[(\S+)\]//g; }	
				
				#---- check the presence of GAA in the last 6 bases ------------------------------------------------------------
				
				my $str1_0=substr($repeat,0,4);
				my $str1=substr($repeat,-5,3);
				
				my $str2_0=substr($reverse_comp_repeat,0,4);
				my $str2=substr($reverse_comp_repeat,-5,3);
				

				#------------- make sure only one prediction is used for each unique repeat-----------------------------------
				
				if(not defined $hash_of_mfe_records->{$key})
					{
						$hash_of_nr_repeats{$arr_tmp[$#arr_tmp]}=1;
						$hash_of_mfe_records->{$key}=$arr_tmp[9];					
					}
				else{
						#next;
					}
					
				
				#if(defined $hash_of_unique_repeats{$repeat} or defined $hash_of_unique_repeats{$reverse_comp_repeat})
				#	{
				#		next;
				#	}	
				#else{					
				#		$hash_of_unique_repeats{$repeat}=1;
				#	}	
						
				
				
				#--------------------------------- check reference and result --------------
				#print "\t\t$arr_tmp[3] \t $arr_tmp[$i]\n";
							
				if($arr_tmp[3] eq $arr_tmp[$i])
					{
						$count2++;
					}
				elsif($arr_tmp[$i] ne "NA")
					{
						$count3++;
					}
				else{
						$count4++;
					}
					
				#$line=~s/\s+\[/_\[/g;$line=~s/\s+/\t/g;$line=~s/_/ /g;	
				#print "$line\n";			
				$count1++;	
			}
		#print "$count1,$count2,$count3\n";	
		if($count2==0 or ($count2+$count3)==0)
			{
				$mfe_ppv_tmp=0;
			}
		else{
				$mfe_ppv_tmp=($count2/($count2+$count3))*1;
				$mfe_ppv_tmp=sprintf("%.2f",$mfe_ppv_tmp);
			}
		
		#if($j==6)
		#	{
		my $mfe_total=$count1;
		my $mfe_ppv=$mfe_ppv_tmp;
		my $mfe_tp=$count2;
		my $mfe_fp=$count3;
		my $mfe_np=$count4;
		#	}
		#print qq~Excluded $j\t$mfe_ppv_tmp\t$count2\t$count3\t$count4\n~;	
		
		#}
		#print "Total:$count1\n\n";
		
		#---------------------------------------------------------------------------
		
		return($mfe_total,$mfe_ppv,$mfe_tp,$mfe_fp,$mfe_np);
	}


sub analysis_of_longer_leader()
	{
		my($longer_leader_total,$longer_leader_ppv,$longer_leader_tp,$longer_leader_fp,$longer_leader_np,$hash_of_current_sample,$hash_of_longer_leader_records)=@_;
		
		#print "\nAnalysing longer leader:\n";
	
		my $file="longer_leader_prediction.txt";

		open(RD,"$file") or print "$!";
		my @arr_file=<RD>;
		close(RD);
	
		my $count2_tab5=0;my $count3_tab5=0;my $count4_tab5=0;
		my $count2_tab6=0;my $count3_tab6=0;my $count4_tab6=0;
		my $count2_tab7=0;my $count3_tab7=0;my $count4_tab7=0;
		my $count2_tab8=0;my $count3_tab8=0;my $count4_tab8=0;
		my $count2_tab9=0;my $count3_tab9=0;my $count4_tab9=0;
		my $count2_tab10=0;my $count3_tab10=0;my $count4_tab10=0;
	
		my $count=0;
		
		foreach my $line(@arr_file)
			{
				chomp $line; $line=~s/\r+//g;$line=~s/>//g;#$line=~s/ \[\S+\]//g;
			
			
				my @arr_tmp=split('\t',$line);
			
				my $accession=$arr_tmp[0];
				my $key="$arr_tmp[0]-$arr_tmp[1]-$arr_tmp[2]";
			
				#---- print the LB and UB ---------------------
				my $dist=$arr_tmp[4];$dist=~s/^Dist: //;
				if($dist=~/\[(\d+)-(\d+)\]/)
					{
						#if($1<$lb){$lb=$1;}
						#if($2<$lb){$lb=$2;}
						
						#if($1>$ub){$ub=$1;}
						#if($2>$ub){$ub=$2;}				
					}
				#----------------------------------------------
				
				
				if(not $hash_of_current_sample->{$key}){next;}
				
				
				#----- store this record in hash, will be used to optimize scoring parameters --------------------------------
				if(not defined $hash_of_longer_leader_records->{$key})
					{
						$hash_of_longer_leader_records->{$key}=$arr_tmp[10];
					}
				else{
						#next;
					}		
				#-------------------------------------------------------------------------------------------------------------	
						
				#print "$line\n";
				
				#$hash_of_ref_repeat_direction{$key}=$arr_tmp[3];
				
				
			
				if($arr_tmp[3]=~/ \[\S+\]/ ){$arr_tmp[3]=~s/ \[(\S+)\]//g; }
				
				$arr_tmp[5]=~s/ \[1\]//g;
				$arr_tmp[6]=~s/ \[1\.2\]//g;
				$arr_tmp[7]=~s/ \[1\.4\]//g;
				$arr_tmp[8]=~s/ \[1\.6\]//g;
				$arr_tmp[9]=~s/ \[1\.8\]//g;
				$arr_tmp[10]=~s/ \[2]//g;
				
				if($arr_tmp[3] eq $arr_tmp[5]){	$count2_tab5++;	}
				elsif($arr_tmp[5] ne "NA")
					{
						$count3_tab5++;
						#print "$line\n";
					}
				else{
						$count4_tab5++;
						#print "$line\n";
					}	
						
				if($arr_tmp[3] eq $arr_tmp[6]){	$count2_tab6++;	}elsif($arr_tmp[6] ne "NA"){$count3_tab6++;}else{$count4_tab6++;}		
				if($arr_tmp[3] eq $arr_tmp[7]){	$count2_tab7++;	}elsif($arr_tmp[7] ne "NA"){$count3_tab7++;}else{$count4_tab7++;}		
				if($arr_tmp[3] eq $arr_tmp[8]){	$count2_tab8++;	}elsif($arr_tmp[8] ne "NA"){$count3_tab8++;}else{$count4_tab8++;}		
				if($arr_tmp[3] eq $arr_tmp[9]){	$count2_tab9++;	}elsif($arr_tmp[9] ne "NA"){$count3_tab9++;}else{$count4_tab9++;}		
				if($arr_tmp[3] eq $arr_tmp[10]){$count2_tab10++;}elsif($arr_tmp[10] ne "NA"){$count3_tab10++;}else{$count4_tab10++;}		
	
			
				#$count1_tab5++;	
				$count++;
			}
		
		my $ppv_tab5=($count2_tab5/($count2_tab5+$count3_tab5))*1;$ppv_tab5=sprintf("%.2f",$ppv_tab5);
		my $ppv_tab6=($count2_tab6/($count2_tab6+$count3_tab6))*1;$ppv_tab6=sprintf("%.2f",$ppv_tab6);
		my $ppv_tab7=($count2_tab7/($count2_tab7+$count3_tab7))*1;$ppv_tab7=sprintf("%.2f",$ppv_tab7);
		my $ppv_tab8=($count2_tab8/($count2_tab8+$count3_tab8))*1;$ppv_tab8=sprintf("%.2f",$ppv_tab8);
		my $ppv_tab9=($count2_tab9/($count2_tab9+$count3_tab9))*1;$ppv_tab9=sprintf("%.2f",$ppv_tab9);
		my $ppv_tab10=($count2_tab10/($count2_tab10+$count3_tab10))*1;$ppv_tab10=sprintf("%.2f",$ppv_tab10);
		
		#print "Total: $count\tPPV\tTP\tFP\tNP\nMinimum 1 base longer\t$ppv_tab5\t$count2_tab5\t$count3_tab5\t$count4_tab5\nMinimum 20% longer\t$ppv_tab6\t$count2_tab6\t$count3_tab6\t$count4_tab6\nMinimum 40% longer\t$ppv_tab7\t$count2_tab7\t$count3_tab7\t$count4_tab7\nMinimum 60% longer\t$ppv_tab8\t$count2_tab8\t$count3_tab8\t$count4_tab8\nMinimum 80% longer\t$ppv_tab9\t$count2_tab9\t$count3_tab9\t$count4_tab9\nMinimum double long\t$ppv_tab10\t$count2_tab10\t$count3_tab10\t$count4_tab10\n";	
		
		
		$longer_leader_total=$count;
		$longer_leader_ppv=$ppv_tab10;
		$longer_leader_tp=$count2_tab10;
		$longer_leader_fp=$count3_tab10;
		$longer_leader_np=$count4_tab10;
		
		my $ll_total1=$count;my $ll_ppv1=$ppv_tab5; my $ll_tp1=$count2_tab5; my $ll_fp1=$count3_tab5; my $ll_np1=$count4_tab5;
		my $ll_total2=$count;my $ll_ppv2=$ppv_tab6; my $ll_tp2=$count2_tab6; my $ll_fp2=$count3_tab6; my $ll_np2=$count4_tab6;
		my $ll_total3=$count;my $ll_ppv3=$ppv_tab7; my $ll_tp3=$count2_tab7; my $ll_fp3=$count3_tab7; my $ll_np3=$count4_tab7;
		my $ll_total4=$count;my $ll_ppv4=$ppv_tab8; my $ll_tp4=$count2_tab8; my $ll_fp4=$count3_tab8; my $ll_np4=$count4_tab8;
		my $ll_total5=$count;my $ll_ppv5=$ppv_tab9; my $ll_tp5=$count2_tab9; my $ll_fp5=$count3_tab9; my $ll_np5=$count4_tab9;
		my $ll_total6=$count;my $ll_ppv6=$ppv_tab10; my $ll_tp6=$count2_tab10; my $ll_fp6=$count3_tab10; my $ll_np6=$count4_tab10;
		
		return($ll_total1,$ll_ppv1,$ll_tp1,$ll_fp1,$ll_np1,$ll_total2,$ll_ppv2,$ll_tp2,$ll_fp2,$ll_np2,$ll_total3,$ll_ppv3,$ll_tp3,$ll_fp3,$ll_np3,$ll_total4,$ll_ppv4,$ll_tp4,$ll_fp4,$ll_np4,$ll_total5,$ll_ppv5,$ll_tp5,$ll_fp5,$ll_np5,$ll_total6,$ll_ppv6,$ll_tp6,$ll_fp6,$ll_np6);
	}






sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        
        $average=sprintf("%.2f",$average);
        return $average;
}

sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        
        $std=sprintf("%.2f",$std);
        
        return $std;
}
