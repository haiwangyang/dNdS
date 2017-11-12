#!/usr/bin/perl
use warnings;
use strict;
use Bio::Perl;
use lib "/home/haiwangyang/program/perl_lib/alignDB/lib";
use YAML qw(Dump Load DumpFile LoadFile);
use AlignDB::IntSpan;
use AlignDB::Util;
use Data::Dumper;;
use List::MoreUtils qw(firstidx all any uniq );
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Math::Combinatorics;
use Statistics::Lite qw(mean);
$|++;

######################################
### get a hash that links codon and aa
######################################
our %genetic_code = get_genetic_code();

##############################################################
###  get a hash describing syn-site and non-site of each codon
###  The values were retrieved from MEGA4
##############################################################
my ($codon2syn_site, $codon2non_site) = get_codon2syn_site_AND_codon2non_site();
my %codon2syn_site = %$codon2syn_site;
my %codon2non_site = %$codon2non_site;

################################################################
### count syn-num and non-num of each path when codon1 => codon2
################################################################
my %Codon1_Codon2_syn_non_num_info = get_Codon1_Codon2_syn_non_num_info();

########################################################
### calculate dNdS for CDS fasta files in current folder
########################################################
open W, ">dNdS_Results.csv";
print W "filename\tdNdS\tdN\tdS\n";
foreach my $filename (glob "*.fa"){
    print "now in $filename\n";
    my %name2seq = get_hash_in_fasta($filename);
    my @seqs = values %name2seq;
    my ($dn_ds,$dn,$ds) = get_dn_ds_avg(\@seqs);
    print W "$filename\t$dn_ds\t$dn\t$ds\n";
  }
close W;

#######################################
### all subroutines are below this line
#######################################
# get average dN/dS dN dS
sub get_dn_ds_avg {
    my $array_now = shift;
    my @seqs = @$array_now;
    my $combinat_t = Math::Combinatorics->new(
        count => 2,
        data  => [ 0 .. @seqs - 1 ],
    );
    my $dn_ds_avg = 0;  my $dn_ds_count=0;
    my $ds_avg = 0; my $ds_count=0;
    my $dn_avg = 0; my $dn_count=0;
    while ( my ( $idx1, $idx2 ) = $combinat_t->next_combination ) {
        my ($pN_num, $pN_site, $pS_num, $pS_site) = get_pn_ps($seqs[$idx1], $seqs[$idx2] );
        my $ps = $pS_num /$pS_site;
	my $pn = $pN_num /$pN_site;
	my ($dS,$dN);
	if ((1 - 4/3 * $pN_num/$pN_site)>0 && (1 - 4/3 * $pS_num/$pS_site)>0){
	    $dS = jc_correct($ps);
	    $dN = jc_correct($pn);
	}
	else{ 
	    $dS="N/C";
	    $dN="N/C";
	}
	my $dN_dS;
	if ($dN eq "N/C" or $dS eq "N/C"){
	    $dN_dS="N/C";
	}elsif($dN==0 and $dS==0){
	   # $dNdS=-2;
	    $ds_count++;
	    $dn_count ++;
	}elsif($dS ==0 and $dN!=0){
	   # $dNdS=-1 ;
	    $ds_count++;
	    $dn_count ++;
	    $dn_avg+=$dN;
	}else{	    
	    $dN_dS = $dN /$dS;
	    $dn_ds_avg +=$dN_dS;
	    $dn_avg +=$dN;
	    $ds_avg +=$dS;
	    $dn_ds_count ++;
	    $ds_count++;
	    $dn_count ++;
	}
    }
    $dn_ds_avg = "NA" if $dn_ds_count == 0;
    $dn_avg = "NA" if $dn_count == 0;
    $ds_avg = "NA" if $ds_count == 0;
    $dn_ds_avg = $dn_ds_avg / $dn_ds_count if $dn_ds_count > 0;
    $dn_avg = $dn_avg / $dn_count if $dn_count > 0;
    $ds_avg = $ds_avg /$ds_count if $ds_count > 0;
    return ($dn_ds_avg,$dn_avg,$ds_avg);
}

# get synonymous number and site, non-synonymous number and site
sub get_pn_ps {
    my ($seq1, $seq2) = @_;
    my ($pN_num, $pN_site1, $pN_site2, $pS_num, $pS_site1, $pS_site2) = (0,0,0,0);
    #print ">Seq1\n", $seq1, "\n", ">Seq2\n", $seq2, "\n\n";
    
    my $align_len_seq1 = length $seq1;
    ###print "Length: $align_len_seq1\n";
    my (@base1, @base2) = ();
    foreach my $index ( 0 .. $align_len_seq1 - 1){
        $base1[$index] = substr($seq1,$index,1);
        $base2[$index] = substr($seq2,$index,1);
    }
    my $triple_switch = 1;
    my ($codon1, $codon2) = ($base1[0], $base2[0]);
    #print "Codon1: $codon1; Codon2: $codon2\n";
    foreach my $index ( 1 .. $align_len_seq1 - 1){
        my $b1 = $base1[$index];
        my $b2 = $base2[$index];
        if ($b1 ne '-'){
            if ($triple_switch < 3){
                #print $codon1, '.=', $b1, "\n"; sleep 1;
                $codon1 .= $b1;
                $codon2 .= $b2;
                $triple_switch ++;
            }
            if ($triple_switch == 3){
                if ($genetic_code{$codon1} and $genetic_code{$codon2}
                    and $genetic_code{$codon1} ne '_'
                    and $genetic_code{$codon2} ne '_'){
                    my $codon1_non = $codon2non_site{$codon1};
                    my $codon2_non = $codon2non_site{$codon2};
                    my $codon1_syn = $codon2syn_site{$codon1};
                    my $codon2_syn = $codon2syn_site{$codon2};
                    $pN_site1 += $codon1_non;
                    $pN_site2 += $codon2_non;
                    $pS_site1 += $codon1_syn;
                    $pS_site2 += $codon2_syn;

                    my $pS_num_this_codon = 0;
                    my $pN_num_this_codon = 0;
                    if (get_diff_base_number($codon1,$codon2) == 0){ # 0 base diff in codon
                    }
                    elsif (get_diff_base_number($codon1,$codon2) >= 1){
                        $pS_num += $Codon1_Codon2_syn_non_num_info{$codon1}{$codon2}->[0];          
                        $pN_num += $Codon1_Codon2_syn_non_num_info{$codon1}{$codon2}->[1];
                    }   
                }
                $codon1 = '';
                $codon2 = '';
                $triple_switch = 0;
            }
        }
    }
    my $pN_site = ($pN_site1 + $pN_site2) / 2;
    my $pS_site = ($pS_site1 + $pS_site2) / 2;
    return ($pN_num, $pN_site, $pS_num, $pS_site);
}

# get intermediate codons for codon pairs
sub get_intermediate_codon {
    my ($codon1, $codon2) = @_;
    my $diff_num = 0;
    my %hash = ();
    for my $index (0..2){
        my $base1 = substr($codon1,$index,1);
        my $base2 = substr($codon2,$index,1);
        $hash{$index} = [$base1, $base2] if $base1 ne $base2;
    }
    
    my @intermediate_codons = ();
    foreach my $ind (keys %hash){
        my $new_codon1 = $codon1;
        substr($new_codon1,$ind,1,$hash{$ind}->[1]);
        if ($new_codon1 ne 'TAA' and
            $new_codon1 ne 'TAG' and
            $new_codon1 ne 'TGA'){
            push @intermediate_codons,$new_codon1;
        }
    }
    return @intermediate_codons;
}

# get number of different bases for codon pairs
sub get_diff_base_number {
    my ($codon1, $codon2) = @_;
    my $diff_num = 0;
    for (0..2){
        $diff_num++ if substr($codon1,$_,1) ne substr($codon2,$_,1);
    }
    return $diff_num;
}

# a hash connecting codon and amino acid
sub get_genetic_code {
    my (%genetic_code) = (
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '_',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );
    return %genetic_code;
}

sub codon2aa {
    my ($codon) = @_;
    $codon = uc $codon;
    my (%genetic_code) = get_genetic_code();
    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{
        print STDERR "Bad codon \"$codon\"!!\n";
        exit;
    }
}

# calculate mutational information and GC-content
sub cal_pi_stuff {
    my ($human_CDS, $chimp_CDS) = @_;
    $human_CDS = uc $human_CDS;
    $chimp_CDS = uc $chimp_CDS;
    my  ($seq_length,           $number_of_comparable_bases,
        $number_of_identities,  $number_of_differences,
        $number_of_gaps,        $number_of_n,
        $number_of_align_error, $pi,
        $first_seq_gc,          $second_seq_gc) = @{&AlignDB::Util::pair_seq_stat(\$human_CDS,\$chimp_CDS)};
    my @indel_sites = @{&AlignDB::Util::pair_indel_sites(\$human_CDS, \$chimp_CDS)};
    my $indel_num = scalar @indel_sites;
    return [$seq_length,           $number_of_comparable_bases,
            $number_of_identities,  $number_of_differences,
            $number_of_gaps,        $number_of_n,
            $number_of_align_error, $pi,
            $first_seq_gc,          $second_seq_gc,
            $indel_num];
}

# get indel number
sub get_indel_num {
    my ($seq1, $seq2) = @_;
    my $indel_num = 0;
    while ($seq1 =~ /-+/g){
        my $len_gap = length $&;
        $indel_num++;
    }
    while ($seq2 =~ /-+/g){
        my $len_gap = length $&;
        $indel_num++;
    }
    return $indel_num;
}     

# put fasta into a hash
sub get_hash_in_fasta {
    my $filename = shift;
    open R, "$filename";
    my @lines = <R>;
    close R;
    my %name_seq = ();
    my $name     = ();
    my $switch   = 0;
    foreach my $this_line (@lines) {
        if ( $this_line =~ />(.*)/ ) {
            $name   = $1;
            $switch = 1;
        }
        elsif ( $switch == 1 ) {
            chomp $this_line;
            $name_seq{$name} = $this_line;
            $switch = 0;
        }
    }
    return %name_seq;
}


# get synonymous and non-synonymous infor for Codon1 => Codon2 
sub get_Codon1_Codon2_syn_non_num_info {
    foreach my $codon1 (grep !/TAA/, grep !/TAG/, grep !/TGA/, keys %genetic_code){   # stop codon does not count, followed the rules of MEGA4
        foreach my $codon2 (grep !/TAA/, grep !/TAG/, grep !/TGA/, keys %genetic_code){   # stop codon does not count, followed the rules of MEGA4
            if ($codon1 ne $codon2){
                if (get_diff_base_number($codon1, $codon2) == 1){  # 2 Codon Diff 1
                    if (codon2aa($codon1) eq codon2aa($codon2)){
                        $Codon1_Codon2_syn_non_num_info{$codon1}{$codon2} = [1,0];
                    }
                    else {
                        $Codon1_Codon2_syn_non_num_info{$codon1}{$codon2} = [0,1];
                    }
                }
                elsif (get_diff_base_number($codon1, $codon2) == 2){  # 2 Codon Diff 1
                    #print "$codon1 $codon2\n"; 
                    my @intermediate_codons = get_intermediate_codon($codon1,$codon2);
                    my $syn_way_number = 0;
                    my $non_way_number = 0;
                    foreach my $intermediate_codon (@intermediate_codons){
                        #next if codon2aa($intermediate_codon) eq '_';
                        if (codon2aa($codon1) eq codon2aa($intermediate_codon)){ 
                            $syn_way_number += 1;
                        }
                        elsif (codon2aa($codon1) ne codon2aa($intermediate_codon)){
                            $non_way_number += 1;
                        }
                        if (codon2aa($intermediate_codon) eq codon2aa($codon2)){ 
                            $syn_way_number += 1;
                        }
                        elsif (codon2aa($intermediate_codon) ne codon2aa($codon2)){
                            $non_way_number += 1;
                        }                        
                    }
                    my $syn_way_number_final = 0;
                    my $non_way_number_final = 0;
                    if (($syn_way_number + $non_way_number) != 0){
                        $syn_way_number_final = 2 * $syn_way_number / ($syn_way_number + $non_way_number);
                        $non_way_number_final = 2 * $non_way_number / ($syn_way_number + $non_way_number);
                    }
                    $Codon1_Codon2_syn_non_num_info{$codon1}{$codon2} = [$syn_way_number_final,$non_way_number_final];
                }
                elsif (get_diff_base_number($codon1, $codon2) == 3){  # 2 Codon Diff 1
                    my @intermediate_codons = get_intermediate_codon($codon1,$codon2);
                    my $syn_way_number = 0;
                    my $non_way_number = 0;
                    my $empty_syn_and_non_number = 0;
                    foreach my $intermediate_codon (@intermediate_codons){  #AAA - 'BBB' - CCC - DDD
                        #next if codon2aa($intermediate_codon) eq '_';                   
                        my @intermediate_codons_further = get_intermediate_codon($intermediate_codon,$codon2);
                        foreach my $intermediate_codon_further (@intermediate_codons_further){ #AAA - BBB - 'CCC' - DDD
                            #next if codon2aa($intermediate_codon_further) eq '_';
                            if (codon2aa($codon1) eq codon2aa($intermediate_codon)){ 
                                $syn_way_number += 1;
                            }
                            elsif (codon2aa($codon1) ne codon2aa($intermediate_codon)){
                                $non_way_number += 1;
                                #print "+\n";
                            }     
                            
                            #print "$codon1 ", codon2aa($codon1), " $intermediate_codon ", codon2aa($intermediate_codon), " $intermediate_codon_further ", codon2aa($intermediate_codon_further), " $codon2 ", codon2aa($codon2), "\n";
                            if (codon2aa($intermediate_codon) eq codon2aa($intermediate_codon_further)){ 
                                $syn_way_number += 1;
                            }
                            elsif (codon2aa($intermediate_codon) ne codon2aa($intermediate_codon_further)){
                                $non_way_number += 1;
                            }
                            if (codon2aa($intermediate_codon_further) eq codon2aa($codon2)){ 
                                $syn_way_number += 1;
                            }
                            elsif (codon2aa($intermediate_codon_further) ne codon2aa($codon2)){
                                $non_way_number += 1;
                            }
                        }
                    }
                    my $syn_way_number_final = 0;
                    my $non_way_number_final = 0;
                    #print "Syn: $syn_way_number\n", "Non: $non_way_number\n\n";
                    if (($syn_way_number + $non_way_number) != 0){
                        $syn_way_number_final = 3 * $syn_way_number / ($syn_way_number + $non_way_number);
                        $non_way_number_final = 3 * $non_way_number / ($syn_way_number + $non_way_number);
                    }
                    $Codon1_Codon2_syn_non_num_info{$codon1}{$codon2} = [$syn_way_number_final,$non_way_number_final];
                }
            }
        }
    }
    return %Codon1_Codon2_syn_non_num_info;
}

# get synonymous and non-synonymous site for codons
sub get_codon2syn_site_AND_codon2non_site {
    my %codon2syn_site = ();
    my %codon2non_site = ();
    %codon2syn_site = (
        'AAA' => 1/3,
        'AAC' => 1/3,
        'AAG' => 1/3,
        'AAT' => 1/3,
        'ACA' => 1,
        'ACC' => 1,
        'ACG' => 1,
        'ACT' => 1,
        'AGA' => 5/6,
        'AGC' => 1/3,
        'AGG' => 2/3,
        'AGT' => 1/3,
        'ATA' => 2/3,
        'ATC' => 2/3,
        'ATG' => 0,
        'ATT' => 2/3,
        'CAA' => 1/3,
        'CAC' => 1/3,
        'CAG' => 1/3,
        'CAT' => 1/3,
        'CCA' => 1,
        'CCC' => 1,
        'CCG' => 1,
        'CCT' => 1,
        'CGA' => 3/2,
        'CGC' => 1,
        'CGG' => 4/3,
        'CGT' => 1,
        'CTA' => 4/3,
        'CTC' => 1,
        'CTG' => 4/3,
        'CTT' => 1,
        'GAA' => 1/3,
        'GAC' => 1/3,
        'GAG' => 1/3,
        'GAT' => 1/3,
        'GCA' => 1,
        'GCC' => 1,
        'GCG' => 1,
        'GCT' => 1,
        'GGA' => 1,
        'GGC' => 1,
        'GGG' => 1,
        'GGT' => 1,
        'GTA' => 1,
        'GTC' => 1,
        'GTG' => 1,
        'GTT' => 1,
        'TAC' => 1,
        'TAT' => 1,
        'TCA' => 1,
        'TCC' => 1,
        'TCG' => 1,
        'TCT' => 1,
        'TGC' => 1/2,
        'TGG' => 0,
        'TGT' => 1/2,
        'TTA' => 2/3,
        'TTC' => 1/3,
        'TTG' => 2/3,
        'TTT' => 1/3
    );
    foreach my $codon (keys %codon2syn_site){
        $codon2non_site{$codon} = 3 - $codon2syn_site{$codon};
    }
    return (\%codon2syn_site, \%codon2non_site);
}

# jc correct pi value
sub jc_correct {
    my ($pi) = @_;
    $pi = -0.75 * log(1 - (4.0 / 3.0) * $pi);
    return $pi;
}
