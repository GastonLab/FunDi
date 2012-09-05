#Copyright 2009, 2010 Daniel Gaston, andrew Roger Lab
#This code is copyrighted under the GNU General Public License Version 3.0
#
#	This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


package Extras;

use 5.010000;
use strict;
use warnings;
use QmmRAxML;
use RAxML;
use Seq;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use Extras ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw() ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw(formatDate numeric_sort getRates compareRates siteLikelihoodCompare siteRateCompare maskAlignment);

our $VERSION = '0.1';


# Preloaded methods go here.
sub maskAlignment{
    my ($inalign, $cutoff) = @_;
    my %seqs = Seq::readPhylip($inalign);
    my @seqs = Seq::convertToArray(\%seqs);
    
    my @masked_sites;
    my %out;
    $inalign = "masked.phy";
    my $length = scalar(@seqs);
    
    foreach(keys %seqs){
        $out{$_}{seq} = "";
    }
    
    for(my $i = 2; $i < scalar(@{$seqs[0]}); $i++){
        my $num_gaps = 0;
        for(my $j = 0; $j < scalar(@seqs); $j++){
            if($seqs[$j][$i] eq '-'){
                $num_gaps++;
            }
        }
        my $frac = $num_gaps / $length;
        if($frac > $cutoff){
            push @masked_sites, ($i - 2);
        }else{
            for(my $j = 0; $j < scalar(@seqs); $j++){
                $out{$seqs[$j][0]}{seq} = $out{$seqs[$j][0]}{seq} . "$seqs[$j][$i]";
                $out{$seqs[$j][0]}{desc} = $seqs[$j][1];
            }
        }
    }
    
    Seq::printPhylip($inalign, \%out);
    return ($inalign, \@masked_sites);
}

sub formatDate{
    my @time = localtime(time());
    
    my $mday = $time[3];
    my $year = $time[5] + 1900;
    my $month = $time[4] + 1;
    
    if(length $mday == 1){
        $mday = "0" . "$mday";
    }
    
    if(length $month == 1){
        $month = "0" . "$month";
    }
    
    my $datestring = "$year" . "$month" . "$mday";
    
    return $datestring;
}

sub getRates{
    my $file = shift;
    my @rates;
    
    open(IN, "$file"), or die "Could not read $file!\n";
    my @data = <IN>;
    close(IN);
    
    @rates = split " ", $data[1];
    shift @rates;
    
    return @rates;
}

sub compareRates{
    my ($rate1, $rate2, $norm) = @_;
    my $rate_diff;
    
    if($norm eq "sum_divide"){
        $rate_diff = ($rate1 - $rate2) / ($rate1 + $rate2);   
    }elsif($norm eq "log_ratio"){
        $rate_diff = log($rate1) / log($rate2);
    }elsif($norm eq "mean_divide"){
        
    }else{
        $rate_diff = $rate1 - $rate2;   
    }
    
    return $rate_diff;
}

sub siteRateCompare{
    my ($data_ref, $sim_ref, $norm) = @_;
    
    my $output = "siterate_comparison.out";
    
    my @data_all_files = @{$data_ref};
    my @sim_all_files = @{$sim_ref};
    
    my @data_files;
    my @sim_files;
    
    foreach(@data_all_files){
        if($_ =~ /subtree/){
            push @data_files, $_;
        }
    }
    
    foreach(@sim_all_files){
        if($_ =~ /subtree/){
            push @sim_files, $_;
        }
    }
    
    if(scalar(@data_files) != scalar(@sim_files)){
        print "Subtree file lists between data and simulated data are not equal!\n";
        exit();
    }
    
    open(OUT, ">$output") or die "Could not create $output!\n";
    
    for(my $i = 0; $i < scalar(@data_files); $i++){
        my @data_rates_tree0 = getRates($data_files[$i]);
        my @sim_rates_tree0 = getRates($sim_files[$i]);
        
        for(my $y = $i + 1; $y < scalar(@data_files); $y++){
    
            my @data_rates_tree1 = getRates($data_files[$y]);
            my @sim_rates_tree1 = getRates($sim_files[$y]);
                
            if(scalar(@data_rates_tree0) != scalar(@data_rates_tree1)){
                print "Number of sites in $data_files[$i] and $data_files[$y] unequal!\n";
                exit();
            }elsif(scalar(@sim_rates_tree0) != scalar(@sim_rates_tree1)){
                print "Number of sites in $sim_files[$i] and $sim_files[$y] unequal!\n";
                exit();
            }
                
            my @sim_rate_diffs;
            for(my $w = 0; $w < scalar(@sim_rates_tree0); $w++){
                my $sim_rate_diff = compareRates($sim_rates_tree0[$w], $sim_rates_tree1[$w], $norm);
                push @sim_rate_diffs, $sim_rate_diff;
            }
                
            my @sorted_sim_rate_diffs = sort numeric_sort @sim_rate_diffs;
                
            my $tail_005 = scalar(@sorted_sim_rate_diffs) * 0.05;
            my $tail_010 = scalar(@sorted_sim_rate_diffs) * 0.10;
                
            my $lower_005_thresh = $sorted_sim_rate_diffs[$tail_005];
            my $lower_010_thresh = $sorted_sim_rate_diffs[$tail_010];
            my $upper_005_thresh = $sorted_sim_rate_diffs[(scalar(@sorted_sim_rate_diffs) - $tail_005)];
            my $upper_010_thresh = $sorted_sim_rate_diffs[(scalar(@sorted_sim_rate_diffs) - $tail_010)];
                
            my @tail_005;
            my @tail_010;
                
            for(my $x = 0; $x < scalar(@data_rates_tree0); $x++){
                my $data_rate_diff = compareRates($data_rates_tree0[$x], $data_rates_tree1[$x], $norm);
                my $site_info;
                    
                if($data_rate_diff < $lower_005_thresh || $data_rate_diff > $upper_005_thresh){
                    $site_info = $x + 1 . "\t" . "$data_rate_diff";
                    push @tail_005, $site_info;   
                }
                    
                if($data_rate_diff < $lower_010_thresh || $data_rate_diff > $upper_010_thresh){
                    $site_info = $x + 1 . "\t" . "$data_rate_diff";
                    push @tail_010, $site_info;
                }
            }
                
            print OUT "0.05 lower and upper bounds are: $lower_005_thresh and $upper_005_thresh\n";
            print OUT "#" x 50;
            print OUT "\n";
            foreach(@tail_005){
                print OUT "$_\n";   
            }
                
            print OUT "0.10 lower and upper bounds are: $lower_010_thresh and $upper_010_thresh\n";
            print OUT "#" x 50;
            print OUT "\n";
            foreach(@tail_010){
                print OUT "$_\n";   
            } 
        }
    }


    close(OUT);
}

sub siteLikelihoodCompare{
    my @files = @_;
    #seperate the files for the log likelihoods of the subtrees into an array
    my @subtrees;
    my $whole;
    foreach(@files){
        if($_ =~ /subtree/){
            push @subtrees, $_;
        }else{
            $whole = $_;   
        }
    }
    
    #get the site likelihoods of the whole tree
    open(WHOLE, "$whole") or die "Could not open $whole for reading...\n";
    my @whole;
    while(<WHOLE>){
    	push @whole, $_;
    }
    close(WHOLE);
    my @whole_likelihoods = split " ", $whole[1];
    #remove the first element which is an identifier
    shift @whole_likelihoods;
    
    #Get the proper site likelihood lines from each of the subtree site likelihood files
    #and construct a matrix of site likelihoods
    my @sub_likelihoods;
    my @likelihood_matrix;
    my $count = 0;
    my $innerarray_length = 0;
    foreach(@subtrees){
	    my @temp;
	    open(SUB, "$_") or die "Could not open $_ for reading...\n";
	    while(my $temp_line = <SUB>){
	    	push @temp, $temp_line;
	    }
	    close(SUB);
	    my @temp_lh = split " ", $temp[1];
	    shift @temp_lh;
	    my $temp_len = scalar(@temp_lh);
	    if($temp_len > $innerarray_length){
	    	$innerarray_length = $temp_len;
	    	$count++;
	    }
	    push @likelihood_matrix, [@temp_lh];
    }
	
    if($count != 1){
    	print "Not all arrays are the same length or they are empty!\n";
    }
	
    #loop through the site likelihood matrix constructed from the subtrees and add
    #the log likelihoods for each subtree together on a sitewise basis. This will construct
    #an array that consists of the sum of the log likelihoods for each site across the subtrees
    #which will then be compared to the site likelihoods for the whole tree
    #taken as a unit.
	
    my @summed_likelihoods;
	
    for(my $i = 0; $i < $innerarray_length; $i++){
	    my $sum = 0;
	    for(my $x = 0; $x < scalar(@likelihood_matrix); $x++){
	    	$sum = $sum + $likelihood_matrix[$x][$i];
	    }
	    push @summed_likelihoods, $sum;
    }
	
    #Compare the summed likelihoods with the likelihoods from the whole tree
    my @delta_lh;
    if(scalar(@whole_likelihoods) != scalar(@summed_likelihoods)){
    	print "Whole likelihoods and Summed likelihoods arrays are not the same size!\n";
    }

    for(my $i = 0; $i < scalar(@whole_likelihoods); $i++){
    	my $diff = $whole_likelihoods[$i] - $summed_likelihoods[$i];
    	push @delta_lh, $diff;
    }
	
    return @delta_lh;
}

sub numeric_sort { $a <=> $b}

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Extras - Perl functions for use with FunDi

=head1 SYNOPSIS

  use Extras;

=head1 DESCRIPTION

This Extras Module stores functions used in FunDi

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Daniel Gaston (daniel.gaston@dal.ca)
For Andrew Roger Lab, Dalhousie University
Permanent Contact: Andrew.Roger@dal.ca

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009 by Daniel Gaston

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.


=cut
