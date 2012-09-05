#Copyright 2009, 2010 Daniel Gaston, andrew Roger Lab
#This code is copyrighted under the GNU General Public License Version 3.0
#
#       This program is free software: you can redistribute it and/or modify
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


package QmmRAxML;

use 5.010000;
use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use QmmRAxML ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw() ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw(getSiteCombLH getSiteClassPost callQmmRAxML calcFreqClassData getEmpFreqs);

our $VERSION = '0.1';


# Preloaded methods go here.
sub getSiteCombLH{
    my $file = shift;
    print "Retrieveing site likelihoods from $file...\n";
    open(IN, "$file") or die "Could not open $file...\n";
    my @raw = <IN>;
    close(IN);
    
    $raw[0] =~ s/[\r\n]+//g;
    
    my @lhs = split / +/, $raw[0];
    
    return @lhs;
}

sub getEmpFreqs{
    my $file = shift;
    open(IN, "$file") or die "Could not open $file.\n\n";
    my @raw = <IN>;
    close(IN);
    
    my $vector = "";
    foreach(@raw){
        if($_ =~ /pi\([A-Z]\):/){
            $vector = "$vector" . "$_";   
        }
    }
    $vector =~ s/pi\([A-Z]\)://g;
    $vector =~ s/[\r\n]+//g;
    my @vector = split / +/, $vector;
    shift @vector;
 
    return @vector;
}

#Retrieve the Site Class Posterior probabilities for components for all treefiles.
#Order: Wholetree then subtrees
sub getSiteClassPost{
    my ($file, $subref, $num_classes) = @_;
    my @files = @{$subref};
    unshift @files, $file;
    print "Retrieving site class posterior probabilities from $file...\n";
    open(IN, "$file") or die "Could not open $file...\n";
    my @posts;
    my $file_count = 0;
    foreach my $file(@files){
        open(IN, "$file") or die "Could not open $file...\n";
        my $count = 0;
        while (<IN>){
            my $line  = $_;
            $line =~ s/[\r\n]+//g;
            my @temp = split / +/, $line;
            if($temp[0] =~ /^-?\d+$/){
                shift @temp;
                if(scalar(@temp) != $num_classes){
                    print "ERROR: Site class posterior probabilities have incorrect number of entries.\n\n";
                    exit();
                }
                for(my $i = 0; $i < scalar(@temp); $i++){
                    my %temp;
                    $temp{class} = ($i + 1);
                    $temp{post}  = $temp[$i];
                    $temp[$i] = \%temp;
                }
                $posts[$count][$file_count] = [@temp];
                $count++;
            }
        }
        $file_count++;
    }
    
    return @posts;
}

#Calculate the various maesures of FD based on Site Frequency Class Preferences
#Tree Order: Wholetree then subtrees
sub calcFreqClassData{
    my ($model, $ref1, $ref2) = @_;
    my @class_posteriors = @{$ref1};
    my @raxml_info_files = @{$ref2};
    
    my $class_file = "9componentsDirichlet_$model.dat";
    my $sep = $/;
    $/ = undef;
    open(CIN, "$class_file") or die "Could not open $class_file.\n\n";
    my $raw = <CIN>;
    close(CIN);
    $/ = $sep;
    
    #Components of the Dirichlet mixture model (Not including the base model)
    my @temp_components = split />/, $raw;
    shift @temp_components;
    my @components;
    foreach(@temp_components){
        my @temp = split /[\r\n]+/, $_;
        my @vector = split / /, $temp[-1];
        push @components, [@vector];
    }
    
    #my $n = scalar(@components) + 1;
    #my $mean = 1.0 / $n;
    
    #Order Wholetree then subtrees. Standard model (JTT/WAG) then dirichlet components
    my @freqs;
    foreach(@raxml_info_files){
        my @freq_vector = getEmpFreqs($_);
        my @temp;
        push @temp, [@freq_vector];
        push @temp, @components;
        push @freqs, [@temp];
    }
    
    open(CLASS, ">Frequency_Class_Eval.txt") or die "Could not create Frequency_Class_Eval.txt\n\n";
    open(CLOG, ">Class.log") or die "Could not create Class.log\n\n";
    print CLASS "Site\tWholetree Class\tSubtree 1 Class\tSubtree 2 Class\tFlag\t95% Set\tSubtree Distance\tWhole-Subtree1 Dist\tWhole-Subtree 2 Dist\tPost Sub Dist\t Post Whole-Sub 1 Dist\tPost Whole-Sub 2 Dist\n";
    
    my @sorted_class_posteriors;
    my @class_freq_sets;
    my @output;
    #loop over sites
    for(my $i = 0; $i < scalar(@class_posteriors); $i++){
        my $site = $i + 1;
        push @output, "$site\t";
        my $flag = 0;
        my $wtc;
        my @previous;
        my @sets;
        my @comp_vectors;
        #Loop over trees (Wholetree then subtrees)
        for(my $y = 0; $y < scalar(@{$class_posteriors[$i]}); $y++){
            my @temp_sorted = sort { $b -> {post} <=> $a -> {post} } @{$class_posteriors[$i][$y]};
            $sorted_class_posteriors[$i][$y] = [@temp_sorted];
            push @output, "$temp_sorted[0]{class} ($temp_sorted[0]{post})\t";
            my $top = $temp_sorted[0];
            my $class_id = $top -> {class} - 1;
            
            #Calculate the 95% set of class frequency components as well as the mean posterior
            my $total = 0;
            my %set;
            $set{top} = $top -> {class};
            $set{top_post} = $top -> {post};
            TEMP: foreach(@temp_sorted){
                $total = $total + $_ -> {post};
                $set{$_ -> {class}} = 1;
                if($total >= 0.9500000000000000000000000){
                    $set{total} = $total;
                    last TEMP;
                }
            }
                
            #$set{sfactor} = $set{top_post} - $mean;
            $set{sfactor} = $set{top_post};
            
            my @weighted_vector;
            foreach my $v (@{$freqs[$y][$class_id]}){
                $v = $v * $set{sfactor};
                push @weighted_vector, $v;
            }
            push @sets, \%set;
            push @comp_vectors, [@weighted_vector];
            
            #Check for switch of class frequency component of highest posterior
            #When $y == 0 we are looking at the wholetree data file    
            if($y == 0){
                $wtc = $set{top};
            }else{
                if($set{top} != $wtc){
                    if($previous[0]){
                        CLASSES: foreach(@previous){
                            if($set{top} != $_){
                                $flag++;
                                last(CLASSES);
                            }
                        }
                    }else{
                        $flag++;
                    }
                }
                push @previous, $set{top};
            }
        }
        push @output, "$flag\t";
        
        #currently this is coded to only work with two subtrees
        #Fix this at some point, perhaps by doing an average of the
        #Euclidean distances?
        if(scalar(@comp_vectors) > 3){
            print "Currently FunDi is only designed to operate with two subtrees. Exiting...\n\n";
            exit();
        }
        
        #######################Check if Top Component in 95% sets##################################
        if(exists $sets[1]{$sets[0]{top}}){
            push @output, "0\t";
            print CLOG "$site\t $sets[0]{top} of subtree 1 exists in 95% set of subtree 2 ( ";
            while( my($k,$v) = each %{$sets[1]} ){
                if($k ne 'total' && $k ne 'top'){
                    print CLOG "$k ";   
                }
            }
            print CLOG ")\n";
        }elsif(exists $sets[0]{$sets[1]{top}}){
            push @output, "0\t";
            print CLOG "$site\t $sets[1]{top} of subtree 2 exists in 95% set of subtree 1 ( ";
            while( my($k,$v) = each %{$sets[0]} ){
                if($k ne 'total' && $k ne 'top'){
                    print CLOG "$k ";   
                }
            }
            print CLOG ")\n";
        }else{
            push @output, "1\t";
            print CLOG "$site\t Neither $sets[0]{top} of subtree 1 or $sets[1]{top} " .
                        "of subtree 2 exist in either set ( ";
            while( my($k,$v) = each %{$sets[0]} ){
                if($k ne 'total' && $k ne 'top'){
                    print CLOG "$k ";   
                }
            }
            print CLOG ") ( exit";
            
            while( my($k,$v) = each %{$sets[1]} ){
                if($k ne 'total' && $k ne 'top'){
                    print CLOG "$k ";   
                }
            }
            print CLOG ")\n";
        }
        ###########################################################################
        
        #Calculate the euclidean distance between frequency vectors
        my @totals;
        my $total_subs = 0;
        my $total1 = 0;
        my $total2 = 0;
        if(scalar(@{$comp_vectors[0]}) != scalar(@{$comp_vectors[1]}) && scalar(@{$comp_vectors[1]}) != scalar(@{$comp_vectors[2]})){
            print "Components vectors not equal length. Exiting...\n\n";
            exit();
        }
        for(my $x = 0; $x < scalar(@{$comp_vectors[0]}); $x++){
            my $diff_subs = $comp_vectors[1][$x] - $comp_vectors[2][$x];
            $diff_subs = $diff_subs * $diff_subs;
            $total_subs = $total_subs + $diff_subs;
            
            my $diff1 = $comp_vectors[0][$x] - $comp_vectors[1][$x];
            $diff1 = $diff1 * $diff1;
            $total1 = $total1 + $diff1;
            
            my $diff2 = $comp_vectors[0][$x] - $comp_vectors[2][$x];
            $diff2 = $diff2 * $diff2;
            $total2 = $total2 + $diff2;
        }
        $total_subs = sqrt($total_subs);
        $total1 = sqrt($total1);
        $total2 = sqrt($total2);
        
        push @output, "$total_subs\t$total1\t$total2\t";
        
        #Calculate euclidean difference between class posterior vectors
        $total_subs = 0;
        $total1 = 0;
        $total2 = 0;
        for(my $j = 0; $j < scalar(@{$class_posteriors[$i][0]}); $j++){
            my $diff_subs = $class_posteriors[$i][0][$j]{post} - $class_posteriors[$i][1][$j]{post};
            $diff_subs = $diff_subs * $diff_subs;
            $total_subs = $total_subs + $diff_subs;
            
            my $diff1 = $class_posteriors[$i][-1][$j]{post} - $class_posteriors[$i][0][$j]{post};;
            $diff1 = $diff1 * $diff1;
            $total1 = $total1 + $diff1;
            
            my $diff2 = $class_posteriors[$i][-1][$j]{post} - $class_posteriors[$i][1][$j]{post};;
            $diff2 = $diff2 * $diff2;
            $total2 = $total2 + $diff2;
        }
        $total_subs = sqrt($total_subs);
        $total1 = sqrt($total1);
        $total2 = sqrt($total2);
        
        push @output, "$total_subs\t$total1\t$total2\n";
    }
    #End loop through sites
    
    foreach my $outline (@output){
        print CLASS $outline;    
    }
    close(CLASS);
    close(CLOG);   
}

sub callQmmRAxML{
    my $string = shift;
    `qmmraxmlHPC $string`;
}

1;
__END__

=head1 NAME

QmmRAxML - Perl extension for calling qmmRAxML and parsing output files

=head1 SYNOPSIS

  use QmmRAxML;
  blah blah blah

=head1 DESCRIPTION

This set of functions are useful for calling the phylogenetic program qmmRAxML as
well as parsing the output files.

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
for the Andrew Roger Lab, Dalhousie University
Permanent Contact: Andrew.Roger@dal.ca

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009 by Daniel Gaston

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.


=cut
