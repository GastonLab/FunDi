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


package Seq;
use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(readFasta
                    readPhylip
                    readPhylipN
                    readStockholm
                    printFasta
                    printFastaNoDesc
                    printPhylip
                    printPhylipN
                    findIdenticalSeqs
                    convertToArray
                    calcSiteConservation);

sub _getData{
    my $file = shift;
    open(IN, "$file") or die "Could not open $file...\n";
    my @data = <IN>;
    close(IN);
    
    return @data;
}

#Sequence Format Parsers
sub readFasta{
    my $file = shift;
    my %fasta_seqs;
    
    my $sep = $/;
    $/ = undef;
    open(IN, "$file") or die "Could not open $file...\n";
    my $raw = <IN>;
    $/ = $sep;
    close(IN);
    
    my @data = split />/, $raw;
    shift @data;
    my $num_sequences = scalar(@data);
    foreach(@data){
        if($_ =~ /(.+?)[\r\n]+/){
            my $info = $1;
            my $id;
            my $desc;
            if($info =~ /^(gi\|.+ *)/){
                $id = $1;
                $desc = $info;
                $desc =~ s/$id//;
            }elsif($info =~ /(.+?)\|(.+)/){
                $id = $1;
                $desc = $2;
            }elsif($info =~ /(.+?) (.+)/){
                $id = $1;
                $desc = $2;
            }else{
                $id = $info;
                $desc = $info;
            }
            my @temp = split /[\n\r]+/, $_;
            shift @temp;
            my $seq = "";
            foreach my $section (@temp){
                $seq = "$seq" . "$section";
            }
            $seq =~ s/ //g;
            $id =~ s/[\r\n]//g;
            $id =~ s/ //g;
            $desc =~ s/[\r\n]//g;
            
            $fasta_seqs{$id}{desc} = $desc;
            $fasta_seqs{$id}{seq}  = $seq;
        }else{
            print "Having trouble parsing file, check format and try again.\n";
            exit();
        }
    }
    
    return %fasta_seqs;
}

#Reads both interleaved and non-interleaved Phylip format
sub readPhylip{
    my $file = shift;
    
    my $verbose;
    $_[0] ? $verbose = $_[0] : $verbose = 0;
    
    my @raw_data = _getData($file);
    my %phylip_seqs;
    
    print "Parsing sequence file $file...\n";
    
    my $num_seqs;
    my $len;
    my $info = shift @raw_data;
    if($info =~ /^ ([0-9]+) ([0-9]+)[\r\n]+$/){
        $num_seqs = $1;
        $len = $2;
    }else{
        print "Could not parse header line of phylip file, number of sequences not found.\n";
        print "Please check the format of your input phylip file and try again. Exiting...\n";
        exit();
    }
    for(my $i = 0; $i < $num_seqs; $i++){
        my $id;
        my $seq;
        
        if($raw_data[$i] =~ /^[\r\n]+/){
            next();
        }
        
        if($raw_data[$i] =~ /\s/){
            if($i == 0){
                if($verbose){
                    print "Whitespace detected in line. Possibly relaxed phylip format. Splitting on whitespace\n";
                }
            }
            my @temp = split /\s+/, $raw_data[$i];
            $id = $temp[0];
            $seq = $temp[1];
        }else{
            if($verbose){
                print "No whitespace detected in line. First 10 characters being used as id.\n";
            }
            $id = substr($raw_data[$i], 0, 10);
            $seq = substr($raw_data[$i], 10);
        }
        $id =~ s/ //g;
        $seq =~ s/ //g;
        $seq =~ s/[\n\r]+//g;
        if(length($seq) != $len){
            print "ERROR: Error in sequence $id, sequence did not have expected number of sites ($len)\n";
            print "Exiting...\n";
            exit();
        }
        if($phylip_seqs{$id}){
            $phylip_seqs{$id}{seq}  = $phylip_seqs{$id}{seq} . $seq;
        }else{
            $phylip_seqs{$id}{desc} = $id;
            $phylip_seqs{$id}{seq}  = $seq;
        }
    }
    return %phylip_seqs;
}

#Parser for Stockholm sequence format. Includes support for reading FSA support output line
#and incorporates it directly in to the sequence.
#Reads both interleaved and non-interleaved Phylip format
sub readStockholm{
    my $file = shift;
    my @raw_data = _getData($file);
    my %alignment;
    
    my $count = 0;
    foreach my $line (@raw_data){
        $line =~ s/[\r\n]+//;
        if($line =~ /^\#/){
            #Will read both FSA and HMMER3 support values
            if($line =~ /^\#=GC Accuracy\s+(\d+)/ || $line =~ /^\#=GC PP_cons\s+(.+)/) {
                my $seq = $1;
                if($alignment{"Confidence"}){
                    $alignment{"Confidence"}{seq} = $alignment{"Confidence"}{seq} . $seq;
                }else{
                    $alignment{"Confidence"}{desc} = "Conf";
                    $alignment{"Confidence"}{seq} = $seq;
                }
            #read HMMER3 individual sequence support values
            }elsif($line =~ /^\#=GR\s+(.+)\s+PP (.+)/){
                my $id = $1;
                my $seq = $2;
                $id =~ s/ +$//;
                $seq =~ s/[\r\n]+//g;
                $id = $id . "_conf";
                if($alignment{$id}){
                    $alignment{$id}{seq} = $alignment{$id}{seq} . $seq;
                }else{
                    $alignment{$id}{desc} = "Conf";
                    $alignment{$id}{seq} = $seq;
                }
            }
        }elsif($line =~ /^\/\/\s*$/){
            #end of alignment
            print "Finished parsing sequences in file $file\n";
            last();
        }elsif($line =~ /^[A-Z0-9]/i){
            if($line =~ /^([A-Z0-9]+?)\/(.+)\s+([A-Z\-\.]+)/i){
                my $id = $1;
                my $desc = $2;
                my $seq = $3;
                $id =~ s/\s//g;
                if($alignment{$id}){
                    $alignment{$id}{seq} = $alignment{$id}{seq} . $seq;
                }else{
                    $alignment{$id}{desc} = $desc;
                    $alignment{$id}{seq} = $seq;
                }
            }elsif($line =~ /^(.+)\s+([A-Z\-\.]+)/i){
                my $id = $1;
                my $seq = $2;
                my $desc = $1;
                $id =~ s/\s//g;
                if($alignment{$id}){
                    $alignment{$id}{seq} = $alignment{$id}{seq} . $seq;
                }else{
                    $alignment{$id}{desc} = $desc;
                    $alignment{$id}{seq} = $seq;
                }
            }else{
                print "Error in parsing putative sequence line $count: $line\n";
            }
        }elsif($line eq ''){
            #print "Break found. Sequence file probably in interleaved stockholm format\n";
        }else{
            print "Error in parsing line $count: $line\n";    
        }
        $count++;
    }
    
    return %alignment;
}

sub readPhylipN{
    my $file = shift;
    my @raw_data = _getData($file);
    
    print "Phylip Interleaved format not currently supported\n";
    exit();
}

#Sequence printing formaters
sub printFasta{
    my ($file, $seq_ref) = @_;
    my %seq = %{$seq_ref};
    
    my @guide;
    if($seq{'Guide'}){
        @guide = @{$seq{'Guide'}};   
    }else{
        @guide = sort(keys %seq);
    }
    
    my @unprinted;
    open(OUT, ">$file") or die "Could not create $file for writing...\n";
    foreach my $id (@guide){
        if($seq{$id}){
            if($seq{$id}{desc}){
                print OUT ">$id | $seq{$id}{desc}\n";
            }else{
                print OUT ">$id\n";
            }
            my $sequence = $seq{$id}{seq};
            $sequence =~ s/(.{60})/$1\n/g;
            print OUT "$sequence\n";
        }else{
            push @unprinted, $id;
        }
    }
    close(OUT);
    
    if(scalar(@unprinted) > 0){
        print "Guide elements not found in alignment:\n";
        foreach(@unprinted){
            print "$_\n";
        }
    }
}

sub printFastaNoDesc{
    my ($file, $seq_ref) = @_;
    my %seq = %{$seq_ref};
    
    my @guide;
    if($seq{'Guide'}){
        @guide = @{$seq{'Guide'}};   
    }else{
        @guide = sort(keys %seq);
    }
    
    my @unprinted;
    open(OUT, ">$file") or die "Could not create $file for writing...\n";
    foreach my $id (@guide){
        if($seq{$id}){
            print OUT ">$id\n";
            my $sequence = $seq{$id}{seq};
            $sequence =~ s/(.{60})/$1\n/g;
            print OUT "$sequence\n";
        }else{
            push @unprinted, $id;
        }
    }
    close(OUT);
    
    if(scalar(@unprinted) > 0){
        print "Guide elements not found in alignment:\n";
        foreach(@unprinted){
            print "$_\n";
        }
    }
}

sub printPhylip{
    my ($file, $seq_ref) = @_;
    my %seq = %{$seq_ref};
    
    my @guide;
    if($seq{'Guide'}){
        @guide = @{$seq{'Guide'}};   
    }else{
        @guide = sort(keys %seq);
    }
    
    my @unprinted;
    my @temp = keys %seq;
    my $num = scalar(@temp);
    my $length = length($seq{$temp[0]}{seq});
    open(OUT, ">$file") or die "Could not create $file for writing...\n";
    print OUT " $num $length\n";
    foreach my $id (@guide){
        if($seq{$id}){
            print OUT "$id $seq{$id}{seq}\n";
        }else{
            push @unprinted, $id;
        }
    }
    close(OUT);
    
    if(scalar(@unprinted) > 0){
        print "Guide elemnts not found in alignment:\n";
        foreach(@unprinted){
            print "$_\n";
        }
    }
}

sub printStockholm{
    my ($file, $seq_ref) = @_;
    my %seq = %{$seq_ref};
    
    my @guide;
    if($seq{'Guide'}){
        @guide = @{$seq{'Guide'}};   
    }else{
        @guide = sort(keys %seq);
    }
    
    my @unprinted;
    open(OUT, ">$file") or die "Could not create $file for writing...\n";
    print OUT "# STOCKHOLM 1.0\n";
    foreach my $id (@guide){
        if($seq{$id}){
            print OUT "$id\t$seq{$id}{seq}\n";
        }else{
            push @unprinted, $id;
        }
    }
    if($seq{'Confidence'}){
        print OUT "#=GC PP_cons\t$seq{'Confidence'}{seq}\n";
    }
    if(scalar(@unprinted) > 0){
        print "Guide elements not found in alignment:\n";
        foreach(@unprinted){
            print "$_\n";
        }
    }
    print OUT "//\n";
    close(OUT);
}

sub printGDE{
    my ($file, $seq_ref) = @_;
    my %seq = %{$seq_ref};
    
    my @guide;
    if($seq{'Guide'}){
        @guide = @{$seq{'Guide'}};   
    }else{
        @guide = sort(keys %seq);
    }
    my @unprinted;
    open(OUT, ">$file") or die "Could not create $file for writing...\n";
    foreach my $id (@guide){
        if($seq{$id}){
            print OUT "%$id $seq{$id}{desc}\n$seq{$id}{seq}\n";
        }else{
            push @unprinted, $id;
        }
    }
    if(scalar(@unprinted) > 0){
        print "Guide elements not found in alignment:\n";
        foreach(@unprinted){
            print "$_\n";
        }
    }
    close(OUT);
}

#sub printPhylipN{
#    my ($file, $seq_ref) = @_;
#    my %seq = %{$seq_ref};
#    
#    print "Phylip Interleaved format not supported. Exiting...\n\n";
#    
#    open(OUT, ">$file") or die "Could not create $file for writing...\n";
#    foreach(sort(keys %seq)){
#        
#    }
#    close(OUT);
#}

#Accessory Methods
sub removeConf{
    my $ref = shift;
    my %seqs = %{$ref};
    
    delete($seqs{"Confidence"});
    
    return %seqs;
}

sub removeCons{
    my $ref = shift;
    my %seqs = %{$ref};
    
    delete($seqs{"Conservation"});
    
    return %seqs;
}

sub readGuide{
    my ($ref, $file) = @_;
    my %seqs = %{$ref};
    
    my %copy = %seqs;
    delete($copy{"Confidence"});
    my @align_taxa = keys (%copy);
    
    open(IN, "$file") or die "Could not open $file for reading...\n\n";
    my @data;
    while(<IN>){
        if($_ !~ /^[\r\n]+$/){
            my $line = $_;
            $line =~ s/[\r\n]+//;
            push @data, $line;
        }
    }
    close(IN);
    
    $seqs{'Guide'} = \@data;
    
    my (@union, @isect, @diff) = ();
    my (%union, %isect) = ();
    my %count = ();
    
    foreach my $e (@align_taxa, @data) { $count{$e}++ }

    foreach my $e (keys %count){
        push(@union, $e);
        if ($count{$e} == 2) {
            push @isect, $e;
        }else{
            push @diff, $e;
        }
    }
    
    if(scalar(@diff) > 0){
        if(scalar(@diff) == 1 && $diff[0] ne 'Confidence'){
            print "Taxa found in alignment or guide but not the other:\n";
            foreach(@diff){
                print "$_\n";
            }
            print "Check guide and alignment and try again...\n\n";
            exit();
        }
    }
    
    return %seqs;
}

sub readGuideTree{
    
}

sub cleanSeq{
    my $ref = shift;
    my %seqs = %{$ref};
    
    foreach(keys %seqs){
        if($_ ne 'Confidence' && $_ ne 'Guide'){
            my $seq = $seqs{$_}{seq};
            $seq =~ s/\./\-/g;
            $seq = uc($seq);
            $seqs{$_}{seq} = $seq;
        }
    }
    
    return %seqs;
}

#Adds several X's as boundary markers at the beginning and end of a sequence.
sub addXMarkers{
    my $ref = shift;
    my %seqs = %{$ref};
    
    foreach(keys %seqs){
        my $seq = $seqs{$_}{seq};
        $seqs{$_}{seq} = "XXX" . $seq . "XXX";
    }
    
    return %seqs;
}

sub trimMarkers{
    my $ref = shift;
    my %seqs = %{$ref};
    
    foreach(keys %seqs){
        my $seq = $seqs{$_}{seq};
        if($seq =~ /^X{3}(.+)X{3}$/){
            $seq = $1;
        }else{
            print "Cannot trim marker ends of $_\n";
        }
        $seqs{$_}{seq} = $seq;
    }
    
    return %seqs;
}

sub findIdenticalSeqs{
    my $file = shift;
    my %seqs = readPhylip($file);
    my @identical;
    
    my @sorted = sort(keys(%seqs));
    for(my $i = 0; $i < scalar(@sorted); $i++){
        my $seq1 = $seqs{$sorted[$i]}{seq};
        for(my $x = $i + 1; $x < scalar(@sorted); $x++){
            my $seq2 = $seqs{$sorted[$x]}{seq};
            if($seq1 eq $seq2){
                push @identical, "$sorted[$i] $sorted[$x]";
            }
        }
    }
    
    if(scalar(@identical) > 0){
        my $ok = 0;
        while(!$ok){
            print "The following sets of identical sequences were found:\n";
            foreach(@identical){
                print "$_\n";   
            }
            print "Do you wish to continue?(y/n): ";
            my $response = <STDIN>;
            chomp $response;
            if($response eq 'N' || $response eq 'n'){
                print "Exiting...\n";
                $ok = 1;
                exit();
            }elsif($response eq 'Y' || $response eq 'y'){
                print "continuing...\n";
                $ok = 1;
            }else{
                print "Invalid selection, try again...\n";   
            }
        }
    }
}

sub convertToArray{
    my $ref = shift;
    my %seqs = %{$ref};
    my @seqs_array;
    foreach(sort (keys (%seqs))){
        my $seq = $seqs{$_}{seq};
        my @temp_seq = split //, $seq;
        unshift @temp_seq, $seqs{$_}{desc};
        unshift @temp_seq, $_;
        push @seqs_array, [@temp_seq];
    }
    
    #Check lengths
    my $length = scalar(@{$seqs_array[0]});
    my $i = 1;
    foreach(@seqs_array){
        if(scalar(@{$_}) != $length){
            print "Array of unequal length found at entry $i. Exiting.\n\n";
            exit();
        }
        $i++;
    }
    
    return @seqs_array;
}

#Needs a two dimensional sequence array
sub calcSiteConservation{
    my ($ref, $cutoff) = @_;
    my @seqs = @{$ref};
    my @conservation = ("Conservation", "Score");
    my $length =  scalar(@{$seqs[0]});
    my $num_col = scalar(@seqs);
    my @conserved;
    
    for(my $i = 2; $i < $length; $i++){
        $conservation[$i] = 0.0000000000;   
    }
    
    for(my $i = 2; $i < $length; $i++){
        my %aas = ('A' => 0.0, 'R' => 0.0, 'N' => 0.0, 'D' => 0.0, 'C' => 0.0,
                   'Q' => 0.0, 'E' => 0.0, 'G' => 0.0, 'H' => 0.0, 'I' => 0.0,
                   'L' => 0.0, 'K' => 0.0, 'M' => 0.0, 'F' => 0.0, 'P' => 0.0,
                   'S' => 0.0, 'T' => 0.0, 'W' => 0.0, 'Y' => 0.0, 'V' => 0.0,
                   '-' => 0.0);
        
        for(my $j = 0; $j < $num_col; $j++){
            $aas{$seqs[$j][$i]}++;   
        }
        foreach my $aa (keys %aas){
            my $cons_score = $aas{$aa} / $num_col;
            if($cons_score > $conservation[$i]){
                $conservation[$i] = $cons_score;
            }
        }
        if($conservation[$i] > $cutoff){
            push @conserved, "$i $conservation[$i]";
        }
    }
    push @seqs, [@conservation];
    my $num_highly_conserved = scalar(@conserved);
    print "There are a total of $num_highly_conserved sites with a conservation score over $cutoff\n";
    return @seqs;
}

1;
