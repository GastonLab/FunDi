#Copyright 2016 KR, andrew Roger Lab
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


package iqtree;

use 5.010000;
use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead
# Do not simply export all your public functions/methods/constants

# This allows declaration use QmmRAxML ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw() ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw(calliqtree getiqtreeSiteLH calliqtreeomp);

our $VERSION = '0.1';

sub calliqtree{
    my $string = shift;
    print "Calling iqtree-omp with the following: $string\n";
    `iqtree-omp $string`;
}


sub calliqtreeomp{
    my $string = shift;
    print "Calling iqtree-omp-avx with the following: $string\n";
    `iqtree-omp-avx $string`;
}

sub getiqtreeSiteLH{
    my $file = shift;
    open(IN, "$file") or die "Could not open $file. Exiting...\n\n";
    my @data = <IN>;
    close(IN);

    my @lhs = split /\s+/, $data[1];
    shift @lhs;

    print "Retrieved " . scalar(@lhs) . " site-likelihoods\n";

    return @lhs;
}
1;
