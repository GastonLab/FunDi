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


package RAxML;

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
our @EXPORT = qw(callRAxML callParallelRAxML getRAxMLSiteLH);

our $VERSION = '0.1';


# Preloaded methods go here.
sub callRAxML{
    my $string = shift;
    print "Calling raxmlHPC with the following: $string\n";
    `raxmlHPC $string`;
}

sub callParallelRAxML{
    my $string = shift;
    print "Calling raxmlHPC-PTHREADS-SSE3 with the following: $string\n";
    `raxmlHPC-PTHREADS-SSE3 $string`
}

sub getRAxMLSiteLH{
    my $file = shift;
    open(IN, "$file") or die "Could not open $file. Exiting...\n\n";
    my @data = <IN>;
    close(IN);
    
    my @lhs = split /\s+/, $data[1];
    shift @lhs;
    
    return @lhs;
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
