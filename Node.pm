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


=head1 Node

  Node: Object to represent nodes of a phylogenetic tree

=head1 Synopsis

  use Node;
  
  Put object methods here
  
=head1 Author
  
  Daniel Gaston
  
=head1 Copyright

  This Module is published under the GNU Public License and may be used freely
  under the terms of that license so long as proper attribution is given to the
  original author.

=cut

package Node;

use strict;
use warnings;

#Class Data and Methods, refer to all instances of class
{ 
  my $_count = 0; 
  sub get_count{
    return $_count;
  }
  sub _incr_count{
    ++$_count;
  }
  sub _decr_count{
    --$_count;
  }
}

#Constructor
sub new {
    my ($class, %arg) = @_;
    my $self = bless {
        _node_id       => $arg{node_id} || die "No Node ID given",
        _name          => $arg{name}, 
        _is_leaf       => $arg{is_leaf},
        _is_internal   => $arg{is_internal},
        _is_root       => $arg{is_root},
        _descendents   => $arg{descendents},
        _branchlength  => $arg{branchlength},
        _bootstrap     => $arg{bootstrap},
        _parent        => $arg{parent},      
    }, $class;

    $class -> _incr_count();
    return $self;
}

#Clone Method
sub clone {
  my ($caller, %arg) = @_;
  my $class = ref($caller);
  my $self = bless{}, $class;
  
  foreach my $attribute ($self -> _all_attr() ){
    my($argument) = ($attribute =~ /^_(.*)/);
    if(exists $arg{$argument}){
      $self -> {$attribute} = $arg{$argument};
    }else{
      $self -> {$attribute} = $caller->{$attribute};
    }
  }
  
  $self -> _incr_count();
  return $self;
}

#Simple Accessor Methods
sub get_node_id         {$_[0]   ->  {_node_id} }
sub get_name            {$_[0]   ->  {_name} }
sub is_leaf             {$_[0]   ->  {_is_leaf} }
sub is_internal         {$_[0]   ->  {_is_internal} }
sub is_root             {$_[0]   ->  {_is_root} }
sub _get_desc_arr       {@{$_[0] ->  {_descendents}} }
sub _get_desc_ref       {$_[0]   ->  {_descendents} }
sub get_branchlength    {$_[0]   ->  {_branchlength} }
sub get_parent          {$_[0]   ->  {_parent} }
sub get_bootstrap       {$_[0]   ->  {_bootstrap} }

sub get_descendents{
    my $self = shift;
    my @descendents;
    
    if(defined($self -> _get_desc_ref())){
        @descendents = $self -> _get_desc_arr();   
    }
    
    return @descendents;
}

#Other Methods
=head2 Name

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 
 Note    : 

"Desc"

=cut
sub get_all_descendents{
  my($self) = @_;
  my @all_descendents;
  
  $self -> _get_all_descendents(\@all_descendents);
  
  return @all_descendents;
}

=head2 Name

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 
 Note    : 

"Desc"

=cut
sub _get_all_descendents{
  my ($self, $arr_ref) = @_;
  my @descendents = $self -> get_descendents();
  push @{$arr_ref}, @descendents;
  foreach(@descendents){
    my $leaf = $_ -> is_leaf();
    if(!$leaf){
      $_ -> _get_all_descendents(\@{$arr_ref});
    }
  }
}

=head2 get_all_ancestors

 Title   : get_all_ancestors
 Usage   : $node_obj -> get_all_ancestors();
 Function: Compile a list of ancestors for a given node back to the root node
 Returns : Array of Node Objects
 Args    : None
 Note    : 

"Desc"

=cut

sub get_all_ancestors{
    my ($self) = @_;
    my @ancestors;
    
    $self -> _get_all_ancestors(\@ancestors);
    
    return @ancestors;
}
sub _get_all_ancestors{
    my ($self, $ances_ref) = @_;
    my $parent = $self -> get_parent();
    if($parent){
        push @{$ances_ref}, $parent;
        $parent -> _get_all_ancestors($ances_ref);
    }
}
#Mutator Methods
sub set_name {
  my($self, $name) = @_;
  $self -> {_name} = $name if $name;
}

sub set_descendents {
  my($self, $descendents) = @_;
  my @desc = @{$descendents};
  $self -> {_descendents} = \@desc if @desc;
}

sub clear_descendents {
    my ($self) = @_;
    $self -> {_descendents} = \[];
}

sub set_bootstrap {
  my($self, $bootstrap) = @_;
  $self -> {_bootstrap} = $bootstrap if $bootstrap;
}

sub set_parent {
  my($self, $parent) = @_;
  $self -> {_parent} = $parent if $parent;
}

sub set_branchlength {
  my($self, $branchlength) = @_;
  $self -> {_branchlength} = $branchlength if $branchlength;
}

sub clear_branchlength{
    my $self = shift;
    $self -> {_branchlength} = '';
}

sub add_descendent {
    my($self, $desc) = @_;
    my @descendents = $self -> get_descendents();
  
    push @descendents, $desc;
    $self -> set_descendents(\@descendents);
}

sub set_root{
  my($self, $desc) = @_;
  $self -> {_is_root} = 1;
}

sub set_internal{
  my($self, $desc) = @_;
  $self -> {_is_internal} = 1;
}

sub set_leaf{
  my($self, $desc) = @_;
  $self -> {_is_leaf} = 1;
}

#Destroy Object Method
sub DESTROY {
  my($self) = @_;
  $self -> _decr_count();
}

1;
