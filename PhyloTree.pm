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


=head1 PhyloTree

  PhyloTree: An object to represent a phylogenetic tree. This moduel serves to
  collect together the nodes of a tree and provides methods for handling trees
  
  This moduel currently will only read and handle a single tree, and is not meant
  for use with multiple trees unless they are all in seperate files and called
  individually. If there are multiple tress in a single file only the first tree
  will be read from the file.

=head1 Synopsis

  use PhyloTree;
  
  
  
  Put object methods here
  
=head1 Author
  
  Daniel Gaston
  
=head1 Copyright

  This Module is published under the GNU Public License and may be used freely
  under the terms of that license so long as proper attribution is given to the
  original author.

=cut

package PhyloTree;
use Node;
use strict;
use warnings;

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
        _file          => $arg{file} || "",
        _tree_id       => $arg{tree_id} || $class -> get_count(),
        _tree_name     => $arg{tree_name} || "",
        _total_length  => $arg{total_length} || "",
        _nodes         => $arg{nodes} || [],      
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

#Accessor Methods
sub get_nodes     {@{$_[0] -> {_nodes}} }
sub get_file      {$_[0]   -> {_file}   }

sub set_nodes{
    my ($self, $nodes_ref) = @_;
    my @nodes = @{$nodes_ref};
    
     $self -> {_nodes} = \@nodes if @nodes;
}

sub set_file{
    my ($self, $file) = @_;
    $self -> {_file} = $file if $file;    
}

sub get_root_node{
    my $self = shift;
    my @nodes = $self -> get_nodes();
    
    foreach(@nodes){
        if($_ -> is_root()){
           return $_;
        }
    }
}

#Methods
=head2 read_tree

 Title   : read_tree
 Usage   : 
 Function: 
 Returns : 
 Args    : 
 Note    : 

"Desc"

=cut

sub read_tree{
    #Set attributes
    my ($self) = @_;
    
    #Read the file data
    my $infile = $self -> {_file};
    open(INTREE, "$infile") or die "Could not open" . $self -> {_file} . "\n";
    my $sep = $/;
    $/ = undef;
    my $raw_data = <INTREE>;
    $/ = $sep;
    close(INTREE);
    
    $raw_data =~ s/\n//g;
    my @trees = split /;/, $raw_data;
    foreach(@trees){
        $_ = $_ . ";";
    }
    
    $self -> {_nodes} = set_all_nodes($self, $trees[0]);
}

sub read_tree_string{
    my ($self, $string) = @_;
    
    $self -> {_nodes} = set_all_nodes($self, $string);
}

=head2 set_all_nodes

 Title   : set_all_nodes
 Usage   : $tree_obj -> set_all_nodes($newick_string)
 Function: Parses a newick tree format string into a tree object
 Returns : List of node objects
 Args    : Newick Format String
 Note    : 

"Desc"

=cut
sub set_all_nodes{
    
    my ($self, $tree_string) = @_;
    #turn the newick format tree string in to an array of 'tokens'
    my $token = '';
    my @tokens = [];
    my $DELIMITER = qr/[():,;]/;

    for(my $i = 0; $i < length($tree_string); $i++){
        $token .= substr($tree_string, $i, 1);
        if($token =~ $DELIMITER){
            if(length($token) > 1){
                my $end = substr($token, (length($token)-1), 1);
                $token =~ s/[():,;]//;
                push @tokens, $token;
                push @tokens, $end;
            }else{
                push @tokens, $token;
            }
            $token = '';
        }
    }
    
    shift @tokens;
    
    
    my $branch = 0;
    if($tree_string =~ /:/g){
        $branch = 1;   
    }
    
    my @open_nodes;
    my $last = "";
    my @nodes;
    my $string;
    my $label;
    my $label_set = 0;
    my $num_nodes = 1;
    
    for(my $i = 0; $i < scalar(@tokens); $i++){
        if($tokens[$i] eq '('){
        	#Create new node
            my $node = Node -> new(node_id => "$num_nodes");
        	push @open_nodes, $node;
        	if($num_nodes == 1){
        		#set the node as the root
                $open_nodes[$#open_nodes] -> set_root();
        	}else{
        		#set as descendent of previous node in stack
        		$open_nodes[($#open_nodes - 1)] -> add_descendent($open_nodes[$#open_nodes]);
        		#set previous node in stack as parent
        		$open_nodes[$#open_nodes] -> set_parent($open_nodes[($#open_nodes - 1)]);
        	}
        	$num_nodes++;
        }elsif($tokens[$i] eq ','){
            #Do nothing, closing nodes has already been incorporated when adding the branchlengths
            #of leaf nodes and when ')' is encountered marking the end of an internal node
        }elsif($tokens[$i] eq ')'){
        	#close open node
        	my $closed_node = pop @open_nodes;
        	push @nodes, $closed_node;
        }elsif($tokens[$i] eq ':'){
        	#Do nothing this is merely inserted as a placeholder
            #so that text such as node labels and branchlengths fall in to
            #the else statement
        }elsif($tokens[$i] eq ';'){
        	#done tree
        }else{
        	if($tokens[($i - 1)] eq '(' || $tokens[($i - 1)] eq ','){
                #found label for a leaf node, create a new node object and set
                #initial values
        		my $node = Node -> new(node_id => "$num_nodes");
                push @open_nodes, $node;
                #set as descendent of previous node in stack
        		$open_nodes[($#open_nodes - 1)] -> add_descendent($open_nodes[$#open_nodes]);
        		#set previous node in stack as parent
        		$open_nodes[$#open_nodes] -> set_parent($open_nodes[($#open_nodes - 1)]);
        		my $label = $tokens[$i];
                $open_nodes[$#open_nodes] -> set_name($label);
        		$num_nodes++;
        	}elsif($tokens[($i + 1)] eq ')' || $tokens[($i + 1)] eq ',' ){
        		if(($tokens[($i - 3)] eq ')') or ($tokens[($i - 2)] eq ')') ){
                    #found the branchlength of an internal node which was just
                    #closed and moved in to the completed node array
        			my $branchlength = $tokens[$i];
        			$nodes[$#nodes] -> set_branchlength($branchlength);
        		}else{
                    #Branchlength of open leaf node
        			my $branchlength = $tokens[$i];
        			$open_nodes[$#open_nodes] -> set_branchlength($branchlength);
        			#close the node
        			my $closed_node = pop @open_nodes;
        			push @nodes, $closed_node;
        		}
        	}elsif($tokens[($i - 1)] eq ')' && $tokens[($i + 1)] eq ':' ){
                #Label/bootstrap for internal node
        		my $label = $tokens[$i];
        		$nodes[$#nodes] -> set_name($label);
                $nodes[$#nodes] -> set_bootstrap($label);
        	}else{
                #branchlength for internal node
        		my $branchlength = $tokens[$i];
        		$nodes[$#nodes] -> set_branchlength($branchlength);
        	}
        }
    }
    foreach(@nodes){
        my @desc = $_ -> get_descendents();
        if(@desc){
            $_ -> set_internal();
        }else{
            $_ -> set_leaf();
        }
    }
    
    $self -> {_nodes} = \@nodes if @nodes;
}

=head2 get_root

 Title   : get_root
 Usage   : $tree_obj -> get_root()
 Function: Find the root node
 Returns : Node Object
 Args    : None
 Note    : 

"Desc"

=cut
sub get_root{
    my $self = shift;
    my @nodes = $self -> get_nodes();
    
    my $root = 0;
    foreach(@nodes){
        my $parent = $_ -> get_parent();
        if(!$parent){
            return $_;
        }
    }
}

=head2 get_leaf_nodes

 Title   : get_leaf_nodes
 Usage   : $tree_obj -> get_leaf_nodes()
 Function: Find all leaf nodes on a tree
 Returns : Array of node objects
 Args    : None
 Note    : 

"Desc"

=cut
sub get_leaf_nodes{
    my $self = shift;
    my @nodes = $self -> get_nodes();
    
    my @leaf_nodes;
    foreach(@nodes){
        if($_ -> is_leaf()){
            push @leaf_nodes, $_;   
        }
    }
    
    return @leaf_nodes;
}


=head2 prune

 Title   : prune
 Usage   : $node_obj -> prune($node)
 Function: Prunes a subtree from the tree at a given node
 Returns : Nothing
 Args    : Node Object
 Note    : 

=cut

sub prune{
    my ($self, $node) = @_;
    my $parent_node = $node -> get_parent();
    my @parent_desc = $parent_node -> get_descendents();
    my $i = 0;
    foreach(@parent_desc){
       if($_ == $node){
            splice @parent_desc, $i, 1; 
       }
       $i++;
    }
    
    $parent_node -> set_descendents(\@parent_desc);
    $self -> collapse_internal($parent_node);
    $node -> set_parent('');
    $node -> set_branchlength('');
}

=head2 get_subtree

 Title   : get_subtree
 Usage   : $tree_obj -> get_subtree($node_obj)
 Function: Creates a new tree object with the root as the node_obj
 Returns : Tree Object
 Args    : Node Object
 Note    : 

"Desc"

=cut

sub get_subtree{
    my ($self, $node) = @_;
    
    my @subtree_nodes = $node -> get_all_descendants();
    my $subtree = new PhyloTree(-nodes => @subtree_nodes);
    
    return $subtree;
}

=head2 get_lca

 Title   : get_lca
 Usage   : $tree_obj -> get_lca($list_of_node_obj)
 Function: Finds the Last Common Ancestor of a group of leaf nodes
 Returns : Node Obj
 Args    : Array of Node Objects
 Note    : 

"Desc"

=cut

sub get_lca{
    my ($self, $nodes_ref) = @_;
    
    my @parents;
    my @nodes = @{$nodes_ref};
    my @ancestor_2d;
    
    foreach(@nodes){
        my @ancestors = $_ -> get_all_ancestors();
        push @ancestor_2d, \@ancestors;
    }
    
   my $lca = 0;
   
   my @ancestor_2d_rev;
   
   for(my $x = 0; $x < scalar(@ancestor_2d); $x++){
        my @temp_arr = reverse @{$ancestor_2d[$x]};
        push @ancestor_2d_rev, \@temp_arr;
   }

    my $i = 0;
    while(!$lca){
        my $last_id = "start";
        for(my $y = 0; $y < scalar(@ancestor_2d_rev); $y++){
            my $id = $ancestor_2d_rev[$y][$i];
            if($last_id eq "start"){
                $last_id = $id;
            }else{
                if($id != $last_id){
                    $lca = $ancestor_2d_rev[$y][$i-1];
                }else{
                    $last_id = $id;
                }
            }
        }
        $i++;
    }
    
    return $lca;
}

=head2 print_tree

 Title   : print_tree
 Usage   : $tree_obj -> print_tree($file);
 Function: Formats a tree object into a newick string and prints to a file
 Returns : None
 Args    : File Name
 Note    : 

"Desc"

=cut
sub print_tree{
    my ($self, $outfile, $subtree) = @_;
    
    my @data = _print_tree_Helper($self -> get_root_node());
    if($data[-1] !~ /\)$/){
        $data[0] = "(" . $data[0];
        $data[-1] .= ")";
    }
    
    my $tree_string = join(',', @data);
    $tree_string .= ";\n";
    
    if($outfile && $subtree){
        $tree_string =~ s/^\((.+\)).*?\:.*\);$/$1/;
        $tree_string =~ s/[\n\r]//;
        $tree_string = "$tree_string" . ";\n";
        open(OUT, ">$outfile") or die "Could not open $outfile for writing!\n";
        print OUT $tree_string;
        close(OUT);
    }elsif($outfile){
        open(OUT, ">$outfile") or die "Could not open $outfile for writing!\n";
        print OUT $tree_string;
        close(OUT);
    }else{
        return $tree_string;
    }
}

sub _print_tree_Helper{
    my $node = shift;
    
    return () if (!defined $node);
    
    my @data;
    
    foreach($node -> get_descendents()){
        push @data, _print_tree_Helper($_);   
    }
    
    if(@data > 1){
        $data[0] = "(" . $data[0];
        $data[-1] .= ")";
        
        my $bootstrap;
        my $name;
        
        if( defined($bootstrap = $node -> get_bootstrap())){
            $data[-1] .= $bootstrap;
        }elsif( defined($name = $node -> get_name())){
            $data[-1] .= $name;
        }
        
        $data[-1] .= ":" . $node -> get_branchlength() if (defined $node -> get_branchlength());
    }else{
        if(defined($node -> get_name() || defined $node -> get_branchlength())){
            push @data, sprintf("%s%s",
				defined $node -> get_name() ? $node -> get_name() : '', 
				defined $node -> get_branchlength() ? ":" .
				$node -> get_branchlength() : '');
        }
        
    }
    
    return @data;
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
sub get_mca{
    
}
=head2 collapse_internal

 Title   : collapse_internal
 Usage   : $tree_obj -> collapse_internal($node_obj)
 Function: "Collapses" an internal node that has only one descendant
 Returns : None
 Args    : Node Obj
 Note    : 

"Desc"

=cut

sub collapse_internal{
    my ($self, $node) = @_;
    my @descendents = $node -> get_descendents();
    if(scalar (@descendents) > 1){
        return 0;   
    }else{
        my $descendent = $descendents[0];
        my $new_branchlen = ($descendent -> get_branchlength()) + ($node -> get_branchlength());
        my $parent = $node -> get_parent();
        $node -> clear_descendents();
        $node -> set_parent('');
        $descendent -> set_parent($parent);
        $descendent -> set_branchlength($new_branchlen);
        my @parent_desc = $parent -> get_descendents();
        my $i = 0;
        foreach(@parent_desc){
           if($_ == $node){
                splice @parent_desc, $i, 1; 
           }
           $i++;
        }
        push @parent_desc, $descendent;
        $parent -> set_descendents(\@parent_desc);
        $self -> cleanup_node($node);
    }
    
    return 1;
}

=head2 cleanup_node

 Title   : cleanup_node
 Usage   : $tree_obj -> cleanup_node($node_ref)
           $tree_obj -> cleanup_node($array_ref)
 Function: Removes a node or list of nodes from a tree objects list of nodes
 Returns : None
 Args    : Reference to Node Object or Reference to Array of Node Objects
 Note    : 

"Desc"

=cut
sub cleanup_node{
    my ($self, $node_ref) = @_;
    my @rem_nodes;
    
    if(ref($node_ref) eq "ARRAY"){
        @rem_nodes = @{$node_ref};
    }else{
        @rem_nodes = $node_ref;
    }
    
    my @nodes = $self -> get_nodes();
    foreach(@rem_nodes){
        my $node_id = $_ -> get_node_id();
        for(my $i = 0; $i < scalar(@nodes); $i++){
            my $id = $nodes[$i] -> get_node_id();
            if($id == $node_id){
                splice(@nodes, $i, 1);   
            }
        }
    }
}


#Mutator Methods

#Destroy Object Method
sub DESTROY {
  my($self) = @_;
  $self -> _decr_count();
}
1;
