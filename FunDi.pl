#Copyright 2009, 2010 Daniel Gaston, Andrew Roger Lab
#This code is copyrighted under the GNU General Public License Version 3.0
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
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


#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use PhyloTree_new;
use Seq;
use Node;
use Extras;
use QmmRAxML;
use RAxML;
use iqtree;
use Math::BigFloat;
use List::Util qw(min max);

################################################################################
#Author:              	Daniel Gaston
#Last Modified:            July 2016
#Version:                       1.1
################################################################################

################################################################################
#                               Constants                                      #
################################################################################

my $OPTDIR = "Optimization/";
my $NUM_CLASSES = 10;
my $phy_file_name;
my $optimal_lh;
my $branch_opt;
my $max_opt;
my $step = 0.02;#####

################################################################################
#                               Main                                           #
################################################################################

if(scalar(@ARGV) == 0){
    runMenuMode(@_);
}else{
    runCLIMode(@_);
}


################################################################################
#                            Subroutines                                       #
################################################################################

sub runCLIMode{
    if(scalar(@ARGV) == 0){
        @ARGV = @_;
    }

    open(LOG, ">>FunDi.log") or die "Could not open FunDi.log...\n\n";
    print LOG "FunDi running with specified options @ARGV\n";

    my $intree = "";
    my $inalign = "";
    my $outroot = "subtree";
    my $subtree = "subtree.def";
    my $midroot;
    my $model = "LG+C20+F+G";  ####
    my $num_rates = 4;
    my $parallel = 0;
    my $help = 0;
    my $verbose = 0;
    my $program = "iqtree";   ####
    my $length = 10000;
    my $debug = 0;
    my $mask_cut = 1.0;
    my $num_cpus = 1;
    my $retree_command = "retree";
    my @masked_sites;
    my $seed = 20120821;
	my $branch_len_opt = 1;

    GetOptions(
               'tree=s'     =>   \$intree,   'align=s'   =>   \$inalign,
               'outroot=s'  =>   \$outroot,  'Midroot'   =>   \$midroot,
               'subtree=s'  =>   \$subtree,  'parallel'  =>   \$parallel,
               'debug'      =>   \$debug,    'help'      =>   \$help,
               'verbose'    =>   \$verbose,  'length=i'  =>   \$length,
               'model=s'    =>   \$model,    'rates=i'   =>   \$num_rates,
               'Program=s'  =>   \$program,  'cutoff=f'  =>   \$mask_cut,
               'Num_cpus=s' =>   \$num_cpus, 'Opt_bl' 	 =>   \$branch_len_opt 
              );

    if($inalign eq "" || $help){
        printHelp();
        exit();
    }

    $program = lc $program;

    if($program ne 'qmmraxml' && $program ne 'raxml' && $program ne 'iqtree'){
        print "Unknown phylogeny program selected. Must be either qmmraxml, iqtree or raxml.\nExiting...\n\n";
        exit();
    }

    if($parallel && $num_cpus == 1){
        $num_cpus = 2;
    }

    # Not being used
    if($mask_cut < 1.0){
        my $ref;
        ($inalign, $ref) = maskAlignment($inalign, $mask_cut);
        @masked_sites = @{$ref};
    }

    # Convert FASTA file to PHYLIP
    if($inalign =~ /\.fasta$/ || $inalign =~ /\.fst$/ || $inalign =~ /\.faa$/){
        print "Detected by file extension that $inalign is a FASTA formated file. Exporting Phylip\n";
        my %seq = Seq::readFasta($inalign);
        my $out = $inalign . ".phy";
        print "Printing as $out\n";
        Seq::printPhylip($out, \%seq);
        $inalign = $out;
    }

  $phy_file_name = $inalign; # This gets used for output files but $inalign is already .phy file name by this point

  # Retrieving from .def and creating a 2D array and subtrees
  print "Retrieving subtree definitions from $subtree...\n";
  my @subtrees = getSubtreesFile($subtree);

  # Optimize branch lengths of the input tree using the chosen ML program and model
  if($program eq 'qmmraxml'){
      print "Optimizing branch lengths of user input tree using QmmRAxML\n";
      my $raxml_out = "$inalign" . "_optimized.out";
      my $raxml_model = "PROTGAMMA" . "$model";
      my $dirichlet = "9componentsDirichlet_" . "$model" . ".dat";
      my $raxml_string = "-s $inalign -t $intree -n $raxml_out -m $raxml_model -u $dirichlet -f e";
      print "Executing QmmRAxML with following command-line: $raxml_string\n";
      callQmmRAxML($raxml_string);
  }elsif($program eq 'iqtree'){
      print "Optimizing branch lengths of user input tree using IQ-Tree\n";
      my $iqtree_out = "$inalign" . ".treefile";  ####
      my $iqtree_string = "-s $inalign -te $intree -m $model -nt $num_cpus -wsl";
      calliqtree($iqtree_string);
  }else{
      print "Optimizing branch lengths of input user tree using RAxML\n";
      my $raxml_out = "$inalign" . "_optimized.out";
      my $raxml_model = "PROTGAMMA" . "$model" . "F";
      my $raxml_string = "-s $inalign -t $intree -n $raxml_out -m $raxml_model -f e";
      if($parallel){
          $raxml_string = $raxml_string . " -T $num_cpus";
          callParallelRAxML($raxml_string);
      }else{
          callRAxML($raxml_string);
      }
  }
  print "Finished optimizing\n";

  # Get the optimized likelihood value of the tree from the appropriate info file
  # depending on ML program used
  my $lh;
  my $info_file;
  if ($program eq 'iqtree'){
      $info_file = "$inalign" . ".log";   ####
      open(IN, "$info_file") or die "Could not open $info_file. Exiting...\n\n";
      my @info = <IN>;
      close(IN);
      foreach(@info){
          if($_ =~ /Optimal log-likelihood:/){
              my @temp = split /: /, $_;
              $lh = $temp[1];
              $lh =~ s/[\r\n]+//g;
              last();
          }
      }
  }else{
      $info_file = "RAxML_info." . "$inalign" . "_optimized.out";
      open(IN, "$info_file") or die "Could not open $info_file . Exiting...\n\n";
      my @info = <IN>;
      close(IN);
      foreach(@info){
          if($_ =~ /Final GAMMA  likelihood:/){
              my @temp = split /: /, $_;
              $lh = $temp[1];
              $lh =~ s/[\r\n]+//g;
              last();
          }
      }
  }


  if($lh eq ""){
      print "WARNING: Could not locate likelihood value in $info_file . Exiting...\n\n";
      exit();
  }
	$optimal_lh = $lh;
  print "Finished Likelihood calculation: $optimal_lh\n";

  # Renaming tree file
  if($program eq 'iqtree'){
      $intree = "$inalign" . ".treefile";
  }else{
      $intree = "$inalign" . "_optimized.tre";
  }

  # Reading optimized newick tree file and printing out to a new file name (treefile or optimized)
  my $file;
  if($program eq 'iqtree'){
      $file = "$inalign" . ".treefile";
  }else{
      $file = "RAxML_result." . "$inalign" . "_optimized.out";
  }

  print "Reading input tree file $file\n";
  open(IN, "$file") or die "Could not find $file . Exiting...\n";
  my $sep = $/;
  $/ = undef;
  my $data = <IN>;
  close(IN);
  $/ = $sep;

  print "Tree with branchlengths re-optimized by $program had lh of $lh...\n";
  print LOG "Tree with branchlengths re-optimized by $program had lh of $lh...\n";
  open(OUT, ">$intree") or die "could not create $intree!\n";
  print OUT $data;
  close(OUT);


  # Code for handling some special cases of running in midroot mode (Not recommended)
  # to re-roote the tree
  if($midroot && $program ne 'puzzle'){
      print "Checking to see if file outtree already exists...\n";

      if(-f "outtree"){
            print "File outtree already exists, moving...\n";
            `mv outree old_outree.tre-bak`;
      }

      print "Creating Phylip Retree config file...\n";
      open(CONFIG, ">retree.config") or die "Could not create retree.config!\n";
      print CONFIG "y\n$intree\nm\nw\nr\nq\n";
      close(CONFIG);

      print "Passing $intree to retree to midroot tree...\n";
      `$retree_command < retree.config`;

      print "Renaming outtree to outtree_midroot.tre\n";
      $intree = "outtree_midroot.tre";
      `mv outtree $intree`;
  }

  print "Parsing $intree and $inalign to create subtree files...\n";    # $intree = "$inalign" . ".treefile"
  my @files = treeParse($debug, $intree, $inalign, $outroot, $midroot, \@subtrees);
  print "files" . @files . "\n";  ####
  my @subfiles = @{$files[0]};
  print @subfiles;  ####
  my @lcas = @{$files[1]};
  print @lcas; ####
  my $tree = $files[2];
  print $tree;  ####
  print "\n=============================================\n\n";

  my $test_file = runMixtureModel($inalign, \@subtrees, \@subfiles, \@lcas, $tree, $midroot, $program, $num_rates, $model, $parallel, $num_cpus, \@masked_sites, $branch_len_opt);
	#my $command = "rm *_subtree*";
	#system($command);
  print "\nFinished running!\n";
  
  close(LOG);
}

sub runMixtureModel{
    my $startTime = time;
    my ($alignment, $subref, $subfilesref, $lcas_ref, $tree, $midroot, $program, $num_rates, $model, $parallel, $num_cpus, $masked_ref, $branch_len_opt) = @_;###added $branch_len_opt
    my @subtrees = @{$subref};
    my @subfiles = @{$subfilesref};
    my @old_lcas = @{$lcas_ref};
    my @masked_sites = @{$masked_ref};
    my $branchlength = 0;

    print "Running Mixture Model\n";

    # Using the LCA Node objects to read the input tree file and get the starting internal
    # branch length between the subtrees
    my $tree_file = $tree -> get_file();
    my $new_tree = new PhyloTree_new(file => "$tree_file");
    $new_tree -> read_tree();
    my @nodes = $new_tree -> get_nodes();
    my @ids;

    foreach(@old_lcas){
        my $id = $_ -> get_node_id();
        push @ids, $id;
    }

    my @lcas;
    foreach my $node (@nodes){
        my $temp_id = $node -> get_node_id;
        foreach my $id (@ids){
            if($temp_id == $id){
                if(! $node -> is_root()){
                    push @lcas, $node;
                }
            }
        }
    }

    if($midroot){
        print "Running in Midroot Mode...\n";
        if(scalar(@lcas) != 2){
            print "Running in midroot mode but there are " + scalar(@lcas) + " lcas\n";
            exit();
        }
    }

    print "Retrieving branch length information...\n";
    if(scalar(@lcas) == 2 && $midroot){
        foreach(@lcas){
            my $temp_branch = $_ -> get_branchlength();
            $branchlength = $branchlength + $temp_branch;
        }
        print "Internal branch length of midrooted tree is $branchlength\n";
        print LOG "Internal branch length of midrooted tree is $branchlength\n";
    }elsif(scalar(@lcas) == 1){
        $branchlength = $lcas[0] -> get_branchlength();
        print "Internal branch length is $branchlength\n";
        print LOG "Internal branch length is $branchlength\n";
    }else{
        print "More than one LCA?\n";
        print "Currently only support use of two clusters in midroot or unrooted mode\n";
        exit();
    }

    print "\n=============================================\n\n";
	
    # Subtree files should already be sorted but this guarantees that subgroups are together and that for each
    # subgroup the alignment file comes first [i] and the treefile comes next [j]
    my @sorted_out = sort @subfiles;
    my @subtree_lh_fh;
    print "Evaluating subtrees...(@sorted_out)\n";

    # These don't appear to be used anywhere. Commenting out for now
    # my @raxml_info_files;
    # my @iqtree_info_files;
    # Run ML program of choice on defined subtrees and place the resulting sitelh
    # files into the subtree_lh_fh array
    for(my $i = 0; $i < scalar(@sorted_out); $i = $i + 2){
        my $lh_fh;
        if($program eq 'qmmraxml'){
            my $j = $i + 1;
            print "Running QmmRAxML step for $sorted_out[$i] and $sorted_out[$j]\n";
            my $raxml_out = "$sorted_out[$j]" . "_optimized.out";
            # my $info = "RAxML_info." . "$raxml_out";
            # push @raxml_info_files, $info;
            my $raxml_model = "PROTGAMMA" . "$model";
            my $dirichlet = "9componentsDirichlet_" . "$model" . ".dat";
            $lh_fh = "$sorted_out[$j]" . ".sitelh";
            my $raxml_string = "-s $sorted_out[$i] -t $sorted_out[$j] -n $raxml_out -m $raxml_model -u $dirichlet -L $lh_fh -f e";
            callQmmRAxML($raxml_string);
        }elsif($program eq 'iqtree'){
            my $j = $i + 1;
            print "Running IQ-Tree step for $sorted_out[$i] and $sorted_out[$j]\n";
            my $iqtree_out = "$sorted_out[$i]" . ".treefile";   ####
            $lh_fh = "$sorted_out[$i]" . ".sitelh";
            #push @iqtree_info_files, $info;
            my $iqtree_model = "$model";
            my $iqtree_string = "-s $sorted_out[$i] -te $sorted_out[$j] -m $model -nt $num_cpus -wsl -fixbr";
            calliqtree($iqtree_string);
        }elsif($program eq 'raxml'){
            my $j = $i + 1;
            print "Running RAxML step for $sorted_out[$i] and $sorted_out[$j]\n";   # _subtree0.phy and _subtree0.tre
            my $raxml_out = "$sorted_out[$j]" . "_optimized.out";
            # my $info = "RAxML_info." . "$raxml_out";
            # push @raxml_info_files, $info;
            my $raxml_model = "PROTGAMMA" . "$model";
            my $raxml_string = "-s $sorted_out[$i] -z $sorted_out[$j] -n $raxml_out -m $raxml_model -f g";
            $lh_fh = "RAxML_perSiteLLs." . $raxml_out;
            if($parallel){
                $raxml_string = $raxml_string . " -T $num_cpus";
                callParallelRAxML($raxml_string);
            }else{
                callRAxML($raxml_string);
            }
        }else{
            print "Unknown phylogeny program. Exiting...\n\n";
            exit();
        }
        print ("*******************************************\n");
        push @subtree_lh_fh, $lh_fh;
    }

    # Calculate the constant value for the subtree site likelihoods
    print "Calculating subtree site likelihoods (independence model)\n";
    my @subtree_lhs;
    my $last;

    # Retrieving subtree site likelihoods and pushing into a 2D array
    foreach(@subtree_lh_fh){
        my @lhs ;
        if($program eq 'puzzle'){
            @lhs = getSiteLH($_);
        }elsif($program eq 'qmmraxml'){
            @lhs = getSiteCombLH($_);
        }elsif($program eq 'iqtree'){
            @lhs = getiqtreeSiteLH($_);
        }else{
            @lhs = getRAxMLSiteLH($_);
        }

        my $length = scalar(@lhs);
        if($last){
            if($length != $last){
                print "Subtree likelihood files did not contain the same number of entries.\n";
                print "Exiting...\n";
                exit();
            }
        }
        $last = $length;
        push @subtree_lhs, [@lhs];
    }

    # Sum the likelihoods of the subtrees at each site
    # Modified by kr to calculate site log-likelihoods for FD on May 13, 2016
    # This modification prevents underflow issues with very small likelihoods
    my @subtree_site_lnls;
    # Loop through the 2D array, doing each subtree independently. The outer for Loop
    # loops through the total number of sites which were counted previously into $last
    for(my $i = 0; $i < $last; $i++){
        my $total_lh = 0;
        # This foreach loop loops over each subgroup/subtree
        foreach(@subtree_lhs){
            my @temp = @{$_};   # Temp holds the column
            $total_lh = $total_lh + $temp[$i]; # Sum of the subtree likelihoods at the site
            # print ("$total_lh\n");
        }
        push @subtree_site_lnls, $total_lh; # 1D array
    }
    print "\n=============================================\n\n";

    # Coarse grain grid search optimization of internal branchlength and rho
    print "Performing coarse grain optimization...\n";
    my $p_lower = 0.00000000001;
    my $p_upper = 1;
    my $p_step = 0.01;
    my $upper = $branchlength;
    my $lower = 0.0;
    my @branch_files;
	my @add_tree_files;

	if ($branch_len_opt eq 1){
		@add_tree_files = create_new_tree_file($tree_file, $branchlength, $step);
	} else {
		$tree_file =~ s/[\r\n]+//g;
		push @add_tree_files, $tree_file;
		print "No internal optmization of internal branchlength. Evaluating site likelihoods...\n";
	}
	my $array_size = scalar(@add_tree_files);
	print "\n $array_size \n";
	foreach my $temp_tree_file(@add_tree_files){
	
		# In the case of QmmRAxML insert the wholetree RAxML info file onto the first of the data files array to
		# preserve the internal order. Wholetree followed by subtrees.	
		print "Using Tree File: $temp_tree_file\n"; 
		if($program eq 'qmmraxml'){
			print "Running QmmRAxML for $alignment and $temp_tree_file...\n";
			my $raxml_out = "$temp_tree_file" . ".out";
			#my $info = "RAxML_info." . "$raxml_out";
			#unshift @raxml_info_files, $info;
			my $raxml_model = "PROTGAMMA" . "$model";
			my $site_file = "$temp_tree_file" . ".sitelh";
			my $dirichlet = "9componentsDirichlet_" . "$model" . ".dat";
			my $raxml_string = "-s $alignment -t $temp_tree_file -n $raxml_out -m $raxml_model -u $dirichlet -L $site_file -f e";
			callQmmRAxML($raxml_string);
		}elsif($program eq 'iqtree'){
			print "Running IQ-Tree for $alignment and $temp_tree_file...\n";  
			my $iqtree_out = "$temp_tree_file" . ".treefile"; 
			my $iqtree_model = "$model";
			my $iqtree_string;
			if ($branch_len_opt eq 1){
				$iqtree_string = "-s $alignment -te $temp_tree_file -m $model -nt $num_cpus -wsl -fixbr";
			} else {
				$iqtree_string = "-s $alignment -te $temp_tree_file -m $model -nt $num_cpus -wsl";
			}
        
			# Rename the sitelh file to match up to assumptions later in the code
			# which are based on the treefile name as opposed to the alignment name
			# This is because of pre-existing code that was for dealing with the internal
			# branchlength optimization code that previously existed
			my $site_file = "$alignment" . ".sitelh";
			my $lh_fh = "$temp_tree_file" . ".sitelh";
			`cp $site_file $lh_fh`;
			my $command = "rm *ckp.gz";
			system($command);
			calliqtree($iqtree_string);   ####
		}elsif($program eq 'raxml'){
			print "Running RAxML for $alignment and $temp_tree_file...\n";
			my $raxml_out = "$temp_tree_file" . ".out";
			#my $info = "RAxML_info." . "$raxml_out";
			#unshift @raxml_info_files, $info;
			my $raxml_model = "PROTGAMMA" . "$model";
			my $raxml_string = "-s $alignment -z $temp_tree_file -n $raxml_out -m $raxml_model -f g";
			my $site_file = "RAxML_perSiteLLs." . "$raxml_out";
			if($parallel){
				$raxml_string = $raxml_string . " -T $num_cpus";
				callParallelRAxML($raxml_string);
			}else{
				callRAxML($raxml_string);
			}

        # Rename the sitelh file to match up to assumptions later in the code
        # which are based on the treefile name as opposed to the alignment name
        # This is because of pre-existing code that was for dealing with the internal
        # branchlength optimization code that previously existed
        my $lh_fh = "$temp_tree_file" . ".sitelh";
        `cp $site_file $lh_fh`;
		}else{
			print "Unknown phylogeny program. Exiting...\n\n";
			exit();
		}

		# Branchfiles is an array, left over from dealing with optimizing the internal
		# branch length between subgroups as well as the Rho parameter.
		push @branch_files, $temp_tree_file;
	}
 
    my @values = optimizeRho("Rough", 0, 0.000000000001, $branch_files[1], \@branch_files, \@subtree_site_lnls, $p_lower, $p_upper, $p_step, $program);
    print "\n=============================================\n\n";
	#print @branch_files;  # so far contains only one treefile, after adding bl optimization, should contain more than one treefile files.  -added
	#print @subtree_site_lnls;   # contains site_likelihoods -the result of each site is the sum of the 2 subtrees
	#print "\n=============================================\n\n";

    # Fine scale optimization of branchlength and rho
    my $rough_p    = $values[0];
    my $rough_b_fh = $values[1];
    my $prior_max  = $values[2];

    my @fine_branch_files;
    my $branch_upper;
    my $branch_lower;

    $tree_file =~ s/[\r\n]+//g;
    push @fine_branch_files, $tree_file;


    my $new_p_lower = $rough_p - (2 * $p_step);
    my $new_p_upper = $rough_p + (2 * $p_step);
    my $new_p_step = $p_step / 100;
    if($new_p_upper > 1){
        $new_p_upper = 1.0;
    }

    my @new_values = optimizeRho("Fine", $prior_max, $rough_p, $rough_b_fh, \@branch_files, \@subtree_site_lnls, $new_p_lower, $new_p_upper, $new_p_step, $program);
    print "\n=============================================\n\n";

    my @new_temp = split /-/, $new_values[1];
    my $bfh = $new_temp[-1];
    $bfh =~ s/[\n\r]+//g;
    $bfh =~ s/\.sitelh//;

    my @mm_site_loglhs = @{$new_values[3]};   
	# print "\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n";
	# print @mm_site_loglhs;   
	# print "\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n";
	
	
    # Class posteriors tree file order: Wholetree then subtrees
    my @site_likelihoods;
    my $rho = $new_values[0];
    my @class_posteriors;
    if($program eq 'puzzle'){
        @site_likelihoods = getSiteLH($new_values[1]);
    }elsif($program eq 'qmmraxml'){
        @site_likelihoods = getSiteCombLH($new_values[1]);
        @class_posteriors = getSiteClassPost($new_values[1], \@subtree_lh_fh, $NUM_CLASSES);
    }elsif($program eq 'iqtree'){    
        @site_likelihoods = getiqtreeSiteLH($new_values[1]);
    }else{
        @site_likelihoods = getRAxMLSiteLH($new_values[1]);
    }

    # Calculate the posterior probability of functional divergence
    open(POST, ">FunDi_Posterior_Scores." .$phy_file_name .".txt") or die "Could not create FunDi_Posterior_Scores.txt\n\n";
    #print POST "Site\tP(FD)\tMixture Model Site Log-Likelihood\tMixture Model Likelihood\tSubtree Combined Likelihood\tSWeighted Subtree Likelihood\n";
    print POST "Site\tPosterior Prob of FD\tMixture Model Site Log-Likelihood\tPosterior Prob Denominator\tSubtree Combined Log-Likelihood\n";
    my $site = 0;
    my @sites;
    for(my $i = 0; $i < scalar(@site_likelihoods); $i++){
        if(scalar(@masked_sites) > 0 && $site == $masked_sites[0]){
            print POST "$site\tNONE\tNONE\n";
            shift @masked_sites;
            $site++;
            $i--;
            next();
        }

        my $lnjntnfd = log($rho) + $site_likelihoods[$i];
        my $lnjntfd = log(1 - $rho) + $subtree_site_lnls[$i];
        my @lnjnts = ($lnjntfd, $lnjntnfd);
        my $lnnfd_denominator = max(@lnjnts) + log(1 + exp(min(@lnjnts) - max(@lnjnts)));
        my $lnpost = log($rho) + $site_likelihoods[$i] - $lnnfd_denominator;
        my $nfd_posterior = exp($lnpost);
        my $posterior = 1 - $nfd_posterior;

        print POST "$site\t$posterior\t$mm_site_loglhs[$i]\t$lnnfd_denominator\t$subtree_site_lnls[$i]\n";

        if($posterior > 0.5){
            push @sites, $posterior;
        }
        $site++;
    }

    close(POST);

    my $fraction = scalar(@sites) / scalar(@site_likelihoods);
    print "Iterated over $site sites\n";
    print "Posterior probabilities written to Posteriors.txt\n";
    print "Fraction of sites with posterior probability of functional divergence above 0.5 is: $fraction\n";
    print LOG "Fraction of sites with posterior probability of functional divergence above 0.5 is: $fraction\n";
    print "Rho (FD Weight): " . (1 - $rho) . "\n";   # Rho is actually the weight associated with Non FD in the code
    my $endTime = time;
    print $endTime - $startTime, " seconds to run " . "\n";
    close(LOG);
    my $test_file = "rho." . $phy_file_name . ".txt";
    print $test_file;
    open (TEST, ">$test_file") or die "file could not be opened";
    print TEST "Fraction of sites with posterior probability of FD above 0.5: " . $fraction . "\n";
    print TEST "Rho (FD weight): " . (1 - $rho) . "\n";
    print TEST "Internal branch length: " . $branchlength ."\n";
	print TEST "Optimal Log-Likelihood (IQ-Tree): " . $optimal_lh ."\n";
	print TEST "Optimal Log-Likelihood (After Branch Length Optimization): " . $max_opt ."\n";
	print TEST "Optimal Internal branch length: " . ($branch_opt) ."\n";  #####
    close(TEST);
	return $test_file;
}

# Optimization of the p parameter (mixture model weight) for each branchlength
sub optimizeRho{
    my ($method, $prior_max, $opt_p, $opt_b, $branch_ref, $subtree_ref, $lower, $upper, $in_function_step, $program) = @_;
    my @branch_files = @{$branch_ref};
    #my @subtree_site_probs = @{$subtree_ref};
    my @subtree_site_lnls = @{$subtree_ref};
    my @values;

    print "Optimizing rho parameter... $method\n";
    my $max;

    if($prior_max){
        $max = $prior_max;
    }

    my @site_loglhs;
    # Loop over all provided branch files (Leftover from branch length optimization)
    foreach my $file_name (@branch_files){
        $file_name = "$file_name" . ".sitelh";

        print "Retrieving site likelihoods for $file_name...\n";
        my @lhs;
        if($program eq 'puzzle'){
            @lhs = getSiteLH($file_name);
        }
        elsif($program eq 'qmmraxml'){
            @lhs = getSiteCombLH($file_name);
        }
        elsif($program eq 'iqtree'){
            #$file_name =~ s/treefile.//ig;
            @lhs = getiqtreeSiteLH($file_name);
        }
        else{
            @lhs = getRAxMLSiteLH($file_name);
        }

        # Calculation of mixture model equation with a given value of p
        print ("Lower bound: $lower Upper Bound: $upper Step Size: $in_function_step\n");
        for(my $p = $lower; $p <= $upper; $p = $p + $in_function_step){
            if(scalar(@lhs) != scalar(@subtree_site_lnls)){
                print "ERROR: During optimization of rho with value of $p arrays for probabilities of whole tree and subtrees not same length\n";
                print "Exiting...\n";
                exit();
            }
            my $total = 0;
            my @temp_site_loglhs;
            for(my $i = 0; $i < scalar(@lhs); $i++){
                my $lnjntnfd = log($p) + $lhs[$i];
                my $lnjntfd = log(1-$p) + $subtree_site_lnls[$i];
                my @lnjnts = ($lnjntfd, $lnjntnfd);
                my $lnprob = max(@lnjnts) + log(1 + exp(min(@lnjnts) - max(@lnjnts)));
                push @temp_site_loglhs, $lnprob;
                $total = $total + $lnprob;
            }

            my @temp = split /-/, $file_name;
            my $branchlength = $temp[-1];
            $branchlength =~ s/\.sitelh//;
            if($max && $total > $max){
                $opt_p = $p;
                $opt_b = $file_name;
                $max = $total;
                @site_loglhs = @temp_site_loglhs;   # temp_site_loglhs contains # of columns in the alignment
            }
            elsif($max && $total < $max){
                next();
            }
            else{
				 $opt_p = $p;
                 $opt_b = $file_name;
                 $max = $total;
                 @site_loglhs = @temp_site_loglhs;
            }
      }
    }
	
	
    #print "Opt P: $opt_p Opt b: $opt_b Max: $max\n";
	$max_opt = $max;
    push @values, $opt_p;
    push @values, $opt_b;
    push @values, $max;
    push @values, \@site_loglhs; 
	my @tab = split( /\./ ,$opt_b);
	if (scalar(@tab) >4 ){
		# $branch_opt = $tab[6].".".$tab[7];
		# $branch_opt = $tab[3].".".$tab[4];
		$branch_opt = $tab[3].".".$tab[4];
		#Taking the branch_opt minus the last step as we always do a step forward
		if ($branch_opt - $step > 0){
			$branch_opt = $branch_opt - $step;
		}
		#print $branch_opt . "\n";
		#print $max_opt . "\n";
	}
	print "Opt P: $opt_p branch_Opt: $branch_opt Max: $max\n";
	#push @values, $branch_opt; 
    return @values;
}

sub getSubtreesFile{
    my $file = shift;
    open(IN, "$file") or die "Could not open $file\n";
    my @subtrees;
    my @temp;
    while(<IN>){
        if($_ !~ m/^#/ && $_ ne "\n"){
            @temp = split;
            foreach(@temp){
                $_ =~ s/\n//g;
                $_ =~ s/ //g;
            }
            push(@subtrees, [ @temp ]);
        }
    }
    close(IN);

   return @subtrees;
}

sub treeParse{
    my ($debug, $intree, $align, $out_root, $midroot, $subtree_ref) = @_;
    my @subtrees = @{$subtree_ref};
    my @outfiles;
    my @lcas;
	print "creating phylotree using tree $intree\n";
	
    my $tree = new PhyloTree_new(file => "$intree");  # $intree = "$inalign" . "_optimized" . ".treefile"
    $tree -> read_tree();

    my @taxa = $tree -> get_leaf_nodes();
    print "Parsed tree with " . scalar(@taxa) . " leaf nodes\n";
	print "entering loop ===================\n";  ####
	
    my $done = 0;
    my $count = 0;
	my $num_of_times_i_ran = 0;
    while(!$done){
		
        for(my $i = 0; $i < scalar(@subtrees); $i++){
			
            my @temp_subtree = @{$subtrees[$i]};
            my @subtree_taxa;
			#print "this is the value of temp_subtree:" . @temp_subtree . "\n";  ####
			#print "this is the value of taxa:" . @taxa . "\n";  ####
            foreach my $node (@temp_subtree){
				print "this is the value of node:" . $node . "\n";  ####
                foreach my $taxon (@taxa){
                    my $taxon_name = $taxon -> get_name();
                    if($taxon_name eq $node){
						print "this is the value of taxon_name:" . $taxon_name . "\n";  ####
                        push @subtree_taxa, $taxon;
                    }
                }
            }
			#print "finished that loop ====\n";####
            my $subtree_out = "$out_root" . "$count" . ".tre";
            my $lca = $tree -> get_lca(\@subtree_taxa);
			#print "got lca" . $lca . "\n"; ####
            my $subtree = create_subtree($lca, $subtree_out);
			#print "got subtree" . $subtree . "\n";
            my $ok = check_subtree($subtree, \@subtrees, $i);
			print "checking if ok\n";
			
            if($ok){
                if(!$debug){
                    # The 1 indicates that this is a subtree, needed for print_tree
                    $subtree -> print_tree($subtree_out, 1);
                }
				print "Pushing lca to lcas \n";
                push @lcas, $lca;

                if($align){
                    my $subtree_align_out = "$out_root" . "$count" . ".phy";
                    push @outfiles, $subtree_align_out;

                    if(!$debug){
                        format_subtree_align($align, $subtree, $subtree_align_out);
                    }

                }else{
                    die "Alignment file was not provided!\n";
                }

                push @outfiles, $subtree_out; # Pushes into an array and at the end creates subtreex.tre file

                if(scalar(@subtrees) > 1){
                    if(!$midroot){
                        $tree -> prune($lca);
                    }
                    splice @subtrees, $i, 1;
                    $i--;
                    $count++;
                }else{
                    $done = 1;
                }
            }
        }
        if(scalar(@subtrees) == 0){
            $done = 1;
        }
		$num_of_times_i_ran++;
		if ($num_of_times_i_ran > 5 ){
			die "Cannot parse the subtrees, exiting...\n";
		}

    }

    return (\@outfiles, \@lcas, $tree);
}

sub create_subtree{
    my $lca = shift;
    my $subtree_out = shift;

    $lca -> set_root();
    $lca -> clear_branchlength('');
    my @subtree_nodes = $lca -> get_all_descendents();
    unshift @subtree_nodes, $lca;
    my $subtree = new PhyloTree_new(file => "$subtree_out");
    $subtree -> set_nodes(\@subtree_nodes);

    return $subtree;
}

sub check_subtree{
    my ($subtree, $list_ref, $element) = @_;
    my @subtree_list = @{$list_ref};
    splice @subtree_list, $element, 1;
    my @other_taxa;
    foreach(@subtree_list){
        my @temp = @{$_};
        foreach my $taxa (@temp){
            push @other_taxa, $taxa;
        }
    }

    my $ok = 1;
    my @subtree_taxa = $subtree -> get_leaf_nodes();
    foreach(@subtree_taxa){
        my $id = $_ -> get_name();
        foreach my $name (@other_taxa){
            if($name eq $id){
                $ok = 0;
            }
        }
    }
    if($ok){
        return 1;
    }else{
        return 0;
    }
}


sub format_subtree_align{
    my ($align, $subtree, $out) = @_;

    open(ALIGN, "$align") or die "Could not read $align!\n";
    my @align = <ALIGN>;
    close(ALIGN);
    open(OUT, ">$out") or die "Could not open $out for writing!\n";

    my @taxa = $subtree -> get_leaf_nodes();
    my @lines;
    shift @align;
    foreach my $line (@align){
        my $name;
        if($line =~ /\s/){
            my @temp = split /\s+/, $line;
            $name = $temp[0];
        }else{
            $name = substr($line, 0, 10);
        }
        $name =~ s/ //g;
        foreach my $taxon (@taxa){
            my $id = $taxon -> get_name();
            if($name eq $id){
                push @lines, $line;
            }
        }
    }

    my $num_taxa = scalar(@lines);
    my $sequence = $lines[0];
    $sequence =~ s/^(.+) +(.+)[\n\t]$/$2/;
    $sequence =~ s/[\n\r]+//;
    my $length = length $sequence;

    print OUT " $num_taxa $length\n";
    foreach(@lines){
        print OUT $_;
    }

    close(OUT);
}

sub runMenuMode{
    my $set = 0;
    my $intree;
    my $inalign;
    my $outroot = "subtree";
    my $subtree = "subtree.def";
    my $parallel = 0;
    my $mask_cut = 1.0;
    my $program = "iqtree";  #qmmraxml
    my $model = "LG+C20+F+G";  #LG
    my $num_rates = 4;
    my $help = 0;
    my $verbose = 0;
    my $midroot = 0;
    my $Num_cpus = 1;
	my $branch_len_opt = 1; 

    print "Now running in Menu Mode...\n";
    print "\nPlease enter an alignment file: ";
    $inalign = <STDIN>;
    chomp $inalign;
    print "\n";

    my @trees0 = <RAxML_result.*>;
    my @tree0 = <*.treefile>;  ###
    my @trees1 = <*.tre>;

	# Not relevant
    my @trees2 = <*.tree>;
    my @trees3 = <RAxML_bipartitions.*>;

    my @trees;
    push @trees, @trees0;
    push @trees, @trees1;
    push @trees, @trees2;
    push @trees, @trees3;

    if(@trees){
        for(my $i = 0; $i < scalar(@trees); $i++){
            print "$i \t$trees[$i]\n";
        }
        print "\nThe following tree files were found. Press the corresponding number to select one or n for none:";
        my $input = <STDIN>;
        chomp $input;

        if($input ne 'n' || $input ne 'N'){
            $intree = $trees[$input];
        }
    }

    do{
        print "FILE OPTIONS\n";
        print "t \tUser input tree: \t\t\t$intree\n";
        print "a \tUser input alignment: \t\t\t$inalign\n";
        print "c \tMask Alignment Column Gap Cut-off: \t$mask_cut\n";
        print "s \tSubtree definition file: \t\t$subtree\n";
        print "o \tSubtree outfile root name: \t\t$outroot\n\n";
        print "GENERAL OPTIONS\n";
        print "p \tParallel Mode: \t\t\t\t" , ($parallel ? "yes" : "no"), "\n";
        if($parallel){
            print "N \t Number of CPUs: \t\t\t\t $Num_cpus\n";
        }
        print "P \tPhylogenetic Program: \t\t\t$program\n";
        print "m \tModel of Evolution: \t\t\t$model\n";
        print "M \tMidroot Mode: \t\t\t\t" , ($midroot ? "yes" : "no") , "\n";
        print "Quit [q], Help [h], Run [Y] or Selection (case-sensitive): ";

        my $sel = <STDIN>;
        chomp $sel;

        if($sel eq 'q'){
            exit();
        }elsif($sel eq 'Y'){
            $set = 1;
        }elsif($sel eq 't'){
            print "\nPlease enter a tree file: ";
            $intree = <STDIN>;
            chomp $intree;
            print "\n";
        }elsif($sel eq 'a'){
            print "\nPlease enter an alignment file: ";
            $inalign = <STDIN>;
            chomp $inalign;
            print "\n";
        }elsif($sel eq 'c'){
            print "Please enter new column gap score cut-off. 1 leaves alignment as is: ";
            $mask_cut = <STDIN>;
            chomp $mask_cut;
            if($mask_cut == 1){
                $mask_cut = 1.0;
            }
            if($mask_cut < 0.0 || $mask_cut > 1.0){
                print "Cutoff must be expressed as a number betwen 0 and 1\n";
                $mask_cut = 1.0;
            }
            print "\n";
        }elsif($sel eq 's'){
            print "\nPlease enter subtree definition file: ";
            $subtree = <STDIN>;
            chomp $subtree;
            print "\n";
        }elsif($sel eq 'o'){
            print "\nPlease subtree outfile root name: ";
            $outroot = <STDIN>;
            chomp $outroot;
            print "\n";
        }elsif($sel eq 'p'){
            print "Would you like to run in parallel mode (y/n)? ";
            my $choice = <STDIN>;
            chomp $choice;
            if($choice eq 'y'){
                $parallel = 1;
            }elsif($choice eq 'n'){
                $parallel = 0;
            }else{
                print "$choice is not a valid selection, keeping previous\n";
            }
		}elsif($sel eq 'O'){ 
            print "Would you like to perform Branch Length Optimization (y/n)? "; 
            my $choice = <STDIN>;
            chomp $choice;
            if($choice eq 'y'){
                $branch_len_opt = 'true';   
            }elsif($choice eq 'n'){
                $branch_len_opt = 'false';   
			}else{
				print "$choice is not a valid selection, keeping previous\n";  
            }
        }elsif($sel eq 'N'){
            print "Enter number of processors to use: ";
            my $choice = <STDIN>;
            chomp $choice;
            $Num_cpus = $choice;
            print "\n";
        }elsif($sel eq 'P'){
            print "Please choose program (raxml/qmmraxml/iqtree): ";
            my $choice = <STDIN>;
            chomp($choice);
            $choice = lc $choice;
            if($choice ne 'qmmraxml' && $choice ne 'raxml' && $choice ne 'iqtree'){
                print "Choice invalid. Leaving as default...\n";
            }else{
                $program = $choice;
                if($program eq 'qmmraxml'){
                    print "Setting number of gamme rates to 4\n";
                    $num_rates = 4;
        }elsif($program eq 'iqtree'){    ##########
          print "Setting relevant model\n";
                    $model = 'LG+C20+F+G';
                }else{
                    if($model eq 'LG'){
                        print "LG model of amino acid substitution only available with qmmraxml. Switching to JTT\n";   ####fix 
                        $model = 'JTT';
                    }
                }
            }
        }elsif($sel eq 'M'){
            print "Would you like to run in midpoint rooted mode (y/n)? ";
            my $choice = <STDIN>;
            chomp $choice;
            if($choice eq 'y'){
                $midroot = 1;
            }elsif($choice eq 'n'){
                $midroot = 0;
            }else{
                print "$choice is not a valid selection, keeping previous\n";
            }
        }elsif($sel eq 'h'){
            printHelp();
            print "Press any key to continue: ";
            my $key = <STDIN>;
        }elsif($sel eq 'r'){
            if($program eq 'qmmraxml'){
                print "Cannot change number of rate categories from 4 when using qmmraxml.\n";
            }else{
                print "Enter number of rate categories: \n";
                my $choice = <STDIN>;
                chomp $choice;
                if($choice !~ /[0-9]+/){
                    print "Invalid rate choice. Must enter a positive integer.\n";
                }elsif($choice == 0){
                    print "Invalid rate choice. Must enter a positive integer.\n";
                }else{
                    $num_rates = $choice;
                }
            }
        }elsif($sel eq 'm'){
            print "Model of protein sequence evolution JTT/WAG/LG/LG+C20+F+G: ";
            my $choice = <STDIN>;
            chomp $choice;
            $choice = uc $choice;
            if($choice ne 'JTT' && $choice ne 'WAG' && $choice ne 'LG' && $choice ne 'LG+C20+F+G'){   #####
                print "Invalid selection. Keeping previous choice\n";
            }elsif($choice eq 'LG' && $program eq 'puzzle'){
                print "LG model of sequence evolution is only usable with qmmraxml. Keeping previous\n";
            }else{
                $model = $choice;
            }
        }else{
            print "I'm sorry but $sel is not a valid selection, please try again\n\n";
        }

    }while(!$set);

    my $command = "-a $inalign -o $outroot -m $model -s $subtree -P $program -r $num_rates ";

    if($intree){
        $command = $command . "-t $intree ";
    }
    if($midroot){
        $command = $command . "-M ";
    }
    if($verbose){
        $command = $command . "-v ";
    }
    if($mask_cut < 1.0){
        $command = $command . "-c $mask_cut ";
    }
    if($parallel){
        $command = $command . "-p -N $Num_cpus "
    }
	if($branch_len_opt){
        $command = $command . "-O ";
    }

    open("LOG", ">>FunDi.log") or die "Could not open FunDi.log...\n\n";
    print "Running FunDi with specified options $command\n";
    print LOG "Running FunDi with specified options $command\n";
    close(LOG);

    @_ = split / /, $command;
    runCLIMode(@_);
}

sub printHelp{
    print <<HELP;

    Usage (Command Line Mode):        perl FunDi.pl [options] -t intree -a inalignment
    Usage (Phylip Style Menu Mode):   perl FunDi.pl

    Options:

    -h  Help:             Print this Help message and exit
    -o  Output root:      Sets the root filename for subtree files. Default: subtree
    -n  Normalization:    Defines the normalization method for site rate comparisons.
                          Default: Simple division by sum

                          Other options:
                          -sum_divide:  Divide by sum of rates
                          -log_ratio:   Ratio of logs of rates
                          -mean_divide: Divide by mean of rates

    -v  Verbose:          Run in verbose mode (Not yet Implemented. Always verbose)
    -p  Parallel:         Use Parallel Version of RAxML
    -s  Subtree:          Name of subtree definition file (Default: subtree.def)
    -M  Midpoint:         Run in midpoint mode rooted mode. (Default: unrooted tree)
    -l  Length:           Set length of simulated sequences (Only in non mixture mode. Default: 1000)
    -x  Mixture Model:    Run mixture model mode (Default)
    -m  Model:            Model of sequence evolution (JTT, WAG, LG, LG+C20+F+G) LG in QMMRAxML only (Default: LG+C20+F+G)
    -O  Optimization:     Run internal branchlength optimization? (yes/no). (Available in PUZZLE and IQ-Tree only)
    -r  Gamma rates:      Number of gamma rates. Selectable in PUZZLE only. (QmmRAxML always 4)
    -P  Program:          Phylogenetic program (puzzle/qmmraxml). QmmRAxML contains class frequency models (Default: qmmraxml)

HELP
}

sub create_new_tree_file{
	my ($tree_file, $branch_len, $step) = @_;
	
	#my $count;
	my @list_of_new_trees;
	if ($branch_len){
		my $branch_len_num_of_digits = length($branch_len);
		my $temp_branch_length = $branch_len - 1 > rem($branch_len,$step) ? $branch_len - 1 : rem($branch_len,$step);
		while ($temp_branch_length <= $branch_len + 1.001){
			if (length($temp_branch_length) > $branch_len_num_of_digits){
					$temp_branch_length = substr($temp_branch_length, 0, $branch_len_num_of_digits);
				}
			# Creating a class variable of PhyloTree_new from $tree_file
			my $new_tree = new PhyloTree_new(file => $tree_file);
			$new_tree -> read_tree();
			my @nodes = $new_tree -> get_nodes();
			my $replaced='false';
			# Searching for the node that holds the branch length and then replacing it with the modified branchlength
			foreach(@nodes){
				my $bl = $_ -> get_branchlength();
				if ($bl == $branch_len){
					$_ -> set_branchlength($temp_branch_length);
					$replaced='true';
				}
			}
			if ($replaced eq 'true'){ 
				print "Successfully replaced old branch length $branch_len with $temp_branch_length in $tree_file \n";
				my $subtree_out = "$tree_file" . "." . "$temp_branch_length" . ".nw"; 
				$new_tree -> print_tree($subtree_out, 0);
				push @list_of_new_trees, $subtree_out;
				$temp_branch_length = $temp_branch_length + $step;
			}
			else {
				die "Failed to replace branch length $branch_len in file $tree_file";
			}		
		}
	} else {
		die "No branch length given";
	}
	return @list_of_new_trees;
}

sub rem { $_[1]*frac($_[0]/$_[1]) }
sub frac { $_[0]-int($_[0]) }