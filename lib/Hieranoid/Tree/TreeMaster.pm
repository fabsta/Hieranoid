package TreeMaster;
=head1 InnerNode

hieranoid::Tree::InnerNode - Container of taxon objects

=head1 SYNOPSIS

 Abstract class

=head1 DESCRIPTION

Abstract class to handle how orthologous group are summarized

=head1 METHODS

=cut
use Moose;
use strict;
use warnings;
use Bio::Phylo::IO qw(parse unparse);
#use Bio::Phylo::IO;
use Data::Dumper;
use Log::Log4perl qw(:easy);
use Hieranoid::Tree::InnerNode;
use Hieranoid::Tree::LeafNode;
use Hieranoid::Tree::GroupFile;
use Hieranoid::FileInformation;


## ATTRIBUTES
has 'currentTreeNode', is => 'ro', isa => 'Str', default => '';
has 'configuration', is => 'ro', isa => 'Object';
# DIRECTORIES
has 'resultsDirectory', is => 'ro', isa => 'Str';
has 'speciesFolder', is => 'ro', isa => 'Str';
has 'speciesTreeFile', is => 'rw', isa => 'Str';
has 'speciesTreeObject', is => 'rw', isa => 'Object'; 
has 'comparisons', is => 'rw', isa => 'ArrayRef[Object]' ; 
has 'comparisonsNames', is => 'rw', isa => 'ArrayRef[Str]' ; 
has 'terminalNodes', is => 'rw', isa => 'ArrayRef[Str]'; 
## ATTRIBUTES - END

=head2 CONSTRUCTOR

=over

=item new()

Initialises a summarizer object

        Reads a tree
        Checks taxa
        Gets order of comparisons

=cut
sub BUILD {
      my $self = shift;
      $self->read_species_tree($self->speciesTreeFile);
      $self->check_existence_of_taxa();
      
      # Make tree bifurcating
      #if($self->speciesTreeObject->is_binary){
       #           print "\t\tSpecies tree is binary\n";
        #  }
        #  else{
        #          $self->resolveTree();               
        #  }
      #exit;
      $self->get_comparison_order();
}
=item getNextComparison()

to be written...

 Title   : getNextComparison
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub getNextComparison{
        print "Not yet implemented\n";
        
}
=item getNextComputableNode()

to be written...

 Title   : getNextComputableNode
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub getNextComputableNode{}
=item read_species_tree()

to be written...

 Title   : read_species_tree
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut

sub read_species_tree{
        my ($self,$tree_file) = (@_);
        if(!-e $tree_file || !-s $tree_file){
                print "\tError reading tree file ($tree_file)\n";
                exit;
        }
        # Read tree
        my $tree_object = Bio::Phylo::IO->parse(
                '-file' => $tree_file,
                '-format' => 'newick')->first;
                $self->speciesTreeObject($tree_object);
                my @terminals_nodes = @{ $tree_object->get_terminals };
                my @taxa_array;
                foreach(@terminals_nodes){
                        # $_;
                        next if ! defined $_;
                        push(@taxa_array,$_->get_name());
                }
                $self->terminalNodes(\@taxa_array);

                return 1;
}
=item check_existence_of_taxa()

to be written...

 Title   : check_existence_of_taxa
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut

sub check_existence_of_taxa{
                my $self = shift;
                my %missing_species_files = ();
                foreach my $species_name(@{$self->terminalNodes}){
                        #my $species_name = $_->get_name;
                        my $species_file_to_check = $self->configuration->speciesFilesDirectory."/".$species_name;

                        # my $species_file_to_check = ();
                        # add corresponding suffix
                        $species_file_to_check .= ($self->configuration->sequenceInputFormat eq "xml")? ".xml" : ".fa";
                        if(!-e $species_file_to_check || ! -s $species_file_to_check){
                                $missing_species_files{$species_name} = 1;
                        }
                }	
                if(keys(%missing_species_files)){
                        print "\tThe following taxa are missing in folder ".$self->configuration->speciesFilesDirectory.": ".join(",",keys(%missing_species_files))."\n";
                        return 0;
                }
                else{
                        return 1;
                }
}
=item get_comparison_order()

to be written...

 Title   : get_comparison_order
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
# list of inner nodes with daughters
sub get_comparison_order{
        my $self = shift;
        my $counter = 0;
        my (@order_of_comparisons,@order_of_comparisons_names, @parallel_computation_nodes_array);
        ($self->speciesTreeObject)->visit_depth_first(
                -post_daughter => sub {
                        my $current_node = shift;
                        return if $current_node eq "";
                        return if !defined $current_node;

                        my $current_node_name = $current_node->get_name;
                        return if !defined $current_node_name;
                        my @terminals = @{ $current_node->get_terminals };
                        my @children = @{ $current_node->get_children };
                        if(@children){
                                #make innernode object
                                #DEBUG("checking node $current_node_name\n");
                            if(!$current_node_name){
                                    print "\tCould not find label for node. Exiting\n";
                             #       exit;
                            }
                            
                            
                            if(!$current_node->get_children()){
                                    print "InnerNode $current_node_name has no children\n";
                                    return 0;
                            }
                    # left
                            my $outputDirectory4Node = $self->configuration->hieranoidResultsDirectory."/".$current_node_name;
                # Set node name
                            if( $current_node_name eq ''){
                                    print "setting node to ".$current_node->get_first_daughter->get_name."_".$current_node->get_last_daughter->get_name."\n";
                                    $self->name($current_node->get_first_daughter->get_name."_".$current_node->get_last_daughter->get_name)
                            }
                            my $innerNode = InnerNode->new(
                                    nodeObject => $current_node,
                                    name => $current_node_name,
                                    configuration => $self->configuration);
                            
                            
                            if(!$current_node->get_first_daughter()){
                                    ERROR("InnerNode $current_node_name has no left daughter\n");
                                    return 0;
                            }
                            else{
                                    #LEFT DAUGHTER  
                                    my $leftDaughterObject = $current_node->get_first_daughter;
                                    #  my $leftDaughterObjectName = $leftDaughterObject->get_name;
                                       #DEBUG("\tleft: ".$leftDaughterObject->get_name."\n");
                                    if($leftDaughterObject->is_terminal){
                                            $innerNode->leftDaughter(LeafNode->new(
                                                    nodeObject => $leftDaughterObject, 
                                                    name => $leftDaughterObject->get_name,
                                                    fileInformation => (FileInformation->new(
                                                            sequenceSearchInputFile => $outputDirectory4Node."/".$leftDaughterObject->get_name.".".$self->configuration->sequenceInputFormat,
                                                            outputDirectory => $outputDirectory4Node."/".$leftDaughterObject->get_name,
                                                            speciesFilesDirectory => $self->configuration->speciesFilesDirectory,
                                                            profileDirectory => $self->configuration->hieranoidProfilesDirectory,
                                                            profileFile => $self->configuration->hieranoidProfilesDirectory."/".$leftDaughterObject->get_name,
                                                            hmmFile => $self->configuration->hieranoidProfilesDirectory."/".$leftDaughterObject->get_name."_hhm_db",
                                                            csdbFile => $self->configuration->hieranoidProfilesDirectory."/".$leftDaughterObject->get_name.".cs219",
                                                            consensusDirectory => $self->configuration->hieranoidConsensusDirectory,
                                                            consensusFile => $self->configuration->hieranoidConsensusDirectory."/".$leftDaughterObject->get_name.".cons.".$self->configuration->sequenceInputFormat,
                                                            alignmentFile => $self->configuration->hieranoidConsensusDirectory."/".$leftDaughterObject->get_name.".aln.".$self->configuration->sequenceInputFormat,
                                                            orthologGroupsFormat => $self->configuration->orthologGroupsFormat,
                                                            sequenceInputFormat => $self->configuration->sequenceInputFormat,
                                                            sourceFile => $self->configuration->speciesFilesDirectory."/".$leftDaughterObject->get_name.".".$self->configuration->sequenceInputFormat))
                                                            ));
                                                            $innerNode->leftDaughterType("leaf");
                                                    }
                                                    else{
                                                            $innerNode->leftDaughter(InnerNode->new(
                                                                    nodeObject => $leftDaughterObject, 
                                                                    name => $leftDaughterObject->get_name,
                                                                    configuration => $self->configuration));
                                                                    $innerNode->leftDaughterType("innerNode");
                                                            }
                                    }        

                # RIGHT DAUGHTER
                            if(!$current_node->get_last_daughter()){
                                    ERROR("InnerNode $current_node_name has no right daughter\n");
                                    return 0;
                            }
                            else{
                                    my $rightDaughterObject = $current_node->get_last_daughter;
                                    #my $rightDaughterObjectName = $rightDaughterObject->get_name;
                                    #DEBUG("\tright: ".$rightDaughterObject->get_name."\n");
                                    if($rightDaughterObject->is_terminal){
                                            $innerNode->rightDaughter(LeafNode->new(
                                                    nodeObject => $rightDaughterObject, 
                                                    name => $rightDaughterObject->get_name,
                                                    fileInformation => (FileInformation->new(
                                                            sequenceSearchInputFile => $outputDirectory4Node."/".$rightDaughterObject->get_name.".".$self->configuration->sequenceInputFormat,
                                                            outputDirectory => $outputDirectory4Node."/".$rightDaughterObject->get_name,
                                                            speciesFilesDirectory => $self->configuration->speciesFilesDirectory,
                                                            profileDirectory => $self->configuration->hieranoidProfilesDirectory,
                                                            profileFile => $self->configuration->hieranoidProfilesDirectory."/".$rightDaughterObject->get_name,
                                                            hmmFile => $self->configuration->hieranoidProfilesDirectory."/".$rightDaughterObject->get_name."_hhm_db",
                                                            csdbFile => $self->configuration->hieranoidProfilesDirectory."/".$rightDaughterObject->get_name.".cs219",
                                                            consensusDirectory => $self->configuration->hieranoidConsensusDirectory,
                                                            consensusFile => $self->configuration->hieranoidConsensusDirectory."/".$rightDaughterObject->get_name.".cons.".$self->configuration->sequenceInputFormat,
                                                            alignmentFile => $self->configuration->hieranoidConsensusDirectory."/".$rightDaughterObject->get_name.".aln.".$self->configuration->sequenceInputFormat,
                                                            orthologGroupsFormat => $self->configuration->orthologGroupsFormat,
                                                            sequenceInputFormat => $self->configuration->sequenceInputFormat,
                                                            sourceFile => $self->configuration->speciesFilesDirectory."/".$rightDaughterObject->get_name.".".$self->configuration->sequenceInputFormat
                                                            ))
                                                                    ));

                                            $innerNode->rightDaughterType("leaf");
                                    }
                                    else{
                                            $innerNode->rightDaughter(InnerNode->new(
                                                    nodeObject => $rightDaughterObject, 
                                                        name => $rightDaughterObject->get_name,
                                                        configuration => $self->configuration));
                                            $innerNode->rightDaughterType("innerNode");
                                    }
                            }
                            $innerNode->childrenType($innerNode->leftDaughterType."_".$innerNode->rightDaughterType);
                   #    print Dumper $innerNode;
                   #    exit;
                                        unshift(@order_of_comparisons, $innerNode);
                                        unshift(@order_of_comparisons_names, $current_node_name);
                        }
                                if(scalar(@children) == scalar(@terminals)){
                                        push(@parallel_computation_nodes_array, $current_node_name);
                                }
                        }
                        );
                        $self->comparisons(\@order_of_comparisons);
                        $self->comparisonsNames(\@order_of_comparisons_names);
}

=item resolveTree()

# Dealing with multi-furcations
# e.g. Multi-furcation with A,B,C
# then do
# A-
 Title   : resolveTree
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut

sub resolveTree{
        my $self = shift;
        my $counter = 0;

        # Workflow
        # Identify all multifurcating nodes
        my @multifurcatingNodes;
        my @internals = @{ $self->speciesTreeObject->get_internals };
        if(!@internals){
                print "\tCould not get internal nodes from species tree. Exiting\n";
                exit;
        }
        foreach(@internals){
                my @children = @{ $_->get_children };
                push(@multifurcatingNodes,$_) if(scalar(@children) > 2);
        }
        if(!scalar(@multifurcatingNodes)){
                print "\tDid not find multifurcating nodes. Continue...\n";
                return 1;
        }
        # Now try to resolve the single nodes
        foreach(@multifurcatingNodes){
                $self->resolveNode($_);
        }
        
        # Iterate over nodes
        


}


sub resolveNode{
      my ($self,$Node2resolve) = (@_);
      my @children = @{ $Node2resolve->get_children };
      print "\tTrying to resolve node: ".$Node2resolve->get_name." with children ".join(",",map{$_->get_name}@children)."\n";
      my ($firstDaughter,$lastDaughter);
      for(my $i=0;$i<scalar(@children);$i=$i+2){
              # test if second Node exists
              # e.g. multifurcation is A,B,C
              # i = 0 => combine A,B 
              # i = 1 => C alone -> can keep this node
              if(!defined $children[$i+1]){
                      print "\tSingle node ".$children[$i]->get_name."\n";
                      $lastDaughter = $children[$i];
                      next;
              }
              if($i==2){
                      ;
              }
              
              my ($firstNode,$secondNode) = ($children[$i],$children[$i+1]);
              print "\tnew node with daughters ".$firstNode->get_name." and ".$secondNode->get_name."\n";
              
              # Make a new node for the two daughters
              my $newInnerNode = Bio::Phylo::Forest::Node->new(
                                -parent => $Node2resolve,
                                #-taxon => $firstNode."_".$secondNode,
                                -first_daughter  => $firstNode,
                                -last_daughter   => $secondNode,
                                -name            => $firstNode->get_name."_".$secondNode->get_name,
                                );
                                #bless $newInnerNode,'Bio::Phylo::Forest::Node';
                                print "name: ".$newInnerNode->get_name."\n";
                                print "parent: ".$newInnerNode->get_parent."\n";
                                print "first D: ".$newInnerNode->get_first_daughter."\n";
                                print "last D: ".$newInnerNode->get_last_daughter."\n";
                                
                                #exit;
                                if($i == 0){
                                        $firstDaughter = $newInnerNode;
                                }
                                else{
                                        $lastDaughter = $newInnerNode;
                                }
                # set first_daughter parameter from Node2resolve
      }
      print "\t\tset first to ".$firstDaughter->get_name."\n";
              $Node2resolve->set_first_daughter($firstDaughter);
      print "\t\tset last to ".$lastDaughter->get_name."\n";
              $Node2resolve->set_last_daughter($lastDaughter);
      
      return 1;   
}

sub get_outgroup_node{
        my ($self, $innerNode) = (@_);
        my $previous_sister;
        if(defined($previous_sister = $innerNode->nodeObject->get_previous_sister)){
                #print "\thave previous sister\n";
                if(defined(my $previous_sisters_left_daughter = $previous_sister->get_first_daughter)){
                        return $previous_sisters_left_daughter;
                }
                elsif(defined(my $previous_sisters_right_daughter = $previous_sister->get_last_daughter)){
                        return $previous_sisters_right_daughter;
                }
                else{
                        return $previous_sister;
                }
        }
        elsif(defined($previous_sister = $innerNode->nodeObject->get_next_sister)){
                #print "\thave next sister\n";
                if(defined(my $previous_sisters_left_daughter = $previous_sister->get_first_daughter)){
                        return $previous_sisters_left_daughter;
                }
                elsif(defined(my $previous_sisters_right_daughter = $previous_sister->get_last_daughter)){
                        return $previous_sisters_right_daughter;
                }
                else{
                        return $previous_sister;
                }
        }
        else{
                #print "\tCould not find outgroup for this node. Skipping\n";
                return 0;
        }
        return 1;
}

1;