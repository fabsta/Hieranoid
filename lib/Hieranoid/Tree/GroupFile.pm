package GroupFile;
=head1 GroupFile

hieranoid::Tree::GroupFile - Container of taxon objects

=head1 SYNOPSIS

 Abstract class

=head1 DESCRIPTION

Abstract class to handle how orthologous group are summarized

=head1 METHODS

=cut
use Moose;
use Carp;
use English;
use Log::Log4perl qw(:easy);
use Data::Dumper;
use List::MoreUtils qw(any all);
use File::Temp qw/ tempfile tempdir /;
use File::Copy;
use Bio::Phylo::IO;
use IO::String;
use Bio::TreeIO;
use Hieranoid::Tree::InnerNode;

## ATTRIBUTES


## SQL files: Inparanoid output
has 'orthologyPredictionFile', is => 'rw', isa => 'Str', required => 1; # InparanoidOutput
has 'originalOrthologyPredictionFile', is => 'rw', isa => 'Str';
## Prediction file for node : Homininae4: A_HUMAN, A_PAN
has 'groupsFile', is => 'rw', isa => 'Str'; # .groups.txt
has 'expandedGroupsFile', is => 'rw', isa => 'Str'; # .expandedGroups.txt
has 'OGTreeFile', is => 'rw', isa => 'Str'; # .OGTree.txt
has 'mappingsFile', is => 'rw', isa => 'Str'; # .mappings.txt 
has 'groupPredictionHash', is => 'rw', isa => 'HashRef[Str]';
has 'species', is => 'rw', isa => 'Str';


=head2 CONSTRUCTOR

=over

=item new()

Initialises a summarizer object

=cut
sub BUILD {
      my $self = shift;
      $self->originalOrthologyPredictionFile($self->orthologyPredictionFile."_original");
}

# how to add species information
#  save as 
# Homininae1:ENSP00000364178:Homo_sapiens:score,ENSPTRP00000021693:Pan_troglodytes:score

=item get_group2IDsHash()

to be written...

 Title   : get_group2IDsHash
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub get_group2IDsHash{
     #   print Dumper @_;
       my ($self,$innerNode,$og_assignments_hashref) = (@_);
        my $species = $self->species;
        my $og_counter = 0;
        

       DEBUG("\tget_group2IDsHash: reading from ".$self->groupsFile."\n");
        ## CHECK DEFINEDNESS
        croak q{One/several of the parameters for 'read_og_assignments' was/were not defined.\n}
        if any {!defined $_} $self->groupsFile;
        
        if(!-e $self->groupsFile || ! -s $self->groupsFile){
                WARN("\tThere is no prediction file ".$self->groupsFile." (but there should be one). Exiting\n");
                exit;
        }
        ### OPENING FILE
        open my $OG_FH, '<', $self->groupsFile or croak "Couldn't open ".$self->groupsFile.": $OS_ERROR";
        my $counter = 0;
        ### READING FROM FILE
        while (my $line = <$OG_FH>) {
                next if $line =~ /^$/;
                $line =~ /(\w+):(.*)/;
                my ($og,$members) = ($1,$2);
                if(!defined $og || !defined $members){
                        ERROR("Could not parse line $line in file ".$self->groupsFile.". Values are $og,$members\n");
                        #exit;
                }
                #print "$og,$members\n";
              #  if(defined $species && $og =~ /$species/){
                        #				print "og matches $species";
                        $og_assignments_hashref->{$og} = $members;
              #  }
              #  else{
              #          $og_assignments_hashref->{$og} = $members;
              #  }
                $og_counter++;
        }
        if(defined $species){
                DEBUG("\tread $og_counter groups for $species") ;
                #DEBUG("\tread $og_counter groups for $species") ;
        }
        else{
                #print("\tread $og_counter groups") ;
                #DEBUG("\tread $og_counter groups") ;
        }
        ### CLOSING FILE
        #print "\t\tRead og assignment file\n";
        #print "found ".keys(%{$og_assignments_hashref})." ids\n";
        close $OG_FH or croak "Couldn't close '".$self->groupsFile."': $OS_ERROR";
        if(!keys(%{$og_assignments_hashref})){
                ERROR("Could not get groupPredictions from ".$self->groupsFile." \n");
                exit;
        };
}

=item get_group2expandedIDsHash()

to be written...

 Title   : get_group2expandedIDsHash
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
# returns 
       # Hash:
               # Catarrhini3 => "CCNK_MACMU,Homininae2"
               # Catarrhini2 => MmuSTS.4752.1.S1_at_MACMU,Homininae1"

sub get_group2expandedIDsHash{
      my ($self,$innerNode,$og_assignments_hashref) = (@_);
      my $species = $self->species;
      my $og_counter = 0;
      if(! -e $self->expandedGroupsFile || ! -s $self->expandedGroupsFile){
                ERROR("\tCould not read resolved mappings ".$self->expandedGroupsFile."\n");
                return 0;
        }

        DEBUG("\tget_group2expandedIDsHash: reading from ".$self->expandedGroupsFile."\n");
        ## CHECK DEFINEDNESS
        croak q{One/several of the parameters for 'read_og_assignments' was/were not defined.\n}
        if any {!defined $_} $self->expandedGroupsFile;

        ### OPENING FILE
        open my $OG_FH, '<', $self->expandedGroupsFile or croak "Couldn't open ".$self->expandedGroupsFile.": $OS_ERROR";
        my $counter = 0;
        ### READING FROM FILE
        while (my $line = <$OG_FH>) {
                next if $line =~ /^$/;
                $line =~ /(\w+):(.*)/;
                my ($og,$members) = ($1,$2);
                if(!defined $og || !defined $members){
                        ERROR("Could not parse line $line in file ".$self->expandedGroupsFile.". Values are $og,$members\n");
                        #exit;
                }
                        $og_assignments_hashref->{$og} = $members;
                $og_counter++;
        }
        close $OG_FH or croak "Couldn't close '".$self->expandedGroupsFile."': $OS_ERROR";
        if(!keys(%{$og_assignments_hashref})){
                ERROR("Could not get groupPredictions from ".$self->expandedGroupsFile." \n");
                exit;
        };
}


sub get_treeString4ID{
      my ($self,$innerNode,$og_assignments_hashref) = (@_);
      my $species = $self->species;
      my $og_counter = 0;
      if(! -e $self->OGTreeFile || ! -s $self->OGTreeFile){
                ERROR("\tCould not read resolved mappings ".$self->OGTreeFile."\n");
                return 0;
        }
		#print "\tstarted reading ".$self->OGTreeFile."\n";
        DEBUG("\tget_group2expandedIDsHash: reading from ".$self->OGTreeFile."\n");
        ## CHECK DEFINEDNESS
        croak q{One/several of the parameters for 'read_og_assignments' was/were not defined.\n}
        if any {!defined $_} $self->OGTreeFile;

        ### OPENING FILE
        open my $OG_FH, '<', $self->OGTreeFile or croak "Couldn't open ".$self->OGTreeFile.": $OS_ERROR";
        my $counter = 0;
        ### READING FROM FILE
        while (my $line = <$OG_FH>) {
                next if $line =~ /^$/;
                chomp($line);
                # Using bitscores on inner nodes
                # ((ENSPTRP00000031469:1.000,ENSP00000358400:1.000)Homininae1:11100)root;
                
                # Using inparalog scores on inner nodes
                # ((ENSPTRP00000031469:1.000,ENSP00000358400:1.000)Homininae1:1.000)root;
                $line =~ /HieranoidOG-\w+:\((.*\))(\w+):(\d+\.?\d+)\)root;/;
                my ($treeString,$ogName,$bitscore) = ($1,$2,$3);
                #print "\tline: $line  -> $treeString,$ogName\n";
                if(!defined $treeString || !defined $treeString){
                        ERROR("Could not parse line $line in file ".$self->OGTreeFile.". Values are $treeString,$ogName\n");
                        print "Could not parse line $line in file ".$self->OGTreeFile.". Values are $treeString,$ogName\n";
                        exit;
                }
                		#print "\t$ogName -> treestring = $treeString\n";
                		#print "\t$ogName -> bitscore = $bitscore\n";
                        $og_assignments_hashref->{$ogName}{'treestring'} = $treeString;
                        $og_assignments_hashref->{$ogName}{'bitscore'} = $bitscore;
                $og_counter++;
        }
        close $OG_FH or croak "Couldn't close '".$self->OGTreeFile."': $OS_ERROR";
        if(!keys(%{$og_assignments_hashref})){
                ERROR("Could not get groupPredictions from ".$self->OGTreeFile." \n");
                exit;
        };
}

sub get_treeleafNodes
{
	my ($self,$innerNode,$og_assignments_hashref) = (@_);
    my $species = $self->species;
    
    my $cat_cmd = "cat ".$innerNode->orthologGroups->OGTreeFile;
	my @all_current_ogs = `$cat_cmd`;
	foreach my $string(@all_current_ogs){
		$string =~ /HieranoidOG-\w*:(.*)/;
        $string = $1;
		my $io = IO::String->new($string);
		my $treeio = Bio::TreeIO->new(-fh => $io,
                              -format => 'newick');
    	my $tree = $treeio->next_tree; # we'll assume it worked for demo purposes
                                   # you might want to test that it was defined

    	my $rootnode = $tree->get_root_node;
		#print "root note: ".$rootnode->id."\n";
    	# process just the next generation
    	my @allDescendents = $rootnode->each_Descendent();
    	my $realroot = $allDescendents[0];
    	 
    	foreach my $node ( $realroot->each_Descendent() ) {
        	#print "branch len is ", $node->id, "\n";
        	$og_assignments_hashref->{$realroot->id} .= $node->id.","; 
    	}
	}
}

sub get_alltreeleafNodes
{
	my ($self,$innerNode,$og_assignments_hashref) = (@_);
    my $species = $self->species;
    #print "\tidentifying treegroup IDs\n";
    my $cat_cmd = "cat ".$innerNode->orthologGroups->OGTreeFile;
	my @all_current_ogs = `$cat_cmd`;
	foreach my $string(@all_current_ogs){
		#print "treestring: $string\n";
		$string =~ /HieranoidOG-\w*:(.*)/;
        $string = $1;
		my $io = IO::String->new($string);
		my $treeio = Bio::TreeIO->new(-fh => $io,
                              -format => 'newick');
    	my $tree = $treeio->next_tree; # we'll assume it worked for demo purposes
                                   # you might want to test that it was defined
		#my @leaves = $tree->get_leaf_nodes;
    	my $rootnode = $tree->get_root_node;
    	my $sortby = 'height';
		my @nodes = $rootnode->get_all_Descendents($sortby);
		my $realroot = $nodes[0];
    	my $realroot2 = $nodes[1];
    	my $realroot3 = $nodes[2];
    	$og_assignments_hashref->{$realroot->id} .= $realroot2->id.","; 
    	$og_assignments_hashref->{$realroot->id} .= $realroot3->id.","; 
    	
    	#print "realroot: ".$realroot->id.":".$realroot2->id."\t".$realroot3->id." \n";
	}
}


=item getIDmappingsHash()

to be written...

 Title   : getIDmappingsHash
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
# returns 
       # Hash:
               # Catarrhini3 => "CCNK_MACMU,Homininae2"
               # Catarrhini2 => MmuSTS.4752.1.S1_at_MACMU,Homininae1"

sub getIDmappingsHash{
      my ($self,$innerNode,$og_assignments_hashref) = (@_);
      my $species = $self->species;
      my $og_counter = 0;
      if(! -e $self->mappingsFile || ! -s $self->mappingsFile){
                ERROR("\tCould not read resolved mappings ".$self->mappingsFile."\n");
                return 0;
        }

        DEBUG("\tget_group2expandedIDsHash: reading from ".$self->mappingsFile."\n");
        ## CHECK DEFINEDNESS
        croak q{One/several of the parameters for 'read_og_assignments' was/were not defined.\n}
        if any {!defined $_} $self->mappingsFile;

        ### OPENING FILE
        open my $OG_FH, '<', $self->mappingsFile or croak "Couldn't open ".$self->mappingsFile.": $OS_ERROR";
        my $counter = 0;
        ### READING FROM FILE
        while (my $line = <$OG_FH>) {
                next if $line =~ /^$/;
                $line =~ /(\w+):(.*)/;
                my ($og,$members) = ($1,$2);
                if(!defined $og || !defined $members){
                        ERROR("Could not parse line $line in file ".$self->mappingsFile.". Values are $og,$members\n");
                        exit;
                }
                        $og_assignments_hashref->{$og} = $members;
                $og_counter++;
        }
        close $OG_FH or croak "Couldn't close '".$self->mappingsFile."': $OS_ERROR";
        if(!keys(%{$og_assignments_hashref})){
                ERROR("Could not get groupPredictions from ".$self->mappingsFile." \n");
                exit;
        };
}

=item getAllUsedIDS()

to be written...

 Title   : getAllUsedIDS
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
        # returns 
        # Hash:
                # Catarrhini2 => "MmuSTS.4752.1.S1_at_MACMU,VPS45_PANTR,VPS45_HUMAN"
                # Catarrhini3 => "CCNK_MACMU,CCNK_F3_PANTR,CCNK_HUMAN"

sub getAllUsedIDS{
       my ($self,$innerNode,$og_assignments_hashref) = (@_);
       my $species = $self->species;
       my $og_counter = 0;
        DEBUG("\tgetAllUsedIDS: reading from ".$self->groupsFile."\n");
                ## CHECK DEFINEDNESS
        croak q{One/several of the parameters for 'read_og_assignments' was/were not defined.\n}
        if any {!defined $_} $self->groupsFile;
                if(!-e $self->groupsFile || ! -s $self->groupsFile){
                        ERROR("\tThere is no prediction file ".$self->groupsFile." (but there should be one). Exiting\n");
                        exit;
        }
                ### OPENING FILE
        open my $OG_FH, '<', $self->groupsFile or croak "Couldn't open ".$self->groupsFile.": $OS_ERROR";
        my $counter = 0;
        ### READING FROM FILE
        while (my $line = <$OG_FH>) {
                next if $line =~ /^$/;
                $line =~ /(\w+):(.*)/;
                my ($og,$members) = ($1,$2);
                if(!defined $og || !defined $members){
                        ERROR("Could not parse line $line in file ".$self->groupsFile.". Values are $og,$members\n");
                        exit;
                }
                my @members_array = split(",",$members);
                if(scalar(@members_array) < 2){
                        ERROR("\t\tWarning: line $line could not be parsed ($og:$members)\n");
                        next;
                }
                $og_assignments_hashref->{$og} = 1;
                foreach(split(",",$members)){
                        $og_assignments_hashref->{$_} = 1;
                }
                $og_counter++;
        }
        if(defined $species){
                #print("\tread $og_counter groups for $species") ;
                #DEBUG("\tread $og_counter groups for $species") ;
        }
        else{
                #print("\tread $og_counter groups") ;
                #DEBUG("\tread $og_counter groups") ;
        }
        ### CLOSING FILE
        #print "\t\tRead og assignment file\n";
        #print "found ".keys(%{$og_assignments_hashref})." ids\n";
        close $OG_FH or croak "Couldn't close '".$self->groupsFile."': $OS_ERROR";
        if(!keys(%{$og_assignments_hashref})){
                ERROR("Could not get groupPredictions from ".$self->groupsFile." \n");
                exit;
        };
}
=item convert2expandedGroupFile()

to be written...

 Title   : convert2expandedGroupFile
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut

sub convert2expandedGroupFile{
        my ($self,$innerNode) = (@_);
        # workflow
        DEBUG("Converting 2 resolved group file");
    # read all group prediction
        my %allGroupPredictions;
        $innerNode->orthologGroups->get_group2IDsHash($innerNode,\%allGroupPredictions);
    # read mapping file
        my %allGroupIDMappings;
        if($innerNode->leftDaughterType eq 'innerNode'){
                $innerNode->leftDaughter->orthologGroups->get_group2expandedIDsHash($innerNode->leftDaughter,\%allGroupIDMappings);
        }
        if($innerNode->rightDaughterType eq 'innerNode'){
                $innerNode->rightDaughter->orthologGroups->get_group2expandedIDsHash($innerNode->rightDaughter,\%allGroupIDMappings);
        }
        #print Dumper %allGroupIDMappings;
        # Could not find orthologous groups, but expected to find (one daughter is an inner node)
        if(!keys(%allGroupIDMappings) &&
                ($innerNode->leftDaughterType eq 'innerNode' || $innerNode->rightDaughterType eq 'innerNode')){
                        ERROR("Could not find orthologous groups, but expected to find (one daughter is an inner node)");
                        exit;
                }
        #print Dumper %allGroupIDMappings;
        # iterate over group predictions
        my $allResolvedGroupsString;
        foreach my $og(sort keys(%allGroupPredictions)){
                #print "og: $og: ".$allGroupPredictions{$og}."\n";
                my $copy_og = 1;
                my @resolved_ids;
                #       resolve them
                if(exists $allGroupIDMappings{$og}){
                        push(@resolved_ids, $allGroupIDMappings{$og});
                }
                else{
                        foreach my $curr_id(split(/,/,$allGroupPredictions{$og})){
                                #print "\tcurr_id: $curr_id: ".$allGroupPredictions{$og}."\n";
                                if(exists $allGroupIDMappings{$curr_id}){
                                        push(@resolved_ids, $allGroupIDMappings{$curr_id});
                                }
                                else{
                                        #print "\tError: Could not find ID $curr_id\n";
                                        push(@resolved_ids, $curr_id);
                                }
                                #print "\t\tresolved: ".join(",",@resolved_ids)."\n";
                        }
                }
                $allResolvedGroupsString .= "$og:".join(",",@resolved_ids)."\n";
        }
         #    # Converted Orthology predictions        
                     write_to_file({file_name => $innerNode->orthologGroups->expandedGroupsFile, text => $allResolvedGroupsString});
#exit;
        return 1;
}
=item write_to_file()

to be written...

 Title   : write_to_file
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub write_to_file{
	#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $file_name = $arg_ref->{file_name};
	my $text      = $arg_ref->{text};
	### OPENING FILE
	open my $out, '>', $file_name
	  or croak "Couldn't open '$file_name': $OS_ERROR";
	### Writing file
	print {$out} $text;
	### CLOSING FILE
	close $out or croak "Couldn't close '$file_name': $OS_ERROR";
	return 1;
}


1;
