package OrthoXML;
=head1 OrthoXML

hieranoid::Tree::OrthoXML - Container of taxon objects

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
use Hieranoid::Tree::InnerNode;

## ATTRIBUTES


## SQL files: Inparanoid output
has 'orthologyPredictionFile', is => 'rw', isa => 'Str', required => 1; # InparanoidOutput
has 'originalOrthologyPredictionFile', is => 'rw', isa => 'Str';
## Prediction file for node : Homininae4: A_HUMAN, A_PAN
has 'groupsFile', is => 'rw', isa => 'Str'; # .groups.txt
has 'expandedGroupsFile', is => 'rw', isa => 'Str'; # .expandedGroups.txt
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
                        exit;
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
                        exit;
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