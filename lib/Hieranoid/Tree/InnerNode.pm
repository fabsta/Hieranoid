package InnerNode;
=head1 InnerNode

hieranoid::Tree::InnerNode - Container of taxon objects

=head1 SYNOPSIS

 Abstract class

=head1 DESCRIPTION

Abstract class to handle how orthologous group are summarized

=head1 METHODS

=cut
use Moose;
with 'Daughter';
use strict;
use warnings;
use Carp;
use English;
use List::MoreUtils qw(any all);
use Data::Dumper;
use Log::Log4perl qw(:easy);
use File::Basename;
use Hieranoid::Tree::Daughter;
use Hieranoid::Tree::GroupFile;
use Hieranoid::FileInformation;



#extends 'Daughter';

## ATTRIBUTES
has 'leftDaughter', is => 'rw', isa => 'Object';
has 'rightDaughter', is => 'rw', isa => 'Object';
##### Outgroup
has 'outgroupDaughter', is => 'rw', isa => 'Object';
### can be removed
has 'rightDaughterType', is => 'rw', isa => 'Str';
has 'leftDaughterType', is => 'rw', isa => 'Str'; # leaf || innerNode
has 'outgroupDaughterType', is => 'rw', isa => 'Str'; # leaf || innerNode
has 'childrenType', is => 'rw', isa => 'Str'; 
#Objects 
has 'orthologGroups', is => 'rw', isa => 'Object';
has 'configuration', is => 'rw', isa => 'Object';
has 'outputDirectory4Node', is => 'rw', isa => 'Str';
has 'sequencesByID', is => 'rw', isa => 'HashRef[Str]';  # ID2 => "HCYAHDTSKKSPF";
=head2 CONSTRUCTOR

=over

=item new()

Initialises a summarizer object

=cut
sub BUILD {
        my $self = shift;
        #print Dumper $self;
        # SET ORTHOLOGY GROUPS
        my $outputDirectory4Node = $self->configuration->hieranoidResultsDirectory."/".$self->name;
        if($self->configuration->orthologGroupsFormat eq 'groupFile'){
                $self->orthologGroups(GroupFile->new(
                        orthologyPredictionFile => $outputDirectory4Node."/sqltable.".$self->name,
                        originalOrthologyPredictionFile => $outputDirectory4Node."/sqltable.".$self->name."_original",
                        expandedGroupsFile => $outputDirectory4Node."/".$self->name.".expandedGroups.txt",
                        OGTreeFile => $outputDirectory4Node."/".$self->name.".OGTree.txt",
                        groupsFile => $outputDirectory4Node."/".$self->name.".groups.txt",
                        mappingsFile => $outputDirectory4Node."/".$self->name.".mappings.txt"));
                        
        }
        else{
        	$self->orthologGroups();
                # OrthoXML Format - not implemented yet
                print "\tOrthoXML format not implemented yet\n";
                exit;
                ERROR("\tOrthoXML format not implemented yet\n");
        }
        my $sequenceSearchInputFile;
        if($self->nodeObject->get_parent){
                $sequenceSearchInputFile = $self->configuration->hieranoidResultsDirectory."/".$self->nodeObject->get_parent->get_name."/".$self->name.".fa";
              #  print "sequence search: $sequenceSearchInputFile\n";
        }
        else{
                $sequenceSearchInputFile = $self->configuration->hieranoidResultsDirectory."/".$self->name."/".$self->name.".fa";
        }
                #print Dumper $self;
                # SET FILE INFORMATIONS
        $self->fileInformation(FileInformation->new(
                        sourceFile => $self->configuration->hieranoidResultsDirectory."/".$self->name."/".$self->name.".fa",
                        sequenceSearchInputFile => $sequenceSearchInputFile,
                        outputDirectory => $self->configuration->hieranoidResultsDirectory."/".$self->name,
                        speciesFilesDirectory => $self->configuration->speciesFilesDirectory,
                        profileDirectory => $self->configuration->hieranoidProfilesDirectory,
                        profileFile => $self->configuration->hieranoidProfilesDirectory."/".$self->name,
                        hmmFile => $self->configuration->hieranoidProfilesDirectory."/".$self->name."_hhm_db",
                        csdbFile => $self->configuration->hieranoidProfilesDirectory."/".$self->name.".cs219",
                        consensusDirectory => $self->configuration->hieranoidConsensusDirectory,
                        consensusFile => $self->configuration->hieranoidConsensusDirectory."/".$self->name.".cons.fa",
                        alignmentFile => $self->configuration->hieranoidConsensusDirectory."/".$self->name.".aln.fa",
                        orthologGroupsFormat => $self->configuration->orthologGroupsFormat,
                        sequenceInputFormat => $self->configuration->sequenceInputFormat));
}
=item get_terminalSpeciesBelonging2Node()

to be written...

 Title   : get_terminalSpeciesBelonging2Node
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut


sub setOutgroup{
          my ($self,$outgroupObject) = (@_);
          my $outputDirectory4Node = $self->fileInformation->outputDirectory;
          if($outgroupObject->is_terminal){
                      $self->outgroupDaughter(LeafNode->new(
                              nodeObject => $outgroupObject, 
                              name => $outgroupObject->get_name,
                              fileInformation => (FileInformation->new(
                                      sequenceSearchInputFile => $self->configuration->speciesFilesDirectory."/".$outgroupObject->get_name.".".$self->configuration->sequenceInputFormat,
                                      #sequenceSearchInputFile => $outputDirectory4Node."/".$outgroupObject->get_name.".".$self->configuration->sequenceInputFormat,
                                      outputDirectory => $outputDirectory4Node."/".$outgroupObject->get_name,
                                      speciesFilesDirectory => $self->configuration->speciesFilesDirectory,
                                      profileDirectory => $self->configuration->hieranoidProfilesDirectory,
                                      profileFile => $self->configuration->hieranoidProfilesDirectory."/".$outgroupObject->get_name,
                                      hmmFile => $self->configuration->hieranoidProfilesDirectory."/".$outgroupObject->get_name."_hhm_db",
                                      csdbFile => $self->configuration->hieranoidProfilesDirectory."/".$outgroupObject->get_name.".cs219",
                                      consensusDirectory => $self->configuration->hieranoidConsensusDirectory,
                                      consensusFile => $self->configuration->hieranoidConsensusDirectory."/".$outgroupObject->get_name.".cons.".$self->configuration->sequenceInputFormat,
                                      alignmentFile => $self->configuration->hieranoidConsensusDirectory."/".$outgroupObject->get_name.".aln.".$self->configuration->sequenceInputFormat,
                                      orthologGroupsFormat => $self->configuration->orthologGroupsFormat,
                                      sequenceInputFormat => $self->configuration->sequenceInputFormat,
                                      sourceFile => $self->configuration->speciesFilesDirectory."/".$outgroupObject->get_name.".".$self->configuration->sequenceInputFormat
                                      ))
                                              ));

                      $self->outgroupDaughterType("leaf");
           }
           else{
                      $self->outgroupDaughter(InnerNode->new(
                              nodeObject => $outgroupObject, 
                                  name => $outgroupObject->get_name,
                                  configuration => $self->configuration));
                      $self->outgroupDaughterType("innerNode");
          }
              

          return 1;
              
}

sub get_terminalSpeciesBelonging2Node {
        my $self = shift;
        my @terminals;
        foreach(@{ $self->get_terminals }){
            push(@terminals,$_->get_name);    
        }
        return @terminals;
}
=item get_alignmentsByID()

to be written...

 Title   : get_alignmentsByID
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub get_alignmentsByID{
        my $self = shift;
        my $ID2alignmentsHashRef = shift;
      #  print "alignments files in  ".$self->fileInformation->alignmentFile."\n";
        if(! -e $self->fileInformation->alignmentFile && ! -s $self->fileInformation->alignmentFile){
              WARN("Missing file (".$self->fileInformation->alignmentFile."). No alignments to read\n");
              return 0;
        }
        open my $ALIGNMENT_FILE, '<', $self->fileInformation->alignmentFile or croak "Couldn't open '".$self->fileInformation->alignmentFile."': $OS_ERROR";
	
	### READING FROM FILE
	while ( my $line = <$ALIGNMENT_FILE> ) {
		chomp($line);
	        my ($id,$alignment) = split(/\t/,$line);
	        $alignment =~ s/\#/\n/g;
	        $ID2alignmentsHashRef->{$id} = $alignment;
	}
        close $ALIGNMENT_FILE or croak "Couldn't close '".$self->fileInformation->alignmentFile."': $OS_ERROR";
}
=item get_alignmentsByID()

to be written...

 Title   : get_alignmentsByID
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub get_HmmsByID{
        my $self = shift;
        my $ID2HMMHashRef = shift;
      #  print "alignments files in  ".$self->fileInformation->alignmentFile."\n";
        if(! -e $self->fileInformation->hmmFile && ! -s $self->fileInformation->hmmFile){
              WARN("Missing file (".$self->fileInformation->hmmFile."). No HMMs to read\n");
              return 0;
        }
        open my $HMM_FILE, '<', $self->fileInformation->hmmFile or croak "Couldn't open '".$self->fileInformation->hmmFile."': $OS_ERROR";
	
	my $HMMBody;
	my $currentHMMName;
	### READING FROM FILE
	while ( my $line = <$HMM_FILE> ) {
		chomp($line);
		next if $line =~ /^$/;
		$HMMBody .= $line."\n";
		#print "\tsaving $line: \n";
		if($line =~ /^NAME\s+(.*)/){
		        $currentHMMName = $1;
		        #print "\tset name to $currentHMMName\n";
		        if(!defined($currentHMMName) || $currentHMMName eq ''){
		                print "\tProblem reading Hmm. Could not find name \n"; 
		                exit;
		        }
		}
		if($line =~ /^\/\//){
		        #print "\treached \/\/\n";
		        #exit;
		        # new entry
		        if($currentHMMName eq '' || $HMMBody eq ''){
		                print "\tEmpty Hmmname or Body (Hmmname is $currentHMMName)\n";
		                exit;
		        }
		        $ID2HMMHashRef->{$currentHMMName} = $HMMBody;
		        $currentHMMName = '';
		        $HMMBody = '';
		        
		}
	}
        close $HMM_FILE or croak "Couldn't close '".$self->fileInformation->hmmFile."': $OS_ERROR";
}
=item getDaughterSequencesById()

        e.g. Node Homininae
                gets all sequences by Id in files *.fa in directory Homininae
                Pan_troglodytes.fa
                Homo_sapiens.fa
                
 Title   : getDaughterSequencesById
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub getDaughterSequencesById{
          my $self = shift;
          my $ID2sequencesHashRef = shift;
          my %iD2SequenceHash;

          my @sourceFileToRead = glob($self->fileInformation->outputDirectory."/*.fa");
          #my @sourceFileToRead = ($self->leftDaughter->fileInformation->sequenceSearchInputFile,
          #              $self->rightDaughter->fileInformation->sequenceSearchInputFile);

         # print "\t read sequence files in directory ".$self->fileInformation->consensusDirectory."/*.fa : ".join(",",@sourceFileToRead)."\n";

          #exit;
          foreach my $species_sequence_file(@sourceFileToRead){
                  next if $species_sequence_file =~ /_masked_query/;
                  #print "\t\t\tcladesequences from ".$self->name." file: $species_sequence_file\n";
                  # SEQXML FILE
                  if ( $self->fileInformation->sequenceInputFormat eq q{xml} ) {
                       #   print("\treading sequence information from ".basename($species_sequence_file)." (file: $species_sequence_file");

                          if(! -e $species_sequence_file && ! -s $species_sequence_file){
                                  WARN("sequence search input file ($species_sequence_file) does not exist\n");
                                  return 0;
                          }
                          readSeqxmlFileIntoHash({
                                  sequence_file => $species_sequence_file,
                                  sequence_href => \%$ID2sequencesHashRef});
                          } 
                          else {
                                  # FASTA FILE
                                  if(! -e $species_sequence_file && ! -s $species_sequence_file){
                                          WARN(print "sequence search input file ($species_sequence_file) does not exist\n");
                                          return 0;
                                  }
                         #         print("\treading sequence information from ".$self->name." (file: $species_sequence_file");
                                  readSequenceFileIntoHash({
                                          sequence_file => $species_sequence_file,
                                          sequence_href => \%$ID2sequencesHashRef});
                                  }
                            #      print Dumper $ID2sequencesHashRef;
                            #      exit;
                                  WARN("Could not read sequence information for node ".$self->name."\n")  if !( keys(%$ID2sequencesHashRef) );
                                  #$self->sequencesByID(%iD2SequenceHash);
                          }
                          return 1;
}
=item get_sequencesByID()

to be written...

 Title   : get_sequencesByID
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub get_sequencesByID{
        my $self = shift;
        my $ID2sequencesHashRef = shift;
        my @terminalSpeciesBelonging2Node = get_terminalSpeciesBelonging2Node($self->nodeObject);
        #print Dumper @terminalSpeciesBelonging2Node;
        my %iD2SequenceHash;
        my $species_sequence_file;
        SPECIES_TO_ADD:
        foreach my $current_species (@terminalSpeciesBelonging2Node) {
                # Have to skip missing sequence files. Just ignore them, we don not need them
                #print "$current_species ".$self->fileInformation->sequenceInputFormat."\n";
                # SEQXML FILE
                if ( $self->fileInformation->sequenceInputFormat eq q{xml} ) {
                        $species_sequence_file = $self->fileInformation->speciesFilesDirectory."/".$current_species.".xml";
                     #   print("\treading sequence information from ".basename($species_sequence_file)." (file: $species_sequence_file");

                        if(! -e $species_sequence_file && ! -s $species_sequence_file){
                                next SPECIES_TO_ADD;
                        }
                        readSeqxmlFileIntoHash({
                                sequence_file => $species_sequence_file,
                                sequence_href => \%$ID2sequencesHashRef});
                } 
                else {
                                # FASTA FILE
                        $species_sequence_file = $self->fileInformation->speciesFilesDirectory."/".$current_species.".fa";
                     #   print "$current_species ".$self->fileInformation->sequenceInputFormat." ($species_sequence_file)\n";
                        if(! -e $species_sequence_file && ! -s $species_sequence_file){
                                next SPECIES_TO_ADD;
                        }
                        #print("\treading sequence information from ".$current_species." (file: $species_sequence_file");
                        readSequenceFileIntoHash({
                                sequence_file => $species_sequence_file,
                                sequence_href => \%$ID2sequencesHashRef});
                        }
                }
                        print("Could not read sequence information for node ".$self->name."\n")  if !( keys(%$ID2sequencesHashRef) );
#                        $self->sequencesByID(%iD2SequenceHash);
}
=item get_cladeSequencesByID()

to be written...

 Title   : get_cladeSequencesByID
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub get_cladeSequencesByID{
        my $self = shift;
        my $ID2sequencesHashRef = shift;
        my %iD2SequenceHash;
        DEBUG("consensus files in dir  ".$self->fileInformation->consensusDirectory."\n");
        
        my @sourceFileToRead = glob($self->fileInformation->consensusDirectory."/*.fa");
        
       # print "\t read sequence files in directory ".$self->fileInformation->consensusDirectory."/*.fa : ".join(",",@sourceFileToRead)."\n";
        
        #exit;
        foreach my $species_sequence_file(@sourceFileToRead){
                next if $species_sequence_file =~ /_masked_query/;
                #print "\t\t\tcladesequences from ".$self->name." file: $species_sequence_file\n";
                # SEQXML FILE
                if ( $self->fileInformation->sequenceInputFormat eq q{xml} ) {
                     #   print("\treading sequence information from ".basename($species_sequence_file)." (file: $species_sequence_file");

                        if(! -e $species_sequence_file && ! -s $species_sequence_file){
                                WARN("sequence search input file ($species_sequence_file) does not exist\n");
                                return 0;
                        }
                        readSeqxmlFileIntoHash({
                                sequence_file => $species_sequence_file,
                                sequence_href => \%$ID2sequencesHashRef});
                        } 
                        else {
                                # FASTA FILE
                                if(! -e $species_sequence_file && ! -s $species_sequence_file){
                                        WARN(print "sequence search input file ($species_sequence_file) does not exist\n");
                                        return 0;
                                }
                       #         print("\treading sequence information from ".$self->name." (file: $species_sequence_file");
                                readSequenceFileIntoHash({
                                        sequence_file => $species_sequence_file,
                                        sequence_href => \%$ID2sequencesHashRef});
                                }
                          #      print Dumper $ID2sequencesHashRef;
                          #      exit;
                                WARN("Could not read sequence information for node ".$self->name."\n")  if !( keys(%$ID2sequencesHashRef) );
                                #$self->sequencesByID(%iD2SequenceHash);
                        }
}
=item get_daughterType()

to be written...

 Title   : get_daughterType
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub get_daughterType{
        my $self = shift;
        my @species_belonging_to_pseudospecies = ();
        foreach(@{ $self->get_terminals }){
                push(@species_belonging_to_pseudospecies, $_->get_name);
        }
        # Determine whether left daughter is a single species or a pseudospecies
        return (scalar(@species_belonging_to_pseudospecies) > 1)? q{profile} : q(sequence);
}
=item readSequenceFileIntoHash()

to be written...

 Title   : readSequenceFileIntoHash
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub readSequenceFileIntoHash {
	#### PARAMETER VARIABLES
	#print Dumper @_;
	my ($arg_ref)     = @_;
	my $sequence_file = $arg_ref->{sequence_file};
	my $sequence_href = $arg_ref->{sequence_href};
	#### TESTING ################################################################################
	my $seq_id = q();
	my $seq    = q();
	## CHECK DEFINEDNESS
	#croak q{One/several of the parameters for 'read_sequence_file_into_hash' was/were not defined.\n}
	#  if any { !defined $_ } $sequence_file;
	
	### CHECK FILE EXISTS AND CAN BE USED
	WARN("Input file $sequence_file is empty or missing\n")   if ( !-e $sequence_file || !-s $sequence_file );
	#print "\t Reading sequences from $sequence_file\n";
	open my $SEQUENCE_FILE, '<', $sequence_file or croak "Couldn't open '$sequence_file': $OS_ERROR";
	
	### READING FROM FILE
	while ( my $line = <$SEQUENCE_FILE> ) {
		chomp($line);
		if ( $line =~ /^>/ ) {
			my @words = split( /\s/, $line );
			$seq_id = $words[0];
			$seq_id =~ s/>//;
             #print "seq_id $seq_id\n";
			$sequence_href->{$seq_id} = q();
		} else {
			$sequence_href->{$seq_id} .= $line;
		}
	}
	
	### CLOSING FILE
	close $SEQUENCE_FILE or croak "Couldn't close '$sequence_file': $OS_ERROR";
#        print $sequence_href->{"LCT_PANTR"};
	
	#my $logger = get_logger("READING FASTA FILES");
        #print("Read ".keys( %{$sequence_href} )." sequences into hash");
    if ( !keys( %{$sequence_href} ) ){
        print("Could not read fasta entries from $sequence_file");
        exit
    }
	return 1;
}
=item readSeqxmlFileIntoHash()

to be written...

 Title   : readSeqxmlFileIntoHash
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub readSeqxmlFileIntoHash {
	#### PARAMETER VARIABLES
	my ($arg_ref)     = @_;
	my $sequence_file = $arg_ref->{sequence_file};
	my $sequence_href = $arg_ref->{sequence_href};
	#### TESTING ################################################################################
	my $seq_id = q();
	my $seq    = q();
	## CHECK DEFINEDNESS
	croak q{One/several of the parameters for 'method' was/were not defined.\n}
	  if any { !defined $_ } $sequence_href, $sequence_file;
	### CHECK FILE EXISTS AND CAN BE USED
	Inparanoid_module::printhelp(
								"Input file $sequence_file is empty or missing")
	  if ( !-e $sequence_file || !-s $sequence_file );
	open my $SEQUENCE_FILE, '<', $sequence_file
	  or croak "Couldn't open '$sequence_file': $OS_ERROR";
	### READING FROM FILE
        # Since we are reading an seqxml file, we just have to care about the following fields
        # <entry id="UniProtKB/Swiss-Prot:Q8N283">
        # <AAseq>MKR...</AAseq>

	my $current_seq = q();
	while ( my $line = <$SEQUENCE_FILE> ) {
		if ( $line =~ /<entry id=\"(.*)\">/ ) {
			$1 =~ /:(\w+)/;
			$seq_id = $1;
			#print $seq_id. "\n" and exit if $seq_id eq "ENSMUSG00000066036";

			#exit;
		}
		if ( $line =~ /<DBRef.*id="(\w+)" \/>/ ) {
			$sequence_href->{$1} = $current_seq;
		}
		if ( $line =~ /<AAseq>(.*)<\/AAseq>/ ) {
			$current_seq = $1;
			$sequence_href->{$seq_id} = $current_seq;
		}
	}
	### CLOSING FILE
	close $SEQUENCE_FILE or croak "Couldn't close '$sequence_file': $OS_ERROR";
	
	my $logger = get_logger("READING SEQXML FILES");
	$logger->debug("Read ".keys( %{$sequence_href} )." sequences into hash");
	
	if ( !keys( %{$sequence_href} ) ){
	  $logger->error("Could not read seqxml entries from $sequence_file");
	}
	return 1;
}

1;