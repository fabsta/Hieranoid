package OrthologySearch;
use Moose;
use Carp;
use English;
use Data::Dumper;
use Benchmark;
use Log::Log4perl qw(:easy);
use Hieranoid::Comparison::OrthologySearch::SimilaritySearch;
use Hieranoid::Comparison::OrthologySearch::SimilaritySearch::Profile;
use Hieranoid::Comparison::OrthologySearch::SimilaritySearch::HierarchicalProfileSearch;
use Hieranoid::Comparison::OrthologySearch::SimilaritySearch::ProfileSearch;
use Hieranoid::Comparison::OrthologySearch;
use Hieranoid::Comparison::OrthologySearch::OrthologyDetection::InparanoidPrediction;
## ATTRIBUTES
has 'nodeObject',    is => 'rw', isa => 'Object';
has 'configuration', is => 'rw', isa => 'Object';

has 'similaritySearch',   is => 'rw', isa => 'Object';
has 'orthologyDetection', is => 'rw', isa => 'Object';

#has 'file', is => 'rw', isa => 'Str', required => 0;

sub BUILD
  {
	my $self = shift;

	# make similarity search

#print "summarizeType: ".$self->configuration->summarizeInformation." childrenType :".$self->nodeObject->childrenType."\n";

	# LEAF - LEAF
	if ( $self->nodeObject->childrenType eq "leaf_leaf" )
	  {

		# Consensus
		if ( $self->configuration->summarizeInformation =~ /consensus/ )
		  {

			#print "matches leaf_leaf consensus\n";
			$self->similaritySearch(
									 SequenceSearch->new(
														  nodeObject    => $self->nodeObject,
														  configuration => $self->configuration
									 )
			);
		  }

		# Profile
		if ( $self->configuration->summarizeInformation eq 'profile' )
		  {

			#print "matches leaf_leaf profile\n";
			$self->similaritySearch(
									 ProfileSearch->new(
														 nodeObject    => $self->nodeObject,
														 configuration => $self->configuration
									 )
			);
		  }
		if ( $self->configuration->summarizeInformation eq 'hierarchical_profile' )
		  {

			#print "matches leaf_leaf profile\n";
			#$self->similaritySearch(HierarchicalProfileSearch->new(nodeObject => $self->nodeObject,
			#						       configuration => $self->configuration));
			$self->similaritySearch(
									 SequenceSearch->new(
														  nodeObject    => $self->nodeObject,
														  configuration => $self->configuration
									 )
			);

			#print "\tsequence mode for leaf leaf\n";
		  }

		#print "Leaf_leaf created\n";
		#exit;
	  }

	# LEAF - InnerNode
	if (    $self->nodeObject->childrenType eq "leaf_innerNode"
		 || $self->nodeObject->childrenType eq "innerNode_leaf" )
	  {

		# Consensus
		if ( $self->configuration->summarizeInformation =~ /consensus/ )
		  {

			#print "matches innerNode_leaf consensus\n";
			$self->similaritySearch(
									 SequenceSearch->new(
														  nodeObject    => $self->nodeObject,
														  configuration => $self->configuration
									 )
			);
		  }

		# Profile
		if ( $self->configuration->summarizeInformation eq 'hierarchical_profile' )
		  {

			#print "matches leaf_innerNode hierarchical_profile\n";
			$self->similaritySearch(
									 HierarchicalProfileSearch->new(
																	 nodeObject    => $self->nodeObject,
																	 configuration => $self->configuration
									 )
			);
		  }
		if ( $self->configuration->summarizeInformation eq 'profile' )
		  {

			#print "matches innerNode_innerNode profile\n";
			$self->similaritySearch(
									 ProfileSearch->new(
														 nodeObject    => $self->nodeObject,
														 configuration => $self->configuration
									 )
			);
		  }
	  }

	# InnerNoder - InnerNode
	if ( $self->nodeObject->childrenType eq "innerNode_innerNode" )
	  {

		# Consensus
		if ( $self->configuration->summarizeInformation =~ /consensus/ )
		  {

			#print "matches leaf_leaf consensus\n";
			$self->similaritySearch(
									 SequenceSearch->new(
														  nodeObject    => $self->nodeObject,
														  configuration => $self->configuration
									 )
			);
		  }

		# Profile
		if ( $self->configuration->summarizeInformation eq 'hierarchical_profile' )
		  {

			#print "matches innerNode_innerNode profile\n";
			$self->similaritySearch(
									 HierarchicalProfileSearch->new(
																	 nodeObject    => $self->nodeObject,
																	 configuration => $self->configuration
									 )
			);
		  }
		if ( $self->configuration->summarizeInformation eq 'profile' )
		  {

			#print "matches innerNode_innerNode profile\n";
			$self->similaritySearch(
									 ProfileSearch->new(
														 nodeObject    => $self->nodeObject,
														 configuration => $self->configuration
									 )
			);
		  }
	  }

	# Orthology Detection
	if ( $self->configuration->orthologyPredictionTool eq 'inparanoid' )
	  {
		$self->orthologyDetection(
								   InparanoidPrediction->new(
														   nodeObject    => $self->nodeObject,
														   resultsObject => $self->similaritySearch->resultsObject,
														   configuration => $self->configuration
								   )
		);
	  }
	else
	  {
		print "\tother orthology predictions not implemented yet\n";
		exit;
	  }

	#print "\t\tOrthologySearch for type ".$self->nodeObject->childrenType." created\n";

	#my $similarity_search = SimilaritySearch->new();

	#make orthology detection
	#      my $orthology_detection =
  }

sub start
  {
	my ( $self, $jobNumber ) = (@_);

	#        print "\t\tOrthologySearch for type ".$self->nodeObject->childrenType." created\n";
	my $allRequiredFilesExist = 0;
	if (    -e $self->similaritySearch->resultsObject->simAA
		 && -s $self->similaritySearch->resultsObject->simAA
		 && -e $self->similaritySearch->resultsObject->simAB
		 && -s $self->similaritySearch->resultsObject->simAB
		 && -e $self->similaritySearch->resultsObject->simBA
		 && -s $self->similaritySearch->resultsObject->simBA
		 && -e $self->similaritySearch->resultsObject->simBB
		 && -s $self->similaritySearch->resultsObject->simBB )
	  {
		$allRequiredFilesExist = 1;
	  }
	if (    $self->nodeObject->outgroupDaughter
		 && -e $self->similaritySearch->resultsObject->simAC
		 && -s $self->similaritySearch->resultsObject->simAC
		 && -e $self->similaritySearch->resultsObject->simBC
		 && -s $self->similaritySearch->resultsObject->simBC )
	  {
		$allRequiredFilesExist = 1;
	  }
	if ( !$allRequiredFilesExist )
	  {

		# Taking time
		# start timer
		my $start_time_ortho = new Benchmark;

		# Similarity Search
		$self->similaritySearch->start($jobNumber);

		# Stop time
		my $end_time_ortho = new Benchmark;
		my $time_diff_ortho = timediff( $end_time_ortho, $start_time_ortho );
		my $debugString =
		    $self->nodeObject->leftDaughter->name . "-"
		  . $self->nodeObject->rightDaughter->name
		  . "\tOrthologyPrediction: "
		  . timestr( $time_diff_ortho, 'all' );
		attach_to_file( { file_name => $self->configuration->timeFile, text => $debugString . "\n" } );

	  }
	else
	  {
		print "\tlooks like all required files already exist. Skipping...\n";
	  }

	# Parallel mode
	# we compute similarities only
	# orthology detection will be done by master node
	# so we can return here
	return 1 if $jobNumber;

	my ( $leftDaughter, $rightDaughter ) = ( $self->nodeObject->leftDaughter, $self->nodeObject->rightDaughter );
	$allRequiredFilesExist = 0;

	# Check results from Similarity Search
	if (    -e $self->similaritySearch->resultsObject->simAA
		 && -s $self->similaritySearch->resultsObject->simAA
		 && -e $self->similaritySearch->resultsObject->simAB
		 && -s $self->similaritySearch->resultsObject->simAB
		 && -e $self->similaritySearch->resultsObject->simBA
		 && -s $self->similaritySearch->resultsObject->simBA
		 && -e $self->similaritySearch->resultsObject->simBB
		 && -s $self->similaritySearch->resultsObject->simBB )
	  {
		$allRequiredFilesExist = 1;
	  }
	if (    $self->nodeObject->outgroupDaughter
		 && -e $self->similaritySearch->resultsObject->simAC
		 && -s $self->similaritySearch->resultsObject->simAC
		 && -e $self->similaritySearch->resultsObject->simBC
		 && -s $self->similaritySearch->resultsObject->simBC )
	  {
		$allRequiredFilesExist = 1;
	  }
	if ( !$allRequiredFilesExist )
	  {
		ERROR("\tNo results from similarity search. Skipping\n");
		return 0;
	  }

	if (    -e $self->nodeObject->orthologGroups->orthologyPredictionFile
		 && -s $self->nodeObject->orthologGroups->orthologyPredictionFile )
	  {
		DEBUG("\tOrthology predictions already exist. Skipping\n");
		return 1;
	  }
	else
	  {

		# Orthology Search

		DEBUG("\tStart Orthology prediction\n");
		$self->orthologyDetection->start;

		#exit;
	  }
	if (    !-e $self->nodeObject->orthologGroups->orthologyPredictionFile
		 && !-s $self->nodeObject->orthologGroups->orthologyPredictionFile )
	  {
		ERROR(   "\tOrthology predictions do not exist (file: "
			   . $self->nodeObject->orthologGroups->orthologyPredictionFile
			   . "). Skipping\n" );
		exit;
		return 0;
	  }

  }

=item attach_to_file()

to be written...

 Title   : attach_to_file
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut

sub attach_to_file()
  {
	#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $file_name = $arg_ref->{file_name};
	my $text      = $arg_ref->{text};
	if ( !defined $text || $text eq '' )
	  {
		print "\tTrying to write nothing to $file_name\n";
		exit;
	  }
	### OPENING FILE
	open my $out, '>>', $file_name
	  or croak "Couldn't open '$file_name': $OS_ERROR";
	### Writing file
	print {$out} $text;
	### CLOSING FILE
	close $out or croak "Couldn't close '$file_name': $OS_ERROR";
	return 1;
  }

1;
