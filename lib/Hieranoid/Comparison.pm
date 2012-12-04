package Comparison;
=head1 Comparison

hieranoid::Comparison - Container of taxon objects

=head1 SYNOPSIS

 use lib::hieranoid::Comparison;
 my $currentComparison = Comparison->new(nodeObject =>$innerNode,
						configuration => $hieranoidConfiguration);

=head1 DESCRIPTION

Compares orthologous groups between daughter nodes in a given phylogenetic tree

=head1 METHODS

=cut

use Moose;
use Carp;
use File::Temp qw/ tempfile tempdir /;
use English;
use Log::Log4perl qw(:easy);
use Data::Dumper;
use Benchmark;
use Hieranoid::Tree::InnerNode;
use Hieranoid::Comparison::OrthologySearch::SimilaritySearch;
use Hieranoid::Comparison::OrthologySearch::SimilaritySearch::SequenceSearch;
use Hieranoid::Comparison::OrthologySearch;
use Hieranoid::Comparison::Summarizer;
use Hieranoid::Comparison::Summarizer::PairwiseConsensusSummarizer;
use Hieranoid::Comparison::Summarizer::CladeConsensusSummarizer;
use Hieranoid::Comparison::Summarizer::ProfileSummarizer;



require 'Configurations/Configuration.pm';
### LOGGING LEVEL
# everything under 'error' will be reported
my $logFile = $Configuration::hieranoid_log;
Log::Log4perl->easy_init({
	level => $DEBUG,
	#level => $WARN,
	layout => '%d %p> %F{1}:%L %M - %m%n',
	file => ">>".$logFile});



#This is just a "factory" that creates objects

## ATTRIBUTES
has 'nodeObject', is => 'rw', isa => 'Object';
has 'configuration', is => 'rw', isa => 'Object';

## Logic
#has 'similaritySearch', is => 'rw', isa => 'Object';
has 'orthologySearch', is => 'rw', isa => 'Object';
has 'orthologySearchDone', is => 'rw', isa => 'Str';
has 'summarizer', is => 'rw', isa => 'Object';





=head2 CONSTRUCTOR

=over

=item new()

Initialises a comparison object for a given inner Node.

 Title   : new
 Usage   : Comparison->new(nodeObject =>$innerNode,
 						configuration => $hieranoidConfiguration);
 Function: Initialises comparison object
 Returns : -
 Args: nodeObject - inner Node of type <InnerNode>

=cut


sub BUILD {
      my $self = shift;
      # defined workflows
      # consensus
      #         similaritySearchTool: Blast | Ublast
      # profile
      #         similaritySearchTool

      if($self->configuration->summarizeInformation eq 'clade_consensus'){
              $self->summarizer(CladeConsensusSummarizer->new(configuration => $self->configuration,
                                                                nodeObject => $self->nodeObject));
      }
              elsif($self->configuration->summarizeInformation eq 'pairwise_consensus'){
                      $self->summarizer(PairwiseConsensusSummarizer->new(configuration => $self->configuration,
                              nodeObject => $self->nodeObject));              
              }
              elsif($self->configuration->summarizeInformation eq 'hierarchical_profile'){
                      $self->summarizer(ProfileSummarizer->new(configuration => $self->configuration,
                                    nodeObject => $self->nodeObject));
              }
              elsif($self->configuration->summarizeInformation eq 'profile'){
                      $self->summarizer(ProfileSummarizer->new(configuration => $self->configuration,
                                    nodeObject => $self->nodeObject));
              }
      else{
              $self->summarizer(ProfileSummarizer->new(configuration => $self->configuration,
                                                        nodeObject => $self->nodeObject));
      }
      # make orthology search
              $self->orthologySearch(OrthologySearch->new(configuration => $self->configuration,
                                                              nodeObject => $self->nodeObject,
                                                              resultsObject => $self->nodeObject->orthologGroups));     
}

=back

=over

=item compareNodes()

Start comparison between two nodes in a tree (either leaves or internal nodes)

 Title   : compareNodes
 Usage   : $currentComparison->compareNodes();
 Function: Finds orthologs between daughter nodes of a phylogeny
 Returns : 1 on success
 Args: -

=cut
sub compareNodes{
        my ($self,$jobNumber,$type_of_analysis) = (@_);
        # Sub-process: write start file
        if($jobNumber){
                my $startFile = $self->nodeObject->fileInformation->outputDirectory;
                if($type_of_analysis eq 'summarize'){
                        $startFile .= "/process_prepareAnalysis_started.$jobNumber";
                }
                if($type_of_analysis eq 'orthologySearch'){
                        $startFile .= "/process_orthologySearch_started.$jobNumber";
                }
                `touch $startFile`;
                if(!-e $startFile){
                        print "\tcould not create start file $startFile for process $jobNumber\n";
                        exit;
                }
                #print " before summarize $jobNumber\n";
        }
        
        # Taking time
        # start timer
        my $start_time_summ =  new Benchmark;
        # Summarize information
        $self->summarizer->summarizeInformation($self->nodeObject,$jobNumber);  # building consensus || profiles
        # Sub-process: finish this step
        if(defined($type_of_analysis) && $type_of_analysis eq 'summarize'){
                print "\tfinishing summarize for $jobNumber\n";
                my $finishFile = $self->nodeObject->fileInformation->outputDirectory."/process_prepareAnalysis_finished.$jobNumber";
                `touch $finishFile`;
                if(!-e $finishFile){
                        print "\tcould not create start file $finishFile for process $jobNumber\n";
                        exit;
                }
                print "\tExiting sub-process $jobNumber\n";
                exit;   
        }
        # Orthology Search
        my $orthologySearchReturn = $self->orthologySearch->start($jobNumber);
        # Stop time
        my $end_time_ortho = new Benchmark;
        my $time_diff_total = timediff($end_time_ortho,$start_time_summ);
        my $debugString = $self->nodeObject->leftDaughter->name."-".$self->nodeObject->rightDaughter->name."\tTotal: ".timestr($time_diff_total, 'all');
        attach_to_file({file_name => $self->configuration->timeFile, text => $debugString."\n"});        

        # Sub-process finishes with error
        # Write error file
        if(!$orthologySearchReturn && 
        (defined($type_of_analysis) && $type_of_analysis eq 'orthologySearch')){
                my $errorFile = $self->nodeObject->fileInformation->outputDirectory."/process_orthologySearch_error.$jobNumber";
                `touch $errorFile`;
                exit;
        }
        # Sub-process: finish this step
        if(defined($type_of_analysis) && $type_of_analysis eq 'orthologySearch'){
                my $finishFile = $self->nodeObject->fileInformation->outputDirectory."/process_orthologySearch_finished.$jobNumber";
                `touch $finishFile`;
                if(!-e $finishFile){
                        print "\tcould not create stop file $finishFile for process $jobNumber\n";
                        exit;
                }  
                print "\tExiting sub-process $jobNumber\n";
                exit; 
        }
        
        if($self->summarizer->saveOrthologyPredictions($self->nodeObject)){
                ;          
        }
        return 1;
}
=item attach_to_file()

to be written...

 Title   : attach_to_file
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub attach_to_file() {
	#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $file_name = $arg_ref->{file_name};
	my $text      = $arg_ref->{text};
	if(!defined $text || $text eq ''){
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
