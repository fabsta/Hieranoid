#!/usr/bin/perl -w
#
#Script name  :    hieranoid.pl
#
#Date created :    Dezember 2010
#
#Author       :    Fabian Schreiber <fschrei@sbc.su.se>
#
#.
#
# Starter file for analysis using Hieranoid
#
#
#INCLUDES
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Carp;
use English;
use Log::Log4perl qw(:easy);
use Data::Dumper;
use Benchmark;
use Hieranoid::Comparison;
use Hieranoid::Tree::TreeMaster;
use Hieranoid::Tree::InnerNode;
use Hieranoid::Config;
#use Config;
use sbc::orthoxml::Database;
use sbc::orthoxml::Gene;
use sbc::orthoxml::Group;
use sbc::orthoxml::Species;
use sbc::orthoxml::Membership;
use sbc::orthoxml::io::OrthoXMLReader;
#use sbc::orthoxml::io::OrthoXMLWriter;

####################
# VARIABLES
####################
### INPUT PARAMETERS #############################################################################

##### SPECIES TREE
my $species_tree_file = q();
##### Folder containing all species
my $species_folder = q();
### OPTION FOR COMPUTATION ON SINGLE CORE, MULTI CORE OR COMPUTER CLUSTER
my $execution_option   = q();
my $use_outgroup       = q();
my $compute_alignments = 0;
my $compute_trees      = 0;
my $jobNumber;
my $compute_node = 0;
my $type_of_analysis;
my $configurationFile = "Configurations/Configuration.pm";
my $options = GetOptions(
						  "mode|m=s"    => \$execution_option,     # numeric
						  "tree|t=s"    => \$species_tree_file,    # string
						  "species|s=s" => \$species_folder,
						  "out|o=s"     => \$use_outgroup,
						  "con|c=s"     => \$configurationFile,
						  "comp_ali"    => \$compute_alignments,
						  "job|j=s"     => \$jobNumber,
						  "comp_tree"   => \$compute_trees,
						  "n=s"         => \$compute_node,
						  "a=s"         => \$type_of_analysis
);

#print "tmp dir is ".$ENV{TMPDIR}."\n";
#$ENV{TMPDIR}= "/scratch/";

### CHECK FILE EXISTS AND CAN BE USED  ###########################################################
BEGIN:
  {
	my $hieranoidConfiguration = Config->new( configurationFile => $configurationFile );
	my $hieranoid_log          = "logs/hieranoid.log";
	my $inparanoid_log         = "logs/inparanoid.log";
	unlink($hieranoid_log)  if -e $hieranoid_log;
	unlink($inparanoid_log) if -e $inparanoid_log;
	print "Removing log files\n";
	require $configurationFile;
	mkdir("log") if !-e "log";
  }

####################
# MAIN
####################
MAIN:
  {

	# make a new Configuration object here
	my $hieranoidConfiguration = Config->new( configurationFile => $configurationFile );
	### LOGGING LEVEL
	#print Dumper $hieranoidConfiguration;
	# everything under 'error' will be reported
	my $logFile = $hieranoidConfiguration->hieranoid_log;
	$logFile .= ".job$jobNumber" if $jobNumber;
	Log::Log4perl->easy_init(
		{
		   level => $DEBUG,

		   #level => $WARN,
		   layout => '%d %p> %F{1}:%L %M - %m%n',
		   file   => ">>" . $logFile
		}
	);

	if ( $hieranoidConfiguration->tmpDir )
	  {
		$ENV{TMPDIR} = $hieranoidConfiguration->tmpDir;
	  }
	if ( defined($jobNumber) )
	  {

		# Test if jobNumber is lower/equal to number of available nodes
		if ( $jobNumber > $hieranoidConfiguration->available_nodes )
		  {
			print "Job number ($jobNumber) higher than number of available nodes "
			  . $hieranoidConfiguration->available_nodes . " \n";
			exit;
		  }
		$hieranoidConfiguration->jobNumber($jobNumber);
	  }
	if ( defined($type_of_analysis) )
	  {
		$hieranoidConfiguration->type_of_analysis($type_of_analysis);
	  }

	# override Configuration file parameters if there are user-defined ones

	# Problematic species tree
	if ($species_tree_file)
	  {
		if ( !-e $species_tree_file || !-s $species_tree_file )
		  {
			print "Tree file provided, but it does not exist\n";
		  }
		$hieranoidConfiguration->treeFile($species_tree_file);
	  }

	# Empty sequences folder
	if ($species_folder)
	  {
		if ( !-e $species_folder || !-s $species_folder )
		  {
			print "Sequence folderprovided, but it does not exist/is empty\n";
		  }
		$hieranoidConfiguration->speciesFilesDirectory($species_folder);
	  }

	my $treemaster = TreeMaster->new( configuration   => $hieranoidConfiguration,
									  speciesTreeFile => $hieranoidConfiguration->treeFile );
	my $counter = 0;

	# start timer
	my $start_time_hieranoid = new Benchmark;

  COMPARISON:
	foreach my $innerNode ( reverse( @{ $treemaster->comparisons } ) )
	  {
		my $currentComparison;
		if ($compute_node)
		  {
			next COMPARISON if $compute_node ne $innerNode->name;
		  }

		#my $currentComparison = Comparison->new(nodeObject =>$innerNode,
		#					configuration => $hieranoidConfiguration);
		print "\tCurrent Node "
		  . $innerNode->name . " ("
		  . $innerNode->leftDaughter->name . " vs. "
		  . $innerNode->rightDaughter->name . ") ";

		#last if $innerNode->name ne 'Euarchontoglires';
		if ( $hieranoidConfiguration->use_outgroup )
		  {

			#print "\ttaking outgroup\n";
			my $outgroupNode = $treemaster->get_outgroup_node($innerNode);

			if ( $outgroupNode eq 0 )
			  {
				print "Root node reached, there exists no outgroup\n";
				print "\n";
			  }
			elsif ( $outgroupNode eq '' )
			  {
				print "Error finding outgroup. Outgroup object is empty\n";
				exit;
			  }
			else
			  {

				#print "\tsetting outgroup to ".$outgroupNode->get_name."\n";
				$innerNode->setOutgroup($outgroupNode);
				print "using " . $outgroupNode->get_name() . " as outgroup\n";
			  }

			# test existence of outgroup object
			if ( !$innerNode->outgroupDaughter )
			  {

				#print "no outgroup object\n"
			  }

			#next;
			$currentComparison =
			  Comparison->new( nodeObject => $innerNode, configuration => $hieranoidConfiguration );
		  }
		else
		  {
			print "using no outgroup\n";
			$currentComparison =
			  Comparison->new( nodeObject => $innerNode, configuration => $hieranoidConfiguration );
			print "\n";
		  }
		$currentComparison->compareNodes( $jobNumber, $type_of_analysis );
	  }

	# Stop time
	my $end_time_hieranoid = new Benchmark;
	my $time_diff_total    = timediff( $end_time_hieranoid, $start_time_hieranoid );
	my $debugString        = "Hieranoid\tTotal: " . timestr( $time_diff_total, 'all' );
	&attach_to_file( { file_name => $hieranoidConfiguration->timeFile, text => $debugString . "\n" } );

  }

####################
# END
####################

END
{
	exit;
	print STDERR <<ENDHELP;
.__   .-".                                  
(o\"\  |  |   Hieranoid:                      
   \_\ |  |                                  
  _.---:_ |     Forget about the rest
 ("-..-" /                                   
  "-.-" /                                    
    /   |                                    
    "--"  FaB                                
ENDHELP

}

sub printhelp
  {
	my $message = shift;
	print "\n------------\n$message\n------------\n";
	print STDERR <<ENDHELP;
    Parameters of $0:
    perl hieranoid_wrapper.pl -m ("statistics"|"analyse"|"prepare"|"split"|translate) -i Fastafile  -a Assignment_file [-nsplits|-pfamA|-pfamT|-pfamH|-pfamF]
	
	
	mode|m=s" => \$execution_option,
						"tree|t=s" => \$species_tree_file,
						"species|s=s" => \$species_folder,
						"out|o=s" => \$use_outgroup

  Required parameters
	-s|species : Folder containing species under study
	-t|tree: File containing the species tree
	-m: modus, available options are: single, multi, cluster
  Example call:
	perl hieranoid_wrapper.pl -s species/orthomcl_fasta/ -m single
ENDHELP
	exit;
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

=cut




