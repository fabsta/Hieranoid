########################################################################
# Script name  :    Configuration.pm
#
# Date created :    December 2011
#
# Author       :    Fabian Schreiber <fab.schreiber@gmail.com>
#  
# This is the configuration file for Hieranoid
# See the User manual for detailed descriptions
# NOTE: Paths must end with "/"
########################################################################
use Cwd; #automatically determins current working directory
package Configuration;
{
#####
# PROJECT PARAMETERS
# 
#####
$rootDirectory = Cwd::cwd;
## attach lib/ to PERL5LIB path
push(@INC,"$rootDirectory/lib");
## Check the following setting and make changes accordingly
	# a directory with 'projectName' will be created in the current directory
	$projectName = "project1";
	# Folder with proteome sequences
	# Note that the file names have to match the species names in the tree
	#$speciesFilesDirectory = "testdata/ensemblv60_12species/longestTranscript_smallSubset";
    #$speciesFilesDirectory = "testdata/inparanoid_data/refGenomes/";
    $speciesFilesDirectory = "data/sequences";
    
    ## Choose sequence input format
	# fasta : ".fa", 
	# seqxml : ".xml"
	$sequenceInputFormat = "fa";
    
    ## Tree File
	# Specify guide tree file
    $treeFile = "data/tree/inparanoid10.tre";
	#$treeFile = "species/trees/inparanoid10RefGe.tre";
	#$treeFile = "species/trees/ensembl60.tre";
    
##############
### BINARIES
#
# Make sure the following binaries are installed
###############

## HMMER BIN DIRECTORY
	# ALIGNMENT PROGRAM
	$muscle = "muscle";
	$kalign = "kalign";
	# SIMILARITY SEARCH
	# Needs one of the following two programs: usearch or blast
	# If using Usearch, you also need segmasker, as usearch does not perform sequence masking
	# Binary of Usearch program http://www.drive5.com/usearch/
	$usearch = "usearch";
	# SEQUENCE MASKING
	$segmasker = "$rootDirectory/extern/segmasker";
	# Binary of Blast program
	$blast = "blastall";
	$formatdb = "formatdb";
	
	# PERL BINARY
	$perl = "perl";
	
	# PROFILE VS PROFILE COMPARISONS
	$hhmake = "/Users/fs9/bin/source/hhsuite-2.0.0-macosx/bin/hhmake";
	$hhblits = "/Users/fs9/bin/source/hhsuite-2.0.0-macosx/bin/hhblits";
	$hhsearch = "/Users/fs9/bin/source/hhsuite-2.0.0-macosx/bin/hhsearch";
	$hmmbuild = "/Users/fs9/bin/source/hmmer-3.0/src/hmmbuild";
	$hmmscan = "/Users/fs9/bin/source/hmmer-3.0/src/hmmscan";
# Parser
	$blast2xml = "blast2xml.pl";
	# Parses similarity search output in xml-format and combines HSP 
	$blastParser = "blast_parser.pl";
	
    
#################################################
## Don't have to change the following
##################################################
	
    $allResultsDirectory = $rootDirectory."/$projectName";
	$hieranoidResultsDirectory = $allResultsDirectory."/nodes";
	$hieranoidProfilesDirectory = $allResultsDirectory."/profiles";
	$hieranoidConsensusDirectory = $allResultsDirectory."/consensus";
	$hieranoidMappingDirectory = $allResultsDirectory."/mapping";


## Similarity Search
        # Summarizing information
        # for consensus: 'clade_consensus'
        # for profiles : 'hierarchical_profile'
        # ...
        $summarizeInformation = 'clade_consensus';
        
## use of outgroup
		# use of outgroup increases computation time, but might detect gene losses
        $use_outgroup = '';
## Profile search
        $profileSearch = "hhsearch";
        # number of hits to look at for profile-profile search
        $noHHsearchHits = "5";
## Orthology prediction        
        # Inparanoid = 'inparanoid'
        # can be expanded later
        $orthologyPredictionTool = 'inparanoid';
        
# Format of orthology predictions 
        # Format of orthology predictions
        # comma-separated genes: 'groupFile'
        # OrthoXML format : 'orthoxml'
        $orthologGroupsFormat = 'groupFile';
        
# Similarity Search        
		# Tool to perform similarity searches
	    # usearch
        # blast
        $similaritySearchTool = 'usearch';
        ## Similarity search specific options
# Blast-specific
        $similaritySearchCutoff = 40;
        # one of the following
        # --maxrejects 0 --maxaccepts 1000
        # --nousort
        $ublastParameters = "--maxlen 100000 --minlen 4 --evalue 0.01 --maxrejects 5 --maxaccepts 5";
        
# Add orphan genes
        $addNonMatchingSequences = 'true';
### Multi-Core specific options
	# Parallel
	# If > 1 node is available, Hieranoid will be started parallel
    	$available_nodes = 1;
    	# Time for the hieranoid master process to wait for child processes to finish
    	$wallTime = 150000;
	## COMPUTATION ON CLUSTER	
		# "single" - single core 
		# "multi" - multi-core
		# "cluster" - (SGE) cluster
		$computationMode = "single";
		$sshCluster = "ssh ferlin.pdc.kth.se";
		## multi-core
		$number_of_cores = "1";       
		## NUMBER OF CLUSTER JOBS
		$number_of_cluster_jobs = "20";

## SEQUENCE TYPE
	# PROTEINS = "p"
	$sequence_type = "p";



####################
## Debugging information
####################
	$log_directory = "$projectName/log";
	$hieranoid_log = "$log_directory/$projectName.hieranoid.log";
	$inparanoid_log = "$log_directory/$projectName.inparanoid.log";
	$timeFile = "$log_directory/$projectName.benchmark.times";
	$tmpDir = "/tmp/";

}

1;
