package Config;
use Moose;
use Carp;
use English;
use Data::Dumper;
use Cwd;

#Reads in Configuration file
## ATTRIBUTES
has 'configurationFile', is => 'rw', isa => 'Str';

has 'treeFile', is => 'rw', isa => 'Str';   


has 'hieranoidResultsDirectory', is => 'rw', isa => 'Str';   
has 'hieranoidProfilesDirectory', is => 'rw', isa => 'Str';   
has 'hieranoidConsensusDirectory', is => 'rw', isa => 'Str';   
has 'allResultsDirectory', is => 'rw', isa => 'Str';   
#has 'hieranoidMappingDirectory', is => 'rw', isa => 'Str';   

has 'profileSearch', is => 'rw', isa => 'Str';   

has 'speciesFilesDirectory', is => 'rw', isa => 'Str';   
has 'rootDirectory', is => 'rw', isa => 'Str';   
has 'log_directory', is => 'rw', isa => 'Str';   

has 'computation_mode', is => 'rw', isa => 'Str';   
has 'typeOfAnalysis', is => 'rw', isa => 'Str';   
has 'orthologyPredictionTool', is => 'rw', isa => 'Str';   
has 'sequenceInputFormat', is => 'rw', isa => 'Str';   
has 'similaritySearchTool', is => 'rw', isa => 'Str';   
has 'summarizeInformation', is => 'rw', isa => 'Str';   
has 'orthologGroupsFormat', is => 'rw', isa => 'Str';   

has 'addNonMatchingSequences', is => 'rw', isa => 'Str';   
has 'computationMode', is => 'rw', isa => 'Str';   
has 'noHHsearchHits', is => 'rw', isa => 'Str';   


has 'use_outgroup', is => 'rw', isa => 'Str';   
# Parallel
has 'available_nodes', is => 'rw', isa => 'Str';   
has 'wallTime', is => 'rw', isa => 'Str';   
has 'jobNumber', is => 'rw', isa => 'Str';   
has 'type_of_analysis', is => 'rw', isa => 'Str';   
##### Tools
# Tree
#has 'fasttree', is => 'rw', isa => 'Str';   
# Alignment
has 'kalign', is => 'rw', isa => 'Str';   
has 'muscle', is => 'rw', isa => 'Str';   
# Profiles
has 'hhmake', is => 'rw', isa => 'Str';   
has 'hhsearch', is => 'rw', isa => 'Str';   
has 'hmmbuild', is => 'rw', isa => 'Str';   
has 'hmmscan', is => 'rw', isa => 'Str';   

has 'create_hhblits_db', is => 'rw', isa => 'Str';   
has 'create_hhblits_csdb', is => 'rw', isa => 'Str';   

has 'hhconsensus', is => 'rw', isa => 'Str';   
has 'profile_comparer', is => 'rw', isa => 'Str';   
# Masking
has 'segmasker', is => 'rw', isa => 'Str';   
# SimilaritySearch
has 'usearch', is => 'rw', isa => 'Str';   
has 'ublastParameters', is => 'rw', isa => 'Str';   

has 'blast', is => 'rw', isa => 'Str';   
has 'formatdb', is => 'rw', isa => 'Str';   
has 'perl', is => 'rw', isa => 'Str';   
has 'timeFile', is => 'rw', isa => 'Str';   



has 'hhblits', is => 'rw', isa => 'Str';   
has 'similaritySearchCutoff', is => 'rw', isa => 'Str';   
has 'blastParser', is => 'rw', isa => 'Str';   
has 'blast2xml', is => 'rw', isa => 'Str';   


#logging
has 'hieranoid_log', is => 'rw', isa => 'Str';   
has 'inparanoid_log', is => 'rw', isa => 'Str';   
has 'timeFile', is => 'rw', isa => 'Str';   
has 'tmpDir', is => 'rw', isa => 'Str';   




sub BUILD {
      my ($self, $configurationFile) = (@_);
      # Reading Configuration file
      require $self->configurationFile;
      #$self->configurationFile($configurationFile);
      $self->sequenceInputFormat($Configuration::sequenceInputFormat);
      $self->treeFile($Configuration::treeFile);
      $self->similaritySearchTool($Configuration::similaritySearchTool);
      $self->orthologyPredictionTool($Configuration::orthologyPredictionTool);
      $self->allResultsDirectory($Configuration::allResultsDirectory);
      $self->profileSearch($Configuration::profileSearch);
      $self->orthologGroupsFormat($Configuration::orthologGroupsFormat);
      $self->computationMode($Configuration::computationMode);
      $self->available_nodes($Configuration::available_nodes);
      $self->wallTime($Configuration::wallTime);
      $self->use_outgroup($Configuration::use_outgroup);
      $self->noHHsearchHits($Configuration::noHHsearchHits);
      $self->addNonMatchingSequences($Configuration::addNonMatchingSequences);
      $self->muscle($Configuration::muscle);
      $self->kalign($Configuration::kalign);
      $self->hhsearch($Configuration::hhsearch);
      $self->hhmake($Configuration::hhmake);
      $self->hmmbuild($Configuration::hmmbuild);
      $self->hmmscan($Configuration::hmmscan);
      $self->segmasker($Configuration::segmasker);
      $self->summarizeInformation($Configuration::summarizeInformation);
      $self->usearch($Configuration::usearch);
      $self->ublastParameters($Configuration::ublastParameters);
      $self->blast($Configuration::blast);
      $self->formatdb($Configuration::formatdb);
      $self->perl($Configuration::perl);
      $self->timeFile($Configuration::timeFile);
      
      #$self->hhblits($Configuration::hhblits);
      #$self->create_hhblits_db($Configuration::create_hhblits_db);
      #$self->create_hhblits_csdb($Configuration::create_hhblits_csdb);
      
      $self->blastParser($Configuration::blastParser);
      $self->blast2xml($Configuration::blast2xml);
      $self->similaritySearchCutoff($Configuration::similaritySearchCutoff);
      $self->rootDirectory($Configuration::rootDirectory);
      $self->log_directory($Configuration::log_directory);
      $self->hieranoidResultsDirectory($Configuration::hieranoidResultsDirectory);
      $self->hieranoidProfilesDirectory($Configuration::hieranoidProfilesDirectory);
      $self->hieranoidConsensusDirectory($Configuration::hieranoidConsensusDirectory);
      $self->speciesFilesDirectory($Configuration::speciesFilesDirectory);
      $self->hieranoid_log($Configuration::hieranoid_log);
      $self->timeFile($Configuration::timeFile);
      $self->inparanoid_log($Configuration::inparanoid_log);
      $self->tmpDir($Configuration::tmpDir);
      #return 1;
      
      
      ## Test if required programs/directories exist
      die 'No root directory specified\n' if(!$self->rootDirectory || !-e $self->rootDirectory);
      die "No treeFile (".$self->treeFile.") specified\n" if(!$self->treeFile || ! -e $self->treeFile);
      die 'No allResultsDirectory specified\n' if(!$self->allResultsDirectory);
      die 'No results directory specified\n' if(!$self->hieranoidResultsDirectory);
      die 'No profiles directory specified\n' if(!$self->hieranoidProfilesDirectory);
      die 'No summarizeInformation specified\n' if(!$self->summarizeInformation);
      die 'No profileSearch specified\n' if(!$self->profileSearch);
      die 'No similaritySearchTool specified\n' if(!$self->similaritySearchTool);
      die 'No orthologyPredictionTool specified\n' if(!$self->orthologyPredictionTool);
      die 'No orthologGroupsFormat specified\n' if(!$self->orthologGroupsFormat);
      die 'No computation_mode specified\n' if(!$self->computationMode);
      die 'No available_nodes specified\n' if(!$self->available_nodes);
      die 'No wallTime specified\n' if(!$self->wallTime);
      die 'No noHHsearchHits specified\n' if(!$self->noHHsearchHits);
      die 'No perl specified\n' if(!$self->perl);
      die 'No timeFile specified\n' if(!$self->timeFile);
      
      die "Directory with sequence files (".$self->speciesFilesDirectory.") does not exists\n" if(!-e $self->speciesFilesDirectory);
      die "Directory with sequence files empty ".$self->speciesFilesDirectory."\n" if(! glob($self->speciesFilesDirectory."/*"));
      #die 'No muscle specified\n' if(!defined($self->muscle) || !-e $self->muscle);
      #die 'No hhmake specified\n' if(!defined($self->hhmake) || !-e $self->hhmake);
      #die 'No hhblits specified\n' if(!$self->hhblits);
      die 'No usearch specified\n' if(!$self->usearch);
      die 'No similaritySearchCutoff specified\n' if(!$self->similaritySearchCutoff);
      die 'No sequenceInputFormat specified\n' if(!$self->sequenceInputFormat);
      die 'No speciesFilesDirectory specified\n' if(!$self->speciesFilesDirectory);

      
      mkdir($self->allResultsDirectory) if ! -e $self->allResultsDirectory;
      mkdir($self->hieranoidResultsDirectory) if ! -e $self->hieranoidResultsDirectory;
      mkdir($self->hieranoidConsensusDirectory) if ! -e $self->hieranoidConsensusDirectory;
      mkdir($self->hieranoidProfilesDirectory) if ! -e $self->hieranoidProfilesDirectory;
      mkdir($self->tmpDir) if ! -e $self->tmpDir;
      #make required directories
      mkdir($self->log_directory) if ! -e $self->log_directory;
      die "No log_directory (".$self->log_directory.") specified\n" if(!$self->log_directory || ! -e $self->log_directory);
      
}



1;