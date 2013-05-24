package Summarizer;
=head1 Summarizer

hieranoid::Comparison::Summarizer - Container of taxon objects

=head1 SYNOPSIS

 Abstract class

=head1 DESCRIPTION

Abstract class to handle how orthologous group are summarized

=head1 METHODS

=cut
use Moose::Role;
use Carp;
use Bio::TreeIO;
use IO::String;
use File::Copy;
use File::Basename;
use English;
use Bio::SeqIO;
use Bio::TreeIO;
use Bio::Phylo::IO;
use Bio::Phylo::Forest;
use Log::Log4perl qw(:easy);
use Data::Dumper;
use Benchmark;
use File::Temp qw/ tempfile tempdir /;
use Hieranoid::Tree::InnerNode;
use Hieranoid::FileInformation;
use Bio::AlignIO;

#require 'Configurations/Configuration.pm';
### LOGGING LEVEL
# everything under 'error' will be reported
#my $logFile = $Configuration::hieranoid_log;
#Log::Log4perl->easy_init({
#	#level => $DEBUG,
#	level => $WARN,
#        layout => '%d %p> %F{1}:%L %M - %m%n',
#	file => ">>".$logFile});
## ATTRIBUTES   
has 'configuration', is => 'rw', isa => 'Object', required => 1;



=head2 CONSTRUCTOR

=over

=item new()

Initialises a summarizer object

=cut

sub BUILD {
      my $self = shift;
      DEBUG("\t\tSummarizer created\n");
}
=item summarizeInformation()

Summarizes results of orthology search of daughter nodes.
two daughter nodes will be merged into one pseudospecies
so that it can be used in an Inparanoid comparison

 Title   : compareNodes
 Usage   : $self->summarizer->summarizeInformation($self->nodeObject);
 Function: Summarizes two daughter nodes 
 Returns : 1 on success
 Args: -

=cut
sub summarizeInformation{
        my ($self,$innerNode) = (@_);
        mkdir($innerNode->fileInformation->outputDirectory) if ! -e $innerNode->fileInformation->outputDirectory;
        # case: Sequence/profile search files already exist
        # e.g. Homininae: Homo_sapiens.fa && Pan_troglodytes.fa already exist
        # this saves a lot of time for deeper inner nodes
        
        if(-e $innerNode->leftDaughter->fileInformation->sequenceSearchInputFile && -s $innerNode->leftDaughter->fileInformation->sequenceSearchInputFile
                && -e $innerNode->rightDaughter->fileInformation->sequenceSearchInputFile && -s $innerNode->rightDaughter->fileInformation->sequenceSearchInputFile
                ){
                        DEBUG("Sequence/Profile search files already exist. Skipping this step\n");
                        return 1; 
                }
   # Left: Leaf
        if(-e $innerNode->leftDaughter->fileInformation->sequenceSearchInputFile && -s $innerNode->leftDaughter->fileInformation->sequenceSearchInputFile){
                DEBUG("\t\t\tAlignments and consensus already exists. Skipping...\n");
        }
        else{
                if($innerNode->leftDaughterType eq 'leaf'){
                        DEBUG("\t\tSummarize information: just copy file for ".$innerNode->leftDaughter->name."\n");
                        copy($innerNode->leftDaughter->fileInformation->sourceFile,$innerNode->leftDaughter->fileInformation->sequenceSearchInputFile) or die "Copy failed: (".$innerNode->leftDaughter->fileInformation->sourceFile.",".$innerNode->leftDaughter->fileInformation->sequenceSearchInputFile.")";
                }
                # Left: InnerNode
                else{
                                DEBUG("\t\t\tSummarize information for left innerNode ".$innerNode->leftDaughter->name." (alignment+profile)\n");
                                $self->alignAndConsense($innerNode,$innerNode->leftDaughter);
                                # Allow addition of orphan genes for master process only 
                                if(!$self->configuration->jobNumber){
                                        DEBUG("Adding orphan sequences\n");
                                        print("\t\tAdding orphan sequences from ".$innerNode->leftDaughter->name."\n");
                                        $self->addDaughterNonGroupSequences($innerNode->leftDaughter);
                                }
                }
        }
        #print "\t\tSummarizing information for right: ".$innerNode->rightDaughter->name."\n";
   # Right: Leaf
        if(-e $innerNode->rightDaughter->fileInformation->sequenceSearchInputFile && -s $innerNode->rightDaughter->fileInformation->sequenceSearchInputFile){
                DEBUG("\t\t\tAlignments and consensus already exists. Skipping...\n");
        }        
        else{
                if($innerNode->rightDaughterType eq 'leaf'){
                        DEBUG("\t\tSummarize information: just copy file for ".$innerNode->rightDaughter->name."\n");
                        copy($innerNode->rightDaughter->fileInformation->sourceFile,$innerNode->rightDaughter->fileInformation->sequenceSearchInputFile) or die "Copy failed: (".$innerNode->rightDaughter->fileInformation->sourceFile.",".$innerNode->rightDaughter->fileInformation->sequenceSearchInputFile.")\n";
                        #         $innerNode->rightDaughter->sequenceSearchInputFile($innerNode->rightDaughter->sourceFile);
                }
                else{
                        DEBUG("\t\t\tSummarize information for right innerNode ".$innerNode->rightDaughter->name." (alignment+profile)\n");
                        $self->alignAndConsense($innerNode,$innerNode->rightDaughter);
                        
                        if(!$self->configuration->jobNumber){
                                DEBUG("Adding orphan sequences\n");
                                print("\t\tAdding orphan sequences from ".$innerNode->rightDaughter->name."\n");
                                $self->addDaughterNonGroupSequences($innerNode->rightDaughter);
                        }
                }
        }
   # Outgroup: Leaf
        if($self->configuration->use_outgroup && $innerNode->outgroupDaughter){
                #print "looking at outgroup\n";
                if(-e $innerNode->outgroupDaughter->fileInformation->sequenceSearchInputFile && -s $innerNode->outgroupDaughter->fileInformation->sequenceSearchInputFile){
                        DEBUG("\t\t\tAlignments and consensus already exists. Skipping...\n");
                }        
                else{
                        if($innerNode->outgroupDaughterType eq 'leaf'){
                                DEBUG("\t\tSummarize information: just copy file for ".$innerNode->outgroupDaughter->name."\n");
                                print "\tcopying outgroup ".$innerNode->outgroupDaughter->fileInformation->sourceFile." - ".$innerNode->outgroupDaughter->fileInformation->sequenceSearchInputFile."\n";
                                copy($innerNode->outgroupDaughter->fileInformation->sourceFile,$innerNode->outgroupDaughter->fileInformation->sequenceSearchInputFile) or die "Copy failed: (".$innerNode->outgroupDaughter->fileInformation->sourceFile.",".$innerNode->outgroupDaughter->fileInformation->sequenceSearchInputFile.")\n";
                                #         $innerNode->rightDaughter->sequenceSearchInputFile($innerNode->rightDaughter->sourceFile);
                        }
                        else{
                                DEBUG("\t\t\tSummarize information for outgroup node ".$innerNode->outgroupDaughter->name." (alignment+profile)\n");
                                $self->alignAndConsense($innerNode,$innerNode->outgroupDaughter);

                                if(!$self->configuration->jobNumber){
                                        DEBUG("Adding orphan sequences\n");
                                        print("\t\tAdding orphan sequences from ".$innerNode->outgroupDaughter->name."\n");
                                        $self->addDaughterNonGroupSequences($innerNode->outgroupDaughter);
                                }
                        }
                }  
        }
    return 1;
}

        # Input
                # groups in format Catarrhini2:MmuSTS.4752.1.S1_at_MACMU,VPS45_PANTR,VPS45_HUMAN
                # completely resolved and containing no Pseudospecies names

=item alignAndConsense()
        Summarizes results of orthology search of daughter nodes.
        two daughter nodes will be merged into one pseudospecies
        so that it can be used in an Inparanoid comparison

=head3 Example: OG1: Euarchontoglires2:Homininae2,Murinae2
                1. Get all sequences from Hominiae2 : CCNK_F3_PANTR,CCNK_HUMAN
                   Get all sequences from Murinae2 : Ccnk_MOUSE,XP_580110.2_RAT
                2. Build an alignment of CCNK_F3_PANTR,CCNK_HUMAN,Ccnk_MOUSE,XP_580110.2_RAT
                3. Build consensus sequence and save it as Euarchontoglires2
                

 Title   : alignAndConsense
 Usage   : $self->summarizer->summarizeInformation($self->nodeObject);
 Function: Summarizes two daughter nodes 
 Returns : 1 on success
 Args: -

=cut
                                                        
sub alignAndConsense{
        my ($self,$parentNode,$innerNode, $jobNumber) = (@_);
        my ($fh2, $alignment_tmp_file) = tempfile(UNLINK => 0);
        my $isMasterProcess = 0;
        $jobNumber = $self->configuration->jobNumber;
        #print "jobnumber $jobNumber\n";
        #exit;
        if($self->configuration->available_nodes > 1 && !$jobNumber){
                $isMasterProcess = 1; 
        }
        #print "\tmaster: $isMasterProcess";
        DEBUG("\t\tAlignment and consensus of node ".$innerNode->name." (OrthologyPredictions in ".$innerNode->orthologGroups->groupsFile.")\n");
        my $predictionsFile = $innerNode->orthologGroups->groupsFile; # can be either in orthoXML or group format
        
        my %groupPredictionHash;
        #$innerNode->orthologGroups->get_group2expandedIDsHash($innerNode,\%groupPredictionHash);
        $innerNode->orthologGroups->get_alltreeleafNodes($innerNode,\%groupPredictionHash);
        #print Dumper %groupPredictionHash;
        #exit;
        my %ID2sequencesHash;
        
        #print "\tTrying to collect sequences from ".$innerNode->fileInformation->sequenceSearchInputFile."\n";
        $innerNode->get_sequencesByID(\%ID2sequencesHash);
        $innerNode->get_cladeSequencesByID(\%ID2sequencesHash);
        # align&consense
        #print "\t\tcompute alignments and consensus\n";
        # numeric sort does not work anymore because clusters might look like "Homininae85"
        my $number_of_ogs_total = keys(%groupPredictionHash);
        my %assigned_sequence_ids_hash;
        my $alignments2FileString;
        ### Important output files
        my $sequenceSearchInputFile = $innerNode->fileInformation->sequenceSearchInputFile;
        my $consensusFile = $innerNode->fileInformation->consensusFile;
        my $alignmentFile = $innerNode->fileInformation->alignmentFile;
        my $print2FileString = '';
                
        if($isMasterProcess){  
                        print "\tusing x nodes\n";
                        # parallel
                        # iterate over %singleSimilaritySearchResultsAB
                        my @jobStarted_files;
                        my @jobFinished_files;
                        my $all_finished = 0;
                   # Start parallel Jobs
                        # no_avaliable nodes
                        # range 1  | range 2 | range 3 | 
                        # 1-20     | 21-40   | 41-60   |
                        for (1 .. $self->configuration->available_nodes){
                                print "\tstart sub process $_\n";
                                my $parallelSimilaritySearch;
                                if($self->configuration->computationMode eq 'cluster'){
                                        $parallelSimilaritySearch = $self->configuration->sshCluster.""
                                }
                                else{
                                        $parallelSimilaritySearch = $self->configuration->perl." $0 -j $_ -n ".$parentNode->name." -a summarize -c ".$self->configuration->configurationFile."";
                                }
                                print "\t$parallelSimilaritySearch\n";
                                #exit;
                                system("$parallelSimilaritySearch &");
                                push(@jobStarted_files,$parentNode->fileInformation->outputDirectory."/process_prepareAnalysis_started.$_");
                                push(@jobFinished_files,$parentNode->fileInformation->outputDirectory."/process_prepareAnalysis_finished.$_");
                        }
                                #exit;
                        if(!scalar(@jobStarted_files) || !scalar(@jobFinished_files)){
                                ERROR("\tCould not initialize start/finish files for parallel similarity search\n");
                                exit;
                        }
                        # initially created array of output files from sub-processes
                        # whenever a sub-process successfully finishes, the corresponding
                        # output file will be deleted from the array
                        # empty array --> finished
                        my $directoryWithResultsFiles = $parentNode->fileInformation->outputDirectory;
                        while(@jobFinished_files){
                                # writes status file in the end
                                for(my $i=0;$i<scalar(@jobFinished_files); $i++){
                                        my $currentFile = $jobFinished_files[$i];
                                        DEBUG("\tchecking $currentFile\n");
                                        if(-e $currentFile){
                                                splice(@jobFinished_files, $i, 1); 
                                                DEBUG("\t job has finished (file does exist)\n");
                                        }
                                        else{
                                                DEBUG("\t job has not finished (file does not exist yet)\n");
                                        }
                                }
                                last if !@jobFinished_files;                                        
                                DEBUG("\tCurrent running jobs: ".join(",",@jobFinished_files)."");
                                # check used time
                                
                                my $wallTime = $self->configuration->wallTime;
                                my $findFiles_cmd = "find $directoryWithResultsFiles/process_prepareAnalysis_started* -mmin \+$wallTime";
                                DEBUG("\t$findFiles_cmd\n");
                                my @findOldFiles = `$findFiles_cmd`;
                                if(@findOldFiles){
                                        print "\t There are files running longer than  minutes\n";
                                        exit;
                                }
                                DEBUG("\twaiting\n");
                                sleep(10);
                                #exit;
                        }
                        # Now combine information
                        # first we save all in a temporary file
                        my $cat_seq_cmd = "cat $sequenceSearchInputFile.job* > $sequenceSearchInputFile";
                        #print "$cat_seq_cmd\n";
                        my $cat_con_cmd = "cat $consensusFile.job* > $consensusFile";
                        #print "$cat_con_cmd\n";
                        my $cat_ali_cmd = "cat $alignmentFile.job* > $alignmentFile";
                        #print "$cat_ali_cmd\n";
                        `$cat_seq_cmd`;
                        `$cat_con_cmd`;
                        `$cat_ali_cmd`;
                        if(! -e $sequenceSearchInputFile || ! -s $sequenceSearchInputFile){
                                ERROR("Could not combine files from subprocesses into $sequenceSearchInputFile\n");
                                exit;
                        }
                        # clean files up    
                        print "Clean up old files\n";
                        
                           foreach(glob("$directoryWithResultsFiles/*job*")){
                                   print "\tTrying to delete $_\n";
                                   unlink $_ or print "Error removing file \"$_\": $!\n";
                           }
                        #unlink glob "$sequenceSearchInputFile.job?.*";
                        #unlink glob "$consensusFile.job?.*";
                        #unlink glob "$alignmentFile.job?.*";
                        # process started | finished  
                        foreach(@jobStarted_files){
                                unlink $_;
                        }
                        foreach(@jobFinished_files){
                                unlink $_;
                        }
                        unlink glob dirname($sequenceSearchInputFile)."/process_prepareAnalysis_started*";
                        unlink glob dirname($sequenceSearchInputFile)."/process_prepareAnalysis_finished*";
                        
                        print "\twe stop here to see if the first stage finishes successfully\n";
                        
        }
        else{
                #print "\tno master process\n";
                  my $current_number_of_og = 0;
                  my %duplicateSequencesTotalHash = ();
                # Change output files in parallel mode
                  if($jobNumber){
                        print "\t\tparallel mode, job number $jobNumber\n";
                        $sequenceSearchInputFile .= ".job".$jobNumber;
			$consensusFile .= ".job".$jobNumber;
                        $alignmentFile .= ".job".$jobNumber;
                  }
                  #print Dumper %groupPredictionHash;
                        OG:
                         foreach my $og(sort keys(%groupPredictionHash)){
                                 $current_number_of_og++;
                                 #print "$og && $current_number_of_og  % ".$self->configuration->available_nodes." ".($current_number_of_og % ($self->configuration->available_nodes+1))." != $jobNumber...";
                                 
                                 if($jobNumber && (($current_number_of_og % $self->configuration->available_nodes)+1) != $jobNumber){
                                         #print "skip\n"; 
                                         next OG;
                                 }
                                 else{  
                                         my $current_sequences = q();
                                         my $number_of_sequences_for_group = 0;
                                         my @ids_for_group = ();
                                         my %duplicateSequencesInGroupHash = ();
					 my $current_seq_string = "";
					 my $tmp_keys = "";
					 my $arbitri = 0;
                                         ### delete old files
                                         ### FILE CONTAINS:
                                         # speciesname04 --> member1, member2,...,memberx
                                         #DEBUG("\tCurrent ortholog group: $og ( $current_number_of_og / $number_of_ogs_total)\n");
                                         foreach my $current_id(split(/,/,$groupPredictionHash{$og})){
					 #foreach my $current_id(@tmp_id){
                                                 next if $current_id eq ',';
                                                 my $sequence = q();
                                                 if(exists $duplicateSequencesInGroupHash{$current_id}){
                                                          WARN("Duplicate sequence $current_id found. Skipping..\n");
                                                          print "Duplicate sequence $current_id found. Skipping..\n";
                                                          next;
                                                  }
                                                 if(exists $duplicateSequencesTotalHash{$current_id}){
                                                           WARN("Multiple occurences of ID $current_id in OGs detected. Exiting..\n");
                                                           print "Multiple occurences of ID $current_id in OGs detected. Exiting..\n";
                                                           exit;
                                                 }
                                                 # CHECK IF THIS ID IS ACTUALLY A PROFILE
						 foreach $tmp_keys (keys %ID2sequencesHash){
                                                 #if(!exists $ID2sequencesHash{$current_id}){
							 if ($tmp_keys =~/$current_id/){
								$arbitri = 1;
								$current_seq_string = $ID2sequencesHash{$tmp_keys};
								last;
							 }
						 }
						 if ($arbitri == 0){
                                                  	 ERROR("\t\tCould not find sequence for id [$current_id] of ".$innerNode->name."\n");
                                                         print "\t\tCould not find sequence for id [$current_id] of ".$innerNode->name."\n";	
                                                         exit;
						 }
                                                 # Avoid duplicate sequences
                                                 $duplicateSequencesInGroupHash{$current_id} = 1;
                                                 $duplicateSequencesTotalHash{$current_id} = 1;
                                                 #my $current_seq_string = $ID2sequencesHash{$tmp_keys};
                                                 my $MATCHPATTERN = 'A-Za-z\-\.\*\?';
						 #my $MATCHPATTERN = 'A-Za-z\';
                                                 #if($current_seq_string !~ /^([$MATCHPATTERN]+)$/){
						 if($current_seq_string !~ /^[$MATCHPATTERN]/){
                                                         WARN("Wrong sequence format for $current_id with sequence $current_seq_string");
                                                         next;
                                                 }
                                                 if($current_seq_string eq '' || length($current_seq_string) < 2 || !defined($current_seq_string)){
                                                         WARN("empty sequence for $current_id with sequence $current_seq_string");
                                                         next;
                                                 }
                                                 $current_sequences .= ">".$current_id."\n".$current_seq_string."\n";
                                                 #print "\tchecking $current_id\n";
                                                 #Bio::Seq->new(-seq   => $current_seq_string, 
                                                 #                        -display_id  => $current_id,
                                                 #                        );
                                                 #next;
                                                 $number_of_sequences_for_group++;
                                                 #DEBUG("\t\tadding $current_id\n");
                                                 push(@ids_for_group, $current_id);
                                                 $assigned_sequence_ids_hash{$current_id} = 1;
                                         }
                                         # TEST seq == ''        
                                         if($current_sequences eq q()){
                                                 ERROR("\tCould not find sequences for $og.");
                                                 exit;
                                         }
                                         # #seqs > 2        
                                         if($number_of_sequences_for_group < 2){
                                                 ERROR("Group $og contains < 2 sequences (".join(",",@ids_for_group).")");
                                                 next OG;
                                         }
                                         # Alignment
                                         
                                         #my $alignment_call = $self->configuration->muscle." -maxiters 1 -quiet -diags -sv -distance1 kbit20_3 > $alignment_tmp_file";
                                         #my $alignment_call = $self->configuration->kalign."  -o $alignment_tmp_file >/dev/null 2>&1";
                                         my %alignmentsMethods = (
                                                 "muscle" => $self->configuration->muscle." -maxiters 1 -quiet -diags -sv -distance1 kbit20_3 > $alignment_tmp_file", # Muscle
                                                 "kalign" => $self->configuration->kalign."  -o $alignment_tmp_file >/dev/null 2>&1" # Kalign
                                                 );
                                         my $alignment_successful = 0;
                                         
                                         ALIGNMENT_METHOD:
                                         foreach(keys %alignmentsMethods){
                                                 my $alignment_call = $alignmentsMethods{$_};
                                                  last if $alignment_successful;
                                                 #print "\t$og has $number_of_sequences_for_group sequences: $current_sequences\n";
                                                 #print "Checking alignment in $alignment_tmp_file : ";
                                           # Alignment         
                                                 if(!doAlignment({ alignmentCall => $alignment_call, alignment_tmp_file => $alignment_tmp_file, current_sequences => $current_sequences })){
                                                 next OG;
                                         }
                                                #print "\tsuccessfull ..";
                                                # Test if alignment was successfull
                                                if(!-e $alignment_tmp_file || !-s $alignment_tmp_file){
                                                        WARN("Alignment for $og with members ".join(",",@ids_for_group)." failed");
                                                        print "Alignment for $og with members ".join(",",@ids_for_group)." failed\n";
                                                        #exit;
                                                        next ALIGNMENT_METHOD;
                                                }
                                                #my $alignmentFile = Bio::SeqIO->new(-file => "$alignment_tmp_file",
                                                #                        -format => 'Fasta');
                                                #while ( my $seq = $alignmentFile->next_seq() ) {
                                                #        my $sequence = $seq->seq();
                                                #        if(($sequence =~ /=|\@|\^/) || !($sequence =~ m/\w+/)){
                                                #                print "Problem with: ".$seq->id."\n";
                                                #                next ALIGNMENT_METHOD;
                                                #        }
                                                #}
                                                $alignment_successful = 1;
                                                #print "$_\n";
                                        }
                                        next OG if !$alignment_successful;
                                # Consensus        
                                         my $consensus_string = doConsensus({ alignment_tmp_file => $alignment_tmp_file });
                                         #print " consensus as well\n";
                                         if(!$consensus_string || $consensus_string eq "" ){
                                                 ERROR("could not compute consensus string for alignment ($alignment_tmp_file)\n");
                                                 print "could not compute consensus string for alignment ($alignment_tmp_file)\n";
                                                 next OG;
                                         }
                                         $print2FileString .= ">$og\n$consensus_string\n";
                                         # Save alignment 
                                         my $alignmentAsString;
                                         my $alignment_in = Bio::AlignIO->new('-file' => "$alignment_tmp_file");
                                         my $alignment = $alignment_in->next_aln();
                                         foreach my $seq ($alignment->each_seq) {
                                                 $alignmentAsString .= ">".$seq->id."#".$seq->seq."#";
						 #print "\n-------------\n";
						 #print $seq->id;
						 #print "\n-------------\n";
                                         }
                                         $alignments2FileString .= "$og\t$alignmentAsString\n";
                                         
                                         #next OG;
                                 }
                        }
                        if($print2FileString eq ''){
                                ERROR("Could not summarize group information. Exiting");
                                exit;
                        }
                                 # SequenceSearchInput        
                                 write_to_file({file_name => $sequenceSearchInputFile, text =>  $print2FileString});
                                 # Consensus        
                                 write_to_file({file_name => $consensusFile, text =>  $print2FileString});
                                 # Alignments
                                 write_to_file({file_name => $alignmentFile, text =>  $alignments2FileString});
                         }
                
                # Successfull?        
                         if(!-e $sequenceSearchInputFile || ! -s $sequenceSearchInputFile){
                                               ERROR("Could not summarize hits to  $sequenceSearchInputFile. Exiting\n");
                                               return 0;
                         }
                        return 1;
                        
}


=item addNonGroupSequences()

 Title   : addNonGroupSequences
 Usage   : $self->summarizer->summarizeInformation($self->nodeObject);
 Function: Summarizes two daughter nodes 
 Returns : 1 on success
 Args: -

=cut
sub addNonGroupSequences{
        my ($self,$innerNode) = (@_);

        # get all current groups : as predicted by Inparanoid
        my %all_used_IDs = ();

        #print "\tadding predictions for node ".$innerNode->leftDaughter->name."\n";
    # left Daughter
        if($innerNode->leftDaughterType eq 'innerNode'){
                #print "\ttrying to add sequences to ".$innerNode->fileInformation->sequenceSearchInputFile."\n";
                $innerNode->orthologGroups->getAllUsedIDS($innerNode->leftDaughter,\%all_used_IDs);
                if(!keys %all_used_IDs){
                        ERROR("\tProblem getting current orthologous groups\n"); 
                        exit;
                }        
                #print Dumper %all_used_IDs;
                addDaughterNonGroupSequences({all_used_IDs => \%all_used_IDs,
                         sequenceFile => $innerNode->fileInformation->sequenceSearchInputFile,
                         daughterNode => $innerNode->leftDaughter});
                }
         #print "\tadding predictions for node ".$innerNode->rightDaughter->name."\n";
    # right Daughter
        if($innerNode->rightDaughterType eq 'innerNode'){
                #print "\ttrying to add sequences to ".$innerNode->fileInformation->sequenceSearchInputFile."\n";
                $innerNode->orthologGroups->getAllUsedIDS($innerNode->rightDaughter,\%all_used_IDs);
                if(!keys %all_used_IDs){
                        ERROR("\tProblem getting current orthologous groups\n"); 
                        exit;
                }
                #print Dumper %all_used_IDs;
                 addDaughterNonGroupSequences({all_used_IDs => \%all_used_IDs,
                         sequenceFile => $innerNode->fileInformation->sequenceSearchInputFile,
                         daughterNode => $innerNode->rightDaughter});
                 }        ### SAVE CLUSTER INFORMATION
         return 1;
}

=item addDaughterNonGroupSequences()
        1. Gets all IDs present in orthologous groups ()

 Title   : addDaughterNonGroupSequences
 Usage   : $self->summarizer->summarizeInformation($self->nodeObject);
 Function: Summarizes two daughter nodes 
 Returns : 1 on success
 Args: -

=cut
sub addDaughterNonGroupSequences{
        my ($self,$daughterNode) = (@_);
        my $sequences2add;
        my %all_used_IDs = ();
        $daughterNode->orthologGroups->getAllUsedIDS($daughterNode,\%all_used_IDs);
        if(!keys(%all_used_IDs)){
                ERROR("Could not get all used IDs for node ".$daughterNode->name."\n");
                exit;
        }
        #print Dumper %all_used_IDs;
        # Get all sequence / cladeSequences
        my %ID2sequencesHash;
	#print keys(%ID2sequencesHash);
	#print "\n-------------------\n";
        $daughterNode->getDaughterSequencesById(\%ID2sequencesHash);
	#print keys(%ID2sequencesHash);
	#print "\n-------------------\n";
        #$daughterNode->get_cladeSequencesByID(\%ID2sequencesHash);
        if(!keys(%ID2sequencesHash)){
                ERROR("Could not get clade sequence for  ".$daughterNode->name."\n");
                exit;
        }
        #print Dumper %ID2sequencesHash;
        foreach my $current_id(keys(%ID2sequencesHash)){
                #print "\t\tchecking $current_id\n";
                if(!exists $all_used_IDs{$current_id}){
                        $sequences2add .= ">$current_id\n".$ID2sequencesHash{$current_id}."\n";
                        #print "\tadd $current_id\n";
                }
        }
        #print "\t\tattaching sequences to ".$daughterNode->fileInformation->sequenceSearchInputFile."\n";
        attach_to_file({file_name => $daughterNode->fileInformation->sequenceSearchInputFile, text =>  $sequences2add}) if $sequences2add ne '';
        #exit;
                return 1;
}

=item addDaughterPredictions()
        to be written

 Title   : addDaughterPredictions
 Usage   : $self->summarizer->summarizeInformation($self->nodeObject);
 Function: Summarizes two daughter nodes 
 Returns : 1 on success
 Args: -

=cut
sub addDaughterPredictions{
        #### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $all_used_IDs_href = $arg_ref->{all_used_IDs};
	my $groupsFile      = $arg_ref->{groupsFile};
	my $daughterNode = $arg_ref->{daughterNode};
	my %og_assignments_of_daughterNode;
        my $groups2add;
        if(! -e $daughterNode->orthologGroups->groupsFile){
                # This is a single species then # we do nothing
                DEBUG("\thas no predition file\n");
        }
        else{
                DEBUG("\t\tget predictions for ".$daughterNode->name." from :".$daughterNode->orthologGroups->groupsFile."\n");
                $daughterNode->orthologGroups->get_group2IDsHash($daughterNode,\%og_assignments_of_daughterNode);
        }
        if(!keys(%og_assignments_of_daughterNode)){
                DEBUG("\t\tNo additional predictions to add\n");
                return 1;
        }
                #print Dumper %og_assignments_of_daughterNode;
                #print "\t\tread ".keys(%og_assignments_of_daughterNode)." assignments from ".$daughterNode->name."\n";

                ALL_OG:
                foreach my $og(sort keys(%og_assignments_of_daughterNode)){
                        my $copy_og = 1;
                        #print "checking ".$og_assignments_of_daughterNode{$og}."\n";
                        foreach my $curr_id(split(/,/,$og_assignments_of_daughterNode{$og})){
                                if(exists $all_used_IDs_href->{$curr_id} || exists $all_used_IDs_href->{$og}){
                                        $copy_og = 0;
                                        #if(exists $all_used_IDs_href->{$curr_id}){print "\t\tnot adding ".$og_assignments_of_daughterNode{$curr_id}."\n";}
                                        #if(exists $all_used_IDs_href->{$og}){print "\t\tnot adding $og: ".$og_assignments_of_daughterNode{$og}."\n";}
                                        next ALL_OG;	
                                }
                        }
                        $groups2add .= "$og:".$og_assignments_of_daughterNode{$og}."\n" if($copy_og);
                        DEBUG("\tAdding $og:".$og_assignments_of_daughterNode{$og});
                }
                attach_to_file({file_name => $groupsFile, text =>  $groups2add}) if (defined($groups2add) && $groups2add ne '');
                return 1;
}

sub addTreeDaughterPredictions{
        #### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $all_used_IDs_href = $arg_ref->{all_used_IDs};
	my $groupsFile      = $arg_ref->{groupsFile};
	my $daughterNode = $arg_ref->{daughterNode};
	my %og_assignments_of_daughterNode;
        my $groups2add;
        if(! -e $daughterNode->orthologGroups->OGTreeFile){
                # This is a single species then # we do nothing
                DEBUG("\thas no predition file\n");
        }
        else{
                DEBUG("\t\tget predictions for ".$daughterNode->name." from :".$daughterNode->orthologGroups->OGTreeFile."\n");
        		# Read all daughter OGs
                #$daughterNode->orthologGroups->get_group2IDsHash($daughterNode,\%og_assignments_of_daughterNode);
                $daughterNode->orthologGroups->get_treeString4ID($daughterNode,\%og_assignments_of_daughterNode);
        }
        if(!keys(%og_assignments_of_daughterNode)){
                DEBUG("\t\tNo additional predictions to add\n");
                return 1;
        }
             ALL_OG:
                foreach my $og(sort keys(%og_assignments_of_daughterNode)){
                        #print "adding $og?..";
                        my $copy_og = 1;
                        if(exists $all_used_IDs_href->{$og})
                        {
                        	#print "\tno\n";
                        	#;
                        }
                        else
                        {
                        	#print "\tyes\n";
                        	
                        	$groups2add .= "HieranoidOG-$og:(".$og_assignments_of_daughterNode{$og}{'treestring'}."$og:".$og_assignments_of_daughterNode{$og}{'bitscore'}.")root;\n";
                        	$all_used_IDs_href->{$og} = 1;
                        }
                }
                attach_to_file({file_name => $groupsFile, text =>  $groups2add});
                return 1;
}

=item addPredictions()
        Retrieves all IDS from current ortholog groups
        # adds allowed groups from left/right daughter

 Title   : addPredictions
 Usage   : $self->summarizer->summarizeInformation($self->nodeObject);
 Function: Summarizes two daughter nodes 
 Returns : 1 on success
 Args: -

=cut
sub addPredictions{
        my ($self,$innerNode) = (@_);

        # get all current groups : as predicted by Inparanoid
        my %all_used_IDs = ();

        $innerNode->orthologGroups->getAllUsedIDS($innerNode,\%all_used_IDs);
        if(!keys %all_used_IDs){
                ERROR("\tProblem getting current orthologous groups\n"); 
                exit;
        }
        #my $cat_copy_cmd = "cp ".$innerNode->orthologGroups->OGTreeFile."  ".$innerNode->orthologGroups->OGTreeFile."_copy";
        #`$cat_copy_cmd`;
        my $cat_cmd = "cat ".$innerNode->orthologGroups->OGTreeFile;
        my @all_current_ogs = `$cat_cmd`;
        my %all_used_IDs_Nodes;
        foreach my $treeString(@all_current_ogs)
        {
        	chomp($treeString);
        	#print "tree string before: $treeString\n";
        	$treeString =~ /HieranoidOG-\w*:(.*)/;
        	$treeString = $1;
        	#print "tree string: $treeString\n";
        	my $tree = Bio::Phylo::IO->parse(
				    -format => 'newick',
    				-string => $treeString)->first;
    		my $root = $tree->get_root;
        	my @terminals = @{ $root->get_terminals };
        	my @internals = @{ $root->get_internals };
        	foreach(@terminals){
        		my $name = $_->id;
        		$all_used_IDs_Nodes{$name} = 1;
        	}
        	foreach(@internals){
        		my $name = $_->id;
        		$all_used_IDs_Nodes{$name} = 1;
        	}
        	#%all_used_IDs_Nodes = map {$_->id => 1} @terminals;
        	#%all_used_IDs_Nodes = map {$_->id => 1} @internals;
        	#print "\tnodes: ".join(",",map {$_->id => 1} @internals)." and ".join(",",map {$_->id => 1} @terminals)."\n";
        }
        #my $saveText;
        #foreach(keys %all_used_IDs_Nodes)
        #{
        #	$saveText .= $_." -> 1\n";
        #}
        #attach_to_file({file_name => $innerNode->orthologGroups->OGTreeFile."_hash", text => $saveText});
			#print Dumper %all_used_IDs_Nodes;
			#exit;
        #print "\tadding predictions for node ".$innerNode->leftDaughter->name."\n";
    # left Daughter
        if($innerNode->leftDaughterType eq 'innerNode')
        {
                 #addDaughterPredictions({all_used_IDs => \%all_used_IDs,
                 #        groupsFile => $innerNode->orthologGroups->groupsFile,
                 #        daughterNode => $innerNode->leftDaughter});
                 addTreeDaughterPredictions({all_used_IDs => \%all_used_IDs_Nodes,
                         groupsFile => $innerNode->orthologGroups->OGTreeFile,
                         daughterNode => $innerNode->leftDaughter});
        }
         #print "\tadding predictions for node ".$innerNode->rightDaughter->name."\n";
    # right Daughter
        if($innerNode->rightDaughterType eq 'innerNode')
        {
                 #addDaughterPredictions({all_used_IDs => \%all_used_IDs,
                 #        groupsFile => $innerNode->orthologGroups->groupsFile,
                 #        daughterNode => $innerNode->rightDaughter});
                 addTreeDaughterPredictions({all_used_IDs => \%all_used_IDs_Nodes,
                         groupsFile => $innerNode->orthologGroups->OGTreeFile,
                         daughterNode => $innerNode->rightDaughter});
         }        ### SAVE CLUSTER INFORMATION
         return 1;
}
=item saveOrthologyPredictions()
        Retrieves all IDS from current ortholog groups
        # adds allowed groups from left/right daughter

 Title   : saveOrthologyPredictions
 Usage   : $self->summarizer->summarizeInformation($self->nodeObject);
 Function: Summarizes two daughter nodes 
 Returns : 1 on success
 Args: -

=cut
sub saveOrthologyPredictions{
        my ($self,$innerNode) = (@_);
        if($self->configuration->orthologGroupsFormat eq 'groupFile'){
        	$self->orthologyPredictions2OGs($innerNode); 
        }
        else
        {
        	$self->orthologyPredictions2OGsOrthoXML($innerNode);
        }
        $self->addPredictions($innerNode) if $self->configuration->addNonMatchingSequences eq 'true';
        
        $innerNode->orthologGroups->convert2expandedGroupFile($innerNode); 
}
=item orthologyPredictions2OGs()
        Retrieves all IDS from current ortholog groups
        # adds allowed groups from left/right daughter

 Title   : orthologyPredictions2OGs
 Usage   : $self->summarizer->summarizeInformation($self->nodeObject);
 Function: Summarizes two daughter nodes 
 Returns : 1 on success
 Args: -

=cut
sub orthologyPredictions2OGs{
        my ($self,$innerNode) =(@_);
        my (%og_assignments1_hash,%og_assignments2_hash);
        my $mapping_cat_cmd = "cat ";
        DEBUG("\t\tSummarizing node information: ".$innerNode->name."\n");
        # e.g. ! Homininae.groups.txt
        if(!-e $innerNode->orthologGroups->orthologyPredictionFile && !-s $innerNode->orthologGroups->orthologyPredictionFile){
                WARN("\t\tNo orthology predictions found\n");
                return 0;
        }
        # e.g. Homininae.groups.txt
        if(-e $innerNode->orthologGroups->originalOrthologyPredictionFile && -s $innerNode->orthologGroups->originalOrthologyPredictionFile){
                DEBUG("\tconverted orthology predictions already exist. Skipping\n");
                return 1;
        }
        # Make a copy of original Inparanoid predictions
        copy($innerNode->orthologGroups->orthologyPredictionFile, $innerNode->orthologGroups->originalOrthologyPredictionFile) || die "Could not copy to ".$innerNode->orthologGroups->originalOrthologyPredictionFile."\n";
        my %all_og2id_hash;
        my %all_id_mappings;
        my %ID2Tree_mappings;
        # Read Mapping files from daughter nodes?
        if($innerNode->leftDaughterType eq 'leaf' && $innerNode->rightDaughterType eq "leaf"){
                #DEBUG("no mappings to read (both are leafes)\n");
                #                exit;
        }
        if($innerNode->leftDaughterType eq 'innerNode'){
                #print "need og assignments from ".$innerNode->leftDaughter->name."\n";
                $innerNode->leftDaughter->orthologGroups->get_group2IDsHash($innerNode->leftDaughter,\%all_og2id_hash);
                $innerNode->leftDaughter->orthologGroups->get_group2IDsHash($innerNode->leftDaughter,\%og_assignments1_hash);
                $innerNode->leftDaughter->orthologGroups->get_group2expandedIDsHash($innerNode->leftDaughter,\%all_id_mappings);
                $innerNode->leftDaughter->orthologGroups->get_treeString4ID($innerNode->leftDaughter,\%ID2Tree_mappings);
              #  print Dumper %all_id_mappings;
                $mapping_cat_cmd .= " ".$innerNode->leftDaughter->orthologGroups->mappingsFile." ";
        }
        if($innerNode->rightDaughterType eq 'innerNode'){
                #print "need og assignments from ".$innerNode->rightDaughter->name." (".$innerNode->rightDaughter->orthologGroups->groupsFile.")\n";
                $innerNode->rightDaughter->orthologGroups->get_group2IDsHash($innerNode->rightDaughter,\%all_og2id_hash);
                $innerNode->rightDaughter->orthologGroups->get_group2IDsHash($innerNode->rightDaughter,\%og_assignments2_hash);
                $innerNode->rightDaughter->orthologGroups->get_group2expandedIDsHash($innerNode->rightDaughter,\%all_id_mappings);
                $innerNode->rightDaughter->orthologGroups->get_treeString4ID($innerNode->rightDaughter,\%ID2Tree_mappings);
              #  print Dumper %all_id_mappings;
              DEBUG("Finished reading expanded groups\n");
                $mapping_cat_cmd .= " ".$innerNode->rightDaughter->orthologGroups->mappingsFile." ";
        }
        #if(!keys(%all_og2id_hash) && !($innerNode->leftDaughterType eq 'leaf' && $innerNode->rightDaughterType eq "leaf")){
        #        ERROR("\tCould not find OG mapping files from daughter nodes\n");
        #        exit;
        #}
       
        # Adding all mappings         
        $mapping_cat_cmd .= " >> ".$innerNode->orthologGroups->mappingsFile;
        if(($innerNode->rightDaughterType eq 'innerNode' || $innerNode->leftDaughterType eq 'innerNode' )
        && (! -e $innerNode->orthologGroups->mappingsFile && ! -s $innerNode->orthologGroups->mappingsFile)){
                #print "concat mappings: $mapping_cat_cmd\n";
                system($mapping_cat_cmd);
        }
       
        my $previous_og = -1;
        my @array_of_ids_for_og;  
        my %og_to_ids;
        my $treeString;
        my $treeFile;
        my ($clade1,$clade2) = (-1,-1);
        my (@clade1species, @clade2species);
        my $previous_bitscore = -1;
        my ($InparanoidSQLFile, $GroupsFile, $GroupsFileResolved, $MappingsFile);
        open my $INPUT_FH, '<', $innerNode->orthologGroups->originalOrthologyPredictionFile or die "Couldn't open '".$innerNode->orthologGroups->originalOrthologyPredictionFile."' $!\n";
        while(<$INPUT_FH>){
                ## e.g.
                # Homininae1	11100	Pan_troglodytes.fa	1.000	ENSPTRP00000031469
				# Homininae1	11100	Homo_sapiens.fa	1.000	ENSP00000358400
				# Homininae2	7694	Pan_troglodytes.fa	1.000	ENSPTRP00000046739
				# Homininae2	7694	Homo_sapiens.fa	1.000	ENSP00000352925
                chomp;
                my ($og,$bit_score,$species,$bootstrap,$id,$bootstrap_percentage) = split(/\t/,$_);
                
                
                
                $clade1 = $species if $clade1 eq '-1';
                $previous_bitscore = $bit_score if $previous_bitscore == -1;
                $previous_bitscore = $bootstrap if $previous_bitscore == -1;
                
                if($clade1 ne $species){
                	$clade2 = $species;
                }
                $previous_og = $og if $previous_og eq '-1';
                if($og ne $previous_og){
                                WARN("NO MEMBERS FOR OG:  $previous_og: ".join(",",@array_of_ids_for_og)."\n") if(!scalar(@array_of_ids_for_og));
                        # convert any pseudospecies in species
                        my $new_cluster = $innerNode->name."$previous_og";
                        my $og_string = join(",",@clade1species).",".join(",",@clade2species);
                        $treeFile .= "HieranoidOG-$new_cluster:((".$og_string.")$new_cluster:$previous_bitscore)root;\n"; 
                        #print "$treeFile\n";
                        # REMEMBER OG2ID
                        foreach(@array_of_ids_for_og){
                                $og_to_ids{$_} = $previous_og;	
                        }
                # $og_to_ids{VP345_HUMAN} = Homininae3;
                        #print "array of ids for ogs ".join(",",@array_of_ids_for_og)."\n";
                # case: VP345_HUMAN, VP345_PANTR
                        if($innerNode->leftDaughterType eq "leaf" && $innerNode->rightDaughterType eq "leaf"){
                                ;#DEBUG("\t\tboth are leafes\n");
                        }
                        else{
                                # e.g. OG2: Homininae4  --> resolve it
                                my @to_resolve = (@array_of_ids_for_og);
                                my @resolved_ids = ();
                                #print "\t\tinner Nodes.trying to resolve ".join(",",@array_of_ids_for_og)."\n";
                                foreach(@to_resolve){
                                        #print "\tresolving single $_\n";
                                        if(exists $all_id_mappings{$_}){
                                                push(@resolved_ids, $all_id_mappings{$_});
                                                shift(@to_resolve);
                                        }
                                        else{
                                                #print "could not find $_ in mappings\n";
                                                push(@resolved_ids, $_);
                                        }
                                }
                                #$all_og2id_hash{$new_cluster} = join(",",@resolved_ids);
                                #print {$OUTPUT_CLUSTER_RESOLVED_FH} "$new_cluster:".join(",",@resolved_ids)."\n";
                                $MappingsFile .= "$new_cluster:".join(",",@resolved_ids)."\n";
                                #$GroupsFileResolved .= "$new_cluster:".join(",",@resolved_ids)."\n";
                                #print "mapping: $new_cluster:".join(",",@resolved_ids)."\n";
                        }
                        #$all_og2id_hash{$new_cluster} = join(",",@array_of_ids_for_og);
                        # save 
                        #print {$OUTPUT_CLUSTER_FH} "$new_cluster:".join(",",@array_of_ids_for_og)."\n";
                        $GroupsFile  .= "$new_cluster:".join(",",@array_of_ids_for_og)."\n";
			#print "****---------\n";
                        #print $GroupsFile;
			#print @array_of_ids_for_og;
                        #print "\n****--------\n";
                        #print "\n\n";
                        @array_of_ids_for_og = ();
                        $previous_og = $og;
                        @clade1species = ();
                        @clade2species = ();
                        $clade1 = $species;
                        
                        # Inner Node annotation: bit score
                        #$previous_bitscore = $bit_score;
                        # Inner Node annotation: bootstrap
                        $previous_bitscore = $bootstrap;
                        
                }
                # replacing ids
                if(exists $ID2Tree_mappings{$id}){
                	my $extended_id = $ID2Tree_mappings{$id}{'treestring'};
                	my $og_bitscore = $ID2Tree_mappings{$id}{'bitscore'};
                	#print "REPL: $id -> $extended_id\n";
                	push(@clade1species, "$extended_id$id:".$og_bitscore) if $clade1 eq $species;
                	push(@clade2species, "$extended_id$id:".$og_bitscore) if $clade2 eq $species;
                }
				else
				{
					push(@clade1species, "$id:$bootstrap") if $clade1 eq $species;
                	push(@clade2species, "$id:$bootstrap") if $clade2 eq $species;
				}                
				
                chomp($id);
                ###### CHECKING	
                #if(!exists $og_assignments1_hash{$id} && $innerNode->leftDaughterType eq "innerNode"){;#print("ID $id of species A consists of several ids, but could not resolve it. Might not be harmful\n");}
                #if(!exists $og_assignments2_hash{$id} && $innerNode->rightDaughterType eq "innerNode"){;#print("ID $id of species B consists of several ids, but could not resolve it. Might not be harmful\n");}
                ###### CHECKING END
                if(exists $og_assignments1_hash{$id} || exists $og_assignments2_hash{$id}){
 #                        print {$OUTPUT_FH}	$innerNode->name."$og\t$bit_score\t$species\t$bootstrap\t$id\n";
                        $InparanoidSQLFile .= $innerNode->name."$og\t$bit_score\t$species\t$bootstrap\t$id\n";						
			#print "id---->",$id,"\n";
                        push(@array_of_ids_for_og,$id);
                }
                else{
                        #IF Id is not present in og2id hash then it is either a single species or an error
                        # og2id hash looks like
                        # OG4 = Homo_sapiens1, Homo_sapiens2, Mus_musculus3
 #                        print {$OUTPUT_FH} $innerNode->name."$og\t$bit_score\t$species\t$bootstrap\t$id\n";	
                        $InparanoidSQLFile .= $innerNode->name."$og\t$bit_score\t$species\t$bootstrap\t$id\n";						
			#print "id---->",$id,"\n";
                        push(@array_of_ids_for_og,$id);
                }
        }

        #### LAST ENTRY - START
        if(!scalar(@array_of_ids_for_og)){
                #ERROR("NO MEMBERS FOR OG:  $previous_og: ".join(",",@array_of_ids_for_og));
                WARN("NO MEMBERS FOR OG:  $previous_og: ".join(",",@array_of_ids_for_og)."\n");
        }
        else{

                # convert any pseudospecies in species
                my $og_string = "(".join(",",@array_of_ids_for_og).")";
                my $new_cluster = $innerNode->name."$previous_og";
                 $og_string = join(",",@clade1species).",".join(",",@clade2species);
                $treeFile .= "HieranoidOG-$new_cluster:((".$og_string.")$new_cluster:$previous_bitscore)root;\n"; 
                        
                # REMEMBER OG2ID
                foreach(@array_of_ids_for_og){
                        $og_to_ids{$_} = $previous_og;	
                }
                # $og_to_ids{VP345_HUMAN} = Homininae3;
                #print "array of ids for ogs ".join(",",@array_of_ids_for_og)."\n";

                # case: VP345_HUMAN, VP345_PANTR
                if($innerNode->leftDaughterType eq "leaf" && $innerNode->rightDaughterType eq "leaf"){
                        ;#DEBUG("\t\tboth are leafes\n");
                }
                else{
                        # e.g. OG2: Homininae4  --> resolve it
                        my @to_resolve = (@array_of_ids_for_og);
                        my @resolved_ids = ();
                        #print "\t\tinner Nodes.trying to resolve ".join(",",@array_of_ids_for_og)."\n";
                        foreach(@to_resolve){
                                #print "\tresolving single $_\n";
                                if(exists $all_id_mappings{$_}){
                                        push(@resolved_ids, $all_id_mappings{$_});
                                        shift(@to_resolve);
                                }
                                else{push(@resolved_ids, $_);}
                        }
                        #        $all_og2id_hash{$new_cluster} = join(",",@resolved_ids);
                        #        print {$OUTPUT_CLUSTER_RESOLVED_FH} "$new_cluster:".join(",",@resolved_ids)."\n";
                        $MappingsFile .= "$new_cluster:".join(",",@resolved_ids)."\n";
                        #$GroupsFileResolved .= "$new_cluster:".join(",",@resolved_ids)."\n";
                        #print "mapping: $new_cluster:".join(",",@resolved_ids)."\n";
                }
                $all_og2id_hash{$new_cluster} = join(",",@array_of_ids_for_og);
                # save 
                #print {$OUTPUT_CLUSTER_FH} "$new_cluster:".join(",",@array_of_ids_for_og)."\n";
                $GroupsFile  .= "$new_cluster:".join(",",@array_of_ids_for_og)."\n";
        }       
        #print "\n\n";
        ##### LAST ENTRY - END
        close $INPUT_FH or die "Couldn't close '".$innerNode->orthologGroups->originalOrthologyPredictionFile."': $!\n";

        # Write output files
        # Converted Orthology predictions        
                write_to_file({file_name => $innerNode->orthologGroups->orthologyPredictionFile, text => $InparanoidSQLFile});
        # Converted Orthology predictions       
		#print  "\ntext--->",text => $GroupsFile,"\n";
		write_to_file({file_name => $innerNode->orthologGroups->groupsFile, text => $GroupsFile});
        # Tree file        
                write_to_file({file_name => $innerNode->orthologGroups->OGTreeFile, text => $treeFile});
    #    # Converted Orthology predictions        
               # write_to_file({file_name => $innerNode->orthologGroups->expandedGroupsFile, text => $GroupsFileResolved});
        # Converted Orthology predictions        
                attach_to_file({file_name => $innerNode->orthologGroups->mappingsFile, text => $MappingsFile});
        # e.g. cp Homininae.groups Homininae.mappings
        if($innerNode->leftDaughterType eq "leaf" && $innerNode->rightDaughterType eq "leaf"){
                copy($innerNode->orthologGroups->groupsFile,$innerNode->orthologGroups->mappingsFile) or ERROR("Could not copy groupsFile to mappingsFile\n");
		}
		unlink($innerNode->orthologGroups->originalOrthologyPredictionFile);
        return 1;
}
=item write_to_file()
        Retrieves all IDS from current ortholog groups
        # adds allowed groups from left/right daughter

 Title   : write_to_file
 Usage   : $self->summarizer->summarizeInformation($self->nodeObject);
 Function: Summarizes two daughter nodes 
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
=item attach_to_file()
        Retrieves all IDS from current ortholog groups
        # adds allowed groups from left/right daughter

 Title   : attach_to_file
 Usage   : $self->summarizer->summarizeInformation($self->nodeObject);
 Function: Summarizes two daughter nodes 
 Returns : 1 on success
 Args: -

=cut
sub attach_to_file{
	#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $file_name = $arg_ref->{file_name};
	my $text      = $arg_ref->{text};
	if(!defined $text){
	        #print "\tTrying to add nothing to file $file_name\n";
	        return 0;
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
=item doAlignment()
        Retrieves all IDS from current ortholog groups
        # adds allowed groups from left/right daughter

 Title   : doAlignment
 Usage   : $self->summarizer->summarizeInformation($self->nodeObject);
 Function: Summarizes two daughter nodes 
 Returns : 1 on success
 Args: -

=cut
sub doAlignment{
        my ($arg_ref) = @_;
        my $alignmentCall    = $arg_ref->{alignmentCall};
        my $current_sequences = $arg_ref->{current_sequences};
        my $alignment_tmp_file    = $arg_ref->{alignment_tmp_file};
        
        my $seq_in_file    = "sequences.tmp";
        write_to_file({file_name => $seq_in_file, text =>  $current_sequences});

        my $alignment_call = "muscle -maxiters 1 -quiet -diags -sv -distance1 kbit20_3 -in $seq_in_file -out $alignment_tmp_file";
        #print "$alignment_call\n";
        system($alignment_call);

        # temorary change to read from source file
        #print "alignment call is : $alignmentCall\n";
        #open(my $MUSCLE_IN, "|  $alignmentCall  ") || die "muscle failed: $!\n";
        #print {$MUSCLE_IN} "$current_sequences\n";
        #close($MUSCLE_IN);
        # SAVE INFORMATION WHICH SEQUENCES BELONGS TO THIS OG
        ##
        # Make alignment
        ##
        if(! -e $alignment_tmp_file || !-s $alignment_tmp_file){
                ERROR("\tCould not build alignment. \nAlignment call: (cat $current_sequences | $alignmentCall)");
                #	die;
                return 0;
        }
        return 1;
}
=item doConsensus()
        Retrieves all IDS from current ortholog groups
        # adds allowed groups from left/right daughter

 Title   : doConsensus
 Usage   : $self->summarizer->summarizeInformation($self->nodeObject);
 Function: Summarizes two daughter nodes 
 Returns : 1 on success
 Args: -

=cut
sub doConsensus{
        my ($arg_ref) = @_;
        my $alignment_tmp_file    = $arg_ref->{alignment_tmp_file};
        my $threshold_percent = 50;
        my $alignment_in = Bio::AlignIO->new('-file' => "$alignment_tmp_file");
        my $alignment = $alignment_in->next_aln();
        if(!defined($alignment) || $alignment eq ''){
                ERROR("could not compute consensus string for alignment ($alignment_tmp_file)\n");
                print "could not compute consensus string for alignment ($alignment_tmp_file)\n";
                return undef;
		}
		my $consensus_string = $alignment->consensus_string($threshold_percent);
        # make sure only allowed symbols are in there
        $consensus_string =~ s/\?|\*//g;
        if($consensus_string eq ""){
                ERROR("could not compute consensus string for alignment ($alignment_tmp_file)\n");
                return undef;
        }
        else{
                return $consensus_string;
        }
}

1;
