package SimilaritySearch;
=head1 SimilaritySearch

hieranoid::Comparison::SimilaritySearch - Container of taxon objects

=head1 SYNOPSIS

 Abstract class

=head1 DESCRIPTION

Abstract class to handle how orthologous group are summarized

=head1 METHODS

=cut
use Moose::Role;
use Carp;
use English;
use File::Copy;
use File::Basename;
use Log::Log4perl qw(:easy);
use File::Temp qw/ tempfile tempdir /;
use Hieranoid::Tree::InnerNode;
use Data::Dumper;
use Benchmark;

#use lib::hieranoid::Comparison::OrthologySearch::SimilaritySearch::SequenceSearch;
#use lib::hieranoid::Comparison::OrthologySearch::SimilaritySearch::HierarchicalProfileSearch;




## ATTRIBUTES
has 'nodeObject', is => 'rw', isa => 'Object';
has 'configuration', is => 'rw', isa => 'Object';
has 'resultsObject', is => 'rw', isa => 'Object';

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
            DEBUG("\t\tInitialising Similarity search object created\n");

            my ($leftDaughter, $rightDaughter) = ($self->nodeObject->leftDaughter->name,$self->nodeObject->rightDaughter->name);
            if($self->configuration->use_outgroup && $self->nodeObject->outgroupDaughter){
                    #print "\tusing outgroup\n";
                    if(!$self->nodeObject->outgroupDaughter){
                            print "\tProblem with outgroup. Object undefined\n";
                    }
                    my $outgroup = $self->nodeObject->outgroupDaughter->name;
                    
                    $self->resultsObject(ResultsSimilaritySearch->new(
                            simAA => $self->nodeObject->fileInformation->outputDirectory."/".$leftDaughter."-".$leftDaughter,
                            simAB => $self->nodeObject->fileInformation->outputDirectory."/".$leftDaughter."-".$rightDaughter,
                            simBA => $self->nodeObject->fileInformation->outputDirectory."/".$rightDaughter."-".$leftDaughter,
                            simBB => $self->nodeObject->fileInformation->outputDirectory."/".$rightDaughter."-".$rightDaughter,
                            simAC => $self->nodeObject->fileInformation->outputDirectory."/".$leftDaughter."-".$outgroup,
                            simBC => $self->nodeObject->fileInformation->outputDirectory."/".$rightDaughter."-".$outgroup,
                    ));    
                    #print Dumper $self->resultsObject;
                    #exit;        
            }
            else{
                    $self->resultsObject(ResultsSimilaritySearch->new(
                            simAA => $self->nodeObject->fileInformation->outputDirectory."/".$leftDaughter."-".$leftDaughter,
                            simAB => $self->nodeObject->fileInformation->outputDirectory."/".$leftDaughter."-".$rightDaughter,
                            simBA => $self->nodeObject->fileInformation->outputDirectory."/".$rightDaughter."-".$leftDaughter,
                            simBB => $self->nodeObject->fileInformation->outputDirectory."/".$rightDaughter."-".$rightDaughter,
                    ));
            }
                    
}

=item start()

Serial mode:
            1. Blast Pass
            2. Pass either Blast or Profile search

Parallel mode
            splitting up the query sequence, i.e. taking every ith sequence
            blasting against whole db


 Title   : start
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub start{
              my ($self,$jobNumber) = (@_);
              DEBUG("\t\tPROFILE-search started\n");
              print "\t\tSimilarity-search started\n";
              #put the code here!
              my $debug_output_file;
              my $isMasterProcess = 0;
              my $addAllowedQueryDatabaseEntities = 0;
             
           # Master process or not?  
              if($self->configuration->available_nodes > 1 && !defined($self->configuration->jobNumber)){
                      $isMasterProcess = 1; 
              }
              #print "is master process: $isMasterProcess for $jobNumber ".$self->configuration->jobNumber."\n";
              #exit;
              
              # A - left Daughter in current Comparison
              # B - right Daughter in current Comparison
              # C - outgroup Daughter in current Comparison
              my %comparisons_hash = (
                      "A-B" => {
                                'outputFile' => $self->resultsObject->simAB,
                                'leftDaughter' => $self->nodeObject->leftDaughter,
                                'leftDaughterType' => $self->nodeObject->leftDaughterType,
                                'rightDaughter' => $self->nodeObject->rightDaughter,
                                'rightDaughterType' => $self->nodeObject->rightDaughterType,
                      },
                      "B-A" => {
                                  'outputFile' => $self->resultsObject->simBA,
                                  'leftDaughter' => $self->nodeObject->rightDaughter,
                                  'leftDaughterType' => $self->nodeObject->rightDaughterType,
                                  'rightDaughter' => $self->nodeObject->leftDaughter,
                                  'rightDaughterType' => $self->nodeObject->leftDaughterType,
                      },
                      "A-A" => {
                                    'outputFile' => $self->resultsObject->simAA,
                                    'leftDaughter' => $self->nodeObject->leftDaughter,
                                    'leftDaughterType' => $self->nodeObject->leftDaughterType,
                                    'rightDaughter' => $self->nodeObject->leftDaughter,
                                    'rightDaughterType' => $self->nodeObject->leftDaughterType,
                      },
                      "B-B" => {
                                      'outputFile' => $self->resultsObject->simBB,
                                      'leftDaughter' => $self->nodeObject->rightDaughter,
                                      'leftDaughterType' => $self->nodeObject->rightDaughterType,
                                      'rightDaughter' => $self->nodeObject->rightDaughter,
                                      'rightDaughterType' => $self->nodeObject->rightDaughterType,
                      },
                      "A-C" => {
                                          'outputFile' => $self->resultsObject->simAC,
                                          'leftDaughter' => $self->nodeObject->leftDaughter,
                                          'leftDaughterType' => $self->nodeObject->leftDaughterType,
                                          'rightDaughter' => $self->nodeObject->outgroupDaughter,
                                          'rightDaughterType' => $self->nodeObject->outgroupDaughterType,
                      },
                      "B-C" => {
                                        'outputFile' => $self->resultsObject->simBC,
                                        'leftDaughter' => $self->nodeObject->rightDaughter,
                                        'leftDaughterType' => $self->nodeObject->rightDaughterType,
                                        'rightDaughter' => $self->nodeObject->outgroupDaughter,
                                        'rightDaughterType' => $self->nodeObject->outgroupDaughterType,
                      }      
                                           
              );
              my @comparison_order = ("A-B","B-A","A-A","B-B");
              if($self->configuration->use_outgroup && $self->nodeObject->outgroupDaughter){
                      push(@comparison_order, ("A-C","B-C"));
              }
              #@comparison_order = ("A-C");
              #print "comparisons ".join(",",@comparison_order)."\n";
              #exit;
              #my @comparison_order = ("B-B");
              
              my %allowedQueryDatabaseEntities = ();
              my $fileWithPrefilterHits;
              my $directoryWithResultsFiles = $self->nodeObject->fileInformation->outputDirectory;

            Comparison:
            foreach my $currentComparison(@comparison_order){
                        # A - B
                        my %singleSimilaritySearchResults = ();
                        # allowedQueryDatabaseEntities = 
                        #         'CNNM2_F2_MACMU' => 1
                        #         'CNNM4_HUMAN' => 1
                        my $similaritySearchOutput = $comparisons_hash{$currentComparison}{outputFile};
                        $fileWithPrefilterHits = $similaritySearchOutput."_prefilter.tab";
                        my ($leftDaughter,$rightDaughter) = ($comparisons_hash{$currentComparison}{leftDaughter},
                                                                  $comparisons_hash{$currentComparison}{rightDaughter});
                        my ($leftDaughterType,$rightDaughterType) = ($comparisons_hash{$currentComparison}{leftDaughterType},
                                                                  $comparisons_hash{$currentComparison}{rightDaughterType});
                                                                  DEBUG("\tchecking existence of $similaritySearchOutput\n");
                        if(!-e $similaritySearchOutput && !-s $similaritySearchOutput){
                                    print "\t current Comparison $currentComparison (".basename($similaritySearchOutput).")\n";
                                    
                                    
                        # client process      
                                    if($jobNumber)
                                    {
                                                if(!keys(%allowedQueryDatabaseEntities))
                                                {
                                                            $addAllowedQueryDatabaseEntities = 1;
                                                }
                                                #print "\tParallel mode, job number $jobNumber. filtered Blast search\n";
                                                INFO("\tParallel mode, job number $jobNumber. filtered Blast search");
                                                if(!$self->singleBlastSearch($similaritySearchOutput,$leftDaughter,
                                                            $rightDaughter,
                                                            \%singleSimilaritySearchResults,
                                                            \%allowedQueryDatabaseEntities,
                                                            $jobNumber))
                                                            {
                                                                    ERROR("Could not find hits ".$leftDaughter->name." - ".$rightDaughter->name.". Exiting\n");
                                                                    return 0;
                                                            }
                                                # Change output file name for sub-process
                                                $similaritySearchOutput .= ".job".$jobNumber;
                                                #print "read ".keys(%singleSimilaritySearchResults)." allowed query sequences \n";
                                                #exit;
                                    }           
                        # Master process  || single node    
                                    else
                                    {
                                            if(!$isMasterProcess)
                                            {
                                                #print "\tMaster mode. 1.Blast pass....\n";
                                                INFO("\tNon-Master mode. 1.Blast pass....");
                                                #next;
                                                # the master process computes all prefilter files A-A,A-B,B-A,B-B
                                                # this way, the sub-processes do not have to do it

                                                if(!$self->singleBlastSearch($similaritySearchOutput,
                                                                $leftDaughter,
                                                            $rightDaughter,
                                                            \%singleSimilaritySearchResults,
                                                            \%allowedQueryDatabaseEntities))
                                                            {
                                                                    ERROR("Could not find hits ".$leftDaughter->name." - ".$rightDaughter->name.". Exiting\n");
                                                                    return 0;
                                                            }
                                            }       
                             
                                    }

                                    #if(!keys(%allowedQueryDatabaseEntities)){DEBUG("\t\tno allowed keys\n");}
                                    #else{DEBUG("\t\thas ".keys(%allowedQueryDatabaseEntities)." allowed keys\n");}
                # 2. pass
                                    if(!$isMasterProcess){
                                                #print "\tusing 1 node\n";
                                                # always do sequence vs sequence search for leaf vs leaf
                                                if($self->nodeObject->childrenType eq "leaf_leaf")
                                                {
                                                        #print "\tleaf - leaf --> sequence\n";
                                                        #exit;
                                                        DEBUG("\tStarting 2.pass and saving in $similaritySearchOutput");
                                                        $self->do_ublastxml_2pass($similaritySearchOutput, 
                                                                $leftDaughter,$rightDaughter,
                                                                $leftDaughterType,$rightDaughterType,
                                                                \%singleSimilaritySearchResults,\%allowedQueryDatabaseEntities);
                                                }
                                                else
                                                {
                                                            if($self->configuration->summarizeInformation =~ /consensus/)
                                                            {
                                                                        #print "\tconsensus --> sequence\n";
                                                                        #exit;
                                                                        $self->do_ublastxml_2pass($similaritySearchOutput, 
                                                                                    $leftDaughter,$rightDaughter,
                                                                                    $leftDaughterType,$rightDaughterType,
                                                                                    \%singleSimilaritySearchResults,\%allowedQueryDatabaseEntities);
                                                            }
                                                            else
                                                            {
                                                                        if($self->configuration->profileSearch eq 'hhsearch')
                                                                        {
                                                                                    #print "\telse --> profile\n";
                                                                                    #exit;
                                                                                    $self->doHHsearchProfileSearch($similaritySearchOutput, 
                                                                                                $leftDaughter,$rightDaughter,
                                                                                                $leftDaughterType,$rightDaughterType,
                                                                                                \%singleSimilaritySearchResults,\%allowedQueryDatabaseEntities);
                                                                        }
                                                                        else
                                                                        {
                                                                                    #print "\telse --> sequence\n";
                                                                                    #exit;
                                                                                    $self->doHHBlitsProfileSearch($similaritySearchOutput, 
                                                                                            $leftDaughter,$rightDaughter,
                                                                                            $leftDaughterType,$rightDaughterType,
                                                                                            \%singleSimilaritySearchResults,\%allowedQueryDatabaseEntities);
                                                                        }
                                                            }
                                                }
                                    }
                                    else{       
                                            # just continue
                                            # sub-processes will read prefilter results and pick what they need
                                            INFO("\tParallel mode, master. do nothing");
                                            print "\t\t\tParallel mode (master). Do nothing\n";
                                            next;
                                    }
                        }
                        else{
                                DEBUG("\t\tOutput file already exists ($similaritySearchOutput)\n");
                        }
            }
            if($isMasterProcess){  
                        #print "\tusing x nodes\n";
                        # parallel
                        # iterate over %singleSimilaritySearchResultsAB
                        my @jobStarted_files;
                        my @jobFinished_files;
                        my $all_finished = 0;
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
                                                $parallelSimilaritySearch = $self->configuration->perl." $0 -j $_ -n ".$self->nodeObject->name." -a orthologySearch -c ".$self->configuration->configurationFile."";
                                    }
                                    print "\t$parallelSimilaritySearch\n";
                                    system("$parallelSimilaritySearch &");
                                    #exit;
                                    
                                    push(@jobStarted_files,$self->nodeObject->fileInformation->outputDirectory."/process_orthologySearch_started.$_");
                                    push(@jobFinished_files,$self->nodeObject->fileInformation->outputDirectory."/process_orthologySearch_finished.$_");
                        }
                        if(!scalar(@jobStarted_files) || !scalar(@jobFinished_files)){
                                    ERROR("\tCould not initialize start/finish files for parallel similarity search\n");
                                    exit;
                        }
                              # initially created array of output files from sub-processes
                              # whenever a sub-process successfully finishes, the corresponding
                              # output file will be deleted from the array
                              # empty array --> finished
                        while(@jobFinished_files){
                                    # writes status file in the end
                                    for(my $i=0;$i<scalar(@jobFinished_files); $i++){
                                                my $currentFile = $jobFinished_files[$i];
                                                #DEBUG("\tchecking $currentFile\n");
                                                if(-e $currentFile){
                                                            splice(@jobFinished_files, $i, 1); 
                                                            DEBUG("\t job has finished (file does exist)\n");
                                                }
                                                else{
                                                            ;#DEBUG("\t job has not finished (file does not exist yet)\n");
                                                }
                                    }
                                    last if !@jobFinished_files;                                        
                                    #DEBUG("\tCurrent running jobs: ".join(",",@jobFinished_files)."");

                                    # Check if errors have occured
                                    #my $findErrorFiles_cmd = "find $directoryWithResultsFiles/process_orthologySearch_error*";
                                    #DEBUG("\t$findErrorFiles_cmd\n");
                                    #my @findErrorFiles = `$findErrorFiles_cmd`;
                                    #if(@findErrorFiles){
                                    #            print "\t At least one of the sub-proccesses reported an error.\n";
                                    #            print "\t\tMaster process will stop here\n";
                                    #            exit;
                                    #}
                                    # check used time
                                    my $wallTime = $self->configuration->wallTime;
                                    my $findFiles_cmd = "find $directoryWithResultsFiles/process_orthologySearch_started* -mmin \+$wallTime";
                                    #DEBUG("\t$findFiles_cmd\n");
                                    my @findOldFiles = `$findFiles_cmd`;
                                    if(@findOldFiles){
                                                print "\t There are files running longer than $wallTime minutes\n";
                                                print "\t\tMaster process will stop here\n";
                                                exit;
                                    }
                                    #DEBUG("\twaiting\n");
                                    sleep(60);
                        }
                        #Summarize results A-A,A-B,B-A,B-B
                        Comparison:
                        foreach my $currentComparison(@comparison_order){
                                    print "\tFinally checking output from  ".$currentComparison."\n";
                                    # have to summarize results first
                                    # Now combine information
                                    # first we save all in a temporary file
                                    print "\tYes, we are the master process\n";
                                    my $temp_file = $comparisons_hash{$currentComparison}{outputFile}."_temp";
                                    my $findResultFiles_cmd = "find ".$comparisons_hash{$currentComparison}{outputFile}.".job*";

                                    #print "\t$findResultFiles_cmd\n";
                                    my @findResultFiles = `$findResultFiles_cmd`;
                                    if(!scalar(@findResultFiles)){
                                                print "\tCould not find result files: ".$comparisons_hash{$currentComparison}{outputFile}.".job* \n";
                                                ERROR("\tCould not find result files: ".$comparisons_hash{$currentComparison}{outputFile}.".job* \n");
                                    }
                                    print "\tBy now, results should be in ".$temp_file."\n";
                                    my $cat_cmd = "cat ".$comparisons_hash{$currentComparison}{outputFile}.".job* > $temp_file";
                                    print $cat_cmd."\n";
                                    system($cat_cmd);
                                    # then we sort it
                                    my $sort_cmd = "sort $temp_file > ".$comparisons_hash{$currentComparison}{outputFile};
                                    print "\tsorted and in: ".$comparisons_hash{$currentComparison}{outputFile}."\n";
                                    #print "\tconcatenating results: $cat_cmd\n$sort_cmd";
                                    system($sort_cmd);
                                    # Results exist?
                                    if(!-e $comparisons_hash{$currentComparison}{outputFile} || ! -s $comparisons_hash{$currentComparison}{outputFile}){
                                                ERROR("Could not find hits ".$comparisons_hash{$currentComparison}{leftDaughter}->name." - ".$comparisons_hash{$currentComparison}{rightDaughter}->name.". Exiting\n");
                                                return 0;
                                    }
                        }
                        print "Clean up old files\n";
                        foreach(glob("$directoryWithResultsFiles/*job*")){
                                    print "\tTrying to delete $_\n";
                                    unlink $_ or print "Error removing file \"$_\": $!\n";
                        }
                        unlink glob $directoryWithResultsFiles."/process_*";
            }
            #unlink glob $directoryWithResultsFiles."/process_finished*";
            # lets make a hash of all query sequences
            #exit;
            return 1;
                                       
            INFO("\t\t Similarity search finished\n");
}

=item singleBlastSearch()

to be written...

 Title   : singleBlastSearch
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub singleBlastSearch{
	#### PARAMETER VARIABLES
        my ($self, $output_file,$leftDaughter,$rightDaughter,$singleSimilaritySearchResults,$allowedQueryDatabaseEntities,$jobNumber) = (@_);
	my $query         = $leftDaughter->fileInformation->sequenceSearchInputFile;
	my $database      = $rightDaughter->fileInformation->sequenceSearchInputFile;
	my ($queryName,$databaseName) = (basename($query),basename($database));
	my $segmasker = $self->configuration->segmasker;
	# using tmp file to avoid heavy writing to afs
	my (undef,$ublast_output_tmp) = tempfile( SUFFIX => '.prefilter.out', UNLINK => 0);
        my (undef,$query_filtered_file) = tempfile( SUFFIX => '.query_filtered.fa', UNLINK => 0);
        my ($noQuerySequences,$noDBSequences);
        $noQuerySequences = `grep ">" $query -c`;
        chomp($noQuerySequences);
        $noDBSequences = `grep ">" $database -c`;
        chomp($noDBSequences);
        
        unlink($query_filtered_file) if -e $query_filtered_file;
        
        my $ublast_output = $output_file.".1.Pass.out";
        if($jobNumber){
            $ublast_output .= ".job$jobNumber";
        }
        my $ublastxml_output = $output_file."_ublast.xml";
        
        if(!-e $query || ! -s $query){
                        ERROR("ERROR. Fasta file empty/missing ($query)\n"); 
                        return 0;# "ERROR. Fasta file empty/missing ($query)\n";
        }
        ### here can we simply split the query database
        # iterate over fasta file and take every i.th sequence
        if($jobNumber){
                DEBUG("\tget filtered file\n");
            # get filtered file
            getFastaFileSplit({
                        fastaFile => $query,
                        fastaFileOut => $query_filtered_file,
                        jobNumber => $jobNumber,
                        availableNodes => $self->configuration->available_nodes
                        });
            # set variables accordingly
            if(!-e $query_filtered_file || ! -s $query_filtered_file){
                        ERROR("ERROR. Could not read assigned part of fasta file ($query_filtered_file)\n"); 
                        return 0;#die "ERROR. Could not read assigned part of fasta file ($query_filtered_file)\n";
            }           
            $query = $query_filtered_file;
        }
        my $addAllowedQueryDatabaseEntities = 1;
        $addAllowedQueryDatabaseEntities = 0 if(keys(%$allowedQueryDatabaseEntities));
        #print "\tchecking $ublast_output (tmp: $ublast_output_tmp)\n";
        #exit;
        # Taking time
        # start timer
        my $start_time =  new Benchmark;
        
         DEBUG("\tsimilarity search tools: ".$self->configuration->similaritySearchTool."\n");
         if($self->configuration->similaritySearchTool eq 'usearch' || $self->configuration->similaritySearchTool eq 'ublast'){
                if(!-e $ublast_output || !-s $ublast_output){
                         DEBUG("\tusing Usearch");
                         my (undef, $masked_query_database_file) = tempfile( SUFFIX => '.masked_query.db',UNLINK => 1);
                         unlink($masked_query_database_file) if -e $masked_query_database_file;
                         
                         my $segmasker_call = "$segmasker -in $query -out $masked_query_database_file -outfmt fasta";
                         #print $segmasker_call."\n";
                         system($segmasker_call);
                         if(!-e $masked_query_database_file || ! -s $masked_query_database_file){
                                 WARN("Could not mask low complexity regions in query database ".basename($query)." empty $masked_query_database_file (call: $segmasker_call)\n");
                                 #print "Could not mask low complexity regions in query database ".basename($query)." empty $masked_query_database_file (call: $segmasker_call)\n"; 
                                 $masked_query_database_file = $query;
                                 #exit;
                         }
                         $noQuerySequences = `grep ">" $masked_query_database_file -c`;
                         chomp($noQuerySequences);
                         if(!defined($noQuerySequences) || $noQuerySequences == 0)
                         {
                         	#ERROR("Empty query file $masked_query_database_file. Blast call: ".$self->configuration->usearch." --quiet --query $query  --db $database  --userout $ublast_output_tmp --userfields query+target+bits+ql+tl+qlo+qhi+tlo+thi+qs+ts ".$self->configuration->ublastParameters."\n");
                         	#return 0;
                         }
                        my $blast_call = $self->configuration->usearch." --quiet --query $query  --db $database  --userout $ublast_output_tmp --userfields query+target+bits+ql+tl+qlo+qhi+tlo+thi+qs+ts ".$self->configuration->ublastParameters;
                        #my $blast_call = $self->configuration->usearch."  --query $query  --db $database  --userout $ublast_output_tmp --userfields query+target+bits+ql+tl+qlo+qhi+tlo+thi+qs+ts ".$self->configuration->ublastParameters;
                        DEBUG("call: $blast_call");
                        #print "call: $blast_call\n";
                        my @blast_error = `$blast_call`;
			system($blast_call);                        
                        if(!-e $ublast_output_tmp || !-s $ublast_output_tmp){
                                print "\tPrefiltering blast resulted in no hits (file: $ublast_output_tmp)\n$blast_call\n";
                                ERROR("\tPrefiltering blast resulted in no hits (file: $ublast_output_tmp)\n$blast_call\n");
                                ERROR("Error was: ".join(",",@blast_error));
                                return 0;
                        }
                        my $blast2xml_call = $self->configuration->perl." blast2xml.pl -i $ublast_output_tmp | ".$self->configuration->perl." ./blast_parser.pl ".$self->configuration->similaritySearchCutoff." > $ublast_output";
                        DEBUG("$blast2xml_call\n");
                        `$blast2xml_call`;
                        #print "check now\n";
                        #exit;
                        DEBUG("\tblast2xml-call: $blast2xml_call");
                        #print "\tblast2xml-call: $blast2xml_call\n";
                        if(!-e $ublast_output || !-s $ublast_output){
                                #print "\tPrefiltering blast resulted in no hits (file: $ublast_output)\n";
                                ERROR("\tPrefiltering blast resulted in no hits (file: $ublast_output)\n");
                                return 0;
                        }
                }
                else
                {
                        print "\t\tOutput file already exist. Reading previous hits...\n";
                }
                ### Reading Hits from previous search 
                # save as 'query' = "hit1 hit2 hit3"
                        my @blast_parser_lines = `cat $ublast_output`;
                        my $noOfQuery = 0;
                        foreach my $aLine (@blast_parser_lines) {
                                chomp($aLine);
                                #	print $aLine."\n";
                                next if $aLine =~ /^query/;
                                my @words = split( /\s+/, $aLine );
                                        if ( exists( $singleSimilaritySearchResults->{ $words[0] } ) ) {
                                                $singleSimilaritySearchResults->{ $words[0] } = $singleSimilaritySearchResults->{ $words[0] } . " " . $words[1];
                                        }
                                        else{
                                                $singleSimilaritySearchResults->{ $words[0] } = $words[1];
                                        }
                                        if($addAllowedQueryDatabaseEntities){
                                                $allowedQueryDatabaseEntities->{$words[0]} = 1;
                                                $allowedQueryDatabaseEntities->{$words[1]} = 1;
                                        }
                                $noOfQuery++;
                        }
                        
        }
                 elsif($self->configuration->similaritySearchTool eq 'blast'){
                         DEBUG("\tusing Blast");
                         
                        system ($self->configuration->formatdb." -i $database") if !-e "$database.pin";
                        
                        #my $noQuerySequences = `grep ">" $query -c`;
                        #chomp($noQuerySequences);
                        #my $noDBSequences = `grep ">" $database -c`;
                        #chomp($noDBSequences);
                        
                        #print STDERR "\nStarting first BLAST pass for $query - $database on ";
                	system("date");
                	#print $self->configuration->blast." -C3 -F\"m S\" -i $query -d $database -p blastp -v $noDBSequences -b $noDBSequences -M BLOSUM62 -z 5000000 -m7 | ./".$self->configuration->blastParser." 40|\n";
                	#exit;
                	open FHR, $self->configuration->blast." -C3 -F\"m S\" -i $query -d $database -p blastp -v $noDBSequences -b $noDBSequences -M BLOSUM62 -z 5000000 -m7 | ./".$self->configuration->blastParser." 40|";
                        while (<FHR>) {
                		my $aLine = $_;
                		chomp ($aLine);
                		my @words = split (/\s+/, $aLine);

                		if ( exists( $singleSimilaritySearchResults->{ $words[0] } ) ) {
                                        $singleSimilaritySearchResults->{ $words[0] } = $singleSimilaritySearchResults->{ $words[0] } . " " . $words[1];
                                }
                                else{
                                        $singleSimilaritySearchResults->{ $words[0] } = $words[1];
                                }
                                if($addAllowedQueryDatabaseEntities){
                                        $allowedQueryDatabaseEntities->{$words[0]} = 1;
                                        $allowedQueryDatabaseEntities->{$words[1]} = 1;
                                }
                	}
                	close (FHR);
                }
                else{
                        die "wrong parameter for similarity search tool\n";
                }
                
        # Stop time
        my $end_time = new Benchmark;
        my $time_diff = timediff($end_time,$start_time);
        my $debugString = "1.Pass\t".$self->configuration->similaritySearchTool."\t$queryName ($noQuerySequences) - $databaseName ($noDBSequences):  ".timestr($time_diff, 'all');
        attach_to_file({file_name => $self->configuration->timeFile, text => $debugString."\n"});        

        if($jobNumber){print "\tfinished job $jobNumber. Output files should be in $ublast_output\n";}
        # Copy tmp results to regular file
        # Mainly for debugging reasons
        copy($ublast_output_tmp,$ublast_output ) || die "\tCould not copy $ublast_output_tmp to $ublast_output";
        if (!keys(%$singleSimilaritySearchResults)){
                ERROR("Could not parse hits from first Blast run (cat $ublastxml_output | ./blastParser ".$self->configuration->similaritySearchCutoff.")\n");
                #print "\tMaybe there were no hits\n";
                return 0;
                #exit;
        }
        unlink($ublast_output_tmp);
        unlink($ublast_output);
        unlink($query_filtered_file);
        return 1;
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
   	ERROR("Input file $sequence_file is empty or missing\n")   if ( !-e $sequence_file || !-s $sequence_file );
   	#print "\t Reading sequences from $sequence_file\n";
   	open my $SEQUENCE_FILE, '<', $sequence_file or croak "Couldn't open '$sequence_file': $OS_ERROR";

   	### READING FROM FILE
   	while ( my $line = <$SEQUENCE_FILE> ) {
   		chomp($line);
   		if ( $line =~ /^>/ ) {
   			my @words = split( /\s/, $line );
   			$seq_id = $words[0];
   			$seq_id =~ s/>//;
   #                        print "seq_id $seq_id\n" if $seq_id =~ /LCT/;
 #   			$sequence_href->{$seq_id} = q();
   			$sequence_href->{$seq_id} = ">$seq_id\n";
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
           ERROR("Could not read fasta entries from $sequence_file");
           exit;
       }
   	return 1;
}    


=item getHits4Process()

to be written...

 Title   : getHits4Process
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut

sub getHits4Process{
                my ($arg_ref) = @_;
                my $hitsFile    = $arg_ref->{hitsFile};
                my $Hits2compute_href    = $arg_ref->{Hits2compute};
                my $ownProcessId    = $arg_ref->{ownProcessId};
                my $NoProcesses = $arg_ref->{NoProcesses};
                my $addAllowedQueryDatabaseEntities = $arg_ref->{addAllowedQueryDatabaseEntities};

                my $allowedQueryDatabaseEntities = $arg_ref->{allowedQueryDatabaseEntities};
                if(!-e $hitsFile || !-s $hitsFile){
                        ERROR("Required prefilter file $hitsFile for process $ownProcessId not found");
                        return 0;
                }
                my %hits4query;
                my $lineCounter = 1;
                my $saveID = 0;
                my $previousQueryID = "NaN";
                my $start = 1;
                open my $hitsFile_FH, '<', $hitsFile or croak "Couldn't open '$hitsFile': ";
                my $noOfQuery = 0;
          	while(defined(my $aLine = <$hitsFile_FH>)) {
          		chomp($aLine);
          	        next if $aLine =~ /^query/;
          		#DEBUG("$ownProcessId: ".$aLine."...");
          		my @words = split( /\s+/, $aLine );
          		my ($query,$hit) = ($words[0],$words[1]);
          		$previousQueryID = $query if $previousQueryID eq "NaN";
          		#print "$previousQueryID = $query\n";
          		if(exists $hits4query{$query}{$hit}){
  		                #print "skipped duplicate $query - $hit\n";
  		                next;
  		        }
                        if($ownProcessId){
                                # print $noOfQuery." % ".$self->configuration->available_nodes." :".$noOfQuery % ($self->configuration->available_nodes +1)." == $jobNumber\n";
                                #print "(".$noOfQuery." % ".($NoProcesses + 1)." :".(($noOfQuery % $NoProcesses) +1)." == $ownProcessId\n";
                                # Note
                                # Increase available nodes by 1 to get the right number
                                # without increase
                                #22 % 3 :1 == 3
                                #23 % 3 :2 == 3
                                #24 % 3 :0 == 3
                                if($previousQueryID ne $query){$noOfQuery++;}
                                if((($noOfQuery % $NoProcesses)+1) == $ownProcessId){
                                        #DEBUG("\t\t$noOfQuery: Process allowed to take $query\n");
                                        if ( exists( $Hits2compute_href->{$query})){
                                                #print "\t\t\tsave: query ".$words[0]." hit: ".$Hits2compute_href->{ $words[0] } . " " . $words[1]."\n";
                                                $Hits2compute_href->{$query} = $Hits2compute_href->{$query}." ".$hit;
                                        }
                                        else{
                                                #print "\t\t\tsave: query ".$words[0]." hit: ".$words[1]."\n";
                                                $Hits2compute_href->{$query} = $hit;
                                        }
                                        if($addAllowedQueryDatabaseEntities){
                                                $allowedQueryDatabaseEntities->{$query} = 1;
                                                $allowedQueryDatabaseEntities->{$hit} = 1;
                                        }
                                        $saveID = 1;
                                }
                                else{
                                        #DEBUG("\t\t$noOfQuery: Process not allowed to take $query\n");
                                        $saveID = 0;
                                }

                        }
  		        else{
          		        if ( exists( $Hits2compute_href->{$query})) {
          			        $Hits2compute_href->{$query} = $Hits2compute_href->{$query}." ".$hit;
          		                }
          		        else{
                                        $Hits2compute_href->{$query} = $hit;
          		        }
          		        if($addAllowedQueryDatabaseEntities){
          		                $allowedQueryDatabaseEntities->{$query} = 1;
          		                $allowedQueryDatabaseEntities->{$hit} = 1;
  		                }
	                }
	                $hits4query{$query}{$hit} = 1;
	                #$noOfQuery++;
	                $previousQueryID = $query;
                        #print Dumper $Hits2compute_href;
                        # if $noOfQuery++ > 3;
                        #exit;
  		       # last if $noOfQuery++ > 50;
          	}
                close $hitsFile_FH or croak "Couldn't close '$hitsFile': ";
                #exit;
                return 1;

}


=item getFastaFileSplit()

to be written...

 Title   : getFastaFileSplit
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub getFastaFileSplit{
            #### PARAMETER VARIABLES
            my ($arg_ref) = @_;
            my $fastaFile = $arg_ref->{fastaFile};
            my $fastaFileOut      = $arg_ref->{fastaFileOut};
            my $jobNumber      = $arg_ref->{jobNumber};
            my $availableNodes      = $arg_ref->{availableNodes};
            # create one SeqIO object to read in,and another to write out
            my $seq_in = Bio::SeqIO->new(
                        -file   => "<$fastaFile",
                        -format => 'fasta'
            );
            my $seq_out = Bio::SeqIO->new(
                        -file   => ">$fastaFileOut",
                        -format => 'fasta'
            );
            my $noOfQuery = 0;
            my $noSequencesCopied = 0;
            #my $noSeqs = `grep ">" $fastaFile -c`;
            #chomp($noSeqs);
            #print "\treading $fastaFile ($noSeqs)\n";
            # write each entry in the input file to the output file
            while (my $inseq = $seq_in->next_seq)
            {
                        #print "$noOfQuery % ($availableNodes)+1  (".($noOfQuery % ($availableNodes)+1).")== $jobNumber ";
                        if(($noOfQuery % ($availableNodes)+1) == $jobNumber)
                        {
                                    $seq_out->write_seq($inseq);
                                    #print "$jobNumber ($fastaFile) -> ".$inseq->id."\t".$inseq->seq."\n";
                                    $noSequencesCopied++;
                 #                   print "take\n";
                        }
                        else
                        {
                #                    print "not take\n";
                        }
                        $noOfQuery++;
            }
            #print "$noSequencesCopied / $noOfQuery copied\n";
            #exit;
            if(!$noSequencesCopied)
            {
                        print "\tCould not filter sequence file\n";
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
