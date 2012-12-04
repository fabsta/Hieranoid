package ProfileSearch;
=head1 HierarchicalProfileSearch

hieranoid::Comparison::OrthologySearch::SimilaritySearch::HierarchicalProfileSearch - Container of taxon objects

=head1 SYNOPSIS

 use lib::hieranoid::Comparison;
 my $currentComparison = Comparison->new(nodeObject =>$innerNode,
						configuration => $hieranoidConfiguration);

=head1 DESCRIPTION

Compares orthologous groups between daughter nodes in a given phylogenetic tree

=head1 METHODS

=over

=cut
use Moose;
with 'SimilaritySearch';
with 'Profile';
use Carp;
use English;
use Data::Dumper;
use File::Copy;
use Log::Log4perl qw(:easy);
use File::Temp qw/ tempfile tempdir /;
use List::MoreUtils qw( any );
use Hieranoid::Tree::InnerNode;
use File::Basename;


=item doHHBlitsProfileSearch()

to be written...

 Title   : doHHBlitsProfileSearch
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut      
sub doHHBlitsProfileSearch(){
        #### PARAMETER VARIABLES
        my ($self, $output_file,$leftDaughter,$rightDaughter,$leftDaughterType,$rightDaughterType,$singleSimilaritySearchResults,$allowedQueryDatabaseEntities) = (@_);
        my $query         = $leftDaughter->fileInformation->sequenceSearchInputFile;
        my $database      = $rightDaughter->fileInformation->sequenceSearchInputFile;
        print "\t$query - $database\n";
        #my $analysis_step = $;
        my $score_cutoff  = $self->configuration->similaritySearchCutoff;
        my $usearch       = $self->configuration->usearch;
        my $hhblits       = $self->configuration->hhblits;
        #my $database_cs_file = $arg_ref->{database_cs_file};
        my $species          = $leftDaughter->name;
        my $species2         = $rightDaughter->name;
        #### TESTING ################################################################################
        ## CHECK DEFINEDNESS
        print "$query, $database, $output_file\n";
        croak q{One/several of the parameters for 'method' was/were not defined.\n}
        if any { !defined $_ } $query, $database, $output_file;
        my %sequencesA;
        my %query2resolvedID_mapping;
        my %database2resolvedID_mapping;
        my %id2SeqAli_mapping;

        my %ID2Alignment_mapping;
        my %sequencesB;
        my $seq;
        my $seqID;
        my $restrictSearch = 0;
        if(!keys(%$allowedQueryDatabaseEntities)){
                $restrictSearch = 1;
        }


        my ($fh2,$databaseHMMFile )       = tempfile( UNLINK => 1);
        my ($fh3, $tmp_alignment_file) = tempfile(UNLINK => 1);
        my ($fh4, $tmp_hmm_output_file) = tempfile(UNLINK => 1);
        my ( $fh5, $parser_output )       = tempfile(UNLINK => 1);
        my $queryHMMFile = $leftDaughter->fileInformation->profileFile;
        $databaseHMMFile = $rightDaughter->fileInformation->profileFile;
        my $hitCSDBFile = $rightDaughter->fileInformation->csdbFile;
        my $hitDBFile = $rightDaughter->fileInformation->hmmFile;

        my $noOfQuery = 1;
        my $totalNoQueries = keys(%$singleSimilaritySearchResults);
        %id2SeqAli_mapping = ();
        readSequenceFileIntoHash({sequence_file => $database, sequence_href => \%id2SeqAli_mapping});
        if($rightDaughterType eq 'innerNode'){
                $rightDaughter->orthologGroups->get_group2expandedIDsHash( $rightDaughter->orthologGroups,
                        $rightDaughter,
                        \%database2resolvedID_mapping);
                        $rightDaughter->get_alignmentsByID(\%id2SeqAli_mapping);
        }
                die "Could not read sequences in $database\n" if(!keys(%id2SeqAli_mapping));

                # Database
        print "\tCreating database...\n";
        if(!-e $databaseHMMFile && ! -s $databaseHMMFile) {

                        createHhblitsProfiles({
                                self => $self,
                                sequenceFile => $database,
                                configuration => $self->configuration,
                                name => basename($database),
                                hhm_file => $databaseHMMFile,
                                id2SeqAli_mapping => \%id2SeqAli_mapping,
                                hhblitsDatabase => 1,
                                hitCSDBFile => $hitCSDBFile,
                                hitDBFile => $hitDBFile,
                                });
           }
                        # create HHBlits database
                        %id2SeqAli_mapping = ();
                        # Query
                        # Read A and B
                        readSequenceFileIntoHash({sequence_file => $query, sequence_href => \%id2SeqAli_mapping});
                        print "read sequences in $database\n";
                        #print Dumper %sequencesB;
                        #exit;

                        if($leftDaughterType eq 'innerNode'){
                                $leftDaughter->orthologGroups->get_group2expandedIDsHash( $leftDaughter->orthologGroups,
                                        $leftDaughter,
                                        \%query2resolvedID_mapping);
                                        $leftDaughter->get_alignmentsByID(\%id2SeqAli_mapping);
                        }
                                die "Could not read sequences in $query\n" if(!keys(%id2SeqAli_mapping));
                               # my ($fh, $tmp_alignment_file) = tempfile();
                               # my ($fh2, $tmp_hmm_output_file) = tempfile();

                                foreach my $queryID(keys %id2SeqAli_mapping){

                                        #    if(!-e $queryHMMFile && ! -s $queryHMMFile) {
                                                print "\tCreating query...\n";
                                                buildHhblitsAlignment({configuration => $self->configuration,
                                                        og_name => $queryID,
                                                        alignment_string => $id2SeqAli_mapping{$queryID},
                                                        tmp_alignment_file => $tmp_alignment_file,
                                                        tmp_hmm_output_file => $tmp_hmm_output_file,
                                                        hhm_file => $queryHMMFile,
                                                        attach => 1});

                                                        #         }
                                                        print "\tNow performing HHblits search\n";
                                                        my $hhblitsCmd = "time $hhblits -i $tmp_hmm_output_file -o out -n 1 -d $databaseHMMFile\n";
                                                        print "$hhblitsCmd\n";
                                                        exit;
                                }
                                                return 1;
}
=item doHierarchicalhhblitsSearch()

to be written...

 Title   : doHierarchicalhhblitsSearch
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut  
sub doHierarchicalhhblitsSearch{
        	#### PARAMETER VARIABLES
                my ($self, $output_file,$leftDaughter,$rightDaughter,$leftDaughterType,$rightDaughterType,$singleSimilaritySearchResults) = (@_);
        	my $query         = $leftDaughter->fileInformation->sequenceSearchInputFile;
        	my $database      = $rightDaughter->fileInformation->sequenceSearchInputFile;
                print "\t$query - $database\n";
        	#my $analysis_step = $;
        	my $score_cutoff  = $self->configuration->similaritySearchCutoff;
        	my $usearch       = $self->configuration->usearch;
        	my $hhblits       = $self->configuration->hhblits;
        	#my $database_cs_file = $arg_ref->{database_cs_file};
        	my $species          = $leftDaughter->name;
        	my $species2         = $rightDaughter->name;
        #### OTHER VARIABLES
        	my ( $fh,          $search_output_file )       = tempfile(UNLINK => 1);
        	my ( $fh2,         $search_domtb_output_file ) = tempfile(UNLINK => 1);
        	my ( $fh_query,    $temp_query_file )          = tempfile(UNLINK => 1);
        	my ( $temp_hmm_fh, $temp_hmm_file )            = tempfile(UNLINK => 1);
        	$search_domtb_output_file = "test_domtb.table";
        	my %id_to_hmm_hash = ();
        #### TESTING ################################################################################
        	## CHECK DEFINEDNESS
        	print "$query, $database, $output_file\n";
        	croak q{One/several of the parameters for 'method' was/were not defined.\n}
        	  if any { !defined $_ } $query, $database, $output_file;
                  my %sequencesA;
                  my %query2resolvedID_mapping;
                  my %database2resolvedID_mapping;
                  my %ID2Alignment_mapping;
                  my %sequencesB;
                  my $seq;
                  my $seqID;

                  # Read A and B
                  readSequenceFileIntoHash({sequence_file => $query, sequence_href => \%sequencesA});
                  readSequenceFileIntoHash({sequence_file => $database, sequence_href => \%sequencesB});
                  print "read sequences in $database\n";
                  #print Dumper %sequencesB;
                  #exit;
                  
                  if($leftDaughterType eq 'innerNode'){
                          $leftDaughter->orthologGroups->get_group2expandedIDsHash( $leftDaughter->orthologGroups,
                                                                                                        $leftDaughter,
                                                                                               \%query2resolvedID_mapping);
                          $leftDaughter->get_alignmentsByID(\%ID2Alignment_mapping);
                                                                                       }
                  if($rightDaughterType eq 'innerNode'){
                       $rightDaughter->orthologGroups->get_group2expandedIDsHash( $rightDaughter->orthologGroups,
                                                                                                     $rightDaughter,
                                                                                                    \%database2resolvedID_mapping);
                        $rightDaughter->get_alignmentsByID(\%ID2Alignment_mapping);
                  }
          	  die "Could not read sequences in $query\n" if(!keys(%sequencesA));
          	  die "Could not read sequences in $database\n" if(!keys(%sequencesB));
                        my ($fh3, $queryHMMFile);
                        my ($fh4, $queryCSDBFile);
                        my ($fh5, $databaseHMMFile);
                        my ($fh6, $databaseCSDBFile);
                        my $type_of_comparison;
                        foreach my $query(keys(%$singleSimilaritySearchResults)){
                                my $results_string;
                                HIT:
                                foreach my $hit(split(" ",$singleSimilaritySearchResults->{$query})){
                                        print "\tCurrent comparison is: $query - $hit\n";
                                        if(exists $query2resolvedID_mapping{$query} && exists $query2resolvedID_mapping{$hit}){
                                                # profile - profile

                                                if(!profile_vs_profile_search({configuration => $self->configuration,
                                                        alignment_stringA => $ID2Alignment_mapping{$query},
                                                        alignment_stringB => $ID2Alignment_mapping{$hit},
                                                        query => $query,
                                                        results_stringref => \$results_string,
                                                        hit => $hit, 
                                                        og_name => $query})){
                                                            #    next HIT;
                                                        }
                                        }
                                                # Profile - Sequence      
                                                if( exists $query2resolvedID_mapping{$query} && not exists $query2resolvedID_mapping{$hit}){
                                                        print "profile vs sequence\n";
                                                        if(!profile_vs_sequence_search({configuration => $self->configuration,
                                                                alignment_stringA => $ID2Alignment_mapping{$query},
                                                                sequence_stringB => $sequencesB{$hit},
                                                                query => $query ,
                                                                results_stringref => \$results_string,
                                                                hit => $hit , 
                                                                og_name => $query})){
                                                                     #   next HIT;
                                                                }
                                                }     
                                                # Sequence - Profile     
                                                if( not exists $query2resolvedID_mapping{$query} && exists $query2resolvedID_mapping{$hit}){
                                                        print "sequence vs profile\n";
                                                                if(!sequence_vs_profile_search({configuration => $self->configuration,
                                                                        sequence_stringA => $sequencesB{$hit},
                                                                        alignment_stringB => $ID2Alignment_mapping{$query},
                                                                        query => $hit ,
                                                                        results_stringref => \$results_string,
                                                                        hit => $query , 
                                                                        og_name => $hit})){
                                                                               # next HIT;
                                                                        }
                                                }     
                                                if(not exists $query2resolvedID_mapping{$query} && not exists $query2resolvedID_mapping{$hit}){
                                                                        if(!sequence_vs_sequence_search({
                                                                                
                                                                                results_stringref => \$results_string
                                                                        })){
                                                                              #  next HIT;
                                                                        }
                                                }
                                                print "\t\tadding line to $output_file\n";
                                                attach_to_file({file_name => $output_file,
                                                        text => $results_string});            
                                        }
                                }  
        	return 1;
}
=item profile_vs_profile_search()

to be written...

 Title   : profile_vs_profile_search
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub profile_vs_profile_search{
        my ($arg_ref) = @_;
        my $configuration    = $arg_ref->{configuration};
        my $alignment_stringA    = $arg_ref->{alignment_stringA};
        my $alignment_stringB    = $arg_ref->{alignment_stringB};
        my $query = $arg_ref->{query};
        my $hit = $arg_ref->{hit};
        my $results_stringref    = $arg_ref->{results_stringref};
        my $og_name    = $arg_ref->{og_name};
        print "profile vs profile\n";
#        my ($fh, $queryHMMFile) = tempfile(SUFFIX => '.hhm');
#        my $database_dir = tempdir( CLEANUP => 1);
        my $database_dir = tempdir(CLEANUP => 1);
        my $template = "$hit";
        my $queryHMMFile = "$database_dir/$query.hhm";
        my $hitHMMFile = "$database_dir/$hit.hhm";
        my $hitDBFile = "$database_dir/$hit.db";
        my $hitCSDBFile = "$database_dir/$hit.cs_db";
        my $hhblitsResultsFile = "$database_dir/$hit.results";
        my $parser_output = "$database_dir/$hit.parsed";

        #Query
        if(!buildHhblitsAlignment({
                configuration => $configuration,
                og_name => $query,
                alignment_string => $alignment_stringA,
                hhm_file => $queryHMMFile,
                search_type => "query"})){
                        print "\t\tCould not build Hmm for $query. Skipping\n";
                        return 0;
                }
                print "\t\tQuery done\n";

                # HIT
                if(!buildHhblitsAlignment({
                        configuration => $configuration,
                        og_name => $hit,
                        alignment_string => $alignment_stringB,
                        hhm_file => $hitHMMFile,
                        search_type => "query"})){
                                print "\t\tCould not build Hmm. Exiting\n";
                                return 0;
                        }

                        if(!buildHhblitsDB4Group({
                                configuration => $configuration,
                                og_name => $hit,
                                hitHMMFile => $hitHMMFile,
                                database_dir => $database_dir,
                                hitDBFile => $hitDBFile,
                                hitCSDBFile => $hitCSDBFile
                                })){
                                        print "\t\tCould not HHblits DB for $hit. Skipping\n";
                                        return 0;
                        }
                                print "database done\n";
                                print "\t\t start hhblits\n";

                                my $hhblits_command = $configuration->hhblits." -i $queryHMMFile -o $hhblitsResultsFile -n 1 -context_data libs/context_data.lib -cs_lib libs/cs219.lib -db $hitCSDBFile -dbhhm $hitDBFile -norealign -nofilter";
                                print "$hhblits_command\n";
                                `$hhblits_command`;
                                my @hhblits_results = `cat $hhblitsResultsFile`;
                                print "\t read ".scalar(@hhblits_results)." lines\n";
                                my $parsed_string;
                                if(!parse_hhblits_output({hhblits_results => \@hhblits_results,
                                                        parsed_string => \$parsed_string,
                                                        query_name => $query})){
                                                return 0;
                                        }
                                write_to_file({file_name => $parser_output,
                                                text => $parsed_string});
                                my $blast2xml_call = $configuration->perl." ".$configuration->blast2xml." -i $parser_output | ".$configuration->perl." ".$configuration->blastParser." ".$configuration->similaritySearchCutoff;

                                $$results_stringref = `$blast2xml_call`;
                                print "\tadding line ".$$results_stringref."\n";
                                
                                print "unlinking: $database_dir/*\n";
                           #     unlink glob "$database_dir/*";
return 1;
}
=item profile_vs_sequence_search()

to be written...

 Title   : profile_vs_sequence_search
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub profile_vs_sequence_search{
        my ($arg_ref) = @_;
        my $configuration    = $arg_ref->{configuration};
        my $sequence_stringB    = $arg_ref->{sequence_stringB};
        my $alignment_stringA    = $arg_ref->{alignment_stringA};
        my $query = $arg_ref->{query};
        my $hit = $arg_ref->{hit};
        my $og_name    = $arg_ref->{og_name};
        my $results_stringref    = $arg_ref->{results_stringref};
        my ($fh, $hmm_file) = tempfile( SUFFIX => '.hmm', UNLINK => 1);
        my ($fh2, $query_file) = tempfile( SUFFIX => '.fa', UNLINK => 1);
        my ($fh3, $hmmer_domtblout_file) = tempfile( SUFFIX => '.tblout', UNLINK => 1);
        my ($fh4, $parser_output) = tempfile( SUFFIX => '.parsed', UNLINK => 1);
        print "profile ($query) vs sequence ($hit)\n";

        # get alignment for A
        if($alignment_stringA eq ''){
                print "\tCould not find alignment for query $query\n";
                exit;
                return 0;
        }
        if($sequence_stringB eq ''){
                print "\tCould not find sequence for hit $hit\n";
                exit;
                return 0;
        }
        # build hmm for A
#        print "\talignment string : $alignment_stringA\n";
        if(!(buildHmmerAlignment({
                configuration => $configuration,
                alignment_string => $alignment_stringA,
                hmm_file => $hmm_file,
                og_name => $query
                }))){
                        return 0;
                }
        # save sequence of B        
#        print "saving sequence $sequence_stringB\n";
                open(my $TMP_ALN_FH, ">$query_file  ") || die "saving tmp alignment failed: $!\n";
                print {$TMP_ALN_FH} "$sequence_stringB\n";
                close($TMP_ALN_FH);
                # start hmmer search
                my $hmmsearch = $configuration->hmmsearch." --domtblout $hmmer_domtblout_file $hmm_file $query_file > /dev/null 2>&1";
                
                print "$hmmsearch\n";
                `$hmmsearch`;
                my $parsed_string;
                if(!parse_hmmer_domtblout({hmmer_domtblout_file => $hmmer_domtblout_file,
                                        parsed_string => \$parsed_string})){
                                return 0;
                        }
                write_to_file({file_name => $parser_output,
                                text => $parsed_string});
                my $blast2xml_call = $configuration->perl." ".$configuration->blast2xml." -i $parser_output | ".$configuration->perl." ".$configuration->blastParser." ".$configuration->similaritySearchCutoff;

                $$results_stringref = `$blast2xml_call`;
                print "\tadding line ".$$results_stringref."\n";
                return 1;
}
=item sequence_vs_profile_search()

to be written...

 Title   : sequence_vs_profile_search
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub sequence_vs_profile_search{
        my ($arg_ref) = @_;
        my $configuration    = $arg_ref->{configuration};
        my $sequence_stringA    = $arg_ref->{sequence_stringA};
        my $alignment_stringB    = $arg_ref->{alignment_stringB};
        my $query = $arg_ref->{query};
        my $hit = $arg_ref->{hit};
        my $results_stringref    = $arg_ref->{results_stringref};
        my $og_name    = $arg_ref->{og_name};
        my ($fh, $hmm_file) = tempfile(SUFFIX => '.hmm');
        my ($fh2, $query_file) = tempfile( SUFFIX => '.fa');
        my ($fh3, $hmmer_domtblout_file) = tempfile( SUFFIX => '.tblout');
        my ($fh4, $parser_output) = tempfile( SUFFIX => '.parsed', UNLINK => 1);
        print "sequence ($query) vs profile ($hit)\n";

        # get alignment for A
        if($alignment_stringB eq ''){
                print "\tCould not find alignment for hit $hit\n";
                return 0;
        }
        if($sequence_stringA eq ''){
                print "\tCould not find sequence for query $query\n";
                return 0;
        }
        # build hmm for A
        if(!(buildHmmerAlignment({
                configuration => $configuration,
                alignment_string => $alignment_stringB,
                hmm_file => $hmm_file,
                og_name => $hit
                }))){
                        return 0;
                }
        # save sequence of B        
                open(my $TMP_ALN_FH, ">$query_file  ") || die "saving tmp alignment failed: $!\n";
                print {$TMP_ALN_FH} "$sequence_stringA\n";
                close($TMP_ALN_FH);
                # start hmmer search
                # hmmsearch || hmmscan??
                my $hmmscan = $configuration->hmmscan." --domtblout $hmmer_domtblout_file $hmm_file $query_file > /dev/null 2>&1";
                #print "$hmmscan\n";
                #exit;
                 my $parsed_string;
                        if(!parse_hmmer_domtblout({hmmer_domtblout_file => $hmmer_domtblout_file,
                                                parsed_string => \$parsed_string})){
                                        return 0;
                                }
                        write_to_file({file_name => $parser_output,
                                        text => $parsed_string});
                my $blast2xml_call = $configuration->perl." ".$configuration->blast2xml." -i $parser_output | ".$configuration->perl." ".$configuration->blastParser." ".$configuration->similaritySearchCutoff;

                $$results_stringref = `$blast2xml_call`;
                print "\tadding line ".$$results_stringref."\n";
                return 1;

}
=item sequence_vs_sequence_search()

to be written...

 Title   : sequence_vs_sequence_search
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub sequence_vs_sequence_search{
        print "sequence vs sequence\n";
}




1;