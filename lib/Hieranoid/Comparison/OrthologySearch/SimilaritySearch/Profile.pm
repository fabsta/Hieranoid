package Profile;
=head1 Profile

hieranoid::Comparison::OrthologySearch::SimilaritySearch::Profile - Container of taxon objects

=head1 SYNOPSIS

 use lib::hieranoid::Comparison;
 my $currentComparison = Comparison->new(nodeObject =>$innerNode,
						configuration => $hieranoidConfiguration);

=head1 DESCRIPTION

Compares orthologous groups between daughter nodes in a given phylogenetic tree

=head1 METHODS

=over

=cut
use Moose::Role;
use Carp;
use English;
use File::Basename;
use Data::Dumper;
use File::Copy;
use Log::Log4perl qw(:easy);
use File::Temp qw/ tempfile tempdir /;
use List::MoreUtils qw( any );
use Hieranoid::Tree::InnerNode;


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
        DEBUG("MODE: profile vs profile\n");
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
                        ERROR("\t\tCould not build Hmm for $query. Skipping\n");
                        return 0;
                }
                DEBUG("\t\tQuery done\n");

                # HIT
                if(!buildHhblitsAlignment({
                        configuration => $configuration,
                        og_name => $hit,
                        alignment_string => $alignment_stringB,
                        hhm_file => $hitHMMFile,
                        search_type => "query"})){
                                ERROR("\t\tCould not build Hmm. Exiting\n");
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
                                        ERROR("\t\tCould not HHblits DB for $hit. Skipping\n");
                                        return 0;
                        }
                                DEBUG("\t\tdatabase done\n");
                                DEBUG("\t\t start hhblits\n");

                                my $hhblits_command = $configuration->hhblits." -i $queryHMMFile -o $hhblitsResultsFile -n 1 -context_data libs/context_data.lib -cs_lib libs/cs219.lib -db $hitCSDBFile -dbhhm $hitDBFile -norealign -nofilter";
                                DEBUG("$hhblits_command\n");
                                `$hhblits_command`;
                                my @hhblits_results = `cat $hhblitsResultsFile`;
                                
                                #print "\t read ".scalar(@hhblits_results)." lines\n";
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
                                DEBUG("\tadding line ".$$results_stringref."\n");
                                #DEBUG("unlinking: $database_dir/*\n");
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
        my ($fh, $hmm_file) = tempfile( SUFFIX => '.hmm');
        my ($fh2, $query_file) = tempfile( SUFFIX => '.fa');
        my ($fh3, $hmmer_domtblout_file) = tempfile( SUFFIX => '.tblout');
        my ($fh4, $parser_output) = tempfile( SUFFIX => '.parsed');
        INFO("\t\tprofile ($query) vs sequence ($hit)\n");

        # get alignment for A
        if($alignment_stringA eq ''){
                ERROR("\tCould not find alignment for query $query\n");
                exit;
                return 0;
        }
        if($sequence_stringB eq ''){
                ERROR("\tCould not find sequence for hit $hit\n");
                exit;
                return 0;
        }
        # build hmm for A
#        ERROR "\talignment string : $alignment_stringA\n";
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
                
                DEBUG("\thmmsearch: $hmmsearch\n");
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
                DEBUG("\tadding line ".$$results_stringref."\n");
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
        my ($fh4, $parser_output) = tempfile( SUFFIX => '.parsed');
        INFO("sequence ($query) vs profile ($hit)\n");

        # get alignment for A
        if($alignment_stringB eq ''){
                ERROR("\tCould not find alignment for hit $hit\n");
                return 0;
        }
        if($sequence_stringA eq ''){
                ERROR("\tCould not find sequence for query $query\n");
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
                DEBUG("\tadding line ".$$results_stringref."\n");
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
=item buildHmmerAlignment()

to be written...

 Title   : buildHmmerAlignment
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub buildHmmerAlignment{
        my ($arg_ref) = @_;
        my $alignment_string    = $arg_ref->{alignment_string};
        my $configuration    = $arg_ref->{configuration};
        my $hmm_file    = $arg_ref->{hmm_file};
        my $og_name    = $arg_ref->{og_name}; 
        my ($fh, $tmp_alignment_file) = tempfile(SUFFIX => '.fas', UNLINK => 1);
               
        open(my $TMP_ALN_FH, ">$tmp_alignment_file  ") || die "saving tmp alignment failed: $!\n";
        print {$TMP_ALN_FH} "$alignment_string\n";
        close($TMP_ALN_FH);
        
        
        my $hmmbuild_call =  $configuration->hmmbuild." --informat afa -n \"$og_name\" --amino  $hmm_file $tmp_alignment_file > /dev/null 2>&1";
        system("$hmmbuild_call");
        DEBUG("$hmmbuild_call\n");
        if(! -e $hmm_file || ! -s $hmm_file){
                ERROR("\tCould not build hmm for $og_name ($hmmbuild_call)\n");
                exit;
        }
        return 1;
}
=item buildHhblitsAlignment()

to be written...

 Title   : buildHhblitsAlignment
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub buildHhblitsAlignment{      
        my ($arg_ref) = @_;
        my $configuration    = $arg_ref->{configuration};
        my $alignment_string    = $arg_ref->{alignment_string};
        my $hhm_file    = $arg_ref->{hhm_file};
        my $tmp_alignment_file = $arg_ref->{tmp_alignment_file};
        my $tmp_hmm_output_file = $arg_ref->{tmp_hmm_output_file};
        my $og_name    = $arg_ref->{og_name};
        my $attach    = $arg_ref->{attach};
        #$hmmbuild_call .= "> /dev/null 2>&1";
        open(my $TMP_ALN_FH, ">$tmp_alignment_file  ") || die "saving tmp alignment failed: $!\n";
        print {$TMP_ALN_FH} "$alignment_string\n";
        close($TMP_ALN_FH);
        my $hmmbuild_call;
        if($attach){
                 $hmmbuild_call = $configuration->hhmake." -i $tmp_alignment_file -name \"$og_name\" -o $tmp_hmm_output_file -M first ";
                 `$hmmbuild_call`;
                 DEBUG("HMMbuild call: $hmmbuild_call\n");
                 if(!-e $tmp_hmm_output_file || ! -s $tmp_hmm_output_file){
                         ERROR("Could not build profile $tmp_hmm_output_file ($hmmbuild_call). Skipping\n");
                         print "Could not build profile $tmp_hmm_output_file ($hmmbuild_call). Skipping\n";
                         exit;
                         return 0;
                         
                 }
                  `cat $tmp_hmm_output_file >> $hhm_file`;
         }
         else{
                 $hmmbuild_call = $configuration->hhmake." -i $tmp_alignment_file -name \"$og_name\" -o $hhm_file -M first ";
                 `$hmmbuild_call`;
                 #DEBUG("HMMbuild call: $hmmbuild_call\n");
                 if(!-e $hhm_file || ! -s $hhm_file){
                         ERROR("Could not build profile $hhm_file ($hmmbuild_call). Skipping\n");
                         print "Could not build profile $hhm_file ($hmmbuild_call). Skipping\n";
                         exit;
                         return 0;
                 }
         }
         
         return 1;
}
=item buildHhblitsDB4Group()

to be written...

 Title   : buildHhblitsDB4Group
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub buildHhblitsDB4Group{
        my ($arg_ref) = @_;
        my $configuration    = $arg_ref->{configuration};
        my $hitHMMFile    = $arg_ref->{hitHMMFile};
        my $hitDBFile    = $arg_ref->{hitDBFile};
        my $database_dir    = $arg_ref->{database_dir};
        my $hitCSDBFile    = $arg_ref->{hitCSDBFile};
        my $og_name    = $arg_ref->{og_name};

        #print "\tLooking at tmp-dir $database_dir\n";
        ## CREATE HHBLITS DATABASE
        my $hhblits_db_command = $configuration->perl." ".$configuration->create_hhblits_db." -hhmdir $database_dir -ohhm $hitDBFile";
        #print $hhblits_db_command."\n";
        system($hhblits_db_command);	
        
        if(!-e $hitDBFile || ! -s $hitDBFile){
                print "Could not create hhblits database (missing hhmfile: $hitDBFile).\n";
                exit;
        }

        # CS-DB
        my $hhblits_cs_db_command = $configuration->perl." ".$configuration->create_hhblits_csdb."  -i $database_dir -o $hitCSDBFile -ext hhm";
        #print "\t\t\tAttaching to CS_DB\n";
        #print "$hhblits_cs_db_command\n";
        system("$hhblits_cs_db_command");
        if(! -e $hitCSDBFile && ! -s $hitCSDBFile){
                ERROR("\tCould not build cs-db for hhblits (missing file: $hitCSDBFile). \n\t\tcommand: $hhblits_cs_db_command\n");
                return 0;
        }
        return 1;
}
=item perform_hhsearch()

to be written...

 Title   : perform_hhsearch
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub perform_hhsearch{
        my ($self, $queryHMMFile,$databaseHMMFile,$results_stringref,$ID2SequenceLengthHref,$debug_hhsearch_output) = (@_);
        
        my $hhsearch_command = $self->configuration->hhsearch." -i $queryHMMFile -d $databaseHMMFile 2> /dev/null";
        # global mode
        #my $hhsearch_command = $self->configuration->hhsearch." -i $queryHMMFile -d $databaseHMMFile -global  2> /dev/null";
        my @hhsearch_results = `$hhsearch_command`;
        DEBUG($hhsearch_command."\n");
        #attach_to_file({file_name => $debug_hhsearch_output, text => join(",",@hhsearch_results)."\n"}); #### DEBUG
        if(!scalar(@hhsearch_results)){
                print "\tdid not find hits\n";
                return 1;
        }
        my $parsed_string;

        if(!parse_hhblits_output({hhblits_results => \@hhsearch_results,
                                parsed_string => \$parsed_string,
                                ID2SequenceLengthHref => $ID2SequenceLengthHref})){
                        print "\tThere was a problem converting the hits from Ublast format\n";
                        return 0;
                }
                
        $$results_stringref .= $parsed_string."\n";
        # clean up
        unlink($queryHMMFile) if -e $queryHMMFile;
        unlink($queryHMMFile.".hhr") if -e $queryHMMFile.".hhr";
        
        return 1;
}
=item parse_hmmer_domtblout()

to be written...

 Title   : parse_hmmer_domtblout
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub parse_hmmer_domtblout {
	#### PARAMETER VARIABLES
	my ($arg_ref)        = @_;
#	print Dumper $arg_ref;
	my $hmmer_domtblout_file = $arg_ref->{hmmer_domtblout_file};
	my $parsed_stringref = $arg_ref->{parsed_string};
	my $statistics_href;

	#print "parsing results\n";
	if(! -e $hmmer_domtblout_file || ! -s $hmmer_domtblout_file){
	        WARN("\tNo hits. Skipping\n");
	        return 1;
	}
	#print "\t reading from $hmmer_domtblout_file\n";
	open my $HMM_TABLE_FH, "<",
	  $hmmer_domtblout_file || die "\tCould not read results from hmmsearch ($hmmer_domtblout_file)\n";
	while ( my $line = <$HMM_TABLE_FH> ) {
		next if $line =~ /^#/;

				#print $line;
		my @splitted_line = split( /\s+/, $line );
		my (
			 $target_name, $accession, $tlen,      $query_name,
			 $accession_q, $qlen,      $e_value,   $score,
			 $bias,        $number,    $number_of, $e_value_d,
			 $i_value_d,   $score_d,   $bias_d,    $hmm_from,
			 $hmm_to,      $ali_from,  $ali_to,    $env_from,
			 $env_to,      $acc
		  )
		  = (
			  $splitted_line[0],  $splitted_line[1],  $splitted_line[2],
			  $splitted_line[3],  $splitted_line[4],  $splitted_line[5],
			  $splitted_line[6],  $splitted_line[7],  $splitted_line[8],
			  $splitted_line[9],  $splitted_line[10], $splitted_line[11],
			  $splitted_line[12], $splitted_line[13], $splitted_line[14],
			  $splitted_line[15], $splitted_line[16], $splitted_line[17],
			  $splitted_line[18], $splitted_line[19], $splitted_line[20],
			  $splitted_line[21]
		  );
                #  query	target	bits	ql	tl	qlo	qhi	tlo	thi	qs	ts
                  $$parsed_stringref .= "$query_name\t$target_name\t$score\t$qlen\t$tlen\t$ali_from\t$ali_to\t$hmm_from\t$hmm_to\t".(($ali_to - $ali_from)+1)."\t".(($hmm_to - $hmm_from)+1);
		      #  print "$$parsed_stringref\n";
		next if ( $query_name eq $target_name );
	}
	close($HMM_TABLE_FH) or die "Could not close file $hmmer_domtblout_file\n";
        if($$parsed_stringref eq ''){
                ERROR("Could not parse from $hmmer_domtblout_file or no hits\n");
                return 0;
        }
	return 1;
}
=item parse_hhblits_output()

to be written...

 Title   : parse_hhblits_output
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub parse_hhblits_output {
        #### PARAMETER VARIABLES
        my ($arg_ref)            = @_;
        my $hhblits_results_aref = $arg_ref->{hhblits_results};
        my $parsed_stringref      = $arg_ref->{parsed_string};
        my $ID2SequenceLengthHref      = $arg_ref->{ID2SequenceLengthHref};

        ## TESTING
        ## CHECK DEFINEDNESS
        croak
        q{One/several of the parameters for 'parse_hhsearch_output' was/were not defined.\n}
        if any { !defined $_ } $hhblits_results_aref;
        my $hits = 0;
        my $start_parsing=0;
        my $query_name;
        my $query_length = 0;
        my %hits_found;
        
        foreach my $line ( @{$hhblits_results_aref} ) {
                next if $line =~ /^#/;
                next if $line  =~ /hits/;
                next if $line  =~ /Done/;
                next if $line =~ /Query file/; 
        #print "line: ".$line."\n";
                if($line =~ /^Query\s+(.*)/){
                        $query_name = $1;
                        if(!exists $ID2SequenceLengthHref->{$query_name}){
                                WARN("Could not find length for query: $query_name. Taking cols\n");
                                #exit;
                        }
                        $query_length  = $ID2SequenceLengthHref->{$query_name};
                        next;
                }
                if($line =~ /^\s+No Hit/){
                        $start_parsing=1;
                        next;
                }
               # Query         Homininae597
                next if !$start_parsing;
                last if $line =~ /^No 1/;
                if ($line =~ /^\s+\d+/){
                        my @splitted_line = split( /\s+/, $line );

                        #	print Dumper @splitted_line;
                        my ($empty_space, $hit_number,  $hit_name, $prob,
                                $e_value,     $p_value,     $score,       $ss,
                                $cols,        $query_coord, $target_coord)
                                = ($splitted_line[0], $splitted_line[1], $splitted_line[2],
                                        $splitted_line[3], $splitted_line[4], $splitted_line[5],
                                        $splitted_line[6], $splitted_line[7], $splitted_line[8],
                                        $splitted_line[9], $splitted_line[10]
                                        );
                                        #if(exists $hits_found{$hit_name}){
                                        #        print "\tneed to check this case\n";
                                                #$line
                                                #exit;
                                        #}
                                        #print "$empty_space, hitnumber:$hit_number,targetname:$hit_name,prob:$prob,evalue:$e_value,pvalue:$p_value,score:$score,ss:$ss,cols:$cols,query:$query_coord,target:$target_coord\n";
                                        #if(!$prob =~ /\d+\.\d+/){print $line."\n";print "problem with cols: $prob\n"; exit;}
                                        #if(!$cols =~ /\d+/){print $line."\n";print "problem with cols: $cols\n"; exit;}
                                        #if(!$e_value =~ /\d+/){print $line."\n";print "problem with cols: $e_value\n"; exit;}
                                        ## parse coordinates
                                        #Query                    Hit id                 BitS     Q_len H_len   Q_lon H_lon  Q_mat H_mat Positions
                                         #         ADH3_YEAST      ADHP_ECOLI      207.2   375     336     316     311     316     311     q:57-372 h:23-333
                                        
                                        $query_coord =~ /(\d+)-(\d+)/;
                                        my ( $query_from, $query_to ) = ( $1, $2 );
                                        croak qq{Could not parse parameter query_from ($query_coord) in line $line\n} if ( !defined $query_from || !defined $query_to );

                                        $target_coord =~ /(\d+)-(\d+)/;
                                        my ( $target_from, $target_to ) = ( $1, $2 );
                                        croak qq{Could not parse parameter target_from in line $line\n} if ( !defined $target_from || !defined $target_to );
                                        if(!$query_length){
                                                $query_length = $cols;
                                        }
                                        my $hit_length;
                                        if(!exists $ID2SequenceLengthHref->{$hit_name}){
                                                WARN("Could not find length for hit $hit_name. Taking cols\n");
                                                $hit_length = $cols;
                                                #exit;
                                        }
                                        else{
                                                $hit_length  = $ID2SequenceLengthHref->{$hit_name};
                                        }
                                    $$parsed_stringref .= "$query_name\t$hit_name\t$score\t$query_length\t$hit_length\t$query_from\t$query_to\t$target_from\t$target_to\t".(($query_to - $query_from)+1)."\t".(($target_to - $target_from)+1)."\n";
                                        #print "$$parsed_stringref\n";
                                }
                        }
                        if($$parsed_stringref eq ''){
                                ERROR("Could not parse from hhblits results or no hits\n");
                                return 0;
                        }
                        return 1;
}
=item createHhblitsProfiles()

to be written...

 Title   : createHhblitsProfiles
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub createHhblitsProfiles{
        my ($arg_ref) = @_;
        my $self    = $arg_ref->{self};
        my $configuration    = $arg_ref->{configuration};
        my $sequenceFile    = $arg_ref->{sequenceFile};
        my $hhm_file    = $arg_ref->{hhm_file};
        my $tmp_alignment_file = $arg_ref->{tmp_alignment_file};
        my $tmp_hmm_output_file = $arg_ref->{tmp_hmm_output_file};
        my $hitDBFile = $arg_ref->{hitDBFile};
        my $hitCSDBFile = $arg_ref->{hitCSDBFile};
        my $id2SeqAli_mapping_href = $arg_ref->{id2SeqAli_mapping};
        my $name    = $arg_ref->{name};
        my $hhblitsDatabase = $arg_ref->{hhblitsDatabase};
        # Iterate over all IDs
        my %id2SeqAli_mapping;
         #my ($fh, $tmp_alignment_file) = tempfile();
         #my ($fh2, $tmp_hmm_output_file) = tempfile();
         my $database_dir = tempdir(  );
         $tmp_hmm_output_file = $database_dir."/$name.hhm";
        my @allHeaders  = `grep ">" $sequenceFile`;
        my @allIDs;
        my $counter = 0;
        foreach my $aLine (@allHeaders){
          		chomp($aLine);
          		my @words = split( /\s+/, $aLine );
          		$words[0] =~ s/>//;
 #         		print "saving ".$words[0]."\n";
#          		exit;
                        push(@allIDs,$words[0]);
        }       
        foreach my $currentHeader(@allIDs){
                DEBUG("\tMaking HHM for $currentHeader ($hhm_file)\n");
                buildHhblitsAlignment({configuration => $self->configuration,
                        og_name => $currentHeader,
                        alignment_string => $id2SeqAli_mapping_href->{$currentHeader},
                        tmp_alignment_file => $tmp_alignment_file,
                        tmp_hmm_output_file => $tmp_hmm_output_file,
                        hhm_file => $hhm_file,
                        attach => 1});
                        last if $counter++ > 200;

                        if($hhblitsDatabase){
                                if(!buildHhblitsDB4Group({
                                        configuration => $configuration,
                                        #og_name => $currentHeader,
                                        #hitHMMFile => $tmp_hmm_output_file,
                                        database_dir => $database_dir,
                                        hitDBFile => $hitDBFile,
                                        hitCSDBFile => $hitCSDBFile
                                        })){
                                                ERROR("\t\tCould not HHblits DB for . Skipping\n");
                                                return 0;
                                        }
                                }
                           #     exit;
        }       
         
        #unlink glob "$database_dir/*";
}
=item write_to_file()

to be written...

 Title   : write_to_file
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub write_to_file() {
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
=item getSequenceLength()

to be written...

 Title   : getSequenceLength
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -

=cut
sub getSequenceLength{
        #### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $sequence = $arg_ref->{sequence};
	my @splitted = split(/\n/,$sequence);
#	print "sequence: $sequence\n\tlength ".length($splitted[1])."\n";
        if(length($splitted[1]) == 0){
                ERROR("\tCould not determine sequence length from ".$splitted[1]."\n");
                exit;
        }
        return length($splitted[1]);
}
     
1;