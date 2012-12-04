package SequenceSearch;
=head1 SequenceSearch

hieranoid::Comparison::OrthologySearch::SimilaritySearch::SequenceSearch - Container of taxon objects

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
use Carp;
use File::Copy;
use English;
use Log::Log4perl qw(:easy);
use File::Temp qw/ tempfile tempdir /;
use File::Basename;
use Hieranoid::Tree::InnerNode;
use Data::Dumper;
use Benchmark;
#use Daughter;
use Hieranoid::Comparison::ResultsSimilaritySearch;





sub do_ublastxml_2pass {
        my ($self, $output_file,$leftDaughter,$rightDaughter,$leftDaughterType,$rightDaughterType,$singleSimilaritySearchResults,$allowedQueryDatabaseEntities) = (@_);
	my $query         = $leftDaughter->fileInformation->sequenceSearchInputFile;
	my $database      = $rightDaughter->fileInformation->sequenceSearchInputFile;
	my ($queryName,$databaseName) = (basename($query),basename($database));
	my $ublast_output = $output_file."_ublast.out";
        INFO("\t$queryName - $databaseName\n");
	#my $analysis_step = $;
	my $score_cutoff  = $self->configuration->similaritySearchCutoff;
	my $usearch       = $self->configuration->usearch;
	my $hhblits       = $self->configuration->hhblits;
	my $blastParser   = $self->configuration->blastParser;
	#my $database_cs_file = $arg_ref->{database_cs_file};
	my $species          = $leftDaughter->name;
	my $species2         = $rightDaughter->name;
   #   print Dumper @Fld;
	# $Fld [0] is query
	# $Fld [1] is database
	# $Fld [2] is query size
	# $Fld [3] is database size
	# $Fld [4] is output name
	# assume the script has already formatted the database
	# we will now do 2-pass approach
	# load sequences
      # MEASURE THE TIME
	my ($prev_time,,, ) = times;
	my $user_time = ();
      # Sequence masking
	my (undef, $tmpi) = tempfile( SUFFIX => '.2query.db',UNLINK => 1);
	my (undef, $tmpd) = tempfile( SUFFIX => '.2hit.db',UNLINK => 1);
	my (undef,$ublast_output_tmp) = tempfile(SUFFIX => '.2ublast.output',UNLINK => 1);
	my ($noQuerySequences,$noDBSequences);

	my %id2Seq_mapping;
        # Read A and B
	readSequenceFileIntoHash({sequence_file => $query, sequence_href => \%id2Seq_mapping});
	readSequenceFileIntoHash({sequence_file => $database, sequence_href => \%id2Seq_mapping});
       
	unlink "$tmpi";
	unlink "$tmpd";

	my $number_of_queries_to_reblast = keys %{$singleSimilaritySearchResults};
	my $number_of_hits_to_reblast = keys %{$singleSimilaritySearchResults};
	DEBUG("Reblasting for $number_of_queries_to_reblast queries");
	DEBUG("Reblasting for $number_of_hits_to_reblast hits");
	#print "Reblasting for $number_of_hits_to_reblast hits\n";
	my $print_to_tmp_BLAST_db_string = ();
	my $print_to_tmp_BLAST_query_string = ();
	my ($noQuery, $noDB);
	my %avoid_duplicates_in_DB_hash = ();
	my %avoid_duplicates_in_query_hash = ();
	
	QUERY:
	foreach my $aQuery ( keys %{$singleSimilaritySearchResults} ) {
                
		# Create single-query file
		if(!exists $id2Seq_mapping{$aQuery}){
		        #print "\tCould not find query $aQuery in Query sequences ($query)\n";
		        WARN("\tCould not find query $aQuery in Query sequences ($query)\n");
                        #exit;
		        next QUERY;
		}
		if(!exists $avoid_duplicates_in_query_hash{$aQuery}){
		        $print_to_tmp_BLAST_query_string .= $id2Seq_mapping{$aQuery} . "\n";
		        $noQuery++;
	        }
	        $avoid_duplicates_in_query_hash{$aQuery} = 1;
		# Create mini-database of hit sequences
		HIT:
		foreach my $aHit ( split( /\s/, $singleSimilaritySearchResults->{$aQuery} ) ) {
			next if exists $avoid_duplicates_in_DB_hash{$aHit};
			if(!exists $id2Seq_mapping{$aHit}){
				#print "\tCould not find hit $aHit in DB ($database)\n";
				WARN("\tCould not find hit $aHit in DB ($database)\n");
                                #exit;
				next HIT;
			}
			$print_to_tmp_BLAST_db_string .= $id2Seq_mapping{$aHit} . "\n";
		        $avoid_duplicates_in_DB_hash{$aHit} = 1;
		        $noDB++;
		}
	}
	if($noQuery == 0 || $noDB == 0){
	        DEBUG("Could not build query ($noQuery) or hits database ($noDB)");
	        exit;
	}
	DEBUG("Trying to write Query database for $noQuery query sequences and length ".length($print_to_tmp_BLAST_db_string)."\n");
        # Writing Query
        open( my $temp_db_FH, ">$tmpd" ) or die "Couldn't open '$tmpd': \n";
        print {$temp_db_FH} "$print_to_tmp_BLAST_db_string" or die "Couldn't write Query db '$tmpd': $!\n";
        close($temp_db_FH) or die "Couldn't close '$tmpd': \n";

	DEBUG("Trying to write hit database for $noDB db sequences and length ".length($print_to_tmp_BLAST_query_string)."\n");
        # Writing DB
        open( my $temp_query_FH, ">$tmpi" ) or die "Couldn't open '$tmpi': \n";
        print {$temp_query_FH} "$print_to_tmp_BLAST_query_string" or die "Couldn't write DB '$tmpi': $!\n";
        close($temp_query_FH) or die "Couldn't close '$tmpi': \n";
        
        
	DEBUG("Starting second UBLAST pass for $query - $database on (".basename($tmpi)." vs. ".basename($tmpd).")");
	#print "\t\t\tStarting second UBLAST...\n";
	# Taking time
        # start timer
        my $start_time = new Benchmark;
        if($self->configuration->similaritySearchTool eq 'usearch' || $self->configuration->similaritySearchTool eq 'ublast'){
                my $blast_call = "$usearch --quiet --query $tmpi  --db $tmpd  --userfields query+target+bits+ql+tl+qlo+qhi+tlo+thi+qs+ts --userout $ublast_output_tmp  ".$self->configuration->ublastParameters;
	        DEBUG("2. blast_call: ".$blast_call."\n");
                `$blast_call`;
                #print "finished second run. copying to $ublast_output\n";
                copy($ublast_output_tmp,$ublast_output) || die "2.Run: Could not copy Ublast output $ublast_output_tmp to $ublast_output\n";
                ( $user_time,,, ) = times;
                $prev_time = $user_time;
                # The second Blast does not have a bitscore cutoff thingy  
                my $blast2xml_call = $self->configuration->perl." blast2xml.pl -i $ublast_output_tmp | ".$self->configuration->perl." ./$blastParser 0 > $output_file ";
                DEBUG("\tConverting to xml ($blast2xml_call)\n");
                `$blast2xml_call`;
                unlink($ublast_output) or die "Could not remove ublast_output: $ublast_output\n";
        }
        if($self->configuration->similaritySearchTool eq 'blast'){
                system ($self->configuration->formatdb." -i $tmpd");
                my $noDBSequences = `grep ">" $tmpd -c`;
                chomp($noDBSequences);
	        system ($self->configuration->blast." -C0 -FF -i $tmpi -d $tmpd -p blastp -v $noDBSequences -b $noDBSequences -M BLOSUM62 -z 5000000 -m7 | ".$self->configuration->perl." ./".$self->configuration->blastParser." 0 > $output_file ");
                #unlink($tmpd);
        }
        
        # Stop time
        my $end_time = new Benchmark;
        my $time_diff = timediff($end_time,$start_time);
        my $debugString = "2.Pass\t".$self->configuration->similaritySearchTool."\t$queryName ($noQuery) - $databaseName ($noDB):  ".timestr($time_diff, 'all');
        attach_to_file({file_name => $self->configuration->timeFile, text => $debugString."\n"});        
   
        unlink($tmpi);
        unlink($tmpd);
        unlink($ublast_output_tmp);
        
        
        return 1;
}

      


1;