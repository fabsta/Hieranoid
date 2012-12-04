package PairwiseConsensusSummarizer;
=head1 PairwiseConsensusSummarizer

hieranoid::Comparison::Summarizer::PairwiseConsensusSummarizer - Container of taxon objects

=head1 SYNOPSIS

 Abstract class

=head1 DESCRIPTION

Abstract class to handle how orthologous group are summarized

=head1 METHODS

=over

=cut
use Moose;
use Carp;
use File::Copy;
use English;
use Log::Log4perl qw(:easy);
use Data::Dumper;
use File::Temp qw/ tempfile tempdir /;
use Hieranoid::Tree::InnerNode;
use Hieranoid::FileInformation;
with 'Summarizer';


=item alignAndConsense()

Summarizes results of orthology search of daughter nodes.
two daughter nodes will be merged into one pseudospecies
so that it can be used in an Inparanoid comparison


=head3 Example: OG1: Euarchontoglires2:Homininae2,Murinae2
                1. Get all sequence from Hominiae2
                   Get all sequence from Murinae2
                2. Build an alignment of Hominiae2,Murinae2
                3. Build consensus sequence and save it as Euarchontoglires2
                
                
 Title   : alignAndConsense
 Usage   : $self->summarizer->summarizeInformation($self->nodeObject);
 Function: Summarizes two daughter nodes 
 Returns : 1 on success
 Args: -



=cut

 # overwrite method because....
 # here we take groups in format
        # Catarrhini3:CCNK_MACMU,Homininae2
        
 
sub alignAndConsense{
      #  print Dumper @_;
        my ($self,$innerNode) = (@_);
        my ($fh2, $alignment_tmp_file) = tempfile(UNLINK => 1);
        print "\t\talignment and consensus of node ".$innerNode->name." (OrthologyPredictions in ".$innerNode->orthologGroups->groupsFile.")\n";
        my $predictionsFile = $innerNode->orthologGroups->groupsFile;
    # read predictions
        my %groupPredictionHash;
        $innerNode->orthologGroups->get_group2IDsHash($innerNode,\%groupPredictionHash);
   # read id2sequences
        my %ID2sequencesHash;
        $innerNode->get_cladeSequencesByID(\%ID2sequencesHash);
        $innerNode->get_sequencesByID(\%ID2sequencesHash);
        # align&consense
        print "\t\tcompute alignments and consensus\n";
   # numeric sort does not work anymore because clusters might look like "Homininae85"
        my $current_number_of_og = 0;
        my $number_of_ogs_total = keys(%groupPredictionHash);
        my %assigned_sequence_ids_hash;
        my $print2FileString;
        # foreach my $og(sort {$a<=>$b} keys(%{$og_hashref})){
                OG:
                foreach my $og(keys(%groupPredictionHash)){
                        my $current_sequences = q();
                        my $number_of_sequences_for_group = 0;
                        my @ids_for_group = ();
                        ### delete old files
                        ### FILE CONTAINS:
                        # speciesname04 --> member1, member2,...,memberx
                        print("\tCurrent ortholog group: $og ( $current_number_of_og / $number_of_ogs_total)");
                                foreach my $current_id(split(/,/,$groupPredictionHash{$og})){
                                        #print "\tsearching for $current_id\n";
                                        my $sequence = q();
                                        # CHECK IF THIS ID IS ACTUALLY A PROFILE
                                        if(! exists $ID2sequencesHash{$current_id}){
                                                print("\t\tCould not find sequence for id '$current_id' of ".$innerNode->name."\n");	
                                                exit;
                                        }	
                                        $current_sequences .= ">".$og."_$number_of_sequences_for_group\n".$ID2sequencesHash{$current_id}."\n";
                                        $number_of_sequences_for_group++;
                                        push(@ids_for_group, $current_id);
                                        $assigned_sequence_ids_hash{$current_id} = 1;
                                }
                                if($current_sequences eq q()){
                                        print("\tERROR: Could not find sequences for $og.");
                                        exit;
                                }
                                if($number_of_sequences_for_group < 2){
                                        print("Group $og contains < 2 sequences (".join(",",@ids_for_group).")");
                                        next OG;
                                }
                        # Alignment
                                my $alignment_call = $self->configuration->muscle." -maxiters 1 -quiet -diags -sv -distance1 kbit20_3 > $alignment_tmp_file";
                                if(!doAlignment({alignmentCall => $alignment_call, alignment_tmp_file => $alignment_tmp_file, current_sequences => $current_sequences})){
                                        next OG;
                                }
                        # Consensus        
                                my $consensus_string = doConsensus({   alignment_tmp_file => $alignment_tmp_file});
                                $print2FileString .= ">$og\n$consensus_string\n";   # Save consensus sequence

                                $current_number_of_og++;
                                next OG;
                }
                        # SequenceSearchInput        
                                write_to_file({file_name => $innerNode->fileInformation->sequenceSearchInputFile, text =>  $print2FileString});
                        # Consensus        
                                write_to_file({file_name => $innerNode->fileInformation->consensusFile, text =>  $print2FileString});
                        # Alignments
                           #     write_to_file({file_name => $innerNode->fileInformation->alignmentFile, text =>  $alignments2FileString});
                        
                        
                        
}
=item orthologyPredictions2OGs()

Summarizes results of orthology search of daughter nodes.
two daughter nodes will be merged into one pseudospecies
so that it can be used in an Inparanoid comparison

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
        print "\t\tSummarizing node information: ".$innerNode->name."\n";
        if(!-e $innerNode->orthologGroups->orthologyPredictionFile && !-s $innerNode->orthologGroups->orthologyPredictionFile){
                print "\t\tNo orthology predictions found\n";
                return 0;
        }
        if(-e $innerNode->orthologGroups->originalOrthologyPredictionFile && -s $innerNode->orthologGroups->originalOrthologyPredictionFile){
                print "\tconverted orthology predictions already exist. Skipping\n";
                return 1;
        }
        # Make a copy of original Inparanoid predictions
        copy($innerNode->orthologGroups->orthologyPredictionFile, $innerNode->orthologGroups->originalOrthologyPredictionFile) || die "Could not copy to ".$innerNode->orthologGroups->originalOrthologyPredictionFile."\n";
        my %all_og2id_hash;
        my %all_id_mappings;
        # Read Mapping files from daughter nodes?
        if($innerNode->leftDaughterType eq 'leaf' && $innerNode->rightDaughterType eq "leaf"){
                ;#print "\t\t\tno mappings to read (both are leafes)\n";
                #                exit;
        }
        if($innerNode->leftDaughterType eq 'innerNode'){
                #print "need og assignments from ".$innerNode->leftDaughter->name."\n";
                $innerNode->leftDaughter->orthologGroups->get_group2IDsHash($innerNode->leftDaughter,\%all_og2id_hash);
                $innerNode->leftDaughter->orthologGroups->get_group2IDsHash($innerNode->leftDaughter,\%og_assignments1_hash);
                $innerNode->leftDaughter->orthologGroups->get_group2expandedIDsHash($innerNode->leftDaughter,\%all_id_mappings);
                $mapping_cat_cmd .= " ".$innerNode->leftDaughter->orthologGroups->mappingsFile." ";
        }
        if($innerNode->rightDaughterType eq 'innerNode'){
                #print "need og assignments from ".$innerNode->rightDaughter->name." (".$innerNode->rightDaughter->orthologGroups->groupsFile.")\n";
                $innerNode->rightDaughter->orthologGroups->get_group2IDsHash($innerNode->rightDaughter,\%all_og2id_hash);
                $innerNode->rightDaughter->orthologGroups->get_group2IDsHash($innerNode->rightDaughter,\%og_assignments2_hash);
                $innerNode->rightDaughter->orthologGroups->get_group2expandedIDsHash($innerNode->rightDaughter,\%all_id_mappings);
                $mapping_cat_cmd .= " ".$innerNode->rightDaughter->orthologGroups->mappingsFile." ";
        }
        if(!keys(%all_og2id_hash) && !($innerNode->leftDaughterType eq 'leaf' && $innerNode->rightDaughterType eq "leaf")){
                print "\tCould not find OG mapping files from daughter nodes\n";
        }
        # Adding all mappings         
        $mapping_cat_cmd .= " >> ".$innerNode->orthologGroups->mappingsFile;
        #if(($innerNode->rightDaughterType eq 'innerNode' || $innerNode->leftDaughterType eq 'innerNode' )
        #&& (! -e $innerNode->orthologGroups->mappingsFile && ! -s $innerNode->orthologGroups->mappingsFile)){
                #print "concat mappings: $mapping_cat_cmd\n";
        #        system($mapping_cat_cmd);
        #}
        my $previous_og = -1;
        my @array_of_ids_for_og;  
        my %og_to_ids;
        my ($InparanoidSQLFile, $GroupsFile, $GroupsFileResolved, $MappingsFile);
        open my $INPUT_FH, '<', $innerNode->orthologGroups->originalOrthologyPredictionFile or die "Couldn't open '".$innerNode->orthologGroups->originalOrthologyPredictionFile."' $!\n";
        while(<$INPUT_FH>){
                ## e.g.
                #	1	13898	M.musculus.fa	1.000	ENSMUSP00000051825	100%
                #	1	13898	R.norvegicus.fa	1.000	ENSRNOP00000050794	100%
                chomp;
                my ($og,$bit_score,$species,$bootstrap,$id,$bootstrap_percentage) = split(/\t/,$_);
                #if($bootstrap_percentage){
                #        print "$og,$bit_score,$species,$bootstrap,$id,$bootstrap_percentage\n";
                #}
                #else{
                #        print "$og,$bit_score,$species,$bootstrap,$id\n";
                #}
                #print "$og,$bit_score,$species,$bootstrap,$id,$bootstrap_percentage\n";
                #	DEBUG("previous set to $og");
                $previous_og = $og if $previous_og == -1;
                if($og != $previous_og){
                        #print "NO MEMBERS FOR OG:  $previous_og: ".join(",",@array_of_ids_for_og)."\n" if(!scalar(@array_of_ids_for_og));
                        # convert any pseudospecies in species
                        my $og_string = "(".join(",",@array_of_ids_for_og).")";
                        my $new_cluster = $innerNode->name."$previous_og";
                        # REMEMBER OG2ID
                        foreach(@array_of_ids_for_og){
                                $og_to_ids{$_} = $previous_og;	
                        }
                        # $og_to_ids{VP345_HUMAN} = Homininae3;
                        #print "array of ids for ogs ".join(",",@array_of_ids_for_og)."\n";
                        # case: VP345_HUMAN, VP345_PANTR
                        if($innerNode->leftDaughterType eq "leaf" && $innerNode->rightDaughterType eq "leaf"){
                                ;#print "\t\tboth are leafes\n";
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
                        #print "\n\n";
                        @array_of_ids_for_og = ();
                        $previous_og = $og;
                }
                chomp($id);

                ###### CHECKING	
                #if(!exists $og_assignments1_hash{$id} && $innerNode->leftDaughterType eq "innerNode"){;#print("ID $id of species A consists of several ids, but could not resolve it. Might not be harmful\n");}
                #if(!exists $og_assignments2_hash{$id} && $innerNode->rightDaughterType eq "innerNode"){;#print("ID $id of species B consists of several ids, but could not resolve it. Might not be harmful\n");}
                ###### CHECKING END
                if(exists $og_assignments1_hash{$id} || exists $og_assignments2_hash{$id}){
                        #                        print {$OUTPUT_FH}	$innerNode->name."$og\t$bit_score\t$species\t$bootstrap\t$id\n";
                        $InparanoidSQLFile .= $innerNode->name."$og\t$bit_score\t$species\t$bootstrap\t$id\n";						
                        push(@array_of_ids_for_og,$id);
                }
                else{
                        #IF Id is not present in og2id hash then it is either a single species or an error
                        # og2id hash looks like
                        # OG4 = Homo_sapiens1, Homo_sapiens2, Mus_musculus3
                        #                        print {$OUTPUT_FH} $innerNode->name."$og\t$bit_score\t$species\t$bootstrap\t$id\n";	
                        $InparanoidSQLFile .= $innerNode->name."$og\t$bit_score\t$species\t$bootstrap\t$id\n";						
                        push(@array_of_ids_for_og,$id);
                }
        }

        #### LAST ENTRY - START
        if(!scalar(@array_of_ids_for_og)){
                #ERROR("NO MEMBERS FOR OG:  $previous_og: ".join(",",@array_of_ids_for_og));
                print("NO MEMBERS FOR OG:  $previous_og: ".join(",",@array_of_ids_for_og)."\n");
        }
        # convert any pseudospecies in species
        my $og_string = "(".join(",",@array_of_ids_for_og).")";
        my $new_cluster = $innerNode->name."$previous_og";
        # REMEMBER OG2ID
        foreach(@array_of_ids_for_og){
                $og_to_ids{$_} = $previous_og;	
        }
        # $og_to_ids{VP345_HUMAN} = Homininae3;
        #print "array of ids for ogs ".join(",",@array_of_ids_for_og)."\n";

        # case: VP345_HUMAN, VP345_PANTR
        if($innerNode->leftDaughterType eq "leaf" && $innerNode->rightDaughterType eq "leaf"){
                ;#print "\t\tboth are leafes\n";
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
        #print "\n\n";
        ##### LAST ENTRY - END
        close $INPUT_FH or die "Couldn't close '".$innerNode->orthologGroups->originalOrthologyPredictionFile."': $!\n";
        # Write output files
        # Converted Orthology predictions        
                write_to_file({file_name => $innerNode->orthologGroups->orthologyPredictionFile, text => $InparanoidSQLFile});
        # Converted Orthology predictions        
                write_to_file({file_name => $innerNode->orthologGroups->groupsFile, text => $GroupsFile});
      #  # Converted Orthology predictions        
      #          attach_to_file({file_name => $innerNode->orthologGroups->expandedGroupsFile, text => $GroupsFileResolved});
                
        # Converted Orthology predictions        
                write_to_file({file_name => $innerNode->orthologGroups->mappingsFile, text => $MappingsFile});

        if($innerNode->leftDaughterType eq "leaf" && $innerNode->rightDaughterType eq "leaf"){
                # e.g. cp Homininae.groups Homininae.mappings
                copy($innerNode->orthologGroups->groupsFile,$innerNode->orthologGroups->mappingsFile) or print "Could not copy groupsFile to mappingsFile\n";
        }
        return 1;
}

1;