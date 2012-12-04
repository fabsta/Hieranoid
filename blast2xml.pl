
use Carp;
use English;
use File::Basename;
use Getopt::Long;

use Data::Dumper;
#use Devel::Size qw(size total_size);
use strict;
use warnings;
#use IO qw(Handle File);
#use XML::Writer;



### Blast 2 xml parser
# 1. read results from UBlast
#    Format is userdefined:
#   
my $verbose = 0;
my $blast_output_file = "blast.out";
my $blast_xmloutput_file = "blast_out.xml";
my $options = GetOptions ("in|i=s" => \$blast_output_file,    # numeric
                        "out|o=s"   => \$blast_xmloutput_file,      # string
			"verbose"  => \$verbose);

#('in|i=s' => \$blast_output_file, 
#                              'out|o=s' => \$blast_xmloutput_file,
#                              );
if(! -e $blast_output_file|| ! -s $blast_output_file){
      print "Could not find Blast output file. Exiting\n";
      exit;
}
my $print_alignment = 1;

my %statistics_hash;
#my @statistics_array;
my $noOfQuery = 0;
my %seq_length_hash = ();
my %seen_hash = ();
my $previous_target_name;
my $previous_query_name;
my @hit_object;
my %duplicate_entries_hash = ();
my $count_query_hits = 1;
my $hit_counter = 1;
my $lastHit = 0;
my @statistics_array;
my $hsp_no = -1;
#print "Reading from $blast_output_file\n";


# print XML start
#&write_xml_start();
print "<BlastOutput>\n";

open my $blast_output_FH, '<', $blast_output_file or croak "Couldn't open '$blast_output_file': ";
my $firstLine = <$blast_output_FH>;
while (<$blast_output_FH>) {
#while (my $line = <$blast_output_FH>) {
        chomp();
        my $line = $_;
	#print "\tline is $line\n";
        #next if /^\s*#/;
        #next if !(/\w+/);
	next if /^query/;
	my @splitted_line = split;
	my @hsp_array;
	if(!scalar(@splitted_line) == 11){
	        print "\tError parsing line $line: has too few entries (should be 11)";
	        print "$line\n";
	        exit;
	}
	# query and target
	my ($query_name, $target_name,
	        $bit_score, $query_length,
	        $hit_length, $query_start,
	        $query_end, $hit_start,
	        $hit_end,$query_ali_length,$hit_ali_length
	        ) = (@splitted_line);
	        
	if(!defined($query_name) || !defined($target_name) || !defined($bit_score) || !defined($query_length) || !defined($hit_length) || !defined($query_start) 
	|| !defined($query_end) || !defined($hit_start) || !defined($hit_end) || !defined($query_ali_length) || !defined($hit_ali_length) ){
	        print "\tproblem parsing line $line\n";
	        #print "queryName: $query_name\n, targetName: $target_name,
        	#        bitScore: $bit_score, QueryLength: $query_length,
        	#        hitLength: $hit_length, QueryStart: $query_start,
        	#        QueryEnd: $query_end, HitStart: $hit_start,
        	#        HitEnd: $hit_end, QueryAliLength: $query_ali_length, HitAliLength: $hit_ali_length\n";
	        #exit;
	        next;
	}
	$previous_target_name = $target_name if !defined($previous_target_name);
	$previous_query_name = $query_name if !defined($previous_query_name);
	
	# write entry
	if($target_name ne $previous_target_name){
		  #print Dumper @statistics_array;
		  #exit;
		  my $hspCounter = 0;
		  next if $previous_query_name =~ /query_length/;
                  my $output_string = q{};
                  #my $writer = new XML::Writer( OUTPUT => \$output_string, DATA_MODE => 1, DATA_INDENT => 1 );
		  
		  
		  #&write_query_start($writer,$count_query_hits,$previous_query_name, $seq_length_hash{$previous_query_name});
		  print "\t<Iteration>\n";
		  print "\t <Iteration_iter-num>$count_query_hits</Iteration_iter-num>\n";
		  print "\t <Iteration_query-def>$previous_query_name</Iteration_query-def>\n";
		  print "\t <Iteration_query-len>$seq_length_hash{$previous_query_name}</Iteration_query-len>\n";

	          # Hit
	          #$writer->startTag( 'Iteration_hits' );
		  print "\t <Iteration_hits>\n";
		  
		  
		  #&write_hit_start($writer,$hit_counter,$previous_target_name, $seq_length_hash{$previous_target_name});
		  print "\t  <Hit>\n";
		  print "\t <Hit_num>$hit_counter</Hit_num>\n";
		  print "\t <Hit_id>$previous_target_name</Hit_id>\n";
		  print "\t <Hit_def>$previous_target_name</Hit_def>\n";
		  print "\t <Hit_len>$seq_length_hash{$previous_target_name}</Hit_len>\n";
		  print "\t   <Hit_hsps>\n";
		  
      
      # Hsp
		  
		  
		  
		  #print "unequal\n";
		  foreach (@statistics_array){
			my @hsp_array = $_;
			my ($bit_score, $query_length, $hit_length, $query_start, $query_end, $hit_start, $hit_end, $query_ali_length, $hit_ali_length
			      ) = ($hsp_array[0][0],$hsp_array[0][1],$hsp_array[0][2],$hsp_array[0][3],$hsp_array[0][4],$hsp_array[0][5],$hsp_array[0][6],$hsp_array[0][7],$hsp_array[0][8]);
			
			      if(!defined($bit_score) || !defined($query_length) || !defined($hit_length) || !defined($query_start) 
                      	|| !defined($query_end) || !defined($hit_start) || !defined($hit_end) || !defined($query_ali_length) || !defined($hit_ali_length) ){
                      	        print "\tproblem parsing line $line\n";
                	        #print "queryName: $query_name\n, targetName: $target_name,
                        	#        bitScore: $bit_score, QueryLength: $query_length,
                        	#        hitLength: $hit_length, QueryStart: $query_start,
                        	#        QueryEnd: $query_end, HitStart: $hit_start,
                        	#        HitEnd: $hit_end, QueryAliLength: $query_ali_length, HitAliLength: $hit_ali_length\n";
                	        #exit;
                	        next;
                	}
			$hspCounter++;
			print "\t    <Hsp>\n";
			print "\t    <Hsp_num>$hspCounter</Hsp_num>\n";
			print "\t    <Hsp_bit-score>$bit_score</Hsp_bit-score>\n";
			print "\t    <Hsp_query-from>$query_start</Hsp_query-from>\n";
			print "\t    <Hsp_query-to>$query_end</Hsp_query-to>\n";
			print "\t    <Hsp_hit-from>$hit_start</Hsp_hit-from>\n";
			print "\t    <Hsp_hit-to>$hit_end</Hsp_hit-to>\n";
			print "\t    <Hsp_align-len>$hit_ali_length</Hsp_align-len>\n";
			if($print_alignment){
	    			print "\t    <Hsp_qseq>A</Hsp_qseq>\n";
	    			print "\t    <Hsp_hseq>A</Hsp_hseq>\n";
			}
			print "\t    </Hsp>\n";
			#print "write hsp with $hspCounter, $bit_score, $query_start, $query_end, $hit_start, $hit_end, $hit_ali_length\n";
			#&write_hsp($writer, $hspCounter, $bit_score, $query_start, $query_end, $hit_start, $hit_end, $hit_ali_length); 
		  }
                        #Finish hit
			#$writer->endTag( 'Hit_hsps' );
			#$writer->endTag( 'Hit' );
			print "\t    </Hit_hsps>\n";
			print "\t    </Hit>\n";
                              #&write_hit_end($writer);
                        #Finish query
			print "  </Iteration_hits>\n";
			print "  </Iteration>\n";
                              #&write_query_end($writer);
#                        #Finish xml
                              #print $output_string;
			      #last if $count_query_hits == 1;
	    $count_query_hits++;
	    @statistics_array = ();
	    $hsp_no = 0;
	    #print "\tclean statistics array\n";
	}
	else{
	    $hsp_no++;
	}
	$previous_target_name = $target_name;
	$previous_query_name = $query_name;
      #Avoid duplicate entries:
      if($duplicate_entries_hash{$line}){
            #print "duplicate line\n";
	    next;
      }
      $duplicate_entries_hash{$line} = 1;
#     query+target+bits+ql+tl+qlo+qhi+tlo+thi+qs+ts+qrow+trow
      # Sequence lenght of query and hit
            $seq_length_hash{$query_name} = $splitted_line[3];
            $seq_length_hash{$target_name} = $splitted_line[4];
            ###
            # 0 -bit_score
            # 1 - query_length
            # 2 - hit_length
            # 3 - query_start
            # 4 - query_end
            # 5 - hit_start
            # 6 - hit_end
            # 7 - query_ali_length
            # 8 - hit_ali_length
            #print "\thsp no: $hsp_no\n";
	    #print "\thsp no: $hsp_no\n";
	    #push @{ $AoA[0] }, "wilma", "betty";
	    $statistics_array[$hsp_no][0] = $bit_score;
	    #print "$statistics_array[$hsp_no][0]\n";
	    $statistics_array[$hsp_no][1] = $query_length;
	    $statistics_array[$hsp_no][2] = $hit_length;
	    $statistics_array[$hsp_no][3] = $query_start;
	    $statistics_array[$hsp_no][4] = $query_end;
	    $statistics_array[$hsp_no][5] = $hit_start;
	    $statistics_array[$hsp_no][6] = $hit_end;
	    $statistics_array[$hsp_no][7] = $query_ali_length;
	    $statistics_array[$hsp_no][8] = $hit_ali_length;
}
close $blast_output_FH or croak "Couldn't close '$blast_output_file': ";


#&write_xml_end();
#exit;

my $output_string;
my $hspCounter = 0;
#print "forgot the last hit\n";
#print Dumper @statistics_array;
      #my $writer = new XML::Writer( OUTPUT => \$output_string, DATA_MODE => 1, DATA_INDENT => 1 );
	          print "\t<Iteration>\n";
		  print "\t <Iteration_iter-num>$count_query_hits</Iteration_iter-num>\n";
		  print "\t <Iteration_query-def>$previous_query_name</Iteration_query-def>\n";
		  print "\t <Iteration_query-len>$seq_length_hash{$previous_query_name}</Iteration_query-len>\n";

		  print "\t <Iteration_hits>\n";
		  
		  
		  #&write_hit_start($writer,$hit_counter,$previous_target_name, $seq_length_hash{$previous_target_name});
		  print "\t  <Hit>\n";
		  print "\t <Hit_num>$hit_counter</Hit_num>\n";
		  print "\t <Hit_id>$previous_target_name</Hit_id>\n";
		  print "\t <Hit_def>$previous_target_name</Hit_def>\n";
		  print "\t <Hit_len>$seq_length_hash{$previous_target_name}</Hit_len>\n";
		  print "\t   <Hit_hsps>\n";
      
		  #&write_query_start($writer,$count_query_hits,$previous_query_name, $seq_length_hash{$previous_query_name});
		  #&write_hit_start($writer,$hit_counter,$previous_target_name, $seq_length_hash{$previous_target_name});
		  
		  foreach (@statistics_array){
			my @hsp_array = $_;
			#print Dumper @hsp_array;
			#exit;
			my ($bit_score, $query_length, $hit_length, $query_start, $query_end, $hit_start, $hit_end, $query_ali_length, $hit_ali_length
			      ) = ($hsp_array[0][0],$hsp_array[0][1],$hsp_array[0][2],$hsp_array[0][3],$hsp_array[0][4],$hsp_array[0][5],$hsp_array[0][6],$hsp_array[0][7],$hsp_array[0][8]);
			$hspCounter++;
		  #$hsp, $bit_score, $query_start, $query_end, $hit_start, $hit_end, $hit_ali_length
      		  #print "write hsp with $hspCounter, $bit_score, $query_start, $query_end, $hit_start, $hit_end, $hit_ali_length\n";
		  #next;
		  print "\t    <Hsp>\n";
			print "\t    <Hsp_num>$hspCounter</Hsp_num>\n";
			print "\t    <Hsp_bit-score>$bit_score</Hsp_bit-score>\n";
			print "\t    <Hsp_query-from>$query_start</Hsp_query-from>\n";
			print "\t    <Hsp_query-to>$query_end</Hsp_query-to>\n";
			print "\t    <Hsp_hit-from>$hit_start</Hsp_hit-from>\n";
			print "\t    <Hsp_hit-to>$hit_end</Hsp_hit-to>\n";
			print "\t    <Hsp_align-len>$hit_ali_length</Hsp_align-len>\n";
			if($print_alignment){
	    			print "\t    <Hsp_qseq>A</Hsp_qseq>\n";
	    			print "\t    <Hsp_hseq>A</Hsp_hseq>\n";
			}
			print "\t    </Hsp>\n";
	    		#&write_hsp($writer, $hspCounter, $bit_score, $query_start, $query_end, $hit_start, $hit_end, $hit_ali_length); 
		  }
                        #Finish hit
			      print "\t    </Hit_hsps>\n";
			      print "\t    </Hit>\n";
                        
                              #&write_hit_end($writer);
                        #Finish query
			print "  </Iteration_hits>\n";
			print "  </Iteration>\n";
                              #&write_query_end($writer);
#                        #Finish xml
                              #print $output_string;
print "\n</BlastOutput>\n";
exit;
&write_xml_end();



exit;
print "\thwaaaat?\n";

my $total_size = int(total_size(\%statistics_hash) / 1048576);

print "\tthe hash to keep all values used ".$total_size." Mbytes\n";
print "\tfinished the program successfully\n";
exit;


my $print_counter = 0;

#print "Writing file\n";
my $last_hit = q{};




sub write_xml_start {
#      my $writer = shift;

      print "<BlastOutput>\n";
  #    print "xml_start\n";
      
 #     $writer->xmlDecl( 'UTF-8' );
 #     $writer->doctype( 'BlastOutput' );
 #     $writer->startTag( 'BlastOutput' );

      ## HEADER
#      $writer->dataElement( 'BlastOutput_program', "blastp" );
 #     $writer->dataElement( "BlastOutput_version", "blastp 2.2.22 [Sep-27-2009]" );
 #     $writer->dataElement( "BlastOutput_reference", "Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), Gapped BLAST and PSI-BLAST: a new generation of protein database search~programs;,  Nucleic Acids Res. 25:3389-3402." );
 #     $writer->dataElement( "BlastOutput_db", "query" );
 #     $writer->dataElement( "BlastOutput_query", "blastp" );
 #     $writer->dataElement( "BlastOutput_query-ID", "lcl|1_0" );
 #     $writer->dataElement( "BlastOutput_query-def", "first_hit" );
 #     $writer->dataElement( "BlastOutput_query-len", "500" );

  #          $writer->startTag( 'BlastOutput_param' );
      # Parameters
  #                $writer->startTag( 'Parameters' );
  #                $writer->dataElement( "Parameters_matrix", "BLOSUM62" );
  #                $writer->dataElement( "Parameters_expect", "blastp" );
  #                $writer->dataElement( "Parameters_gap-open", "11" );
  #                $writer->dataElement( "Parameters_gap-extend", "1" );
  #                $writer->dataElement( "Parameters_filter", "F" );
  #                $writer->endTag( 'Parameters' );
  #                $writer->endTag( 'BlastOutput_param' );
  #                $writer->startTag( 'BlastOutput_iterations' );

      
}


sub write_xml_end {
      
      print "\n</BlastOutput>\n";
#      my $writer = shift;
#      my $output = shift;
  #    print "xml_end\n";
#      $writer->endTag( 'BlastOutput_iterations' );
#      $writer->endTag( 'BlastOutput' );
#      $writer->end();
#      $output->close();
}


sub write_hit_start {
            my ($writer, $hit_counter,$hit, $seq_length)  = (shift,shift,shift,shift);
            $writer->startTag( 'Hit' );
            $writer->dataElement( "Hit_num", "$hit_counter" );                 
            $writer->dataElement( "Hit_id", $hit);                 
            $writer->dataElement( "Hit_def", $hit);                 
            $writer->dataElement( "Hit_accession", "0" );        
            $writer->dataElement( "Hit_len", $seq_length);                 
      
      # Hsp
            $writer->startTag( 'Hit_hsps' );
}

sub write_hit_end {
      my $writer = shift;
      $writer->endTag( 'Hit_hsps' );
      $writer->endTag( 'Hit' );
}




sub write_query_start {
      my ($writer, $count_query_hits,$query, $seq_length)  = (shift,shift,shift,shift);
  #    print "query_start\n";
      $writer->startTag( 'Iteration' );
      # Iteration
      $writer->dataElement( "Iteration_iter-num", "$count_query_hits" );                 
      $writer->dataElement( "Iteration_query-ID", "lcl|1_0" );                 
      $writer->dataElement( "Iteration_query-def", $query );                 
      $writer->dataElement( "Iteration_query-len", $seq_length);                 
      # Hit
      $writer->startTag( 'Iteration_hits' );
}


sub write_query_end {
            my $writer = shift;
  #          print "query_end\n";
            $writer->endTag( 'Iteration_hits' );
            $writer->startTag( 'Iteration_stat' );
            $writer->startTag( 'Statistics' );
            #$writer->dataElement( "Statistics_db-num", "3" );                 
            #$writer->dataElement( "Statistics_db-len", "1777" );                 
            #$writer->dataElement( "Statistics_hsp-len", "0" );                 
            #$writer->dataElement( "Statistics_eff-space", "0" );                 
            #$writer->dataElement( "Statistics_kappa", "0.041" );                 
            #$writer->dataElement( "Statistics_lambda", "0.267" );                 
            #$writer->dataElement( "Statistics_entropy", "0.14" );                 
            $writer->endTag( 'Statistics' );
            $writer->endTag( 'Iteration_stat' );
            $writer->endTag( 'Iteration' );
}



sub write_hsp{
      my ($writer, $hsp, $bit_score, $query_start, $query_end, $hit_start, $hit_end, $hit_ali_length) = @_;
      $writer->startTag( 'Hsp' );
      $writer->dataElement( "Hsp_num", $hsp );                 
      $writer->dataElement( "Hsp_bit-score", $bit_score );                 
      #$writer->dataElement( "Hsp_score", "0" );                 
      #$writer->dataElement( "Hsp_evalue", "0" );                 
      $writer->dataElement( "Hsp_query-from", $query_start );                 
      $writer->dataElement( "Hsp_query-to", $query_end );                 
      $writer->dataElement( "Hsp_hit-from", $hit_start );                 
      $writer->dataElement( "Hsp_hit-to", $hit_end );                 
      #$writer->dataElement( "Hsp_query-frame", "1" );                 
      #$writer->dataElement( "Hsp_hit-frame", "1" );                 
      #$writer->dataElement( "Hsp_identity", "0" );                 
      #$writer->dataElement( "Hsp_positive", "0" );                 
      #$writer->dataElement( "Hsp_gaps", "0" );                 
      $writer->dataElement( "Hsp_align-len", $hit_ali_length );                 
      ## Alignments
      if($print_alignment){
            #      $writer->dataElement( "Hsp_qseq", $statistics_hash{$query}{$hit}{$hsp}{'query_ali'} );                 
            #      $writer->dataElement( "Hsp_hseq", $statistics_hash{$query}{$hit}{$hsp}{'hit_ali'} );                 
            $writer->dataElement( "Hsp_qseq", "A" );                 
            $writer->dataElement( "Hsp_hseq", "A" );                 
      }
      # $writer->dataElement( "Hsp_midline", $statistics_hash{$query}{$hit}{'hit_end'} );                 

      $writer->endTag( 'Hsp' );
}



1;