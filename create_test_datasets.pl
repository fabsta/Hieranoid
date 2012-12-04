use Carp;
use English;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use strict;
use warnings;
use File::Copy;
use File::Temp qw/ tempfile tempdir /;
use Bio::TreeIO;
use List::Util qw(max);
use Bio::SearchIO;
use Math::Combinatorics;
use Bio::SeqIO;
use Bio::AlignIO;
use Benchmark;
require 'Configurations/Configuration.pm';

#use Time::HiRes qw( usleep ualarm gettimeofday tv_interval nanosleep
#		      clock_gettime clock_getres clock_nanosleep clock
#                      stat );
#my $alignment_in = Bio::AlignIO->new('-file' => "/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/alignmentFail/out");
#my $alignment = $alignment_in->next_aln();
#foreach my $seq ($alignment->each_seq) {
#        print $seq->id."\n";
#}
#print "well, looks like we continue\n";
#exit;
#require "lib/hieranoid_module.pm";
#goto GENETREES;
#my %allowedTreeFamGroups = ("TF101012" => 1);
my %allowedTreeFamGroups = (
							 "TF300407" => 1,
							 "TF101001" => 1,
							 "TF101004" => 1,
							 "TF101010" => 1,
							 "TF101012" => 1
);
my %convert_hash = (
	"rnor" => "Rattus_norvegicus",
	"mmus" => "Mus_musculus",
	"hsap" => "Homo_sapiens",
	"ggal" => "Gallus_gallus",
	"mbre" => "Monosiga_brevicollis",
	"atha" => "Arabidopsis_thaliana",
	"spom" => "Schizosaccharomyces_pombe",
	"scer" => "Saccharomyces_cerevisiae",
	## BENCHMARK DATA 'Altenhoff et al.'
	"RATNO" => "Rattus_norvegicus",
	"RAT"   => "Rattus_norvegicus",
	"MOUSE" => "Mus_musculus",
	"HUMAN" => "Homo_sapiens",
	"PANTR" => "Pan_troglodytes",
	"CHICK" => "Gallus_gallus",
	"FUGRU" => "Takifugu_rubripes",
	"TAKRU" => "Takifugu_rubripes",
	"ANOCA" => "Anolis_carolinensis",
	"CULQU" => "Culex_quinquefasciatus",
	"BOTFU" => "Botryotinia fuckeliana",
	"CHAGB" => "Chaetomium_globosum",
	"MAGGR" => "Magnaporthe_grisea",
	"NEUCR" => "Neurospora_crassa",
	"PODAN" => "Podospora_anserina",
	"TAEGU" => "Taeniopygia_guttata",
	"MACEU" => "Macropus_eugenii",
	"BRAFL" => "Branchiostoma_floridae",
	"ASPFU" => "Aspergillus_fumigatus",
	"BOMMO" => "Bombyx mori",
	"CIOIN" => "Ciona_intestinalis",
	"TETNG" => "Tetraodon_nigroviridis",
	"DANRE" => "Danio_rerio",
	"ANOGA" => "Anopheles_gambiae",
	"DROME" => "Drosophila_melanogaster",
	"CAEEL" => "Caenorhabditis_elegans",
	"YEAST" => "Saccharomyces_cerevisiae",
	"SCHPO" => "Schizosaccharomyces_pombe",
	"ARATH" => "Arabidopsis_thaliana",
	"CAERE" => "Caenorhabditis_remanei",
	"BOTAU" => "Bos_taurus",
	"AEDAE" => "Aedes_aegypti",
	"APIME" => "Apis_mellifera",
	"BRARE" => "Danio_rerio",
	"CAEBR" => "Caenorhabditis_briggsae",
	"CAEEL" => "Caenorhabditis_elegans",
	"CANFA" => "Canis_familiaris",
	"DROPS" => "Drosophila_pseudoobscura",
	"FUGRU" => "Fugu_rubripes",
	"GASAC" => "Gasterosteus_aculeatus",
	"MACMU" => "Macaca_mulatta",
	"MONDO" => "Monodelphis_domestica",
	"ORYSA" => "Oryza_sativa",
	"ORYLA" => "Oryzias_latipes",
	"SCHMA" => "Schistosoma_mansoni",
	"TETNI" => "Tetraodon_nigroviridis",
	"XENTR" => "Xenopus_tropicalis",
	"CIOSA" => "Ciona_savignyi",
	"DROGR" => "Drosophila_grimshawi",
	"DROMO" => "Drosophila_mojavensis",
	"DROAN" => "Drosophila_ananassae",
	"DROVI" => "Drosophila_virilis",
	"DROWI" => "Drosophila_willistoni",
	"DROPE" => "Drosophila_persimilis",
	"DROER" => "Drosophila_erecta",
	"DROSE" => "Drosophila_sechellia",
	"DROSI" => "Drosophila_simulans",
	"DROYA" => "Drosophila_yakuba",
	"ORNAN" => "Ornithorhynchus_anatinus",
	"BOVIN" => "Bos_taurus",
	"STRPU" => "Strongylocentrotus_purpuratus",
	"FELCA" => "Felis_catus",
	"CAVPO" => "Cavia_porcellus",
	"LOXAF" => "Loxodonta_africana",
	"SPETR" => "Spermophilus_tridecemlineatus",
	"TUPGB" => "Tupaia_belangeri",
	"SORAR" => "Sorex_araneus",
	"OTOGA" => "Otolemur_garnettii",
	"RABIT" => "Oryctolagus_cuniculus",
	"MYOLU" => "Myotis_lucifugus",
	"OCHPR" => "Ochotona_princeps",
	"ERIEU" => "Erinaceus_europaeus",
	"DASNO" => "Dasypus_novemcinctus",
	"ECHTE" => "Echinops_telfairi",
	"HORSE" => "Equus_caballus",
	"MICMU" => "Microcebus_murinus",
	"PONPY" => "Pongo_pygmaeus",
	"DICDI" => "Dictyostelium_discoideum",
	"BRUMA" => "Brugia_malayi"
);
my $method;
my $type;
my %options = GetOptions(
	"mode|m=s" => \$method,
	"type|t=s" => \$type

	  #			"tree|t=s" => \$species_tree_file,
	  #			"species|s=s" => \$species_folder,
	  #			"out|o=s" => \$use_outgroup
);
if ( !$method )
  {
	die "\tno method selected\n";
  }
my %methodStarter = (
					  clean_fasta_headers             => \&clean_fasta_headers,
					  get_orthobench_dataset          => \&get_orthobench_dataset,
					  blastconvert                    => \&blastconvert,
					  simulations                     => \&simulations,
					  test_score_symmetry             => \&test_score_symmetry,
					  compare2orthobench              => \&compare2orthobench,
					  treefam                         => \&treefam,
					  testing                         => \&testing,
					  inparanoid_testing              => \&inparanoid_testing,
					  inparanoid_testing_pairwise     => \&inparanoid_testing_pairwise,
					  reformat_orthology_predictions  => \&reformat_orthology_predictions,
					  species_overlap                 => \&species_overlap,
					  reduced_inparanoid_dataset      => \&reduced_inparanoid_dataset,
					  testUblastBlastOverlap          => \&testUblastBlastOverlap,
					  BrigitteTrees2species           => \&BrigitteTrees2species,
					  mapOMA2RefOGs                   => \&mapOMA2RefOGs,
					  compare2orthobench              => \&compare2orthobench,
					  compare2orthobench_pairwise     => \&compare2orthobench_pairwise,
					  OrthoMCL2OrthoBench             => \&OrthoMCL2OrthoBench,
					  get_inparanoid_dataset          => \&get_inparanoid_dataset,
					  find_longest_transcript         => \&find_longest_transcript,
					  testInparanoidTime              => \&testInparanoidTime,
					  build_small_extended_orthobench => \&build_small_extended_orthobench
);
if ( not exists $methodStarter{$method} )
  {
	print "\tselected method $method does not exist\n";
	exit;
  }
$methodStarter{$method}->();

sub readSequenceFileIntoHash
  {
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
	die("Input file $sequence_file is empty or missing\n")
	  if ( !-e $sequence_file || !-s $sequence_file );

	#print "\t Reading sequences from $sequence_file\n";
	open my $SEQUENCE_FILE, '<', $sequence_file
	  or croak "Couldn't open '$sequence_file': $OS_ERROR";
	### READING FROM FILE
	while ( my $line = <$SEQUENCE_FILE> )
	  {
		chomp($line);
		if ( $line =~ /^>/ )
		  {
			my @words = split( /\s/, $line );
			$seq_id = $words[0];
			$seq_id =~ s/>//;

			#                        print "seq_id $seq_id\n" if $seq_id =~ /LCT/;
			$sequence_href->{$seq_id} = q();
		  }
		else
		  {
			$sequence_href->{$seq_id} .= $line;
		  }
	  }
	### CLOSING FILE
	close $SEQUENCE_FILE or croak "Couldn't close '$sequence_file': $OS_ERROR";

	#        print $sequence_href->{"LCT_PANTR"};
	#my $logger = get_logger("READING FASTA FILES");
	#print("Read ".keys( %{$sequence_href} )." sequences into hash");
	if ( !keys( %{$sequence_href} ) )
	  {
		print("Could not read fasta entries from $sequence_file");
		exit;
	  }
	return 1;
  }

sub readGtfFile
  {
	my ($arg_ref)         = @_;
	my $gtf_file          = $arg_ref->{gtf_file};
	my $prot2gene_hashref = $arg_ref->{prot2gene_hashref};
	open my $READ_FH, '<', $gtf_file
	  or croak "Couldn't open " . $gtf_file . ": $OS_ERROR";
	my $counter = 0;
	### READING FROM FILE
	while ( my $line = <$READ_FH> )
	  {
		my @splitted_line = split( /\t/, $line );
		my ( $gene_id, $transcript_id, $protein_id );
		foreach (@splitted_line)
		  {
			if (/gene_id\s+\"(.*?)\";/)       { $gene_id       = $1; }
			if (/transcript_id\s+\"(.*?)\";/) { $transcript_id = $1; }
			if (/protein_id\s+\"(.*?)\";/)    { $protein_id    = $1; }
		  }
		if ( $gene_id eq '' || $transcript_id eq '' )
		  {
			die "gene: $gene_id, transcript: $transcript_id, protein: $protein_id\n";
		  }
		next if !$protein_id;

		#print "saving $protein_id - $gene_id\n";
		#exit;
		$prot2gene_hashref->{$protein_id} = $gene_id;
	  }
	close $READ_FH or croak "Couldn't close '" . $gtf_file . "': $OS_ERROR";
	if ( !keys(%$prot2gene_hashref) )
	  {
		print("Empty hash prot2gene_hashref found\n");
		exit;
	  }
	return 1;
  }

sub find_longest_transcript
  {
	my $fastaDir = "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/ensemblv60_12species";
	my $gtfDir =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/ensemblv60_12species/gtf";
	my $outputDir =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/ensemblv60_12species/longestTranscript";
	my %fasta2gtf_mapping = (
							  "Caenorhabditis_elegans.fa"  => "Caenorhabditis_elegans.WS210.60.gtf",
							  "Canis_familiaris.fa"        => "Canis_familiaris.BROADD2.60.gtf",
							  "Ciona_intestinalis.fa"      => "Ciona_intestinalis.JGI2.60.gtf",
							  "Danio_rerio.fa"             => "Danio_rerio.Zv9.60.gtf",
							  "Drosophila_melanogaster.fa" => "Drosophila_melanogaster.BDGP5.25.60.gtf",
							  "Gallus_gallus.fa"           => "Gallus_gallus.WASHUC2.60.gtf",
							  "Homo_sapiens.fa"            => "Homo_sapiens.GRCh37.60.gtf",
							  "Monodelphis_domestica.fa"   => "Monodelphis_domestica.BROADO5.60.gtf",
							  "Mus_musculus.fa"            => "Mus_musculus.NCBIM37.60.gtf",
							  "Pan_troglodytes.fa"         => "Pan_troglodytes.CHIMP2.1.60.gtf",
							  "Rattus_norvegicus.fa"       => "Rattus_norvegicus.RGSC3.4.60.gtf",
							  "Tetraodon_nigroviridis.fa"  => "Tetraodon_nigroviridis.TETRAODON8.60.gtf",
	);
	while ( my ( $fastaFile, $gtfFile ) = each(%fasta2gtf_mapping) )
	  {
		my $output_file = "$outputDir/$fastaFile";
		$fastaFile = "$fastaDir/$fastaFile";
		$gtfFile   = "$gtfDir/$gtfFile";
		my %prot2geneHash;
		my %longestProteinHash;
		my %id2sequences_hash;
		print "$fastaFile - $gtfFile \n";

		# Read sequence file
		print "reading sequences into memory...\n";
		readSequenceFileIntoHash(
								  {
									sequence_file => $fastaFile,
									sequence_href => \%id2sequences_hash
								  }
		);
		print "\t\t\tdone\n";

		# read gtf file
		print "reading gtf into memory...\n";
		readGtfFile( { gtf_file => $gtfFile, prot2gene_hashref => \%prot2geneHash } );
		print "\t\t\tdone\n";
		print "Iterating over " . keys(%id2sequences_hash) . " sequence ids...\n";
		foreach ( keys(%id2sequences_hash) )
		  {
			my $id             = $_;
			my $gene4id        = $prot2geneHash{$id};
			my $sequenceLength = length( $id2sequences_hash{$_} );
			if ( not exists $longestProteinHash{$gene4id} )
			  {
				$longestProteinHash{$gene4id}{"id"}     = $id;
				$longestProteinHash{$gene4id}{"length"} = $sequenceLength;
			  }
			elsif ( $sequenceLength > $longestProteinHash{$gene4id}{"length"} )
			  {
				$longestProteinHash{$gene4id}{"id"}     = $id;
				$longestProteinHash{$gene4id}{"length"} = $sequenceLength;
			  }

			#print "$id -> $gene4id -> $sequenceLength\n";
			#exit;
		  }
		print "\t\t\tdone\n";
		print "\tfound longest protein/transcript for " . keys(%longestProteinHash) . "\n";
		print "\twriting new fasta file\n";
		my $sequences2print;
		foreach ( keys(%longestProteinHash) )
		  {
			my $id          = $longestProteinHash{$_}{"id"};
			my $sequence4id = $id2sequences_hash{$id};
			$sequences2print .= ">$id\n$sequence4id\n";
		  }
		write_to_file( { file_name => $output_file, text => $sequences2print } );
	  }
  }

#        $methodStarter{$method};
sub clean_fasta_headers
  {
	my $fastaDir = "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/ensemblv60_12species/";
	my @fastaFilesToCheck = glob("$fastaDir/*");
	foreach my $fastaFile (@fastaFilesToCheck)
	  {

		#next if!($fastaFile =~ /magnipapillata/);
		my ( undef, $fastatmpfile ) = tempfile();
		my $seq_in  = Bio::SeqIO->new( -format => 'fasta',          -file   => $fastaFile );
		my $seq_out = Bio::SeqIO->new( -file   => ">$fastatmpfile", -format => 'fasta' );
		print "\tLooking at $fastaFile\n";
		while ( my $inseq = $seq_in->next_seq )
		  {

		#print "primary: ".$inseq->primary_id."\nDescription: ".$inseq->desc."\ndisplay: ".$inseq->display_id."\n";
		#exit;
			$inseq->desc('');
			my $sequence = $inseq->seq;

			#print "$sequence before\n" if $inseq->primary_id eq 'ENSP00000421854';
			$sequence =~ s/\*//g;

			#print "$sequence after\n" if $inseq->primary_id eq 'ENSP00000421854';
			#exit if $inseq->primary_id eq 'ENSP00000421854';
			$inseq->seq($sequence);
			$seq_out->write_seq($inseq);
		  }
		if ( -e $fastatmpfile && -s $fastatmpfile )
		  {
			`mv $fastatmpfile $fastaFile`;
		  }
		else
		  {
			print "\tCould not check fasta file $fastaFile\n";
			exit;
		  }
	  }
	exit;
  }

sub get_orthobench_dataset
  {
	my $address   = "http://eggnog.embl.de/orthobench/alignments/";
	my $outputDir = "/Users/fab/Documents/documents/work/current_projects/hieranoid/predictions/orthobench";
	foreach ( 1 .. 70 )
	  {
		my $currentGroup = $_;
		if ( length($currentGroup) == 1 )
		  {
			$currentGroup = "0" . $currentGroup;
		  }
		my $outputFile = "$outputDir/RefOG0$currentGroup.fasta";
		my $group      = "$address/RefOG0$currentGroup.fasta";
		my $wgetCmd    = "wget $group  -O $outputFile";
		print "$wgetCmd\n";
		`$wgetCmd`;

		#exit;
	  }
  }

sub get_inparanoid_dataset
  {
	my $address = "http://inparanoid.sbc.su.se/download/7.0_current/sqltables/";
	my $outputDir =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/predictions/inparanoid/sqltables/";
	my $allSpeciesPairs =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/InparanoidAllPairs.txt";

	#http://inparanoid.sbc.su.se/download/7.0_current/sqltables/sqltable.B.floridae.fa-M.domestica.fa
	my @allowed_species = (
		'C.elegans',
		'D.melanogaster',
		'C.intestinalis',
		'D.rerio',
		'T.nigroviridis',
		'G.gallus',
		'M.domestica',
		'M.musculus',
		'R.norvegicus',
		'C.familiaris',
		'P.troglodytes',
		'H.sapiens',

		#'N.vectensis',
		#'M.brevicollis',
		#'T.adhaerens'
	);

	#my @SpeciesPairs = combine(2,@allowed_species);
	#for (my $i=0; $i< scalar(@SpeciesPairs); $i++){
	#	my @sorted = sort(@{$SpeciesPairs[$i]});
	#	my ($speciesA,$speciesB) = ($sorted[0],$sorted[1]);
	#	my $outputFile = "$outputDir/sqltable.$speciesA.fa-$speciesB.fa";
	#	my $pair = $address."sqltable.$speciesA.fa-$speciesB.fa";
	#	my $wgetCmd = "wget $pair  -O $outputFile";
	#        print "$wgetCmd\n";
	#	#next;
	#        `$wgetCmd`;
	#}
	my $string2save;
	my $og_counter = 1;
	my $previousOG = -1;

	# concatenate all files but change og numbers
	foreach my $pairwiseFile ( glob("$outputDir/*") )
	  {
		print "\tlooking at $pairwiseFile\n";
		open my $OG_FH, '<', $pairwiseFile
		  or croak "Couldn't open " . $pairwiseFile . ": $OS_ERROR";
		### READING FROM FILE
		while ( my $line = <$OG_FH> )
		  {

			#print $line;
			next if $line =~ /^$/;
			next if $line =~ /^#/;
			chomp($line);
			my @splitted_line = split( /\t/, $line );
			my $member        = $splitted_line[4];
			my $og            = $splitted_line[0];
			$previousOG = $og if $previousOG eq '-1';
			if (    !defined $og
				 || $og eq ''
				 || !defined $member
				 || $member eq '' )
			  {
				print( "Could not parse line $line in file " . $pairwiseFile . ". Values are $og,$member\n" );
				exit;
			  }

			# Map IDs to RefOG groups if possible
			# We can later easily check if a ID the e.g. hieranoid prediction is also in the RefOG dataset
			if ( $previousOG ne $og )
			  {
				$og_counter++;
				$previousOG = $og;
			  }

			#print "\treplace group $og with $og_counter\n";
			$splitted_line[0] = $og_counter;
			$string2save .= join( "\t", @splitted_line ) . "\n";

			#print "--> ".join("\t",@splitted_line)."\n";
			#die "enough is enough!\n" if $og_counter > 20;
		  }
		close $OG_FH
		  or croak "Couldn't close '" . $pairwiseFile . "': $OS_ERROR";
	  }
	write_to_file( { file_name => $allSpeciesPairs, text => $string2save } );
  }

sub blastconvert
  {
	my $blast_report = new Bio::SearchIO( -format => 'blastxml',
										  -file   => 'tmp_stuff/blast_parser_test/original_inparanoid/xml' );
	my $result      = $blast_report->next_result;
	my $writertable = new Bio::SearchIO::Writer::ResultTableWriter();
	my $outtab = new Bio::SearchIO( -writer => $writertable,
									-file   => ">searchio.tab" );

	# get a result from Bio::SearchIO parsing or build it up in memory
	$outtab->write_result($result);
	exit;
  }

sub create_orthoxml
  {

#
#		my %predictions_hash = (
#				"MuRa" => {
#							"mmus" => 1,
#							"rnor" => 1
#				},
#				"HoMuRa" => {
#							"hsap" => 1,
#							"mmus" => 1,
#							"rnor" => 1
#				},
#				"HoMuRaGa" => {
#							"ggal" => 1,
#							"hsap" => 1,
#							"mmus" => 1,
#							"rnor" => 1
#				},
#				"HoMuRaGaMo" => {
#							"mbre" => 1,
#							"ggal" => 1,
#							"hsap" => 1,
#							"mmus" => 1,
#							"rnor" => 1
#				},
#				"SiSa" => {
#							"spom" => 1,
#							"scer" => 1
#				},
#				"SiSaAr" => {
#							"atha" => 1,
#							"spom" => 1,
#							"scer" => 1
#				}
#		);
#
#
#		my $current_og;
#		my %species_hash = ();
#		my %og_hash = ();
#		my $current_sequence = q();
#		my $counter = 0;
#
#
#		####
#		# ORTHOMCL
#		####
#		my ($current_species,$current_id);
#		my $input_file = "/Users/fab/Documents/documents/projects/hieranoid/testing/benchmark_datasets/eukaryota.fa";
#		### OPENING FILE
#			open my $IN, '<', $input_file or croak "Couldn't open '$input_file': $!\n";
#			#print "\tbla\n";
#			### READING FROM FILE
#			while (my $line = <$IN>) {
#		#		print $line;
#		#		if($line =~ /^>(\w+)\|(\w+[\.|-]*\w+)/){
#				if($line =~ /^>(([A-Z]+)\d+)/){
#					my ($id,$species) = ($1,$2);
#					croak qq{Could not parse line $line ($species, $id)\n} if (!defined $species || !defined $id);
#					#print "species: $species, id: $id\n";
#					#if($counter++ > 100){
#						#print Dumper %species_hash;
#						#print Dumper %og_hash;
#					#	exit;
#					#}
#					#exit;
#					if($current_sequence){
#					#	print "\tsetting values\n species_hash{$current_species}{$current_id} = $current_sequence;";
#						$species_hash{$current_species}{$current_id} = ">$current_id\n".$current_sequence;
#						#$og_hash{$current_og}{$current_species}{$current_id} = 1;
#						$current_sequence = q();
#					}
#					#print "\tset to $current_species,$current_id\n";
#					($current_species,$current_id) = ($species,$id);
#				}
#				else{
#					$current_sequence .= $line;
#				}
#			}
#		#print Dumper %species_hash;
#		#print Dumper %og_hash;
#		#exit;
#			### CLOSING FILE
#			close $IN or croak "Couldn't close '$input_file': $!\n";
#
#		print "\twriting fasta files";
#		foreach my $species(keys(%species_hash)){
#			next if !exists $convert_hash{$species};
#			my $species_file = "species/benchmark/$convert_hash{$species}.fa";
#			print "\twriting $convert_hash{$species}\n";
#			### OPENING FILE
#			open my $out, '>', $species_file or croak "Couldn't open '$species_file': $!\n";
#
#			foreach my $id(keys(%{$species_hash{$species}})){
#				print {$out} $species_hash{$species}{$id} or croak "Couldn't write '$species_file': $!\n";
#			}
#
#		### CLOSING FILE
#		close $out or croak "Couldn't close '$species_file': $!\n";
#		}
#
#
#
#		exit;
#
#
#		### THIS IS FOR ORTHOMCL
#		my $current_og;
#		my %species_hash = ();
#		my %og_hash = ();
#		my ($current_species,$current_id);
#		my @files =glob("tmp_stuff/orthomcl/*");
#		foreach my $file(@files){
#			$current_og = basename($file);
#			print "\tinspecting $file ($current_og)\n";
#			my $current_sequence = q();
#			### OPENING FILE
#			open my $IN, '<', $file or croak "Couldn't open '$file': $!\n";
#			#print "\tbla\n";
#			### READING FROM FILE
#			while (my $line = <$IN>) {
#		#		print $line;
#				if($line =~ /^>(\w+)\|(\w+[\.|-]*\w+)/){
#					my ($species,$id) = ($1,$2);
#					croak qq{Could not parse line $line ($species, $id)\n} if (!defined $species || !defined $id);
#					#print "species: $species, id: $id\n";
#					if($current_sequence){
#					#	print "\tsetting values\n species_hash{$current_species}{$current_id} = $current_sequence;";
#						$species_hash{$current_species}{$current_id} = ">$current_id\n".$current_sequence;
#						$og_hash{$current_og}{$current_species}{$current_id} = 1;
#						$current_sequence = q();
#					}
#					($current_species,$current_id) = ($species,$id);
#				}
#				else{
#					$current_sequence .= $line;
#				}
#			}
#		#print Dumper %species_hash;
#		#print Dumper %og_hash;
#		#exit;
#			### CLOSING FILE
#			close $IN or croak "Couldn't close '$file': $!\n";
#
#
#		}
#
#		print "\twriting fasta files";
#		foreach my $species(keys(%species_hash)){
#			next if !exists $convert_hash{$species};
#			my $species_file = "species/orthomcl_fasta/$convert_hash{$species}.fa";
#			print "\twriting $convert_hash{$species}\n";
#			### OPENING FILE
#			open my $out, '>', $species_file or croak "Couldn't open '$species_file': $!\n";
#
#			foreach my $id(keys(%{$species_hash{$species}})){
#				print {$out} $species_hash{$species}{$id} or croak "Couldn't write '$species_file': $!\n";
#			}
#
#		### CLOSING FILE
#		close $out or croak "Couldn't close '$species_file': $!\n";
#		}
#
#
#		# predictions in sqltable format
#
#		print "\twriting predictions\n";
#		#print Dumper %og_hash;
#			foreach my $comparison (sort keys(%predictions_hash)){
#				my $comparison_sqltable_file = "predictions/inparanoid/sqltable.$comparison";
#				open my $out_sqltable, '>', $comparison_sqltable_file or croak "Couldn't open '$comparison_sqltable_file': $!\n";
#		print "\twriting predictions ($comparison_sqltable_file)\n";
#
#				foreach my $og(keys(%og_hash)){
#						my $print_string = q{};
#						my $range = 50;
#						my $minimum = 13000;
#						my $bit_score = int(rand($range)) + $minimum;
#					foreach my $species(keys(%{$og_hash{$og}})){
#						next if !exists $predictions_hash{$comparison}{$species};
#							foreach my $id_from_species (keys(%{$og_hash{$og}{$species}})){
#
#								my $og_short = $og;
#								$og_short =~ s/og4_//;
#								$print_string .= qq{$og_short\t$bit_score\t$convert_hash{$species}\t1.000\t$id_from_species\t100\%\n};
#							}
#					}
#				#	print Dumper $og_hash{$og};
#					#print $print_string;
#					print {$out_sqltable} $print_string;
#					#exit;
#				}
#			close $out_sqltable or croak "Couldn't close '$comparison_sqltable_file': $!\n";
#			}
#			exit;
##		 my $orthoxml_start = <<END;
##		<?xml version="1.0" encoding="utf-8"?>\n
##		<orthoXML xmlns="http://orthoXML.org/0.2" version="0.2" origin="inparanoid" originVersion="7.0" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
##		  xsi:schemaLocation="http://orthoXML.org/0.2 http://www.orthoxml.org/0.2/orthoxml.xsd">
##		  <notes>
##		    Example OrthoXML file. Stripped down version of a real InParanoid 7.0 file.
##		  </notes>
##		END
##		my $orthoxml_end = "</orthoXML>\n";
#
#		#print Dumper %{$og_hash{'og4_10034'}};
#		foreach my $comparison (sort keys(%predictions_hash)){
#			my (%comparison_hash,%comparison_genes_hash,%comparison_og) = ();
#			my $comparison_file = "predictions/orthomcl/Inparanoid.$comparison.orthoXML";
#
#			print "\twriting predictions for $comparison ($comparison_file)\n";
#			my $og_id_counter = 1;
#			foreach my $og(keys(%og_hash)){
#				foreach my $species(keys(%{$og_hash{$og}})){
#					next if !exists $predictions_hash{$comparison}{$species};
#					foreach my $id_from_species (keys(%{$og_hash{$og}{$species}})){
#						$comparison_hash{$species}{$og}{$og_id_counter} = $id_from_species;
#						$comparison_genes_hash{$species}{$og_id_counter} = $id_from_species;
#						$comparison_og{$og}{$og_id_counter} = $id_from_species;
#						$og_id_counter++;
#					}
#				}
#			}
#				#print Dumper %og_hash;
#				#exit;
#
#			open my $out, '>', $comparison_file or croak "Couldn't open '$comparison_file': $!\n";
#			print {$out} $orthoxml_start;
#
#		#GENES
#			foreach my $species(keys(%comparison_hash)){
#				print $species."\n";
#				print {$out} " <species name=\"".$convert_hash{$species}."\" NCBITaxId=\"9606\">\n<database name=\"WormBase\" >\n";
#				foreach my $gene(sort keys(%{$comparison_genes_hash{$species}})){
#					print {$out} "<gene id =\"$gene\" geneId=\"".$comparison_genes_hash{$species}{$gene}."\" protId=\"".$comparison_genes_hash{$species}{$gene}."\"/>\n";
#				}
#				print {$out} " </species>\n";
#			}
#			print {$out} " <groups>\n";
#
#			foreach my $og(keys(%og_hash)){
#		    print {$out} "<orthologGroup id=\"$og\">\n";
#				foreach my $id(sort keys(%{$comparison_og{$og}})){
#					print {$out} "<geneRef id=\"$id\">\n";
#					print {$out} "<\geneRef>\n";
#				}
#		    print {$out} "</orthologGroup>\n";
#
#			}
#			print {$out} "</groups>\n";
#
#			print {$out} $orthoxml_end;
#			### CLOSING FILE
#			close $out or croak "Couldn't close '$comparison_file': $!\n";
#		#	die "\twrote to $comparison_file\n";
#
#		}
  }

sub reduced_inparanoid_dataset
  {
	my %file2species_mapping = (
								 "A.mellifera.fa"    => "Apis_mellifera.fa",
								 "A.thaliana.fa"     => "Arabidopsis_thaliana.fa",
								 "C.elegans.fa"      => "Caenorhabditis_elegans.fa",
								 "D.melanogaster.fa" => "Drosophila_melanogaster.fa",
								 "H.sapiens.fa"      => "Homo_sapiens.fa",
								 "M.musculus.fa"     => "Mus_musculus.fa",
								 "N.vectensis.fa"    => "Nematostella_vectensis.fa",
								 "P.troglodytes.fa"  => "Pan_troglodytes.fa",
								 "S.cerevisiae.fa"   => "Saccharomyces_cerevisiae.fa",
								 "X.tropicalis.fa"   => "Xenopus_tropicalis.fa",
	);
	my $inparanoid_predictions_directory =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/inparanoid_data/predictions";
	my $inparanoid_reduced_ids_dir =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/inparanoid_data/reduced_ids/";
	my $inparanoid_small_files =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/inparanoid_data/small/";
	my $inparanoid_small_predictions =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/inparanoid_data/10species_200genes/predictions/10species.pred";
	my $inparanoid_full_files =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/inparanoid_data/fasta/";
	my $reference_inparanoid_predictions = "sqltable.H.sapiens.fa-S.cerevisiae.fa";
	my $file_with_ids                    = "reduced_inparanoid_ids.txt";
	mkdir($inparanoid_reduced_ids_dir) if !-e $inparanoid_reduced_ids_dir;
	my @other_predictions = (
							  "sqltable.A.mellifera.fa-H.sapiens.fa",
							  "sqltable.A.thaliana.fa-H.sapiens.fa",
							  "sqltable.C.elegans.fa-H.sapiens.fa",
							  "sqltable.D.melanogaster.fa-H.sapiens.fa",
							  "sqltable.H.sapiens.fa-M.musculus.fa",
							  "sqltable.H.sapiens.fa-N.vectensis.fa",
							  "sqltable.H.sapiens.fa-P.troglodytes.fa",
							  "sqltable.H.sapiens.fa-S.cerevisiae.fa",
							  "sqltable.H.sapiens.fa-X.tropicalis.fa"
	);
	my %reference_predictions_ids = ();
	my %species_id_hash           = ();
	my %species2id_hash           = ();
	my %humanReferenceIDs;
	my @collectLines;
	my $no_of_groups = 200;
	&read_pairwise_results(
						  {
							profile_file => "$inparanoid_predictions_directory/$reference_inparanoid_predictions",
							id_to_cluster_hashref  => \%reference_predictions_ids,
							no_of_groups           => $no_of_groups,
							humanReferenceIDs_href => \%humanReferenceIDs,
							check_group            => 1,
						  }
	);
	print "\tFound " . keys(%humanReferenceIDs) . " human reference IDS\n";

	# print Dumper %humanReferenceIDs;
	# exit;
	if ( !keys(%reference_predictions_ids) )
	  {
		die
"Could not read $no_of_groups predictions from $inparanoid_predictions_directory/$reference_inparanoid_predictions";
	  }
	foreach my $prediction_file (@other_predictions)
	  {
		print "Reading predictions from $prediction_file\n";
		&read_pairwise_results(
								{
								  profile_file           => "$inparanoid_predictions_directory/$prediction_file",
								  id_to_cluster_hashref  => \%reference_predictions_ids,
								  no_of_groups           => 1000000,
								  check_group            => 1,
								  species2id_href        => \%species2id_hash,
								  humanReferenceIDs_href => \%humanReferenceIDs
								}
		);
	  }

	#print "\tread ".scalar(@collectLines)." lines\n";
	#write_to_file({file_name => $inparanoid_small_predictions, text => join("",@collectLines)});
	#exit;
	print "Read a total of " . keys(%reference_predictions_ids) . " ids\n";
	foreach my $species ( keys(%species2id_hash) )
	  {
		my $species_file = "$inparanoid_reduced_ids_dir/$species";
		my $fasta_file   = "$inparanoid_full_files/" . $file2species_mapping{ basename($species_file) };
		open my $file_with_ids_FH, '>', $species_file
		  or croak "Couldn't open '$species_file': ";
		foreach ( keys( %{ $species2id_hash{$species} } ) )
		  {
			next if /^$/;
			chomp;
			print {$file_with_ids_FH} $_ . "\n"
			  or croak "Couldn't write '$species_file': ";
		  }
		close $file_with_ids_FH or croak "Couldn't close '$species_file': ";

		# avoid duplicates
		`sort $species_file | uniq > lol`;
		`mv lol $species_file`;
		## Now extract all
		my $esl_call = "esl-sfetch -f $fasta_file $species_file > $inparanoid_small_files/"
		  . $file2species_mapping{ basename($species_file) };
		print "Fetching fasta sequences\n$esl_call\n";
		`$esl_call`;
	  }
	print "Everything done. Thanks, it was a pleasure!\n";
  }

sub read_pairwise_results
  {
	my ($arg_ref)              = @_;
	my $profile_file           = $arg_ref->{profile_file};
	my $id_to_cluster_hashref  = $arg_ref->{id_to_cluster_hashref};
	my $cluster_to_id_hashref  = $arg_ref->{cluster_to_id_hashref};
	my $species2id_href        = $arg_ref->{species2id_href};
	my $check_group            = $arg_ref->{check_group};
	my $no_of_groups           = $arg_ref->{no_of_groups};
	my $humanReferenceIDs_href = $arg_ref->{humanReferenceIDs_href};
	my $group_counter          = 0;
	my $collect_refID          = 0;
	my %cluster_to_id_hash;

	if ( !-e $profile_file || !-s $profile_file )
	  {
		print("\t\t[READ PAIRWISE RESULTS] No results file found ('$profile_file')");
		return ();
	  }
	if ( !keys(%$humanReferenceIDs_href) )
	  {
		print "\tReference IDs have to be collected\n";
		$collect_refID = 1;

		#print "\tcollect is $collect_refID\n";
	  }
	my $previous_og = -1;
	my @ids_for_og  = ();
	my %species_ids = ();
	open my $INPUT_FH, '<', $profile_file
	  or die "Couldn't open '$profile_file' $!\n";
	while (<$INPUT_FH>)
	  {
		chomp;

		#push(@$collectLines_aref,$_);
		## e.g.
		#	1	13898	M.musculus.fa	1.000	ENSMUSP00000051825	100%
		#	1	13898	R.norvegicus.fa	1.000	ENSRNOP00000050794	100%
		my ( $og, $bit_score, $species, $bootstrap, $id, $bootstrap_percentage ) = split( /\t/, $_ );

		#print "$og,$bit_score,$species,$bootstrap,$id \n";
		$species_ids{$species}{$id} = 1;
		$previous_og = $og if $previous_og == -1;
		if ( $og != $previous_og )
		  {
			chomp($id);
			my $good_group = 0;
			if ($check_group)
			  {
				foreach (@ids_for_og)
				  {
					if ( !$collect_refID
						 && exists $humanReferenceIDs_href->{$_} )
					  {

						#exists $id_to_cluster_hashref->{$_}){
						print "\t has human reference $_\n";

						#exit;
						$good_group = 1;
						last;
					  }
					if ($collect_refID)
					  {
						$good_group = 1;

						#  print "\t saved human reference id '$_'\n";
						# exit;
						last;
					  }
				  }
			  }

			# Keep this group
			if ( $good_group || !$check_group )
			  {
				foreach (@ids_for_og)
				  {
					$id_to_cluster_hashref->{$_} = $previous_og;
					push( @{ $cluster_to_id_hashref->{$previous_og} }, $_ );

					#print "\tadding cluster2id: $ previous_og - $_\n";
					#exit;
					foreach my $species ( keys(%species_ids) )
					  {
						foreach my $id ( keys( %{ $species_ids{$species} } ) )
						  {
							$species2id_href->{$species}{$id} = 1;
							if ( $species eq 'H.sapiens.fa' && $collect_refID )
							  {
								$humanReferenceIDs_href->{$id} = 1;
							  }
						  }
					  }
				  }

				#print "".join(",", @{$cluster_to_id_hash{$previous_og}} )."\n";
				$group_counter++;
			  }
			$previous_og = $og;
			@ids_for_og  = ();
			%species_ids = ();
		  }
		push( @ids_for_og, $id );
		last if $group_counter >= $no_of_groups;
	  }

#	print "\tRead ".keys(%{$id_to_cluster_hashref})." Ids assigned to ".keys(%{$group_href})." orthologous groups\n";
	close $INPUT_FH or die "Couldn't close '$profile_file': $!\n";
	return 1;
  }

sub simulations
  {
	#####
	# SIMULATIONS
	#####
	# THIS PART IS IMPORTANT FOR GENERATING THE GENE TREES AND SEQUENCE ALIGNMENTS
	my $number_of_gene_trees     = 300;
	my $length_of_alignment      = 100;
	my $folder_to_put_trees      = "./trees";
	my $folder_to_put_alignments = "./alignments";
	my %species_hash             = ();
	my $species_tree             = "species_tree";
	for ( my $count = 1 ; $count <= $number_of_gene_trees ; $count++ )
	  {
		my $temp_length_of_alignment = $length_of_alignment * ( ( $count % 10 ) + 1 );
		my $gene_tree                = $folder_to_put_trees . "/gene_tree$count.tre";
		my $gene_tree_parameter      = $folder_to_put_trees . "/gene_tree$count.param";
		my $gene_alignment           = $folder_to_put_alignments . "/gene_alignment$count.fa";
		my $gene_alignment_parameter = $folder_to_put_alignments . "/gene_alignment$count.param";
		print "gene tree $count\n";
		my $aladen_call = "aladen  -q -o $gene_tree -To $gene_tree_parameter  -Gl 12 -Gh 30 -Gt -Hi $species_tree";

		#print "$aladen_call\n";
		`$aladen_call`;
		my $beep_generateSeqData_call =
"beep_generateSeqData  -o $gene_alignment -To $gene_alignment_parameter -Sm JTT $gene_tree $temp_length_of_alignment";

		#print "$beep_generateSeqData_call\n";
		`$beep_generateSeqData_call`;

		#exit;
	  }

	# CREATE SPECIES FASTA FILES:
	my @files = <$folder_to_put_alignments/*.fa>;
	foreach my $current_alignment (@files)
	  {
		my ( $current_species, $current_sequence, $current_id );
		open my $ALN_FH, '<', $current_alignment
		  or croak "Couldn't open '$current_alignment': $!\n";
		while ( my $line = <$ALN_FH> )
		  {
			$current_alignment =~ /gene_alignment(\d+).fa/;
			my $gene_number = $1;

			#		print $line;
			if ( $line =~ /^>(\w+_\w+)_(\d+)/ )
			  {
				my ( $species, $id ) = ( $1, $2 );
				croak qq{Could not parse line $line ($species, $id)\n}
				  if ( !defined $species || !defined $id );

				#print "species: $species, id: $id\n";
				if ($current_sequence)
				  {

					#	print "\tsetting values\n species_hash{$current_species}{$current_id} = $current_sequence;";
					$species_hash{$current_species}{$current_id} = ">$current_id\n" . $current_sequence;

					#$og_hash{$current_og}{$current_species}{$current_id} = 1;
					$current_sequence = q();
				  }
				( $current_species, $current_id ) = ( $species, $species . "_" . $gene_number . "_" . $id );
			  }
			else
			  {
				$current_sequence .= $line;
			  }
		  }
		close $ALN_FH or croak "Couldn't close '$current_alignment': $!\n";
	  }
	print "\twriting fasta files";
	mkdir("species") if !-e "species";
	foreach my $species ( keys(%species_hash) )
	  {

		#next if !exists $convert_hash{$species};
		my $species_file = "./species/$species.fa";
		print "\twriting $species\n";
		### OPENING FILE
		open my $out, '>', $species_file
		  or croak "Couldn't open '$species_file': $!\n";
		foreach my $id ( keys( %{ $species_hash{$species} } ) )
		  {
			print {$out} $species_hash{$species}{$id}
			  or croak "Couldn't write '$species_file': $!\n";
		  }
		### CLOSING FILE
		close $out or croak "Couldn't close '$species_file': $!\n";
	  }
	exit;
  }

sub treefam
  {

	# required input: alignments, trees
	####
	# TREEFAM
	####
	my %species_hash = ();

#my $input_alignments_dir = "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/treefam/families_work/alignments";
#my $input_trees_dir = "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/treefam/families_work/trees";
	my $input_alignments_dir =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/treefam/families_work/";
	my $input_trees_dir =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/treefam/families_work/";
	my $all_treefam_trees =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/treefam/trees.txt.table";
	my $all_aaSeqs =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/treefam/aa_seq.txt.table";
	my $output_species_dir =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/treefam/small_species";
	my $output_alignment_dir =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/treefam/alignments";
	my $output_trees_dir = "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/treefam/trees";
	mkdir($output_species_dir)   if !-e $output_species_dir;
	mkdir($output_alignment_dir) if !-e $output_alignment_dir;
	mkdir($output_trees_dir)     if !-e $output_trees_dir;
	my %ids_from_clean_trees         = ();
	my %treefamID2transcriptIds_href = ();
	my %saveSpeciesSequences;
	my %ID2species_mapping = ();
	my %ids2aaSeq_hash     = ();

	# FIRST GET ALL IDS FROM CURATED AND CLEANED TREES
	#my @clean_tree_files = <$input_trees_dir/0*/*/clean.nhx>;
	readTreeFamSeedIDs(
						{
						  all_treefam_trees            => $all_treefam_trees,
						  ids_from_clean_trees_href    => \%ids_from_clean_trees,
						  allowedTreeFamGroups_href    => \%allowedTreeFamGroups,
						  treefamID2transcriptIds_href => \%treefamID2transcriptIds_href,
						  ID2species_mapping_href      => \%ID2species_mapping
						}
	);

	#print Dumper $ids_from_clean_trees{'BGIOSIBCA020079_ORYSA'};
	readTreeFamAASeqs(
					   {
						 all_aaSeqs                => $all_aaSeqs,
						 ids2aaSeq_href            => \%ids2aaSeq_hash,
						 allowedTreeFamGroups_href => \%allowedTreeFamGroups
					   }
	);

	#exit;
	foreach my $current_id ( sort keys %ids_from_clean_trees )
	  {
		print "\tsearching for $current_id ...\n";
		if ( exists $ids2aaSeq_hash{$current_id} )
		  {
			print "\tseq found\n";
			my $species = $ID2species_mapping{$current_id};
			$saveSpeciesSequences{$species} .=
			  ">" . $current_id . "_" . $species . "\n" . $ids2aaSeq_hash{$current_id} . "\n";
		  }
		else
		  {
			if ( exists $treefamID2transcriptIds_href{$current_id} )
			  {
				my $transcriptID = $treefamID2transcriptIds_href{$current_id};
				print "\t\tTrying $transcriptID\n";
				if ( exists $ids2aaSeq_hash{$transcriptID} )
				  {
					print "\t\t\tseq found\n";
					print "\t\t\tbelongs to " . $ID2species_mapping{$current_id} . "\n";
					my $species = $ID2species_mapping{$current_id};
					$saveSpeciesSequences{$species} .=
					  ">" . $current_id . "_" . $species . "\n" . $ids2aaSeq_hash{$transcriptID} . "\n";
					next;
				  }
				else
				  {
					print "\t\t\tno seq found\n";
					next;
				  }
			  }
			else
			  {
				print "\t\tno transcript -> no seq\n";
				next;
			  }
			print "\tno seq found\n";
		  }
	  }
	print "found species: " . join( ",", keys(%saveSpeciesSequences) ) . "\n";
	foreach my $curr_species ( keys(%saveSpeciesSequences) )
	  {
		my $longSpeciesName = $convert_hash{$curr_species};
		print "\tsaving $curr_species ($longSpeciesName)\n";
		my $speciesOutputFile = $output_species_dir . "/" . $longSpeciesName . ".fa";
		open my $speciesOutputFile_FH, '>', $speciesOutputFile
		  or croak "Couldn't open '$speciesOutputFile': ";

		#next if /^$/;
		#chomp;
		print {$speciesOutputFile_FH} $saveSpeciesSequences{$curr_species} . "\n"
		  or croak "Couldn't write '$speciesOutputFile': ";
		close $speciesOutputFile_FH
		  or croak "Couldn't close '$speciesOutputFile': ";
	  }
	exit;
  }

sub readTreeFamSeedIDs
  {
	#### PARAMETER VARIABLES
	my ($arg_ref)                    = @_;
	my $all_treefam_trees            = $arg_ref->{all_treefam_trees};
	my $ids_from_trees_href          = $arg_ref->{ids_from_clean_trees_href};
	my $treefamID2transcriptIds_href = $arg_ref->{treefamID2transcriptIds_href};
	my $allowedTreeFamGroups_href    = $arg_ref->{allowedTreeFamGroups_href};
	my $ID2species_mapping_href      = $arg_ref->{ID2species_mapping_href};
	print "\treading tree file ($all_treefam_trees)\n";
	open my $all_treefam_trees_FH, '<', $all_treefam_trees
	  or croak "Couldn't open '$all_treefam_trees': ";

	while ( defined( my $line = <$all_treefam_trees_FH> ) )
	  {
		my $mapped       = 0;
		my %ids_for_tree = ();

		# only take seed
		next if $line =~ /^(TF\d+)\s+[FULL|CLEAN]/;
		$line =~ /^(TF\d+)\s+/;
		my $treefamFamily = $1;
		next if !exists $allowedTreeFamGroups_href->{$treefamFamily};
		my $treefam_tree_id = substr( $line, 0, 8 );

		#print "checking $treefam_tree_id\n";
		#$no_trees++;
		#                  $continue = 0 if $treefam_tree_id eq 'TF101002';
		#     $continue = 0 if $no_trees > 99;
		my @tree_nodes = split( /,|\)/, $line );
		if ( !@tree_nodes )
		  {
			die "problem parsing $treefam_tree_id\n";
		  }

		#                last;
		foreach my $current_tree_node (@tree_nodes)
		  {
			my @treefamId4transcript_array = ();
			my @ids_to_search              = ();
			if ( $current_tree_node =~ /^:/ )
			  {

				#print "SKIPPED: ".$current_tree_node."\n";
				next;
			  }
			if ( !( $current_tree_node =~ /\w:\d+.*\[.*\]/ ) )
			  {

				#print "SKIPPED: ".$current_tree_node."\n";
				next;
			  }

			#print "LINE: ".$current_tree_node."\n";
			#print "\tanalyse!!!!\n";
			my ( $species, $treefam_id, $treefam_id_gene, $treefam_id_transcript );
			$current_tree_node =~ s/TF\d{6}\s+SEED\s+//;
			$current_tree_node =~ /(.*)\[(.*)\]/;
			my ( $id_part, $desc_part ) = ( $1, $2 );

			# CLEAN NEWICK CHARACTERS '(',')'
			$id_part   =~ s/\(|\)//g;
			$id_part   =~ s/\s*//g;
			$desc_part =~ s/\(|\)//g;
			$desc_part =~ s/\s*//g;

			# DISPLAY ID
			my @splitted_id_part = split( /:/, $id_part );
			$treefam_id = $splitted_id_part[0];
			unless ( $treefam_id =~ /_/ )
			  {

				#print "treefam id in wrong format ($treefam_id). Contains no '_'\n";
				#print "SKIPPED: ".$current_tree_node."\n";
				next;
			  }

			#print "\t\tTreeFam Id: $treefam_id\n";
			$ids_from_trees_href->{$treefam_id} = 1;
			$species = ( split( /_/, $treefam_id ) )[1];
			if ( !defined $species )
			  {
				print "\t no species defined: " . $treefam_id . "\n";
				exit;
			  }
			$ID2species_mapping_href->{$treefam_id} = $species;
			$treefam_id = ( split( /_/, $treefam_id ) )[0];
			$ids_from_trees_href->{$treefam_id}     = 1;
			$ID2species_mapping_href->{$treefam_id} = $species;

			# DESCRIPTION
			my @splitted_desc_part = split( /:/, $desc_part );
			foreach (@splitted_desc_part)
			  {

				#              print "$_";
				$species               = $1 if /S=(\w+)/;    # species
				$treefam_id_gene       = $1 if /G=(.*)/;     #gene id
				$treefam_id_transcript = $1 if /O=(.*)/;     #transcript id
			  }
			if ( substr( $treefam_id_transcript, -2, 2 ) =~ /\.\d/ )
			  {
				print "\tcut transcript\n";
				$treefam_id_transcript =~ s/\.\d$//;
				print "\tsaving $treefam_id_transcript\t...";

				#exit;
			  }

			#print "\t\tTranscript: $treefam_id_transcript\n";
			$treefamID2transcriptIds_href->{$treefam_id} = $treefam_id_transcript;
			$ID2species_mapping_href->{$treefam_id}      = $species;
			print "species: $species\ttranscript: $treefam_id_transcript\ttreefam_id: $treefam_id\n";

			#exit;
			#$transcriptIds_from_trees_href->{$treefam_id_transcript} = 1;
		  }

		#exit;
	  }
	close $all_treefam_trees_FH
	  or croak "Couldn't close '$all_treefam_trees': ";
	if ( !keys(%$ids_from_trees_href) )
	  {
		print "\tCould not read IDs from trees from file $all_treefam_trees\n";
		exit;
	  }
	return 1;
  }

sub readTreeFamAASeqs
  {
	#### PARAMETER VARIABLES
	my ($arg_ref)                 = @_;
	my $all_aaSeqs                = $arg_ref->{all_aaSeqs};
	my $ids2aaSeq_href            = $arg_ref->{ids2aaSeq_href};
	my $allowedTreeFamGroups_href = $arg_ref->{allowedTreeFamGroups_href};
	print "\treading aa_seq file ($all_aaSeqs)\n";
	open my $all_aaSeqs_FH, '<', $all_aaSeqs
	  or croak "Couldn't open '$all_aaSeqs': ";
	while ( my $line = <$all_aaSeqs_FH> )
	  {
		my ( $ID, $seq_length, $sequence ) = split( /\t/, $line );
		if ( !defined $ID || !defined $seq_length || !defined $sequence )
		  {
			print "\terror in line: $line ($ID - $seq_length - $sequence)\n";
			exit;
		  }
		print "\tsaving $ID\t...";
		if ( substr( $ID, -2, 2 ) =~ /\.\d/ )
		  {
			print "\tfound a bastard\n";
			my @splitted_ID = split( /\./, $ID );
			$ID =~ s/\.\d$//;
			print "\tsaving $ID\t...";

			#exit;
		  }
		print "\t$ID\n";

		#exit;
		$ids2aaSeq_href->{$ID} = $sequence;
	  }
	close $all_aaSeqs_FH or croak "Couldn't close '$all_aaSeqs': ";
	print "\tread sequences for " . keys(%$ids2aaSeq_href) . " IDs\n";
	if ( !keys(%$ids2aaSeq_href) )
	  {
		print "\tCould not read IDs from trees from file $all_aaSeqs\n";
		exit;
	  }
	return 1;
  }

sub testing
  {
	##### TESTING
	# Iterate over all hieranoid alignments
	my $hieranoid_alignments_to_check =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/nodes/Eukaryota/alignments/";
	my $hieranoid_groupfile =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/results_small/nodes/Bilateria/Bilateria.groups.txt";
	my %id_to_cluster_hieranoid = ();

	#		foreach my $curr_hieranoid_alignment (glob("$hieranoid_alignments_to_check/*")){
	#			print "\tlooking at $curr_hieranoid_alignment\n";
	#			my $ali_name  = basename($curr_hieranoid_alignment);
	#			$ali_name =~ s/\.fa//;
	#			my @headers = `grep ">" $curr_hieranoid_alignment`;
	#			foreach(@headers){
	#				chomp;
	#				/>(.*)/;
	#				$id_to_cluster_hieranoid{$1} = $ali_name;
	#			#	print "\t add $1 to $ali_name\n";
	#			}
	#exit;
	#	}
	my $all_treefam_trees =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/treefam/trees.txt.table";
	print "\treading tree file\n";
	my %ids_from_clean_trees = ();
	open my $all_treefam_trees_FH, '<', $all_treefam_trees
	  or croak "Couldn't open '$all_treefam_trees': ";
	while ( my $line = <$all_treefam_trees_FH> )
	  {
		if ( $line =~ /(TF\d+)\s+SEED\s+(.*)/ )
		  {
			my ( $treefamFamily, $newickTree ) = ( $1, $2 );
			if ( !defined $treefamFamily || !defined $newickTree )
			  {
				print "\terror in line: $line ($treefamFamily - $newickTree)\n";
				exit;
			  }
			next if !exists $allowedTreeFamGroups{$treefamFamily};
			print "\tfound $treefamFamily\n";
			my $input_tree = "test";
			write_to_file( { file_name => $input_tree, text => $newickTree } );
			if ( !-e $input_tree || !-s $input_tree )
			  {
				print "tree file $input_tree does not exists\n";
				exit;
			  }
			my $treeio = Bio::TreeIO->new( -format => 'newick', -file => $input_tree );
			my $tree = $treeio->next_tree;

			# my @nodes  = $tree->get_nodes;
			my @leaves = $tree->get_leaf_nodes;

			# my $root   = $tree->get_root_node;
			foreach (@leaves)
			  {
				print $_->id() . " , ";
				$ids_from_clean_trees{ $_->id() } = 1;
			  }
		  }
	  }
	close $all_treefam_trees_FH
	  or croak "Couldn't close '$all_treefam_trees': ";
	open my $OG_FH, '<', $hieranoid_groupfile
	  or croak "Couldn't open " . $hieranoid_groupfile . ": $OS_ERROR";
	my $counter = 0;
	### READING FROM FILE
	while ( my $line = <$OG_FH> )
	  {
		next if $line =~ /^$/;
		$line =~ /(\w+):(.*)/;
		my ( $og, $members ) = ( $1, $2 );
		if ( !defined $og || !defined $members )
		  {
			print "Could not parse line $line in file " . $hieranoid_groupfile . ". Values are $og,$members\n";
			exit;
		  }
		foreach ( split( ",", $members ) )
		  {
			$id_to_cluster_hieranoid{$_} = $og;
		  }
	  }

	#if(defined $species){
	#        print("\tread $og_counter groups for $species") ;
	#DEBUG("\tread $og_counter groups for $species") ;
	#}
	#else{
	#  print("\tread $og_counter groups") ;
	#DEBUG("\tread $og_counter groups") ;
	#}
	### CLOSING FILE
	#print "\t\tRead og assignment file\n";
	#print Dumper %id_to_cluster_hieranoid;
	#exit;
	print "found " . keys(%id_to_cluster_hieranoid) . " ids\n";
	close $OG_FH
	  or croak "Couldn't close '" . $hieranoid_groupfile . "': $OS_ERROR";

	#		exit;
	my $alignment_counter = 0;

	#exit;
	# save ids to hash: seq_id --> alignment number
	# Iterate over all test alignments
	my $treefam_alignments_to_check =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/treefam/families_work/alignments/";
	my $input_alignments_dir =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/treefam/families_work/";
	my ( $false_negative, $false_negative_total ) = 0;
	my ( $true_positive,  $true_positives_total ) = 0;
	my @files = (
				  "$input_alignments_dir/07/TF300407/full.mfa", "$input_alignments_dir/01/TF101001/full.mfa",
				  "$input_alignments_dir/04/TF101004/full.mfa", "$input_alignments_dir/10/TF101010/full.mfa",
				  "$input_alignments_dir/12/TF101012/full.mfa"
	);

	#foreach my $curr_treefam_alignment (glob("$treefam_alignments_to_check/*/*/seed.mfa")){
	foreach my $curr_treefam_alignment (@files)
	  {
		next if !-s $curr_treefam_alignment;
		print "\tlooking at $curr_treefam_alignment\n";
		$false_negative = 0;
		$true_positive  = 0;
		my %id_to_cluster_treefam = ();
		my $ali_name              = basename($curr_treefam_alignment);
		$ali_name =~ s/\.mfa//;
		my @headers = `grep ">" $curr_treefam_alignment`;

		foreach (@headers)
		  {
			my @split_array = split( /\s/, $_ );
			my $header_to_check = $split_array[0];
			$header_to_check =~ s/>//;
			print "Header to check $header_to_check\n";
			## ADD PREDICTION TO HASH
			if ( !exists $id_to_cluster_hieranoid{$header_to_check} )
			  {

				#	print "warn: $header_to_check not found\n";
				$false_negative++;
				next;
			  }
			else
			  {
				print "\thas " . $id_to_cluster_hieranoid{$header_to_check} . "\n";
				$id_to_cluster_treefam{ $id_to_cluster_hieranoid{$header_to_check} }++;
				print "$header_to_check has " . $id_to_cluster_hieranoid{$header_to_check} . "\n";
				$true_positive++;
			  }
		  }
		## EVALUATE
		if ( keys(%id_to_cluster_treefam) == 0 )
		  {
			print "\t\tNo Hieranoid predictions\n";
			next;
		  }
		if ( keys(%id_to_cluster_treefam) == 1 )
		  {
			foreach ( keys(%id_to_cluster_treefam) )
			  {
				print "\t\tTP: $id_to_cluster_treefam{$_} ($_)\n";
			  }
		  }
		else
		  {
			print "\t\tMultiple Hieranoid predictions: ";

			#print "".join(",",keys(%id_to_cluster_treefam))."\n";
			foreach ( keys(%id_to_cluster_treefam) )
			  {
				print "$_ : $id_to_cluster_treefam{$_},";
			  }
			print "\n";
		  }

		#exit;
		if ( $true_positives_total == 0 )
		  {
			print "\t\tPossible a strange group\n";
		  }
		last if $alignment_counter++ > 40;
		$false_negative_total += $false_negative;
		$true_positives_total += $true_positive;
		print
"\t\ttrue positives: $true_positive ($true_positives_total)| false negatives: $false_negative ($false_negative_total)\n";

		#exit;
	  }

	# iterate over all seq_ids
	# lookup id in hash:
	# if not there: TN
	# if there:
	# after last seq_id to check --> all seq_ids in same?
	exit;
  }

sub mapOMA2RefOGs
  {

# Version nov 2010
#my $omaGroupFile = "/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/oma-groups_nov2010.txt";
#my $omaGroupFileConverted = "/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/oma-groups_nov2010.txt_converted";
	my $omaGroupFile =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/oma-groups_may2011.txt";
	my $omaGroupFileConverted =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/oma-groups_may2011_converted.txt";
	my $OMA2ENSEBMLFIle =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/oma-ensembl.txt";
	my $orthobench_fasta_files =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/predictions/orthobench/";
	my %allowed_species_hash = (
								 HUMAN => 1,
								 MONDO => 1,
								 CANFA => 1,
								 RATNO => 1,
								 MOUSE => 1,
								 PANTR => 1,
								 CHICK => 1,
								 TETNG => 1,
								 DANRE => 1,
								 CIOIN => 1,
								 DROME => 1,
								 CAEEL => 1,
	);

	# Read allowed IDs from OrthoBench
	my %ID2refOG;
	my %refOG_hash;
	my %refOGGenes_hash;

	#  'refOG0001' -> @(ENSTNIP00000005325,ENSDARP00000071531)
	# 2. read refOG predictions
	# Foreach file:
	#               grep fasta headers
	my $lineCounter = 0;
	my $no_saved    = 0;
	my $orthobench_groups =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/IDs_OrthoBench2.txt";
	getOrthoBenchGroups( $orthobench_groups, \%ID2refOG, \%refOG_hash, \%refOGGenes_hash );

	#exit;
	if ( !keys %ID2refOG )
	  {
		print "\tCould not get hash \$ID2refOG. Exiting\n";
		exit;
	  }

	#print "\t\tRead ".keys(%ID2refOG)." IDs ($no_saved)\n";
	#exit;
	my %oma2ensembl_mapping;
	open my $OMA2ENSEBMLFIle_FH, '<', $OMA2ENSEBMLFIle
	  or croak "Couldn't open '$OMA2ENSEBMLFIle': ";
	while (<$OMA2ENSEBMLFIle_FH>)
	  {
		next if /^#/;
		my ( $omaID, $ensemblID ) = split;

		#print " $omaID,$ensemblID\n";
		if ( substr( $ensemblID, 6, 1 ) eq 'P' )
		  {
			$oma2ensembl_mapping{$omaID} = $ensemblID;

			#print " take $ensemblID\n";
		  }
		else
		  {
			next;
		  }

		#exit;
		#print Dumper %oma2ensembl_mapping;
		#exit;
	  }
	close $OMA2ENSEBMLFIle_FH or croak "Couldn't close '$OMA2ENSEBMLFIle': ";
	print "\tFinished reading " . keys(%oma2ensembl_mapping) . " keys\n";
	if ( !keys(%oma2ensembl_mapping) )
	  {
		print "\tno keys from $OMA2ENSEBMLFIle. Exiting\n";
		exit;
	  }
	my $groupCounter = 0;
	open my $omaGroupFile_FH, '<', $omaGroupFile
	  or croak "Couldn't open '$OMA2ENSEBMLFIle': ";
	open my $omaGroupFileConverted_FH, '>', $omaGroupFileConverted
	  or croak "Couldn't open '$omaGroupFileConverted': ";
	while (<$omaGroupFile_FH>)
	  {
		next if /^#/;

		#print $_."\n";
		/^(\d+)\t(\w+)\t(.*)/;
		my ( $omaID, $fingerprint, $groupMember ) = ( $1, $2, $3 );
		my @groupMembers = split "\t", $groupMember;

		#print "\tmembers: ".join(",",@groupMembers)."\n";
		my $foundUsefulID = 0;
		my @group_string_members;
		my $group_string;
		foreach my $groupGene (@groupMembers)
		  {

			# allowed species
			my $species = substr( $groupGene, 0, 5 );

			#print "\t$groupGene --> species is $species \n";
			#exit;
			next if !exists $allowed_species_hash{$species};
			if ( exists $oma2ensembl_mapping{$groupGene} )
			  {

				#print "\tadding $oma2ensembl_mapping{$groupGene}";
				my $mapped_ID = $oma2ensembl_mapping{$groupGene};

				#if(exists $ID2refOG{$mapped_ID}){
				push( @group_string_members, $mapped_ID );

				#$group_string .= $mapped_ID.",";
				$foundUsefulID = 1;

				#}
			  }
			else
			  {
				if ( exists $ID2refOG{$groupGene} )
				  {
					push( @group_string_members, $groupGene );

					#$group_string .= $groupGene.",";
					$foundUsefulID = 1;
				  }
			  }
		  }
		$group_string = "$omaID: " . join( ",", @group_string_members ) . "\n";
		if ($foundUsefulID)
		  {
			print {$omaGroupFileConverted_FH} $group_string
			  or croak "Couldn't write '$omaGroupFileConverted': ";

			#last if $groupCounter++ > 10;
		  }

		#print "$group_string\n";
		#exit;
	  }
	close $omaGroupFile_FH or croak "Couldn't close '$omaGroupFile': ";
	close $omaGroupFileConverted_FH
	  or croak "Couldn't close '$omaGroupFileConverted': ";
	exit;
  }

sub reformat_orthology_predictions
  {
	my $orthologs_table   = "treefam/ortholog.txt.table";
	my $genes_table       = "treefam/genes.txt.table";
	my $new_mapping       = "treefam/new_mapping.txt";
	my %id2geneID_mapping = ();
	### OPENING FILE
	open my $ortholog_FH, '<', $orthologs_table
	  or croak "Couldn't open '$orthologs_table': $OS_ERROR";
	open my $genes_FH, '<', $genes_table
	  or croak "Couldn't open '$genes_table': $OS_ERROR";
	open my $new_mapping_FH, '>', $new_mapping
	  or croak "Couldn't open '$new_mapping': $OS_ERROR";
	### READING FROM FILE
	while ( my $line = <$genes_FH> )
	  {

		#	print $line."\n";
		my @splitted_line = split( /\s/, $line );

		#	print Dumper @splitted_line;
		#	exit;
		my ( $identifier, $gene_id ) = ( $splitted_line[0], $splitted_line[7] );
		if ( !defined $identifier || !defined $gene_id )
		  {
			print "Problem with catching from $line ($identifier, $gene_id)\n";
			exit;
		  }
		$id2geneID_mapping{$identifier} = $gene_id;

		#	exit;
	  }
	if ( !keys(%id2geneID_mapping) )
	  {
		print "Could not read entries from $genes_table\n";
		exit;
	  }
	print "\tread " . keys(%id2geneID_mapping) . " entries\n";
	while ( my $line = <$ortholog_FH> )
	  {
		$line =~ /^(\d+)\s*(\d+)/;
		my ( $geneID_A, $geneID_B ) = ( $1, $2 );
		if ( !defined $geneID_A || !defined $geneID_B )
		  {
			print "Problem with catching from $line ($geneID_A, $geneID_B)\n";
			exit;
		  }
		if ( $geneID_B eq "1223914" || $geneID_A eq "1223914" )
		  {
			print "\t$geneID_A\t$geneID_B\n"
			  . $id2geneID_mapping{$geneID_A} . " "
			  . $id2geneID_mapping{$geneID_B} . "\n";

			#		next;
		  }

		#	next;
		if (    !exists $id2geneID_mapping{$geneID_A}
			 || !exists $id2geneID_mapping{$geneID_B} )
		  {
			print "\tProblem mapping for $line ($geneID_A,$geneID_B)\n";
			exit;
		  }
		print {$new_mapping_FH} $id2geneID_mapping{$geneID_A} . " " . $id2geneID_mapping{$geneID_B} . "\n"
		  or croak "Couldn't write '$new_mapping': $OS_ERROR";
	  }
	### CLOSING FILE
	close $ortholog_FH    or croak "Couldn't close '$orthologs_table': $OS_ERROR";
	close $genes_FH       or croak "Couldn't close '$genes_table': $OS_ERROR";
	close $new_mapping_FH or croak "Couldn't close '$new_mapping': $OS_ERROR";
	exit;
  }
#### ORTHOBENCH
# Read files
sub getPredictionMethodPredictions
  {
	my ( $predictionsFile, $og_assignments_hashref, $ID2Group_href ) = (@_);
	my $og_counter = 0;
	if ( !-e $predictionsFile || !-s $predictionsFile )
	  {
		print( "\tCould not read resolved mappings " . $predictionsFile . "\n" );
		return 0;
	  }

	#DEBUG("\tget_group2expandedIDsHash: reading from ".$predictionsFile."\n");
	## CHECK DEFINEDNESS
	# croak q{One/several of the parameters for 'read_og_assignments' was/were not defined.\n}
	#if any {!defined $_} $predictionsFile;
	### OPENING FILE
	open my $OG_FH, '<', $predictionsFile
	  or croak "Couldn't open " . $predictionsFile . ": $OS_ERROR";
	my $counter = 0;
	### READING FROM FILE
	while ( my $line = <$OG_FH> )
	  {
		next if $line =~ /^$/;
		chomp($line);
		$line =~ /(\w+):\s*(.*),?/;
		my ( $og, $members ) = ( $1, $2 );
		if ( !defined $og || !defined $members )
		  {
			print( "Could not parse line $line in file " . $predictionsFile . ". Values are $og,$members\n" );
			exit;
		  }
		$og_assignments_hashref->{$og} = $members;
		foreach my $id ( split( ",", $members ) )
		  {
			$id =~ s/\s+//g;
			$ID2Group_href->{$id} = $og;
		  }
		$og_counter++;
	  }
	close $OG_FH
	  or croak "Couldn't close '" . $predictionsFile . "': $OS_ERROR";
	if ( !keys( %{$og_assignments_hashref} ) )
	  {
		print( "Could not get groupPredictions from " . $predictionsFile . " \n" );
		exit;
	  }
  }

sub getPredictionMethodPredictionsPairwise
  {
	my ($arg_ref)              = @_;
	my $predictionsFile        = $arg_ref->{'predictionFile'};
	my $og_assignments_hashref = $arg_ref->{'AllPredictionsHash'};
	my $ID2Group_href          = $arg_ref->{'ID2DatabaseOG'};
	my $ID2Gene_href           = $arg_ref->{'ID2GeneHash'};
	my $og_counter             = 0;
	if ( !-e $predictionsFile || !-s $predictionsFile )
	  {
		print( "\tCould not read resolved mappings " . $predictionsFile . "\n" );
		return 0;
	  }
	my %countMissingSpeciesIDs;
	my %missingIDs;

	#DEBUG("\tget_group2expandedIDsHash: reading from ".$predictionsFile."\n");
	## CHECK DEFINEDNESS
	# croak q{One/several of the parameters for 'read_og_assignments' was/were not defined.\n}
	#if any {!defined $_} $predictionsFile;
	### OPENING FILE
	open my $OG_FH, '<', $predictionsFile
	  or croak "Couldn't open " . $predictionsFile . ": $OS_ERROR";
	my $counter = 0;
	my %originalIDs;
	### READING FROM FILE
	while ( my $line = <$OG_FH> )
	  {
		next if $line =~ /^$/;
		chomp($line);

		#print $line. "\n";
		$line =~ /(\w+):\s*(.*),?/;
		my ( $og, $members ) = ( $1, $2 );
		if ( !defined $og || !defined $members )
		  {
			print( "Could not parse line $line in file " . $predictionsFile . ". Values are $og,$members\n" );
			exit;
		  }

		#$og_assignments_hashref->{$og} = $members;
		my @convertedMembers;
		my $no_mapped_member = 0;
		my %avoid_duplicates;
		foreach my $id ( split( ",", $members ) )
		  {
			$id =~ s/\s+//g;

			#print "\t\t$id\n";
			if ( exists $ID2Gene_href->{$id} )
			  {
				my $ID2RefOGGroup = $ID2Gene_href->{$id};
				if ( !exists $avoid_duplicates{$ID2RefOGGroup} )
				  {
					push( @convertedMembers, $ID2RefOGGroup );
				  }
				$avoid_duplicates{$ID2RefOGGroup} = 1;
				$ID2Group_href->{$ID2RefOGGroup} = $og;

				#print "\tmap to " . $ID2Gene_href->{$id} . "\n";
				$no_mapped_member++;
			  }
			else
			  {

				#print "\t\t\tnot present\n";
				my $Species4ID = substr( $id, 0, 7 );
				$Species4ID =~ s/\d+//g;
				$countMissingSpeciesIDs{$Species4ID}++;
				$missingIDs{$og}{$id} = 1;
				push( @convertedMembers, $id );
			  }
		  }
		$originalIDs{$og} = $members;
		if ( $no_mapped_member > 1 && scalar(@convertedMembers) )
		  {

			#print "mapped ".join(",",@convertedMembers)."\n";
			#exit;
			$og_assignments_hashref->{$og} = join( ",", @convertedMembers );
		  }
		$og_counter++;
	  }
	my $printLine;
	foreach my $og ( keys( %{$og_assignments_hashref} ) )
	  {
		$printLine .= "$og (org): $originalIDs{$og}\n";
		$printLine .= "$og:" . $og_assignments_hashref->{$og} . "\n";
		$printLine .=
		  "\tmissing (" . keys( %{ $missingIDs{$og} } ) . "):" . join( ",", keys( %{ $missingIDs{$og} } ) ) . "\n";
	  }
	write_to_file(
				   {
					 file_name => basename($predictionsFile) . "_converted",
					 text      => $printLine
				   }
	);

	#print Dumper %countMissingSpeciesIDs;
	#print "".join("\n\t",keys(%missingIDs))."\n";
	#exit;
	close $OG_FH
	  or croak "Couldn't close '" . $predictionsFile . "': $OS_ERROR";
	if ( !keys( %{$og_assignments_hashref} ) )
	  {
		print( "Could not get groupPredictions from " . $predictionsFile . " \n" );
		exit;
	  }
  }

sub getOrthoBenchGroups
  {
	my ( $predictionsFile, $ID2refOG_href, $refOG_href, $refOGGenes_href ) = (@_);
	if ( !-e $predictionsFile || !-s $predictionsFile )
	  {
		print( "\tCould not read resolved mappings " . $predictionsFile . "\n" );
		return 0;
	  }
### OPENING FILE
	open my $OG_FH, '<', $predictionsFile
	  or croak "Couldn't open " . $predictionsFile . ": $OS_ERROR";
	my $counter           = 0;
	my $gene_counter      = 0;
	my $previousID        = -1;
	my $gene4groupCounter = 0;
	### READING FROM FILE
	while ( my $line = <$OG_FH> )
	  {
		next if $line =~ /^$/;
		chomp($line);
		$line =~ /(\w+):\s*(.*),?/;
		my ( $og, $members ) = ( $1, $2 );
		if ( !defined $og || !defined $members )
		  {
			print( "Could not parse line $line in file " . $predictionsFile . ". Values are $og,$members\n" );
			exit;
		  }
		$previousID = $og if $previousID eq '-1';

		#print "$previousID ne $og\n";
		if ( $previousID ne $og )
		  {
			$gene4groupCounter = 0;
			$previousID        = $og;

			#print "\t\tset geneCounter 0\n";
		  }

		#$og_assignments_hashref->{$og} = $members;
		my @IDarray = split( /\s/, $members );
		foreach my $id (@IDarray)
		  {
			$id =~ s/\s+//g;
			$ID2refOG_href->{$id}                       = $og;
			$refOG_href->{$og}{$gene4groupCounter}{$id} = 1;
			$refOGGenes_href->{$og}{$id}                = $gene4groupCounter;
		  }
		$gene4groupCounter++;
		$gene_counter++;
	  }

	#print Dumper $refOG_href;
	#print Dumper $ID2refOG_href;
	#print "\tfound $gene_counter genes\n";
	#exit;
	close $OG_FH
	  or croak "Couldn't close '" . $predictionsFile . "': $OS_ERROR";
	if ( !keys( %{$ID2refOG_href} ) )
	  {
		print( "Could not get groupPredictions from " . $predictionsFile . " \n" );
		exit;
	  }
  }

sub getOrthoBenchGroupsPairwise
  {
	my ($arg_ref)           = @_;
	my $predictionsFile     = $arg_ref->{'orthobench_groups_file'};
	my $ID2refOG_href       = $arg_ref->{'ID2refOG_href'};
	my $refOG_href          = $arg_ref->{'refOG_href'};
	my $refOGGenes_href     = $arg_ref->{'refOGGenes_href'};
	my $refOG2Members_href  = $arg_ref->{'refOG2Members_href'};
	my $pair2RefOGGene_href = $arg_ref->{'pair2RefOGGene_href'};
	my $ID2Gene_href        = $arg_ref->{'ID2Gene_href'};
	my $og_counter          = 0;
	if ( !-e $predictionsFile || !-s $predictionsFile )
	  {
		print( "\tCould not read resolved mappings " . $predictionsFile . "\n" );
		return 0;
	  }

	#DEBUG("\tget_group2expandedIDsHash: reading from ".$predictionsFile."\n");
	## CHECK DEFINEDNESS
	# croak q{One/several of the parameters for 'read_og_assignments' was/were not defined.\n}
	#if any {!defined $_} $predictionsFile;
	### OPENING FILE
	open my $OG_FH, '<', $predictionsFile
	  or croak "Couldn't open " . $predictionsFile . ": $OS_ERROR";
	my $counter      = 0;
	my $gene_counter = 0;
	my %refOG_array;
	my $previousID        = -1;
	my $gene4groupCounter = 0;
	### READING FROM FILE
	while ( my $line = <$OG_FH> )
	  {
		next if $line =~ /^$/;
		chomp($line);
		$line =~ /(\w+):\s*(.*),?/;
		my ( $og, $members ) = ( $1, $2 );
		if ( !defined $og || !defined $members )
		  {
			print( "Could not parse line $line in file " . $predictionsFile . ". Values are $og,$members\n" );
			exit;
		  }
		$previousID = $og if $previousID eq '-1';
		if ( $previousID ne $og )
		  {
			$gene4groupCounter = 0;
			$previousID        = $og;
		  }
		my @IDarray = split( /\s/, $members );
		foreach my $id (@IDarray)
		  {
			$id =~ s/\s+//g;
			$ID2refOG_href->{$id}                       = $og;
			$refOG_href->{$og}{$gene4groupCounter}{$id} = 1;
			$refOGGenes_href->{$og}{$id}                = $gene4groupCounter;
			$ID2Gene_href->{$id}                        = $og . "_" . $gene4groupCounter;
		  }
		$ID2Gene_href->{$og} = $og . "_" . $gene4groupCounter;
		push( @{ $refOG2Members_href->{$og} }, $og . "_" . $gene4groupCounter );
		$gene4groupCounter++;
		$gene_counter++;
	  }
	close $OG_FH
	  or croak "Couldn't close '" . $predictionsFile . "': $OS_ERROR";
	if ( !keys( %{$refOG2Members_href} ) )
	  {
		print( "Could not get groupPredictions from " . $predictionsFile . " \n" );
		exit;
	  }
  }

sub geteggNOGPredictions
  {
	my ( $predictionsFile, $og_assignments_hashref, $ID2Group_href ) = (@_);
	my $og_counter = 0;
	if ( !-e $predictionsFile || !-s $predictionsFile )
	  {
		print( "\tCould not read resolved mappings " . $predictionsFile . "\n" );
		return 0;
	  }
	my %allowed_species = ();
	my @species =
	  ( "6239", "7227", "7719", "9615", "9606", "9598", "7955", "99883", "13616", "10090", "10116", "9031" );
	foreach my $item (@species) { $allowed_species{$item} = 1; }

	#DEBUG("\tget_group2expandedIDsHash: reading from ".$predictionsFile."\n");
	## CHECK DEFINEDNESS
	# croak q{One/several of the parameters for 'read_og_assignments' was/were not defined.\n}
	#if any {!defined $_} $predictionsFile;
	### OPENING FILE
	open my $OG_FH, '<', $predictionsFile
	  or croak "Couldn't open " . $predictionsFile . ": $OS_ERROR";
	my $counter = 0;
	### READING FROM FILE
	while ( my $line = <$OG_FH> )
	  {
		next if $line =~ /^$/;
		next if $line =~ /^#/;
		chomp($line);
		my @splitted_line = split( /\t/, $line );
		my $member = $splitted_line[0];
		$member =~ /(\d+)\.(.*)/;
		my $species = $1;

		#print "Species is $species\n";
		next if !exists $allowed_species{$species};
		$member = $2;
		my $og = $splitted_line[3];
		if ( !defined $og || $og eq '' || !defined $member || $member eq '' )
		  {
			print( "Could not parse line $line in file " . $predictionsFile . ". Values are $og,$member\n" );
			exit;
		  }
		if ( !$og =~ /meNOG\d{5}/ )
		  {
			print "\tog $og in wrong format\n";
		  }

		#print "\tsaving $og and $member\n";
		$og_assignments_hashref->{$og} .= $member . ",";
		$ID2Group_href->{$member} = $og;
		$og_counter++;
	  }
	close $OG_FH
	  or croak "Couldn't close '" . $predictionsFile . "': $OS_ERROR";
	if ( !keys( %{$og_assignments_hashref} ) )
	  {
		print( "Could not get groupPredictions from " . $predictionsFile . " \n" );
		exit;
	  }

	#exit;
  }

sub geteggNOGPredictionsPairwise
  {
	my ($arg_ref)              = @_;
	my $predictionsFile        = $arg_ref->{'predictionFile'};
	my $og_assignments_hashref = $arg_ref->{'AllPredictionsHash'};
	my $ID2Group_href          = $arg_ref->{'ID2DatabaseOG'};
	my $ID2Gene_href           = $arg_ref->{'ID2GeneHash'};
	my $og_counter             = 0;
	my $previousOG             = -1;
	if ( !-e $predictionsFile || !-s $predictionsFile )
	  {
		print( "\tCould not read resolved mappings " . $predictionsFile . "\n" );
		return 0;
	  }
	my %allowed_species = ();
	my @species =
	  ( "6239", "7227", "7719", "9615", "9606", "9598", "7955", "99883", "13616", "10090", "10116", "9031" );
	foreach my $item (@species) { $allowed_species{$item} = 1; }

	#DEBUG("\tget_group2expandedIDsHash: reading from ".$predictionsFile."\n");
	## CHECK DEFINEDNESS
	# croak q{One/several of the parameters for 'read_og_assignments' was/were not defined.\n}
	#if any {!defined $_} $predictionsFile;
	### OPENING FILE
	open my $OG_FH, '<', $predictionsFile
	  or croak "Couldn't open " . $predictionsFile . ": $OS_ERROR";
	my $counter = 0;
	my @convertedMembers;
	my $has_mapped_member = 0;
	### READING FROM FILE
	while ( my $line = <$OG_FH> )
	  {
		next if $line =~ /^$/;
		next if $line =~ /^#/;
		chomp($line);
		my @splitted_line = split( /\t/, $line );
		my $member = $splitted_line[0];
		$member =~ /(\d+)\.(.*)/;
		my $species = $1;

		#print "Species is $species\n";
		next if !exists $allowed_species{$species};
		$member = $2;
		my $og = $splitted_line[3];

		#next if $og ne 'meNOG06811';
		$previousOG = $og if $previousOG eq '-1';
		if ( !defined $og || $og eq '' || !defined $member || $member eq '' )
		  {
			print( "Could not parse line $line in file " . $predictionsFile . ". Values are $og,$member\n" );
			exit;
		  }
		if ( !$og =~ /meNOG\d{5}/ )
		  {
			print "\tog $og in wrong format\n";
		  }

		#print "$member\n";
		# Map IDs to RefOG groups if possible
		# We can later easily check if a ID the e.g. hieranoid prediction is also in the RefOG dataset
		if ( $previousOG ne $og )
		  {

			#print "\t\t has different og ($previousOG - $og) and mapped_member? $has_mapped_member\n";
			#print "\t\t had theses members: ".join(",",@convertedMembers)."\n";
			if ($has_mapped_member)
			  {

				#print "\t\t\t -> save\n";
				#exit;
				$og_assignments_hashref->{$og} = join( ",", @convertedMembers );
			  }

			#print "\t\t\t -> just continue\n";
			@convertedMembers  = ();
			$has_mapped_member = 0;
			$previousOG        = $og;
		  }

		#next if $og ne 'meNOG06811';
		if ( exists $ID2Gene_href->{$member} )
		  {
			push( @convertedMembers, $ID2Gene_href->{$member} );
			$ID2Group_href->{ $ID2Gene_href->{$member} } = $og;
			$has_mapped_member = 1;

			#print "\t\t\tadded ".$ID2Gene_href->{$member}." (mapped)\n";
		  }
		else
		  {
			push( @convertedMembers, $member );

			#print "\t\t\tadded $member\n";
		  }

		#$ID2Group_href->{$member} = $og;
		$og_counter++;
	  }

	#print Dumper $og_assignments_hashref;
	#exit;
	close $OG_FH
	  or croak "Couldn't close '" . $predictionsFile . "': $OS_ERROR";
	if ( !keys( %{$og_assignments_hashref} ) )
	  {
		print( "Could not get groupPredictions from " . $predictionsFile . " \n" );
		exit;
	  }

	#exit;
  }

sub getInparanoidPredictionsPairwise
  {
	my ($arg_ref)              = @_;
	my $predictionsFile        = $arg_ref->{'predictionFile'};
	my $og_assignments_hashref = $arg_ref->{'AllPredictionsHash'};
	my $ID2Group_href          = $arg_ref->{'ID2DatabaseOG'};
	my $ID2Gene_href           = $arg_ref->{'ID2GeneHash'};
	my $AllUsedIDsFile         = $arg_ref->{'AllUsedIDsFile'};
	my $og_counter             = 0;
	my $previousOG             = -1;
	if ( !-e $predictionsFile || !-s $predictionsFile )
	  {
		print( "\tCould not read resolved mappings " . $predictionsFile . "\n" );
		return 0;
	  }
	### OPENING FILE
	open my $OG_FH, '<', $predictionsFile
	  or croak "Couldn't open " . $predictionsFile . ": $OS_ERROR";
	my $counter = 0;
	my @convertedMembers;
	my $no_mapped_member = 0;
	my %originalIDs;
	my %missingIDs;
	print "\treading all available headers into memory....";
	my @headers = `cat $AllUsedIDsFile`;
	my %allowedIDs;

	foreach (@headers)
	  {
		chomp;
		$allowedIDs{$_} = 1;
	  }
	#print "done\n";
	my $co = 0;

	#exit;
	print "\treading Inparanoid file $predictionsFile\n";

	#exit;
	### READING FROM FILE
  LINE:
	while ( my $line = <$OG_FH> )
	  {
		#print $line;
		next if $line =~ /^$/;
		next if $line =~ /^#/;
		chomp($line);
		my @splitted_line = split( /\t/, $line );
		my $member        = $splitted_line[4];
		my $og            = $splitted_line[0];
		$previousOG = $og if $previousOG eq '-1';

		if ( !defined $og || $og eq '' || !defined $member || $member eq '' )
		  {
			print( "Could not parse line $line in file " . $predictionsFile . ". Values are $og,$member\n" );
			exit;
		  }
		if ( !$og =~ /meNOG\d{5}/ )
		  {
			print "\tog $og in wrong format\n";
		  }

		# Map IDs to RefOG groups if possible
		# We can later easily check if a ID the e.g. hieranoid prediction is also in the RefOG dataset
		if ( $previousOG ne $og )
		  {

			#print "\t\t has different og ($previousOG - $og) and mapped_member? $has_mapped_member\n";
			#print "\t\t had theses members: ".join(",",@convertedMembers)."\n";
			if ( $no_mapped_member > 1 )
			  {

				#print "\t\t\t -> save\n";
				#exit;
				$og_assignments_hashref->{$og} = join( ",", @convertedMembers );
			  }

			#print "\t\t\t -> just continue\n";
			@convertedMembers = ();
			$no_mapped_member = 0;
			$previousOG       = $og;
		  }
		if ( exists $ID2Gene_href->{$member} )
		  {
			push( @convertedMembers, $ID2Gene_href->{$member} );
			$ID2Group_href->{ $ID2Gene_href->{$member} } = $og;
			$no_mapped_member++;

			#print "\t\t\tadded ".$ID2Gene_href->{$member}." (mapped)\n";
		  }
		else
		  {
			if ( !exists $allowedIDs{$member} )
			  {

				#       print "\tnot in dataset\n";
				next LINE;
			  }
			else
			  {
				push( @convertedMembers, $member );
				$missingIDs{$previousOG}{$member} = 1;
			  }

			#print "\t\t\tadded not in dataset? $member\n";
		  }
		push( @{ $originalIDs{$previousOG} }, $member );

		#$ID2Group_href->{$member} = $og;
		$og_counter++;
	  }
	close $OG_FH
	  or croak "Couldn't close '" . $predictionsFile . "': $OS_ERROR";
	my $printLine;
	foreach my $og ( keys( %{$og_assignments_hashref} ) )
	  {
		$printLine .= "$og (org): $originalIDs{$og}\n";
		$printLine .= "$og:" . $og_assignments_hashref->{$og} . "\n";
		$printLine .=
		  "\tmissing (" . keys( %{ $missingIDs{$og} } ) . "):" . join( ",", keys( %{ $missingIDs{$og} } ) ) . "\n";
	  }
	write_to_file(
				   {
					 file_name => basename($predictionsFile) . "_converted",
					 text      => $printLine
				   }
	);
	if ( !keys( %{$og_assignments_hashref} ) )
	  {
		print( "Could not get groupPredictions from " . $predictionsFile . " \n" );
		exit;
	  }

	#exit;
  }

sub OrthoMCL2OrthoBench
  {
	my ( $predictionsFile, $og_assignments_hashref, $ID2Group_href, $ensembl_mappings ) = (@_);

# Version 4
#my $orthomcl_prediction_file = "/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/groups_OrthoMCL-4.txt";
#my $orthomcl_prediction_file_converted = "/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/groups_OrthoMCL-4.txt_converted";
# Version 5
	my $orthomcl_prediction_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/groups_OrthoMCL-5.txt";
	my $orthomcl_prediction_file_converted =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/groups_OrthoMCL-5.txt_converted";
	my %allowedSpecies = (
						   hsap => 1,
						   mdom => 1,
						   clup => 1,    # Canis lupus familiaris
						   rnor => 1,
						   mmus => 1,
						   ptro => 1,
						   ggal => 1,
						   tnig => 1,
						   drer => 1,
						   cint => 1,
						   dmel => 1,
						   cele => 1,
	);
	open my $READ_FH, '<', $orthomcl_prediction_file
	  or croak "Couldn't open " . $orthomcl_prediction_file . ": $OS_ERROR";
	my $og_string;
	### READING FROM FILE
	while ( my $line = <$READ_FH> )
	  {
		my $og;

		#print "line is $line\n";
		my @usuableIDs;
		my @splitted_line = split( /\s/, $line );
		foreach (@splitted_line)
		  {

			#print $_."\n";
			#exit;
			if (/(OG\d{1}_\d+)*:/)
			  {
				$og = $1;
			  }
			else
			  {

				#print $_."\n";
				my ( $species, $id ) = split( /\|/, $_ );

				#print "species: $species,id: $id\n";
				next if ( not exists $allowedSpecies{$species} );
				push( @usuableIDs, $id );
			  }
		  }
		if ( !@usuableIDs )
		  {
			next;
		  }
		else
		  {
			$og_string .= "$og:" . join( ",", @usuableIDs ) . "\n";
		  }

		#print $og_string."\n";
		#exit;
	  }
	close $READ_FH
	  or croak "Couldn't close '" . $orthomcl_prediction_file . "': $OS_ERROR";
### OPENING FILE
	open my $out_FH, '>', $orthomcl_prediction_file_converted
	  or croak "Couldn't open '$orthomcl_prediction_file_converted': $OS_ERROR";
	### Writing file
	print {$out_FH} $og_string;
	### CLOSING FILE
	close $out_FH
	  or croak "Couldn't close '$orthomcl_prediction_file_converted': $OS_ERROR";
  }

sub getOrthoDBPredictionsPairwise
  {
	my ($arg_ref)              = @_;
	my $predictionsFile        = $arg_ref->{'predictionFile'};
	my $og_assignments_hashref = $arg_ref->{'AllPredictionsHash'};
	my $ID2Group_href          = $arg_ref->{'ID2DatabaseOG'};
	my $ID2Gene_href           = $arg_ref->{'ID2GeneHash'};
	my $og_counter             = 0;
	my %originalIDs;
	my %missingIDs;
	my $previousOG = -1;
	my %tmpOG_assignments_hash;
	my %hasConvertedMemberGroups;

	if ( !-e $predictionsFile || !-s $predictionsFile )
	  {
		print( "\tCould not read resolved mappings " . $predictionsFile . "\n" );
		return 0;
	  }
	my %allowed_species = (
		'Caenorhabditis elegans'  => 1,
		'Drosophila melanogaster' => 1,
		'C.intestinalis'          => 1,
		'Danio rerio'             => 1,
		'Tetraodon nigroviridis'  => 1,
		'Gallus gallus'           => 1,
		'Monodelphis domestica'   => 1,
		'Mus musculus'            => 1,
		'Rattus norvegicus'       => 1,
		'Canis familiaris'        => 1,
		'Pan troglodytes'         => 1,
		'Homo sapiens'            => 1,

		#'Nematostella vectensis'=> 1,
		#'M.brevicollis'=> 1,
		#'Hydra magnipapillata'=> 1
	);
	### OPENING FILE
	open my $OG_FH, '<', $predictionsFile
	  or croak "Couldn't open " . $predictionsFile . ": $OS_ERROR";
	my $counter = 0;
	my @convertedMembers;
	my $has_mapped_member = 0;
	### READING FROM FILE
	while ( my $line = <$OG_FH> )
	  {
		next if $line =~ /^$/;
		next if $line =~ /^#/;
		chomp($line);
		my @splitted_line = split( /\t/, $line );
		my $member        = $splitted_line[2];
		my $species       = $splitted_line[3];

		#print "Species is $species\n";
		next if !exists $allowed_species{$species};

		#$member = $2;
		my $og = $splitted_line[1];
		print "og $og, id:$member species: $species\n";
		if ( !defined $og || $og eq '' || !defined $member || $member eq '' )
		  {
			print( "Could not parse line $line in file " . $predictionsFile . ". Values are $og,$member\n" );
			exit;
		  }
		if ( !$og =~ /EOG\w{5}/ )
		  {
			print "\tog $og in wrong format\n";
			exit;
		  }

		#exit;
		# Map IDs to RefOG groups if possible
		# We can later easily check if a ID the e.g. hieranoid prediction is also in the RefOG dataset
		if ( exists $ID2Gene_href->{$member} )
		  {
			$tmpOG_assignments_hash{$og} .= $ID2Gene_href->{$member} . ",";
			$ID2Group_href->{ $ID2Gene_href->{$member} } = $og;

			#push(@convertedMembers, $ID2Gene_href->{$member});
			#$has_mapped_member = 1;
			print "\t\t\tadded converted member: " . $ID2Gene_href->{$member} . " (mapped)\n";
			$hasConvertedMemberGroups{$og} = 1;
		  }
		else
		  {
			$tmpOG_assignments_hash{$og} .= $member . ",";
			$ID2Group_href->{$member} = $og;
			$missingIDs{$og}{$member} = 1;

			#push(@convertedMembers, $member);
			print "\t\t\tadded $member\n";
		  }
		push( @{ $originalIDs{$og} }, $member );
		$og_counter++;
	  }
	close $OG_FH
	  or croak "Couldn't close '" . $predictionsFile . "': $OS_ERROR";
	if ( !keys(%tmpOG_assignments_hash) )
	  {
		print( "Did not find anything to map. Prediction file: " . $predictionsFile . " \n" );
		exit;
	  }
	my $printLine;
	foreach my $og ( keys(%tmpOG_assignments_hash) )
	  {

		#print "checking group $og\n";
		if ( !exists $hasConvertedMemberGroups{$og} )
		  {
			next;
		  }
		else
		  {
			$og_assignments_hashref->{$og} .= $tmpOG_assignments_hash{$og};
			$printLine .= "$og (org): " . join( ",", @{ $originalIDs{$og} } ) . "\n";
			$printLine .= "$og:" . $og_assignments_hashref->{$og} . "\n";
			$printLine .= "\tmissing (" .
			  keys( %{ $missingIDs{$og} } ) . "):" . join( ",", keys( %{ $missingIDs{$og} } ) ) . "\n";
		  }
	  }
	write_to_file(
				   {
					 file_name => basename($predictionsFile) . "_converted",
					 text      => $printLine
				   }
	);
	if ( !keys( %{$og_assignments_hashref} ) )
	  {
		print( "Could not get groupPredictions from " . $predictionsFile . " \n" );
		exit;
	  }

	#exit;
  }

sub convertTreeFamMapping2OrthoBench
  {
	my ( $predictionsFile, $og_assignments_hashref, $ID2Group_href, $ensembl_mappings ) = (@_);
	my $og_counter = 0;
	if ( !-e $predictionsFile || !-s $predictionsFile )
	  {
		print( "\tCould not read resolved mappings " . $predictionsFile . "\n" );
		return 0;
	  }

	#DEBUG("\tget_group2expandedIDsHash: reading from ".$predictionsFile."\n");
	## CHECK DEFINEDNESS
	# croak q{One/several of the parameters for 'read_og_assignments' was/were not defined.\n}
	#if any {!defined $_} $predictionsFile;
	my $treefam_genes =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/treefam/genes.txt.table";

	# READ SPECIES INFORMATION
	open my $treefam_genes_FH, '<', $treefam_genes
	  or croak "Couldn't open " . $treefam_genes . ": $OS_ERROR";
	my %transcript2protein;
	my %allowed_ids;
	my %allowed_species = (
							"9606"  => 1,    # Homo
							"10090" => 1,    # Mouse
							"6239"  => 1,    # Celegans
							"7955"  => 1,    # Danio
							"9615"  => 1,    # Canis
							"9598"  => 1,    # Pan
							"7719"  => 1,    # Ciona
							"99883" => 1,    # Tetrae
							"13616" => 1,    # Monodelphis
							"7227"  => 1,    # Drosophila
							"9031"  => 1,    # Gallus
							"10116" => 1,    # Rat
	);
	### READING FROM FILE
	my $counter = 0;
	while ( my $line = <$treefam_genes_FH> )
	  {
		my @splitted_line = split( /\s/, $line );

		#	print Dumper @splitted_line;
		#	exit;
		#print $line;
		my ( $taxID, $gene_id, $transcript_id, $protein_id ) =
		  ( $splitted_line[8], $splitted_line[1], $splitted_line[2], $splitted_line[4] );

		#rint "$taxID, $gene_id, $transcript_id, $protein_id\n";
		if ( exists $allowed_species{$taxID} )
		  {
			$allowed_ids{$gene_id}       = 1;
			$allowed_ids{$transcript_id} = 1;
			$allowed_ids{$protein_id}    = 1;

			#       $counter++;
		  }
	  }
	close $treefam_genes_FH;
	print "have collected " . keys(%allowed_ids) . " ids\n";
	my $treefam_converted_file = $predictionsFile . "_converted.txt";

	#print "\tsaving output in $treefam_converted_file\n";
	#exit;
	### OPENING FILE
	open my $ENSEMBL_FH, '<', $ensembl_mappings
	  or croak "Couldn't open " . $ensembl_mappings . ": $OS_ERROR";
	my %transcript2protein;
	### READING FROM FILE
	while ( my $line = <$ENSEMBL_FH> )
	  {
		chomp($line);
		my ( $protein, $gene, $protein2, $transcript ) = split( /\t/, $line );
		if ( !$protein || !$transcript )
		  {
			print "Problem parsing line $line (prot: $protein,trans: $transcript) \n";
		  }

		#print "\tsaving $transcript -> $protein\n";
		#exit;
		$transcript2protein{$transcript} = $protein;
	  }

 #open my $OG_OUT_FH, '<', $treefam_converted_file or croak "Couldn't open ".$treefam_converted_file.": $OS_ERROR";
 #close $OG_OUT_FH or die "Could not close $treefam_converted_file\n";
	if ( !keys %transcript2protein )
	  {
		print "\tCould not get hash transcript2protein. Exiting\n";
		exit;
	  }
	close $ENSEMBL_FH or die "Could not close $ensembl_mappings\n";
	### OPENING FILE
	open my $OG_FH, '<', $predictionsFile or croak "Couldn't open " . $predictionsFile . ": $OS_ERROR";
	my $counter = 0;
	### READING FROM FILE
	while ( my $line = <$OG_FH> )
	  {
		next if $line =~ /^$/;
		next if $line =~ /^#/;
		chomp($line);
		my @splitted_line = split( /\t/, $line );
		my $member = $splitted_line[0];

		#$member =~ /\d+\.(.*)/;
		#$member = $1;
		#chop($member);
		#print "\t\t$member -> ";
		if ( $member =~ /\.\d+$/ )
		  {
			$member = substr( $member, 0, ( length($member) - 2 ) );
		  }

		#my @tmp = split(/\./,$member);
		#$member = join("",@tmp[0,-2]);
		#print "$member \n";
		#exit;
		my $og               = $splitted_line[1];
		my $treeFamGroupType = $splitted_line[2];
		if ( $treeFamGroupType ne 'SEED' && $treeFamGroupType ne 'FULL' )
		  {
			print "line $line not readable. Cannot determine treeFamGroupType $treeFamGroupType\n";
			exit;
		  }

		#if($treeFamGroupType eq 'FULL'){
		#    next;
		#}
		if ( !defined $og || $og eq '' || !defined $member || $member eq '' )
		  {
			print( "Could not parse line $line in file " . $predictionsFile . ". Values are $og,$member\n" );
			exit;
		  }
		if ( !$og =~ /TF\d{6}/ )
		  {
			print "\tog $og in wrong format\n";
		  }
		print "\tsaving $og and $member from line:  \n";

		#exit;
		#if(not exists $transcript2protein{$member}){
		#        print "\tmoooooep. Could not map transcript $member\n";
		#        next;
		#}
		#my $protein4transcript = $transcript2protein{$member};
		if ( !exists $allowed_ids{$member} )
		  {
			print "\tId not allowed";
			next;
		  }
		else
		  {
			print "\tyes";
		  }
		print "\n";
		$og_assignments_hashref->{$og} .= $member . ",";
		$ID2Group_href->{$member} = $og;
		$og_counter++;
	  }
	close $OG_FH
	  or croak "Couldn't close '" . $predictionsFile . "': $OS_ERROR";
	open my $OG_OUT_FH, '>', $treefam_converted_file
	  or croak "Couldn't open " . $treefam_converted_file . ": $OS_ERROR";
	foreach my $og ( keys( %{$og_assignments_hashref} ) )
	  {
		print {$OG_OUT_FH} "$og:" . $og_assignments_hashref->{$og} . "\n";
	  }
	close $OG_OUT_FH or die "Could not close $treefam_converted_file\n";
	if ( !keys( %{$og_assignments_hashref} ) )
	  {
		print( "Could not get groupPredictions from " . $predictionsFile . " \n" );
		exit;
	  }
	print "\tdone converting\n";
  }

sub readOrthoBenchPredictions
  {
	#### PARAMETER VARIABLES
	my ($arg_ref)    = @_;
	my $file_name    = $arg_ref->{file};
	my $hashref      = $arg_ref->{hash};
	my $totalHashref = $arg_ref->{totalhashref};
	my $scoreLevel;
	my $has_percentage = 0;
	### OPENING FILE
	#print "\treading reference predictions from $file_name\n";
	open my $in, '<', $file_name
	  or croak "Couldn't open '$file_name': $OS_ERROR";
	### Writing file
	while ( defined( my $line = <$in> ) )
	  {
		chomp($line);
		next if $line =~ /^$/;
		next if $line =~ /^\"/;

		#next if $line =~ /^Total/;
		next if $line =~ /Table/;

		#if($line =~ /RefOG\d+/){
		#print $line."\n";
		#}
		if ( $line =~ /RefOG\d+/ || $line =~ /^Total/ )
		  {
			my ( $RefOG, $RefOG_Size ) = split( "\t", $line );
			$line =~ /\w+\s*\d+\s*(.*)/;
			my $rest                 = $1;
			my @splittedMethodValues = split( "\t", $rest );
			my $method_counter       = 0;
			my $level_counter        = 0;
			my @levels               = ( "Gene-Level", "Gene-Level", "Group-Level", "Group-Level" );
			my @methods = (
							"TreeFam",  "TreeFam",  "TreeFam",  "TreeFam",  "eggNOG",  "eggNOG",
							"eggNOG",   "eggNOG",   "OrthoDB",  "OrthoDB",  "OrthoDB", "OrthoDB",
							"OrthoMCL", "OrthoMCL", "OrthoMCL", "OrthoMCL", "OMA",     "OMA",
							"OMA",      "OMA"
			);
			my @valueTypes = ( "FN", "FP", "Fissions", "Fusion" );
			my $noOfValueInLine = 0;

			foreach my $current_numericalValue (@splittedMethodValues)
			  {
				my $current_method    = $methods[$noOfValueInLine];
				my $current_level     = $levels[ $noOfValueInLine % 4 ];
				my $current_valueType = $valueTypes[ $noOfValueInLine % 4 ];
				$noOfValueInLine++;

#print "\$hashref->{$RefOG}{$current_method}{$current_level}{$current_valueType} = $current_numericalValue;\n" if $line =~ /^Total/;
				$hashref->{$RefOG}{$current_method}{$current_level}{$current_valueType} = $current_numericalValue;

	 #print "\$hashref->{$RefOG}{$current_method}{$current_level}{$current_valueType} = $current_numericalValue\n";
			  }
			next;
		  }
		else
		  {
			my %nopercentageHash = (
									 "Erroneously_assigned_genes" => 1,
									 "Missing_genes"              => 1,
									 "Fusion"                     => 1,
									 "Fissions"                   => 1,
			);
			if ( $line =~ /Accuracy at the gene level/ )
			  {
				$scoreLevel = "Gene-Level";
			  }
			if ( $line =~ /Accuracy at the group level/ )
			  {
				$scoreLevel = "Group-Level";
			  }
			my $statisticsType;
			if ( $line =~ /^Accurately predicted RefOGs/ )
			  {
				$statisticsType = "Accurately_predicted_RefOGs";
			  }
			if ( $line =~ /^Erroneously assigned genes/ )
			  {
				$statisticsType = "Erroneously_assigned_genes";
			  }
			if ( $line =~ /^Missing genes/ )
			  {
				$statisticsType = "Missing_genes";
			  }
			if ( $line =~ /RefOGs affected by erroneously affected genes/ )
			  {
				$statisticsType = "RefOGs_affected_by_erroneously_affected_genes";
			  }
			if ( $line =~ /^RefOGs affected by missing genes/ )
			  {
				$statisticsType = "RefOGs_affected_by_missing_genes";
			  }
			if ( $line =~ /^Fusion/ )
			  {
				$statisticsType = "Fusion";

				#next;
			  }
			if ( $line =~ /^Fissions/ )
			  {
				$statisticsType = "Fissions";
			  }
			if ( $line =~ /^RefOGs affected by fusions/ )
			  {
				$statisticsType = "RefOGs_affected_by_fusions";
			  }
			if ( $line =~ /^RefOGs affected by fissions/ )
			  {
				$statisticsType = "RefOGs_affected_by_fissions";
			  }
			if ($statisticsType)
			  {

				#print "stat type: $statisticsType\n";
				my @splittedMethodValues = split( "\t", $line );
				shift(@splittedMethodValues);

				#print join(",",@splittedMethodValues)."\n";
				my @value = ( "percentage", "counts" );
				my @methods = (
								"TreeFam",  "TreeFam",  "TreeFam",  "TreeFam",  "eggNOG",  "eggNOG",
								"eggNOG",   "eggNOG",   "OrthoDB",  "OrthoDB",  "OrthoDB", "OrthoDB",
								"OrthoMCL", "OrthoMCL", "OrthoMCL", "OrthoMCL", "OMA",     "OMA",
								"OMA",      "OMA",
				);
				my $noOfValueInLine = 0;
				foreach my $current_numericalValue (@splittedMethodValues)
				  {
					if ( !$current_numericalValue =~ /\d/ )
					  {
						$noOfValueInLine++;
						next;
					  }
					my $current_method = $methods[$noOfValueInLine];
					next
					  if ( !defined($current_method) || $current_method eq '' );

					#my $current_level = $levels[$noOfValueInLine%4];
					my $current_valueType = $value[ $noOfValueInLine % 2 ];

#if($statisticsType){
#$hashref->{'Total'}{$current_method}{$scoreLevel}{'Fusion'}{$current_valueType} = $current_numericalValue;
#}
#else{
#print "hashref->{\"Total\"}{$current_method}{$scoreLevel}{$statisticsType}{$current_valueType} = $current_numericalValue\n";
#$hashref->{'Total'}{$current_method}{$scoreLevel}{$statisticsType}{$current_valueType} = $current_numericalValue;
#}
					$totalHashref->{$current_method}{$scoreLevel}{$statisticsType}{$current_valueType} =
					  $current_numericalValue;
					$noOfValueInLine++;
				  }
			  }
		  }

		#exit;
	  }

	#$hashref->{'Total'}{"TreeFam"}{'Group-Level'}{'Fusion'}{'counts'} = 38;
	#$hashref->{'Total'}{"eggNOG"}{'Group-Level'}{'Fusion'}{'counts'} = 24;
	#$hashref->{'Total'}{"OrthoDB"}{'Group-Level'}{'Fusion'}{'counts'} = 12;
	#$hashref->{'Total'}{"OrthoMCL"}{'Group-Level'}{'Fusion'}{'counts'} = 16;
	#$hashref->{'Total'}{"OMA"}{'Group-Level'}{'Fusion'}{'counts'} = 0;
	#exit;	### CLOSING FILE
	close $in or croak "Couldn't close '$file_name': $OS_ERROR";
  }

sub compare2orthobench
  {

# 1. read all hieranoid predictions
#my $hieranoid_prediction_file = "/Users/fab/Documents/documents/work/current_projects/hieranoid/results_ensembl/nodes/Eukaryota/Eukaryota.expandedGroups.txt";
# test
#my $hieranoid_prediction_file = "/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/tmp/Homininae.expandedGroups.txt";
#my $hieranoid_prediction_file = "/Users/fab/Documents/documents/work/current_projects/hieranoid/results_ensembl/nodes/Opisthokonta/Opisthokonta.expandedGroups.txt";
	my $hieranoid_prediction_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/archived_results/results_ensembl/nodes/Bilateria/Bilateria.expandedGroups.txt";
	my $hieranoid_prediction_consensus_outgroup_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/results_orthobench_consensus_outgroup/nodes/Bilateria/Bilateria.expandedGroups.txt";
	my $hieranoid_profile_prediction_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/archived_results/results_ensembl_profile/nodes/Bilateria/Bilateria.expandedGroups.txt";

#my $hieranoid_profile_prediction_file = "/Users/fab/Documents/documents/work/current_projects/hieranoid/results_ensembl_profile/nodes/Opisthokonta/Opisthokonta.expandedGroups.txt";
	my $oma_prediction_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/oma-groups_may2011_converted.txt";
	my $eggNOG_prediction_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/meNOG.mapping.txt";
	my $orthomcl_prediction_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/groups_OrthoMCL-5.txt_converted";
	my $treeFamOldPredictionFile =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/fam_genes.txt.table";
	my $treeFam_prediction_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/fam_genes.txt.table_converted.txt";
	my $ensembl_mapping =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/ensembl_ids.txt";
	my $hieranoid_fasta_file =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/ensemblv60_12species";
	my $orthobench_fasta_files =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/predictions/orthobench/";
	my $orthobench_predictions =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/pictures/R_sources/OrthoBenchReferenceTable.txt";
	my $orthobench_groups =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/IDs_OrthoBench2.txt";
	my $OrthoBenchRawData =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/pictures/R_sources/OrthoBenchRawData.txt";

	#my %hieranoid_prediction;
	my %ID2HieranoidOG;
	my %ID2DatabaseOG;
	my %allPredictions;

#my @predictionMethods_array = qw(Hieranoid_consensus Hieranoid_profile OMA_new eggNOG_new TreeFam_new OrthoMCL_new);
	my @predictionMethods_array =
	  qw(Hieranoid_consensus Hieranoid_profile Hieranoid_consensus_outgroup OMA_new eggNOG_new OrthoMCL_new);

#convertTreeFamMapping2OrthoBench($treeFamOldPredictionFile,\%{$allPredictions{"TreeFam_new"}}, \%{$ID2DatabaseOG{"TreeFam_new"}},$ensembl_mapping);
#exit;
	my %ID2refOG;

	#  'refOG0001' -> @(ENSTNIP00000005325,ENSDARP00000071531)
	my %refOG_hash;
	my %refOGGenes_hash;

	#  {'refOG0001'}{'ENSTNIP00000005325'} = 2
	## OrthoBench groups
	getOrthoBenchGroups( $orthobench_groups, \%ID2refOG, \%refOG_hash, \%refOGGenes_hash );
	if ( !keys %ID2refOG )
	  {
		print "\tCould not get ID2refOGs. Exiting\n";
		exit;
	  }
	if ( !keys %refOG_hash )
	  {
		print "\tCould not get refOG_hash. Exiting\n";
		exit;
	  }
	if ( !keys(%refOGGenes_hash) )
	  {
		print("Empty hash  refOGGenes_hash found\n");
		exit;
	  }
	## Hieranoid predictions (consensus)
	print "\treading Hieranoid predictions (consensus) ....";
	getPredictionMethodPredictions(
									$hieranoid_prediction_file,
									\%{ $allPredictions{"Hieranoid_consensus"} },
									\%{ $ID2DatabaseOG{"Hieranoid_consensus"} }
	);
	if ( !keys %{ $allPredictions{"Hieranoid_consensus"} } )
	  {
		print "\tCould not get Hieranoid_consensus hieranoid_prediction. Exiting\n";
		exit;
	  }
	if ( !keys %{ $ID2DatabaseOG{"Hieranoid_consensus"} } )
	  {
		print "\tCould not get Hieranoid_consensus ID2Group_hash. Exiting\n";
		exit;
	  }
	print "done\n";
	## Hieranoid predictions (profile)
	print "\treading Hieranoid predictions (profile) ....";
	getPredictionMethodPredictions(
									$hieranoid_profile_prediction_file,
									\%{ $allPredictions{"Hieranoid_profile"} },
									\%{ $ID2DatabaseOG{"Hieranoid_profile"} }
	);
	if ( !keys %{ $allPredictions{"Hieranoid_profile"} } )
	  {
		print "\tCould not get Hieranoid_profile hieranoid_prediction. Exiting\n";
		exit;
	  }
	if ( !keys %{ $ID2DatabaseOG{"Hieranoid_profile"} } )
	  {
		print "\tCould not get Hieranoid_profile ID2Group_hash. Exiting\n";
		exit;
	  }
	print "done\n";
	## Hieranoid predictions outgroup (consensus)
	print "\treading Hieranoid predictions outgroup  ....";
	getPredictionMethodPredictions(
									$hieranoid_prediction_consensus_outgroup_file,
									\%{ $allPredictions{"Hieranoid_consensus_outgroup"} },
									\%{ $ID2DatabaseOG{"Hieranoid_consensus_outgroup"} }
	);
	if ( !keys %{ $allPredictions{"Hieranoid_consensus_outgroup"} } )
	  {
		print "\tCould not get Hieranoid_profile hieranoid_prediction. Exiting\n";
		exit;
	  }
	if ( !keys %{ $ID2DatabaseOG{"Hieranoid_consensus_outgroup"} } )
	  {
		print "\tCould not get Hieranoid_profile ID2Group_hash. Exiting\n";
		exit;
	  }
	print "done\n";
	## OMA predictions
	print "\treading OMA predictions....";
	getPredictionMethodPredictions(
									$oma_prediction_file,
									\%{ $allPredictions{"OMA_new"} },
									\%{ $ID2DatabaseOG{"OMA_new"} }
	);
	if ( !keys %{ $allPredictions{"OMA_new"} } )
	  {
		print "\tCould not get oma_prediction. Exiting\n";
		exit;
	  }
	if ( !keys %{ $ID2DatabaseOG{"OMA_new"} } )
	  {
		print "\tCould not get OMA ID2Group_hash. Exiting\n";
		exit;
	  }
	print "done\n";
	## eggNOG predictions
	print "\treading eggNOG predictions....";
	geteggNOGPredictions(
						  $eggNOG_prediction_file,
						  \%{ $allPredictions{"eggNOG_new"} },
						  \%{ $ID2DatabaseOG{"eggNOG_new"} }
	);
	if ( !keys %{ $allPredictions{"eggNOG_new"} } )
	  {
		print "\tCould not get eggNOG. Exiting\n";
		exit;
	  }
	if ( !keys %{ $ID2DatabaseOG{"eggNOG_new"} } )
	  {
		print "\tCould not get eggNOG ID2Group_hash. Exiting\n";
		exit;
	  }
	print "done\n";
	## TreeFam predictions
	print "\treading TreeFam predictions....";
	getPredictionMethodPredictions(
									$treeFam_prediction_file,
									\%{ $allPredictions{"TreeFam_new"} },
									\%{ $ID2DatabaseOG{"TreeFam_new"} }
	);
	if ( !keys %{ $allPredictions{"TreeFam_new"} } )
	  {
		print "\tCould not get TreeFam_new. Exiting\n";
		exit;
	  }
	if ( !keys %{ $ID2DatabaseOG{"TreeFam_new"} } )
	  {
		print "\tCould not get TreeFam_new ID2Group_hash. Exiting\n";
		exit;
	  }
	print "done\n";
	## OrthoMCL predictions
	print "\treading OrthoMCL predictions....";
	getPredictionMethodPredictions(
									$orthomcl_prediction_file,
									\%{ $allPredictions{"OrthoMCL_new"} },
									\%{ $ID2DatabaseOG{"OrthoMCL_new"} }
	);
	if ( !keys %{ $allPredictions{"OrthoMCL_new"} } )
	  {
		print "\tCould not get OrthoMCLnew. Exiting\n";
		exit;
	  }
	if ( !keys %{ $ID2DatabaseOG{"OrthoMCL_new"} } )
	  {
		print "\tCould not get OrthoMCL_new ID2Group_hash. Exiting\n";
		exit;
	  }
	print "done\n";

	#exit;
	#my %sequenceID2species;
	#getspecies4sequence({folder => $hieranoid_fasta_file, hashref => \%sequenceID2species});
	#if(!keys %sequenceID2species){
	#        print "\tCould not get sequence 2 species information. Exiting\n";
	#        exit;
	#}
	my %predictions_hash;
	my %total_prediction_hash;
	&readOrthoBenchPredictions(
								{
								  file         => $orthobench_predictions,
								  hash         => \%predictions_hash,
								  totalhashref => \%total_prediction_hash
								}
	);
	my %overallStatistics;
	my $totalNumberCheckedRefOGGenes       = 0;
	my $totalNumberCheckedHieranoidOGGenes = 0;
	my $totalNumberOfFission               = 0;
	my $totalNumberOfFusion                = 0;
	my $GroupsAffectedByFissions           = 0;
	my $GroupsAffectedByFusion             = 0;
	my %helperHash;
	my $groupCount = 0;

	# For each REFOG
	foreach my $currentRefOG ( sort keys(%refOG_hash) )
	  {
		print "checking $currentRefOG \n";    #with ".join(",",@{$refOG_array{$currentRefOG}})."\n";
		my $fusionDetected4RefOG = 0;
		foreach my $predictionMethod (@predictionMethods_array)
		  {
			my %largestGroupsCount;
			my %predictionMethodOGcount;
			my $noOfMissingGenes     = 0;
			my $noOfFalsePredictions = 0;
			my $noOfFissions         = 0;
			my $noOfFusion           = 0;
			my %speciesSequencesHash;

			#next if $currentRefOG ne 'RefOG055';
			print "$predictionMethod\n";    #with ".join(",",@{$refOG_array{$currentRefOG}})."\n";
			                                # identify Prediction method OG with largest overlap
			foreach my $refOGGene ( keys( %{ $refOG_hash{$currentRefOG} } ) )
			  {
				foreach my $refOGGeneID ( keys( %{ $refOG_hash{$currentRefOG}{$refOGGene} } ) )
				  {
					if ( exists $ID2DatabaseOG{$predictionMethod}{$refOGGeneID} )
					  {
						$predictionMethodOGcount{ $ID2DatabaseOG{$predictionMethod}{$refOGGeneID} }{$refOGGeneID} =
						  1;
					  }
				  }
			  }
			my $numberOfTPs = 0;
			if ( keys(%predictionMethodOGcount) == 0 )
			  {
				print "\tno largest group found\n";
				$noOfMissingGenes += keys( %{ $refOG_hash{$currentRefOG} } );

				#exit;
			  }
			else
			  {

				#print "\tRefGroup had multiple Hieranoid groups\n";
				print "\thas groups (";
				foreach my $og ( keys %predictionMethodOGcount )
				  {
					my $numberOfGenes = keys( %{ $predictionMethodOGcount{$og} } );
					print "$og ($numberOfGenes), ";
				  }
				print ")\n";

				# count as fission
				if ( keys(%predictionMethodOGcount) > 1 )
				  {
					$helperHash{$predictionMethod}{totalNumberOfFission} += ( keys(%predictionMethodOGcount) - 1 );
					$helperHash{$predictionMethod}{GroupsAffectedByFissions}++;
					$noOfFissions = keys(%predictionMethodOGcount);
				  }
				my @largestGroupArray;
				my $max = -1;
				my $largestGroup;
				my %genesFound4RefOG;

				# find group/s with largest overlap
				my $maxNoGenes = -1;
				foreach my $OG ( keys %predictionMethodOGcount )
				  {
					my $noOfGenes = keys( %{ $predictionMethodOGcount{$OG} } );

					#print "\tchecking ".$noOfGenes."\n";
					$maxNoGenes = $noOfGenes if $noOfGenes > $maxNoGenes;
				  }
				my %largestGroups;

				# Add all OGs with max number of Genes
				foreach my $currentPredictionMethodOG ( keys(%predictionMethodOGcount) )
				  {
					my $noOfGenes = keys( %{ $predictionMethodOGcount{$currentPredictionMethodOG} } );
					if ( $noOfGenes == $maxNoGenes )
					  {
						$largestGroups{$currentPredictionMethodOG} = 1;
					  }
				  }
				if ( !keys(%largestGroups) )
				  {
					print("Empty hash found largestGroups\n");
					exit;
				  }
				print "\t\tlargest groups :" . join( ",", keys(%largestGroups) ) . "\n";
				if ( keys(%largestGroups) > 1 )
				  {
					print "\tHave several largest groups\n";

					#exit;
				  }

				#print "\tlargest group: $largestGroup\n";
				# Iterating over groups
				foreach my $currentPredictionMethodOG ( keys(%predictionMethodOGcount) )
				  {
					my $numberOfFPs4group = 0;
					print "\t$currentPredictionMethodOG\n";

					#print "(largest group) \n";
					#print "\tin largest group. Iterating over all Hieranoid predictions in group\n";
					my @predictionMethodOGGroupMember =
					  split( ",", $allPredictions{$predictionMethod}{$currentPredictionMethodOG} );
					if ( !scalar(@predictionMethodOGGroupMember) )
					  {
						print "\tProblem reading hieranoid group member for group $currentPredictionMethodOG\n";
					  }

					# iterate over HieranoidOG
					# {'HieranoidOG1'}{ID1} -> 1
					#        {ID2} -> 2
					#        {ID3} -> 3
					# {'HieranoidOG2'}{ID4}
				  PREDICTIONMETHOD_OG_ID:
					foreach my $currentPredictionMethodOGID (@predictionMethodOGGroupMember)
					  {
						print "\t\t"
						  . $predictionMethod
						  . "_OGID "
						  . $totalNumberCheckedHieranoidOGGenes++
						  . ": $currentPredictionMethodOGID\t";

						# {'HieranoidOG1'}{ID1}
						#        {ID2}
						#        {ID3} -> not in RefOG --> FP
						my $inLargestGroup =
						  ( exists $largestGroups{$currentPredictionMethodOG} )
						  ? 1
						  : 0;
						my $hasRefOGGroup =
						  ( exists $ID2refOG{$currentPredictionMethodOGID} )
						  ? 1
						  : 0;
						my $maps2correctRefOG = 0;
						if ($hasRefOGGroup)
						  {
							$maps2correctRefOG =
							  ( $ID2refOG{$currentPredictionMethodOGID} eq $currentRefOG ) ? 1 : 0;
						  }
						print "(LG:$inLargestGroup\tRefOG:$hasRefOGGroup ("
						  . $ID2refOG{$currentPredictionMethodOGID}
						  . ")\t:CorrRefOG:$maps2correctRefOG) ";
						if ($inLargestGroup)
						  {
							if ($hasRefOGGroup)
							  {
								if ($maps2correctRefOG)
								  {
									my $geneNumber =
									  $refOGGenes_hash{ $ID2refOG{$currentPredictionMethodOGID} }
									  {$currentPredictionMethodOGID};
									if ( !exists $genesFound4RefOG{$geneNumber} )
									  {
										print "(TP (OG: "
										  . $ID2refOG{$currentPredictionMethodOGID}
										  . " - Ref: $currentRefOG), gene number $geneNumber\n";
										$numberOfTPs++;
									  }
									else
									  {
										print " (TP, but gene number already exists: $geneNumber)\n";
									  }
									$genesFound4RefOG{$geneNumber} = 1;
								  }
								else
								  {

									# FP in largest group
									print "-> FP (maps to wrong OG: "
									  . $ID2refOG{$currentPredictionMethodOGID} . ")\n";
									$largestGroupsCount{$currentPredictionMethodOG}{'noOfFalsePredictions'} += 1;
									$noOfFalsePredictions++;
								  }
							  }
							else
							  {

								# FP in largest group
								print "-> FP (largestGroup, not in RefOG)\n";
								$largestGroupsCount{$currentPredictionMethodOG}{'noOfFalsePredictions'} += 1;
								$noOfFalsePredictions++;
							  }
						  }
						else
						  {
							if ($hasRefOGGroup)
							  {
								if ($maps2correctRefOG)
								  {
									print "-> TP (not largest OG group\n";
								  }
								else
								  {
									$numberOfFPs4group++;
									print "-> FP (not largest OG group)\n";
								  }
							  }
							else
							  {
								$numberOfFPs4group++;
								print "-> FP (no RefOG mapping)\n";
							  }
						  }
					  }

					# Count Fusion events
					print "\t\t\tfound $numberOfFPs4group for group)\n";
					if ( $numberOfFPs4group >= 3 )
					  {
						$helperHash{$predictionMethod}{totalNumberOfFusion}++;
						$noOfFusion++;
						$fusionDetected4RefOG = 1;
						$helperHash{$predictionMethod}{noOfFusion}++;
					  }
				  }

				# iterate over RefOG
				# checking for missing genes
				# {'REFOG1'} -> [ID1,
				#                ID2,
				#               ID3]
				foreach my $refOGGene ( sort keys( %{ $refOG_hash{$currentRefOG} } ) )
				  {
					my $foundGene = 0;
					print "\t\tRef_OG Gene $refOGGene\t";
					if ( exists $genesFound4RefOG{$refOGGene} )
					  {
						print "\talready counted in largest group\n";
						next;
					  }
					else
					  {
						foreach my $RefOGID (
											  sort
											  keys( %{ $refOG_hash{$currentRefOG}{$refOGGene} } )
						  )
						  {

							# {'REFOG1'} -> [ID1 -> has Hieranoid prediction in largest/any group,
							#                ID2,
							if ( exists $ID2DatabaseOG{$predictionMethod}{$RefOGID} )
							  {
								my $OG4ID = $ID2DatabaseOG{$predictionMethod}{$RefOGID};
								my $inLargestGroup = ( exists $largestGroups{$OG4ID} ) ? 1 : 0;

								#print "-> ok\n";
								if ($inLargestGroup)
								  {

									#print "-> ok\n";
									$foundGene = 1;
								  }
								else
								  {
									print "found, but not in largest OG";
								  }
							  }
						  }
						if ( !$foundGene )
						  {

							# {'REFOG1'} -> [ID1 -> has no prediction in largest overlapping Hieranoid group,
							#                ID2,
							print "-> FN\n";
							$noOfMissingGenes++;
							$helperHash{$predictionMethod}{noOfMissingGenes}++;
						  }
						else
						  {
							print "ok\n";
						  }
					  }
				  }

				#$helperHash{$predictionMethod}{GroupsAffectedByFusion}++;
				#$GroupsAffectedByFusion++ if $groupHasFusionEvent;
				if ( keys(%largestGroupsCount) > 1 )
				  {

					# now select best overlap group
					my %largestGroupsCounterHash;
					my $smallestNumberErrors = -1;
					print "\t\tFinding best group\n";
					foreach my $largestOG ( keys(%largestGroupsCount) )
					  {
						my $groupScore = (
										   ( defined( $largestGroupsCount{$largestOG}{noOfFalsePredictions} ) )
										   ? $largestGroupsCount{$largestOG}{noOfFalsePredictions}
										   : 0
						);
						$smallestNumberErrors = $groupScore
						  if $smallestNumberErrors == '-1';
						print "\t\t\t\t$groupScore < $smallestNumberErrors : ";
						if ( $groupScore <= $smallestNumberErrors )
						  {
							$smallestNumberErrors = $groupScore;
							$largestGroupsCounterHash{$smallestNumberErrors} = $largestOG;
						  }
					  }

					# set to value from overlapping group with the least number of errors: we want to be friendly.
					my $bestOG = $largestGroupsCounterHash{$smallestNumberErrors};

#$noOfMissingGenes = ((defined($largestGroupsCount{$bestOG}{noOfMissingGenes}))? $largestGroupsCount{$bestOG}{noOfMissingGenes}: 0);
					$noOfFalsePredictions = (
											  ( defined( $largestGroupsCount{$bestOG}{noOfFalsePredictions} ) )
											  ? $largestGroupsCount{$bestOG}{noOfFalsePredictions}
											  : 0
					);
					print map { "\t\t$_ => " . $largestGroupsCount{$_}{noOfFalsePredictions} . " \n" }
					  keys(%largestGroupsCount);
					print "\tbest group is "
					  . $largestGroupsCounterHash{$smallestNumberErrors}
					  . " with $smallestNumberErrors errors\n";
					print map { "\t\t$_ => " . $largestGroupsCounterHash{$_} . "\n" }
					  keys(%largestGroupsCounterHash);
				  }

				#else{
				#    $noOfFalsePredictions = $noOfFalsePredictions;
				#}
				# at least 1 OG with > 3 FP
				if ($fusionDetected4RefOG)
				  {
					$helperHash{$predictionMethod}{GroupsAffectedByFusion}++;
				  }
			  }
			$overallStatistics{$predictionMethod}{$currentRefOG}{"FusionEvents"}     = $noOfFusion;
			$overallStatistics{$predictionMethod}{$currentRefOG}{"FissionEvents"}    = $noOfFissions;
			$overallStatistics{$predictionMethod}{$currentRefOG}{"missingGenes"}     = $noOfMissingGenes;
			$overallStatistics{$predictionMethod}{$currentRefOG}{"falsePredictions"} = $noOfFalsePredictions;
			$overallStatistics{$predictionMethod}{$currentRefOG}{"total"} =
			  keys( %{ $refOG_hash{$currentRefOG} } );
			$overallStatistics{$predictionMethod}{$currentRefOG}{"true_positives"} = $numberOfTPs;

			#$overallStatistics{$predictionMethod}{$currentRefOG}{"true_positives"} = $numberOfTPs;
			print "\t\t\tset no TPs to "
			  . $overallStatistics{$predictionMethod}{$currentRefOG}{"true_positives"} . "\n";
			print
"\t\t$predictionMethod: group statistic: FN: $noOfMissingGenes, FP: $noOfFalsePredictions. TP: $numberOfTPs, Fusion: $noOfFusion, Fission: $noOfFissions (Ref genes "
			  . keys( %{ $refOG_hash{$currentRefOG} } ) . ")\n";

			# Check if we counted correctly
			if ( $noOfMissingGenes + $numberOfTPs != keys( %{ $refOG_hash{$currentRefOG} } ) )
			  {
				print "\tMiscounted something here\n";
				exit;
			  }
		  }    # END Prediction Method
		       #last if $groupCount++ > 5;;
		       #last;
	  }    # END Ref OG
	print Dumper %overallStatistics;

	#
	#     Statistics & printing summary
	#
	my $totalMissingGenes            = 0;
	my $totalfalsePredictions        = 0;
	my $noOfCorrectGroups            = 0;
	my $noOfMissingGeneAffectedGroup = 0;
	my $noOfFalseGeneAffectedGroup   = 0;
	my $AccPredictedRefOGs_group     = 0;
	foreach my $predictionMethod ( sort keys %overallStatistics )
	  {

		foreach my $og ( keys %{ $overallStatistics{$predictionMethod} } )
		  {

#print "og: $og\tmethod: $predictionMethod\n";
#   print $og."\t".$overallStatistics{$predictionMethod}{$og}{"total"}."\t".$overallStatistics{$predictionMethod}{$og}{"missingGenes"}."\t".$overallStatistics{$predictionMethod}{$og}{"falsePredictions"}."\n";
			$helperHash{$predictionMethod}{totalMissingGenes} +=
			  $overallStatistics{$predictionMethod}{$og}{"missingGenes"};
			$helperHash{$predictionMethod}{totalfalsePredictions} +=
			  $overallStatistics{$predictionMethod}{$og}{"falsePredictions"};
			$helperHash{$predictionMethod}{noOfCorrectGroups}++
			  if (    $overallStatistics{$predictionMethod}{$og}{"missingGenes"} == 0
				   && $overallStatistics{$predictionMethod}{$og}{"falsePredictions"} == 0 );
			$helperHash{$predictionMethod}{noOfMissingGeneAffectedGroup}++
			  if ( $overallStatistics{$predictionMethod}{$og}{"missingGenes"} != 0 );
			$helperHash{$predictionMethod}{noOfFalseGeneAffectedGroup}++
			  if ( $overallStatistics{$predictionMethod}{$og}{"falsePredictions"} != 0 );
			$helperHash{$predictionMethod}{AccPredictedRefOGs_group}++
			  if (    $overallStatistics{$predictionMethod}{$og}{"FusionEvents"} == 0
				   && $overallStatistics{$predictionMethod}{$og}{"FissionEvents"} == 0 );
			$helperHash{$predictionMethod}{total_true_positives} +=
			  $overallStatistics{$predictionMethod}{$og}{"true_positives"};
		  }
	  }

	#print Dumper %helperHash;
	#exit;
	#methods on y-axis
	#different values on x-axis
	my @methods = ( "TreeFam", "eggNOG", "OrthoDB", "OrthoMCL", "OMA" );
	my $lines2print;
	$lines2print .=
"Method\tFN\tFP\tFissions\tFusion\tAccPredictedRefOGs\tRefOfByFN\tRefOfByFP\tRefOfByFusion\tRefOfByFissions\tAccPredictedRefOGsGroup\tTotalTPs\tTPPercentage\n";
	foreach my $currentMethod (@methods)
	  {

		#print Dumper $predictions_hash{'Total'}{$currentMethod};
		#print Dumper $total_prediction_hash{$currentMethod};
		#exit;
		my ( $FN, $FP, $Fissions, $Fusion, $AccPredictedRefOGs, $RefOfByFN, $RefOfByFP, $RefOfByFusion,
			 $RefOfByFissions, $AccPredictedRefOGsGroup )
		  = (
			  $predictions_hash{'Total'}{$currentMethod}{"Gene-Level"}{"FN"},
			  $predictions_hash{'Total'}{$currentMethod}{"Gene-Level"}{"FP"},
			  $predictions_hash{'Total'}{$currentMethod}{"Group-Level"}{"Fissions"},
			  $predictions_hash{'Total'}{$currentMethod}{"Group-Level"}{"Fusion"},
			  $total_prediction_hash{$currentMethod}{"Gene-Level"}{"Accurately_predicted_RefOGs"}{"counts"},
			  $total_prediction_hash{$currentMethod}{"Gene-Level"}{"RefOGs_affected_by_missing_genes"}{"counts"},
			  $total_prediction_hash{$currentMethod}{"Gene-Level"}{"RefOGs_affected_by_erroneously_affected_genes"}
				{"counts"},
			  $total_prediction_hash{$currentMethod}{"Group-Level"}{"RefOGs_affected_by_fusions"}{"counts"},
			  $total_prediction_hash{$currentMethod}{"Group-Level"}{"RefOGs_affected_by_fissions"}{"counts"},
			  $total_prediction_hash{$currentMethod}{"Group-Level"}{"Accurately_predicted_RefOGs"}{"counts"},
		  );
		if ( !defined($RefOfByFusion) || $RefOfByFusion eq '' )
		  {
			$RefOfByFusion = 0;
		  }

#$lines2print .= "$currentMethod\t$FN\t$FP\t$Fissions\t$Fusion\t$AccPredictedRefOGs\t$RefOfByFN\t$RefOfByFP\t$RefOfByFusion\t$RefOfByFissions\t$AccPredictedRefOGsGroup\t0\t0\n";
	  }

	# add hieranoid and other
	foreach my $predictionMethod (@predictionMethods_array)
	  {
		my ( $FN, $FP, $Fissions, $Fusion, $AccPredictedRefOGs, $RefOfByFN, $RefOfByFP, $RefOfByFusion,
			 $RefOfByFissions, $AccPredictedRefOGsGroup, $totalTP, $TPPercentage )
		  = (
			  (
				defined( $helperHash{$predictionMethod}{totalMissingGenes} )
				? $helperHash{$predictionMethod}{totalMissingGenes}
				: 0
			  ),
			  (
				defined( $helperHash{$predictionMethod}{totalfalsePredictions} )
				? $helperHash{$predictionMethod}{totalfalsePredictions}
				: 0
			  ),
			  (
				defined( $helperHash{$predictionMethod}{totalNumberOfFission} )
				? $helperHash{$predictionMethod}{totalNumberOfFission}
				: 0
			  ),
			  (
				defined( $helperHash{$predictionMethod}{totalNumberOfFusion} )
				? $helperHash{$predictionMethod}{totalNumberOfFusion}
				: 0
			  ),
			  (
				defined( $helperHash{$predictionMethod}{noOfCorrectGroups} )
				? $helperHash{$predictionMethod}{noOfCorrectGroups}
				: 0
			  ),
			  (
				defined( $helperHash{$predictionMethod}{noOfMissingGeneAffectedGroup} )
				? $helperHash{$predictionMethod}{noOfMissingGeneAffectedGroup}
				: 0
			  ),
			  (
				defined( $helperHash{$predictionMethod}{noOfFalseGeneAffectedGroup} )
				? $helperHash{$predictionMethod}{noOfFalseGeneAffectedGroup}
				: 0
			  ),
			  (
				defined( $helperHash{$predictionMethod}{GroupsAffectedByFusion} )
				? $helperHash{$predictionMethod}{GroupsAffectedByFusion}
				: 0
			  ),
			  (
				defined( $helperHash{$predictionMethod}{GroupsAffectedByFissions} )
				? $helperHash{$predictionMethod}{GroupsAffectedByFissions}
				: 0
			  ),
			  (
				defined( $helperHash{$predictionMethod}{AccPredictedRefOGs_group} )
				? $helperHash{$predictionMethod}{AccPredictedRefOGs_group}
				: 0
			  ),
			  (
				defined( $helperHash{$predictionMethod}{total_true_positives} )
				? $helperHash{$predictionMethod}{total_true_positives}
				: 0
			  ),
			  (
				defined( $helperHash{$predictionMethod}{total_true_positives} )
				? ( ( $helperHash{$predictionMethod}{total_true_positives} * 100 / 1638 ) )
				: 0
			  )
		  );
		$lines2print .=
"new$predictionMethod\t$FN\t$FP\t$Fissions\t$Fusion\t$AccPredictedRefOGs\t$RefOfByFN\t$RefOfByFP\t$RefOfByFusion\t$RefOfByFissions\t$AccPredictedRefOGsGroup\t$totalTP\t$TPPercentage\n";
	  }
	write_to_file( { file_name => $OrthoBenchRawData, text => $lines2print } );
	exit;
	print "Total\t$totalMissingGenes\t$totalfalsePredictions\n";

	#print Dumper %overallStatistics;
	# Accurately predicted RefOGs
	# Erroneously assigned genes
	# Missing genes
	# RefOGs affected by erroneously affected genes
	my $numberOfGroups = keys(%overallStatistics);
	print "Accurately predicted RefOGs\t$noOfCorrectGroups\t"
	  . ( $noOfCorrectGroups * 100 ) / $numberOfGroups . "\n";
	print "Erroneously assigned genes\t$totalfalsePredictions\t\n";
	print "Missing genes\t$totalMissingGenes\t\n";
	print "Affected genes\t$noOfMissingGeneAffectedGroup\t"
	  . ( $noOfMissingGeneAffectedGroup * 100 ) / $numberOfGroups . "\n";
	print "RefOGs affected by erroneously affected genes\t$noOfFalseGeneAffectedGroup\t"
	  . ( $noOfFalseGeneAffectedGroup * 100 ) / $numberOfGroups . "\n";
  }

sub compare2orthobench_pairwise
  {

	# 1. read all hieranoid predictions
	## ALL PREDICTION FILES
#my $hieranoid_prediction_file = "/Users/fab/Documents/documents/work/current_projects/hieranoid/archived_results/results_ensembl/nodes/Bilateria/Bilateria.expandedGroups.txt";
#my $hieranoid_prediction_consensus_outgroup_file = "/Users/fab/Documents/documents/work/current_projects/hieranoid/results_orthobench_consensus_outgroup/nodes/Bilateria/Bilateria.expandedGroups.txt";
#my $hieranoid_prediction_test_file = "/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/test_10genes.txt";
#my $hieranoid_profile_prediction_file = "/Users/fab/Documents/documents/work/current_projects/hieranoid/archived_results/results_ensembl_profile/nodes/Bilateria/Bilateria.expandedGroups.txt";
#my $hieranoid_profile_quickdirty_prediction_file = "/Users/fab/Documents/documents/work/current_projects/hieranoid/archived_results/results_ensembl_profile_quickdirty/nodes/Bilateria/Bilateria.expandedGroups.txt";
	my $hieranoid_prediction_consensus_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/publication_results/results_orthobench_consensus/nodes/Bilateria/Bilateria.expandedGroups.txt";
	my $hieranoid_prediction_consensus_outgroup_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/publication_results/results_orthobench_consensus_outgroup/nodes/Bilateria/Bilateria.expandedGroups.txt";
	my $hieranoid_prediction_profile_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/publication_results/results_orthobench_profile/nodes/Bilateria/Bilateria.expandedGroups.txt";
	my $hieranoid_prediction_profile_outgroup_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/publication_results/results_orthobench_profile_outgroup/nodes/Bilateria/Bilateria.expandedGroups.txt";
	my $hieranoid_Smallprediction_profile_outgroup_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/publication_results/results_orthobenchSmall_profile_outgroup_nousort/nodes/Bilateria/Bilateria.expandedGroups.txt";
	my $hieranoid_profile_prediction_outgroup_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/results_orthobench_profile_outgroup/nodes/Bilateria/Bilateria.expandedGroups.txt";
	my $AllUsedIDsFile =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/Ensembl60_allheaders.txt";
	my $oma_prediction_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/oma-groups_may2011_converted.txt";
	my $eggNOG_prediction_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/meNOG.mapping.txt";
	my $orthomcl_prediction_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/groups_OrthoMCL-5.txt_converted";
	my $treeFamOldPredictionFile =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/fam_genes.txt.table";
	my $treeFam_prediction_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/fam_genes.txt.table_converted.txt";
	my $inparanoid_prediction_file =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/sqltable.all";
	my $orthodb_prediction_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/orthodb4_Metazoa_tabtext.txt";
	my $ensembl_mapping =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/ensembl_ids.txt";
	my $hieranoid_fasta_file =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/ensemblv60_12species";
	my $orthobench_fasta_files =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/predictions/orthobench/";
	my $orthobench_predictions =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/pictures/R_sources/OrthoBenchReferenceTable.txt";

#my $orthobench_groups = "/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/IDs_OrthoBench2.txt";
	my $orthobench_groups =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/IDs_OrthoBench2.txt";
	my $OrthoBenchRawData =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/pictures/R_sources/OrthoBenchRawData.txt";
	my $OrthoBenchRawDataPairwise =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/pictures/R_sources/OrthoBenchRawDataPairwise.csv";
	my %ID2HieranoidOG;
	my %ID2DatabaseOG;
	my %allPredictions;

	my @predictionMethods_array =
	  qw(Hieranoid_consensus Hieranoid_consensus_outgroup Hieranoid_profile Hieranoid_profile_outgroup  Inparanoid OMA eggNOG OrthoMCL TreeFam  OrthoDB);

	#my @predictionMethods_array = qw(TreeFam);

#print Dumper @predictionMethods_array;
#convertTreeFamMapping2OrthoBench($treeFamOldPredictionFile,\%{$allPredictions{"TreeFam_new"}}, \%{$ID2DatabaseOG{"TreeFam_new"}},$ensembl_mapping);
#exit;
	my %ID2refOG;

	#  'refOG0001' -> @(ENSTNIP00000005325,ENSDARP00000071531)
	my %refOG_hash;
	my %refOGGenes_hash;

	#  {'refOG0001'}{'ENSTNIP00000005325'} = 2
	#my %Pair2DatabaseOG;
	# Save all pairs as
	# {group}{gene1}{"ID1"}
	my %refOG2Members;
	my %pair2RefOGGene_hash;
	my %ID2Gene_hash;
	## OrthoBench groups
	print "\tReading OrthoBench data...";
	getOrthoBenchGroupsPairwise(
								 {
								   orthobench_groups_file => $orthobench_groups,
								   ID2refOG_href          => \%ID2refOG,
								   refOG_href             => \%refOG_hash,
								   refOGGenes_href        => \%refOGGenes_hash,
								   refOG2Members_href     => \%refOG2Members,
								   pair2RefOGGene_href    => \%pair2RefOGGene_hash,
								   ID2Gene_href           => \%ID2Gene_hash
								 }
	);
	print "done\n";
	if ( !keys %ID2refOG )
	  {
		print "\tCould not get ID2refOGs. Exiting\n";
		exit;
	  }
	if ( !keys %refOG_hash )
	  {
		print "\tCould not get refOG_hash. Exiting\n";
		exit;
	  }
	if ( !keys(%refOGGenes_hash) )
	  {
		print("Empty hash  refOGGenes_hash found\n");
		exit;
	  }
	if ( !keys(%refOG2Members) )
	  {
		print("Empty hash refOG2Members found\n");
		exit;
	  }
	if ( !keys(%ID2Gene_hash) )
	  {
		print("Empty hash ID2Gene_hash found\n");
		exit;
	  }
	foreach my $predictionMethod (@predictionMethods_array)
	  {
		print "\tReading $predictionMethod ...\n";
		if ( $predictionMethod eq 'OrthoDB' )
		  {
			getOrthoDBPredictionsPairwise(
										   {
											 predictionFile     => $orthodb_prediction_file,
											 AllPredictionsHash => \%{ $allPredictions{"OrthoDB"} },
											 ID2DatabaseOG      => \%{ $ID2DatabaseOG{"OrthoDB"} },
											 ID2GeneHash        => \%ID2Gene_hash
										   }
			  )

			  #next;
		  }
		if ( $predictionMethod eq 'Inparanoid' )
		  {
			getInparanoidPredictionsPairwise(
											  {
												predictionFile     => $inparanoid_prediction_file,
												AllPredictionsHash => \%{ $allPredictions{"Inparanoid"} },
												ID2DatabaseOG      => \%{ $ID2DatabaseOG{"Inparanoid"} },
												ID2GeneHash        => \%ID2Gene_hash,
												AllUsedIDsFile     => $AllUsedIDsFile
											  }
			  )

			  #next;
		  }
		if ( $predictionMethod eq 'Hieranoid_consensus' )
		  {
			getPredictionMethodPredictionsPairwise(
											 {
											   predictionFile     => $hieranoid_prediction_consensus_file,
											   AllPredictionsHash => \%{ $allPredictions{"Hieranoid_consensus"} },
											   ID2DatabaseOG      => \%{ $ID2DatabaseOG{"Hieranoid_consensus"} },
											   ID2GeneHash        => \%ID2Gene_hash
											 }
			  )

			  #next;
		  }
		if ( $predictionMethod eq 'Hieranoid_profile' )
		  {
			getPredictionMethodPredictionsPairwise(
											   {
												 predictionFile     => $hieranoid_prediction_profile_file,
												 AllPredictionsHash => \%{ $allPredictions{"Hieranoid_profile"} },
												 ID2DatabaseOG      => \%{ $ID2DatabaseOG{"Hieranoid_profile"} },
												 ID2GeneHash        => \%ID2Gene_hash
											   }
			  )

			  #next;
		  }
		if ( $predictionMethod eq 'Hieranoid_consensus_outgroup' )
		  {
			getPredictionMethodPredictionsPairwise(
									{
									  predictionFile     => $hieranoid_prediction_consensus_outgroup_file,
									  AllPredictionsHash => \%{ $allPredictions{"Hieranoid_consensus_outgroup"} },
									  ID2DatabaseOG      => \%{ $ID2DatabaseOG{"Hieranoid_consensus_outgroup"} },
									  ID2GeneHash        => \%ID2Gene_hash
									}
			  )

			  #next;
		  }
		if ( $predictionMethod eq 'Hieranoid_profile_outgroup' )
		  {
			getPredictionMethodPredictionsPairwise(
									  {
										predictionFile     => $hieranoid_prediction_profile_outgroup_file,
										AllPredictionsHash => \%{ $allPredictions{"Hieranoid_profile_outgroup"} },
										ID2DatabaseOG      => \%{ $ID2DatabaseOG{"Hieranoid_profile_outgroup"} },
										ID2GeneHash        => \%ID2Gene_hash
									  }
			  )

			  #next;
		  }

		if ( $predictionMethod eq 'Hieranoid_profile_outgroupSmall' )
		  {
			getPredictionMethodPredictionsPairwise(
								 {
								   predictionFile     => $hieranoid_Smallprediction_profile_outgroup_file,
								   AllPredictionsHash => \%{ $allPredictions{"Hieranoid_profile_outgroupSmall"} },
								   ID2DatabaseOG      => \%{ $ID2DatabaseOG{"Hieranoid_profile_outgroupSmall"} },
								   ID2GeneHash        => \%ID2Gene_hash
								 }
			  )

			  #next;
		  }
		if ( $predictionMethod eq 'OMA' )
		  {
			getPredictionMethodPredictionsPairwise(
													{
													  predictionFile     => $oma_prediction_file,
													  AllPredictionsHash => \%{ $allPredictions{"OMA"} },
													  ID2DatabaseOG      => \%{ $ID2DatabaseOG{"OMA"} },
													  ID2GeneHash        => \%ID2Gene_hash
													}
			  )

			  #next;
		  }
		if ( $predictionMethod eq 'eggNOG' )
		  {
			geteggNOGPredictionsPairwise(
										  {
											predictionFile     => $eggNOG_prediction_file,
											AllPredictionsHash => \%{ $allPredictions{"eggNOG"} },
											ID2DatabaseOG      => \%{ $ID2DatabaseOG{"eggNOG"} },
											ID2GeneHash        => \%ID2Gene_hash
										  }
			  )

			  #next;
		  }
		if ( $predictionMethod eq 'TreeFam' )
		  {
			getPredictionMethodPredictionsPairwise(
													{
													  predictionFile     => $treeFam_prediction_file,
													  AllPredictionsHash => \%{ $allPredictions{"TreeFam"} },
													  ID2DatabaseOG      => \%{ $ID2DatabaseOG{"TreeFam"} },
													  ID2GeneHash        => \%ID2Gene_hash
													}
			);
			print "\tallpredictionhash\n";
			print Dumper %{ $allPredictions{"TreeFam"} };
			print "\tid2database\n";
			print Dumper $ID2DatabaseOG{"TreeFam"};

			#next;
		  }
		if ( $predictionMethod eq 'OrthoMCL' )
		  {
			getPredictionMethodPredictionsPairwise(
													{
													  predictionFile     => $orthomcl_prediction_file,
													  AllPredictionsHash => \%{ $allPredictions{"OrthoMCL"} },
													  ID2DatabaseOG      => \%{ $ID2DatabaseOG{"OrthoMCL"} },
													  ID2GeneHash        => \%ID2Gene_hash
													}
			  )

			  #next;
		  }

		#if($predictionMethod eq 'test'){
		#	getPredictionMethodPredictionsPairwise({predictionFile =>$hieranoid_prediction_test_file,
		#					AllPredictionsHash=> \%{$allPredictions{"test"}},
		#					ID2DatabaseOG => \%{$ID2DatabaseOG{"test"}},
		#					ID2GeneHash => \%ID2Gene_hash})
		#next;
		#}
		if ( !keys %{ $allPredictions{$predictionMethod} } )
		  {
			print "\tCould not get $predictionMethod. Exiting\n";
			exit;
		  }
		if ( !keys %{ $ID2DatabaseOG{$predictionMethod} } )
		  {
			print "\tCould not get $predictionMethod ID2Group_hash. Exiting\n";
			exit;
		  }
	  }
	print "\tfinished reading\n";

	#exit;
	my %predictions_hash;
	my %total_prediction_hash;

#&readOrthoBenchPredictions({file => $orthobench_predictions, hash => \%predictions_hash, totalhashref => \%total_prediction_hash});
	my %overallStatistics;
	my %OGStatistics;
	my $totalNumberCheckedRefOGGenes       = 0;
	my $totalNumberCheckedHieranoidOGGenes = 0;
	my $totalNumberOfFission               = 0;
	my $totalNumberOfFusion                = 0;
	my $GroupsAffectedByFissions           = 0;
	my $GroupsAffectedByFusion             = 0;
	my %helperHash;
	my $groupCount = 0;
	my $debugprint = 1;
	my %countMissingSpeciesIDs;

	# For each REFOG
	my $fusionDetected4RefOG = 0;

	#my @RefOGGenes = (1 .. 1638);
	#my @combinations = combine(2,@RefOGGenes);
	my $totalNoOfRefOGPairs = 2681406;
	my $All_P_RefOGPairs    = 0;
	my $All_N_RefOGPairs    = 0;

	# Count all Positives in Orthobench dataset
	foreach my $currentRefOG ( sort keys(%refOG2Members) )
	  {

		#next if $currentRefOG ne 'RefOG041';
		my @GenePairs = combine( 2, @{ $refOG2Members{$currentRefOG} } );
		$All_P_RefOGPairs += scalar(@GenePairs);
	  }

	# Count all true negatives
	my @allOrthoBenchGroups = keys(%refOG2Members);

	my %countedMissingGenes;
	my $wrongIDs;
	my $TN;

	# total number of negatives
	# total number of genes
	my $total_number_of_genes  = 1638;
	my $total_number_genepairs = 1340703;
	my $total_number_positives = 0;
	### REFOGS

	#
	#		#list of other groups
	#		for ( my $i = 0 ; $i < scalar(@allOrthoBenchGroups) ; $i++ )
	#		  {
	#
	#			#next if $allOrthoBenchGroups[$i] ne 'RefOG041';
	#			my @firstgroup_members = @{ $refOG2Members{ $allOrthoBenchGroups[$i] } };
	#			for ( my $j = $i + 1 ; $j < scalar(@allOrthoBenchGroups) ; $j++ )
	#			  {
	#				my @secondGroupMembers = @{ $refOG2Members{ $allOrthoBenchGroups[$j] } };
	#				if ( !@secondGroupMembers )
	#				  {
	#					die "group " . $allOrthoBenchGroups[$j] . " contains no members\n";
	#				  }
	#
	#				# Iterate over all first group members
	#				for ( my $groupA = 0 ; $groupA < scalar(@firstgroup_members) ; $groupA++ )
	#				  {
	#
	#					# Iterate over all first group members
	#					for ( my $groupB = 0 ; $groupB < scalar(@secondGroupMembers) ; $groupB++ )
	#					  {
	#						my $pairDestiny;
	#						my ( $IDA, $IDB ) = ( $firstgroup_members[$groupA], $secondGroupMembers[$groupB] );
	#						print "\t\t$IDA - $IDB: " if $debugprint;
	#						my $group4IDA = $ID2DatabaseOG{$predictionMethod}{$IDA};
	#						my $group4IDB = $ID2DatabaseOG{$predictionMethod}{$IDB};
	#						check_RefOGgroup4TN(
	#											 {
	#											   IDA                    => $IDA,
	#											   IDB                    => $IDB,
	#											   group4IDA              => $group4IDA,
	#											   group4IDB              => $group4IDB,
	#											   pairDestiny            => \$pairDestiny,
	#											   countedMissingGenes    => \%countedMissingGenes,
	#											   countMissingSpeciesIDs => \%countMissingSpeciesIDs,
	#											   predictionMethod       => $predictionMethod,
	#											   TN                     => \$TN,
	#											   wrongIDs               => \$wrongIDs
	#											 }
	#						);
	#						print "-> $pairDestiny\n" if $debugprint;
	#					  }
	#
	#				  }
	#
	#				#my @GenePairs = combine( 2, @secondGroupMembers );
	#				#$All_N_RefOGPairs += scalar(@group_members) * scalar(@GenePairs);
	#			  }
	#		  }
	#$overallStatistics{$predictionMethod}{"TN"} += $TN;
	foreach my $predictionMethod (@predictionMethods_array)
	  {
		### REFOGS
		foreach my $currentRefOG ( sort keys(%refOG2Members) )
		  {

			#next;
			#next if $currentRefOG ne 'RefOG041';
			print "checking $currentRefOG ($predictionMethod)\n" if $debugprint;
			my $noOfMissingGenes = 0;
			my $noOfGenePairs    = 0;
			my $groupMismatches  = 0;
			my $noMatchingPairs  = 0;
			my $wrongIDs         = 0;
			print "genes are " . join( ",", @{ $refOG2Members{$currentRefOG} } ) . "\n";
			my @GenePairs = combine( 2, sort @{ $refOG2Members{$currentRefOG} } );
			my $noOfPairs = scalar(@GenePairs);
			$overallStatistics{$predictionMethod}{"TP"}             += $noOfPairs;
			$overallStatistics{$predictionMethod}{"totalPositives"} += $noOfPairs;

			#$overallStatistics{"All_TP_RefOGPairs"} += $noOfPairs;
			$overallStatistics{$predictionMethod}{"RefOGAllPairs"} = $totalNoOfRefOGPairs;
			$OGStatistics{$predictionMethod}{$currentRefOG}{"RefOGPairs"} = $noOfPairs;
			my %countedMissingGenes;
			my $FN = 0;
			my $TP = 0;
		  PAIR:
			for ( my $i = 0 ; $i < scalar(@GenePairs) ; $i++ )
			  {
				my $pairDestiny = "";
				my ( $IDA, $IDB ) = ( $GenePairs[$i][0], $GenePairs[$i][1] );
				my $group4IDA = $ID2DatabaseOG{$predictionMethod}{$IDA};
				my $group4IDB = $ID2DatabaseOG{$predictionMethod}{$IDB};

				#print "\t\t$IDA ($group4IDA) - $IDB ($group4IDB): " if $debugprint;
				check_RefOGgroup(
								  {
									IDA                    => $IDA,
									IDB                    => $IDB,
									group4IDA              => $group4IDA,
									group4IDB              => $group4IDB,
									pairDestiny            => \$pairDestiny,
									countedMissingGenes    => \%countedMissingGenes,
									countMissingSpeciesIDs => \%countMissingSpeciesIDs,
									predictionMethod       => $predictionMethod,
									FN                     => \$FN,
									wrongIDs               => \$wrongIDs,
									TP                     => \$TP
								  }
				);
				print "-> $pairDestiny\n" if $debugprint;
			  }
			print "There were $noOfMissingGenes missing pairs (" . scalar(@GenePairs) . ")\n" if $debugprint;
			#$overallStatistics{$predictionMethod}{"TP"} += $TP;
			$overallStatistics{$predictionMethod}{"FN"} += $FN;
		  }

		#print Dumper %countMissingSpeciesIDs;
		#exit;
		print "$predictionMethod\n";
		my $ogCounter = 0;

		#$debugprint = 1;
		my $allgroups = keys( %{ $allPredictions{$predictionMethod} } );

		#print Dumper $allPredictions{$predictionMethod};
		my @TN_counter;
		my $TN_groupCounter = 0;
		my $FP              = 0;
		my $TN              = 0;
		my $TP              = 0;
		
		## OTHER METHODS
		foreach my $currentPredictionMethodOG ( keys( %{ $allPredictions{$predictionMethod} } ) )
		  {
			print "\t$currentPredictionMethodOG ($ogCounter / $allgroups)\n";
			my @predictionMethodOGGroupMember =
			  split( ",", $allPredictions{$predictionMethod}{$currentPredictionMethodOG} );
			print "members: " . join( ",", @predictionMethodOGGroupMember ) . "\n";
			if ( !scalar(@predictionMethodOGGroupMember) )
			  {
				print "\tProblem reading hieranoid group member for group $currentPredictionMethodOG\n";
			  }

			# iterate over HieranoidOG
			# {'HieranoidOG1'}{ID1} -> 1
			#        {ID2} -> 2
			#        {ID3} -> 3
			# {'HieranoidOG2'}{ID4}
			my $wrongIDs        = 0;
			my $groupMismatches = 0;

			# Get all pairs for OG
			next if @predictionMethodOGGroupMember < 2;
			my @GenePairs = combine( 2, @predictionMethodOGGroupMember );
			my $noOfPairs = scalar(@GenePairs);
			push( @{ $TN_counter[$ogCounter] }, @predictionMethodOGGroupMember );
			$ogCounter++;
			$overallStatistics{$predictionMethod}{"PreOGPairs"} += $noOfPairs;
			$OGStatistics{$predictionMethod}{$currentPredictionMethodOG}{"PreOGPairs"} = $noOfPairs;
		  PAIR:
			for ( my $i = 0 ; $i < scalar(@GenePairs) ; $i++ )
			  {
				my ( $IDA, $IDB ) = ( $GenePairs[$i][0], $GenePairs[$i][1] );
				print "\t\t$IDA - $IDB: " if $debugprint;
				my $group4A;
				my $group4B;
				my $pairDestiny = "";
				check_group(
							 {
							   IDA         => $IDA,
							   IDB         => $IDB,
							   pairDestiny => \$pairDestiny,
							   FP          => \$FP,
							   wrongIDs    => \$wrongIDs,
							   TP          => \$TP,
							   
							 }
				);
				print "-> $pairDestiny\n" if $debugprint;
			  }

		  }
		$overallStatistics{$predictionMethod}{"TP"} = $TP;
		$overallStatistics{$predictionMethod}{"FP"} = $FP;

		#$overallStatistics{$predictionMethod}{"TN"}    = $TN;
		$overallStatistics{$predictionMethod}{"NoOGs"} = scalar(@TN_counter);
	  }

	#print Dumper %overallStatistics;
	my $lines2print =
"Method\tFalsePositiveRate\tallnegatives\tFP\tFalseNegativeRate\tFN\tall_positives\tfalseDiscoveryRate\tTP\tTotalP\tFP+TP\tOGs\n";

	foreach my $predictionMethod (@predictionMethods_array)
	  {

# sick people = ortholog
#True positive (TP): Sick people correctly diagnosed as sick => orthologs (RefOG) predicted as orthologs (hieranoid)
# e.g. HiG1: OG03_01, OG03_02, OG04_02
# OG03_01 - OG03_02 => TP
#False positive (FP): Healthy people incorrectly identified as sick => non-orthologs (RefOG) predicted as orthologs (hieranoid)
# OG03_01 - OG04_02 => FP
# OG03_02 - OG04_02 => FP
#True negative (TN): Healthy people correctly identified as healthy => non-orthologs (RefOG) predicted as non-orthologs (hieranoid)
# check all Hi1: OG04
#False negative (FN): Sick people incorrectly identified as healthy => orthologs (RefOG) predicted as non-orthologs (hieranoid)
# e.g. RefOG: OG04_01 - OG04_02 -> Hieranoid: group4 - group3

		my $RefOGAllPairs = $overallStatistics{$predictionMethod}{"RefOGAllPairs"};
		## Total P
		my $RefOGtotalPositives = $All_P_RefOGPairs;
		## Total N
		my $RefOGtotalNegatives = $All_N_RefOGPairs;
		my $FP                  = int( $overallStatistics{$predictionMethod}{"FP"} );
		my $FN                  = int( $overallStatistics{$predictionMethod}{"FN"} );
		my $TP                  = int( $overallStatistics{$predictionMethod}{"TP"} );
		my $all_positives       = $overallStatistics{$predictionMethod}{"totalPositives"};
		my $all_negative_pairs  = $total_number_genepairs - $all_positives;

   # all P = all pairs - all negatives
   # all N = all pairs - all positives
   #   False positive rate (α) = type I error = FP / (FP + TN) = 180 / (180 + 1820) = 9% = 1 − specificity
		my $FalsePositiveRate = (
								  ( !defined($all_negative_pairs) || ($all_negative_pairs) == 0 )
								  ? 0
								  : $FP / ($all_negative_pairs)
		);

	#   False negative rate (β) = type II error = FN / (TP + FN) = 10 / (20 + 10) = 33% = 1 − sensitivity
		my $FalseNegativeRate = (
								  ( !defined($RefOGtotalPositives) || ($RefOGtotalPositives) == 0 )
								  ? 0
								  : $FN / ($all_positives)
		);
		my $falseDiscoveryRate = ($FP / ($FP + $TP));
		my $fptp = $FP + $TP;

		#my $fntn = $FN + $TN;
		my $noOGs = $overallStatistics{$predictionMethod}{"NoOGs"};

		# REFOGs
		$lines2print .=
"$predictionMethod\t$FalsePositiveRate\t$all_negative_pairs\t$FP\t$FalseNegativeRate\t$FN\t$all_positives\t$falseDiscoveryRate\t$TP\t$RefOGtotalPositives\t$fptp\t$noOGs\n";

#"$RefOGAllPairs\t$RefOGMatchingPairs\t$RefOGMissingPairs\t$RefOGMissingPairs_WrongIDs\t$RefOGMissingPairs_GroupMismatches\t$RefOGMissingPairsQuote\t";
#$lines2print .= "$PreOGPairs\t$PreOGTruePairs\t$PreOGTruePairsQuote\t$PreOGFalsePairs\t$PreOGFalsePairs_WrongIDs\t$PreOGFalsePairs_GroupMismatches\t$PreOGFalsePairsQuote\n";
	  }
	write_to_file( { file_name => $OrthoBenchRawDataPairwise, text => $lines2print } );
	return 1;
  }

sub check_RefOGgroup
  {

	#print Dumper @_;
	#exit;
	my ($arg_ref)                   = (@_);
	my $IDA                         = $arg_ref->{'IDA'};
	my $IDB                         = $arg_ref->{'IDB'};
	my $group4IDA                   = $arg_ref->{'group4IDA'};
	my $group4IDB                   = $arg_ref->{'group4IDB'};
	my $pairDestiny                 = $arg_ref->{'pairDestiny'};
	my $countedMissingGenes_href    = $arg_ref->{'countedMissingGenes'};
	my $countMissingSpeciesIDs_href = $arg_ref->{'countMissingSpeciesIDs'};
	my $predictionMethod            = $arg_ref->{'predictionMethod'};
	my $FN                          = $arg_ref->{'FN'};
	my $wrongIDs                    = $arg_ref->{'wrongIDs'};

	# Check
	if ( !defined($group4IDA) )
	  {
		$$pairDestiny .= "Missing gene (no prediction for IDA $IDA)";
		my $Species4ID = $IDA;
		$Species4ID =~ s/\d+//g;
		if (    !exists $countedMissingGenes_href->{$IDA}
			 && !exists $countedMissingGenes_href->{$IDB} )
		  {
			$countMissingSpeciesIDs_href->{$predictionMethod}{$Species4ID}++;

			#$noOfMissingGenes++;
			$$wrongIDs++;
			$$FN++;
			return 0;
		  }
		else
		  {
			$$pairDestiny .= " (not counted)";
		  }
	  }
	elsif ( !defined($group4IDB) )
	  {
		$$pairDestiny .= "Missing gene (no prediction for IDB $IDB)";
		if (    !exists $countedMissingGenes_href->{$IDA}
			 && !exists $countedMissingGenes_href->{$IDB} )
		  {

			#$noOfMissingGenes++;
			$$wrongIDs++;
			$$FN++;
			return 0;

			#next PAIR;
		  }
		else
		  {
			$$pairDestiny .= " (not counted)";
			return 0;
		  }
	  }
	elsif ( $group4IDA eq $group4IDB )
	  {
		$$pairDestiny .= "TP (both map to same gene: $group4IDA";

		#$noMatchingPairs++;
		#$TP++;
		return 1;
	  }
	elsif ( $group4IDA ne $group4IDB )
	  {
		$$pairDestiny .= "Missing gene (diff group $group4IDA - $group4IDB)";
		if (    !exists $countedMissingGenes_href->{$IDA}
			 && !exists $countedMissingGenes_href->{$IDB} )
		  {

			#$noOfMissingGenes++;
			#$groupMismatches++;
			$$FN++;
			return 0;

			#next PAIR;
		  }
		else
		  {
			$$pairDestiny .= " (not counted)";
			return 0;
		  }
	  }

  }

sub check_RefOGgroup4TN
  {

	#print Dumper @_;
	#exit;
	my ($arg_ref)                   = (@_);
	my $IDA                         = $arg_ref->{'IDA'};
	my $IDB                         = $arg_ref->{'IDB'};
	my $group4IDA                   = $arg_ref->{'group4IDA'};
	my $group4IDB                   = $arg_ref->{'group4IDB'};
	my $pairDestiny                 = $arg_ref->{'pairDestiny'};
	my $countedMissingGenes_href    = $arg_ref->{'countedMissingGenes'};
	my $countMissingSpeciesIDs_href = $arg_ref->{'countMissingSpeciesIDs'};
	my $predictionMethod            = $arg_ref->{'predictionMethod'};
	my $TN                          = $arg_ref->{'TN'};
	my $wrongIDs                    = $arg_ref->{'wrongIDs'};

	# Check
	if ( !defined($group4IDA) )
	  {
		$$pairDestiny .= "Missing gene (no prediction for IDA $IDA)";
		my $Species4ID = $IDA;
		$Species4ID =~ s/\d+//g;
		if (    !exists $countedMissingGenes_href->{$IDA}
			 && !exists $countedMissingGenes_href->{$IDB} )
		  {
			$countMissingSpeciesIDs_href->{$predictionMethod}{$Species4ID}++;

			#$noOfMissingGenes++;
			$wrongIDs++;
			$$TN++;
			return 0;
		  }
		else
		  {
			$$pairDestiny .= " (not counted)";
		  }
	  }
	elsif ( !defined($group4IDB) )
	  {
		$$pairDestiny .= "Missing gene (no prediction for IDB $IDB)";
		if (    !exists $countedMissingGenes_href->{$IDA}
			 && !exists $countedMissingGenes_href->{$IDB} )
		  {

			#$noOfMissingGenes++;
			$wrongIDs++;
			$$TN++;
			return 0;

			#next PAIR;
		  }
		else
		  {
			$$pairDestiny .= " (not counted)";
			return 0;
		  }
	  }
	elsif ( $group4IDA eq $group4IDB )
	  {
		$$pairDestiny .= "FP (both map to same gene: $group4IDA";

		#$noMatchingPairs++;
		#$TP++;
		#$FP++;
		return 1;
	  }
	elsif ( $group4IDA ne $group4IDB )
	  {
		$$pairDestiny .= "TN(diff group $group4IDA - $group4IDB)";
		if (    !exists $countedMissingGenes_href->{$IDA}
			 && !exists $countedMissingGenes_href->{$IDB} )
		  {

			#$noOfMissingGenes++;
			#$groupMismatches++;
			$$TN++;

			return 0;

			#next PAIR;
		  }
		else
		  {
			$$pairDestiny .= " (not counted)";
			return 0;
		  }
	  }

  }

#===================================================================
sub check_group
  {

	#print Dumper @_;
	#exit;
	my ($arg_ref)   = (@_);
	my $IDA         = $arg_ref->{'IDA'};
	my $IDB         = $arg_ref->{'IDB'};
	my $pairDestiny = $arg_ref->{'pairDestiny'};
	my $FP          = $arg_ref->{'FP'};
	my $wrongIDs    = $arg_ref->{'wrongIDs'};
	my $TP			= $arg_ref->{'TP'};
	my ( $group4A, $group4B );

	# IDA has RefGroup prediction ?
	if ( $IDA =~ /(RefOG\d+)_\d+/ )
	  {
		$group4A = $1;
	  }

	# IDB has RefGroup prediction ?
	if ( $IDB =~ /(RefOG\d+)_\d+/ )
	  {
		$group4B = $1;
	  }

	# Check
	if ( !defined($group4A) )
	  {
		$$pairDestiny .= "FP (no group for IDA $IDA)";
		$$wrongIDs++;
		$$FP++;
		return 0;
	  }
	if ( !defined($group4B) )
	  {
		$$pairDestiny .= "FP (no group for IDB $IDB)";
		$$wrongIDs++;
		$$FP++;
		return 0;
	  }
	if ( $group4A eq $group4B )
	  {
		$$pairDestiny .= "TP (both map to same gene: $group4A";
		$$TP++;
		return 1;
	  }
	if ( $group4A ne $group4B )
	  {
		$$pairDestiny .= "FP (diff group $group4A - $group4B)";
		$$FP++;
		return 0;
	  }

	return 0;

  }

sub testInparanoidTime
  {

	# Goal: Find best Ublast parameters
	# Test Human - 9 other species
	# Compare Hieranoid pairwise predictions with Inparanoid pairwise predictions
	# How: Hieranoid run with only pairs of species
	my $fastaDir = "testdata/ensemblv60_12species/longestTranscript/";
	my $predDir  = "predictions/InparanoidHieranoidComparison/";
	my $timeFile = "$predDir/time";
	mkdir($predDir) if !-e $predDir;

	#my $usearch = "/Users/fab/bin/exe/usearch4.1.93_i86darwin32";
	my $usearch = $Configuration::usearch;
	my @species = (
		"Homo_sapiens",
		"Danio_rerio",
		"Canis_familiaris",
		"Pan_troglodytes",
		"Ciona_intestinalis",
		"Caenorhabditis_elegans",

		#"Mus_musculus",
	);
	my $currSpeciesFile;
	my $start_time_inparanoid = new Benchmark;
	for ( my $i = 0 ; $i < scalar(@species) ; $i++ )
	  {
		my $speciesA         = $species[$i];
		my $speciesA_file    = "$fastaDir/$speciesA.fa";
		my $noQuerySequences = `grep ">" $speciesA_file -c`;
		chomp($noQuerySequences);
		for ( my $j = $i + 1 ; $j < scalar(@species) ; $j++ )
		  {
			my $speciesB      = $species[$j];
			my $speciesB_file = "$fastaDir/$speciesB.fa";
			my $noDBSequences = `grep ">" $speciesB_file -c`;
			chomp($noDBSequences);
			my $inparanoid_call = "perl inparanoid_runtime.pl -q $speciesA_file -d $speciesB_file";
			my $start_time_inparanoid_comparison = new Benchmark;
			print $inparanoid_call. "\n";

			#system($inparanoid_call);
			# Stop time
			my $end_time_hieranoid_comparison = new Benchmark;
			my $time_diff_pair = timediff( $end_time_hieranoid_comparison, $start_time_inparanoid_comparison );
			my $debugString    = "Inparanoid\tPair: $speciesA ($noQuerySequences) - $speciesB ($noDBSequences) "
			  . timestr( $time_diff_pair, 'all' );
			&attach_to_file( { file_name => $timeFile, text => $debugString . "\n" } );
		  }
	  }
	my $end_time_inparanoid = new Benchmark;
	my $time_diff_total     = timediff( $end_time_inparanoid, $start_time_inparanoid );
	my $debugString         = "Inparanoid\tTotal " . timestr( $time_diff_total, 'all' );
	&attach_to_file( { file_name => $timeFile, text => $debugString . "\n" } );
  }
##### Usearch - Blast
sub testUblastBlastOverlap
  {

	# Goal: Find best Ublast parameters
	# Test Human - 9 other species
	# Compare Hieranoid pairwise predictions with Inparanoid pairwise predictions
	# How: Hieranoid run with only pairs of species
	my $fastaDir = "testdata/inparanoid_data/refGenomes/";
	my $timeFile = "predictions/parameterEstimation/time";
	my $predDir  = "predictions/parameterEstimation/";

	#my $usearch = "/Users/fab/bin/exe/usearch4.1.93_i86darwin32";
	my $usearch       = $Configuration::usearch;
	my $humanFile     = "$fastaDir/Homo_sapiens.fa";
	my $humanFile1000 = "$fastaDir/Homo_sapiens1000.test";
	my $referenceFile = "$predDir/Homo-Pan.reference";
	if ( !-e $humanFile1000 || !-s $humanFile1000 )
	  {
		my $noOfHumanGenes = 1000;

		# take first x human genes
		my $seq_in = Bio::SeqIO->new( -file   => "<$humanFile",
									  -format => 'fasta', );
		my $seq_out = Bio::SeqIO->new( -file   => ">$humanFile1000",
									   -format => 'fasta', );

		# write each entry in the input file to the output file
		my $seq_counter = 0;
		while ( my $inseq = $seq_in->next_seq )
		  {
			last if $seq_counter++ >= $noOfHumanGenes;
			$seq_out->write_seq($inseq);
		  }
	  }
	my @otherSpecies = (
		"Pan_troglodytes",

		#"Mus_musculus",
		#"Anopheles_gambiae",
		#"Arabidopsis_thaliana",
		#"Caenorhabditis_elegans",
		#"Drosophila_melanogaster",
		#"Nematostella_vectensis",
		#"Saccharomyces_cerevisiae",
		#"Xenopus_tropicalis"
	);
	my $currSpeciesFile;
	my @maxaccept_params = ( 0, 1, 5, 10, 40, 100, 500 );
	my @maxrejects_params = ( 0, 8, 16, 32, 64, 128, 265, 500, 1000 );
	my @alignment_types = ("");
	if ( $type eq 'test' )
	  {
		foreach my $currSpecies (@otherSpecies)
		  {
			$currSpeciesFile = "$fastaDir/$currSpecies.fa";
			print "checking $referenceFile\n";

			#exit;
			if ( !-e $referenceFile || !-s $referenceFile )
			  {
				print "Building reference set\n";

				#exit;
				$currSpeciesFile = "$fastaDir/$currSpecies.fa";
				my $parameterCombination = "Reference";
				my $output_file          = "$predDir/Homo-Pan.reference";
				my $blast_call =
"$usearch --quiet --query $humanFile1000  --db $currSpeciesFile  --userfields query+target+bits+ql+tl+qlo+qhi+tlo+thi+qs+ts --userout $output_file --evalue 10 --maxlen 100000 --minlen 4 --nousort";

				#print $blast_call."\n";
				#exit;
				my $start_run = time();
				system($blast_call);
				my $end_run  = time();
				my $run_time = $end_run - $start_run;
				print "Job took $run_time seconds\n";

				#system($tmp_inparanoid_call);
				#system($tmpUblastCmd);
				open my $outputFileFH, '>>', $timeFile
				  or croak "Couldn't open '$timeFile': ";
				print {$outputFileFH} "Homo_sapiens - $currSpecies : \"$parameterCombination\" $run_time\n"
				  or croak "Couldn't write '$timeFile': ";
				close $outputFileFH or croak "Couldn't close '$timeFile': ";
			  }
			my $parsed_output_file = "$predDir/Homo-Pan.parsed";
			## Parsing?
			my $output_file = "$predDir/Homo-Pan.reference";
			if ( !-e $parsed_output_file || !-s $parsed_output_file )
			  {
				my $parse_call =
				  "perl blast2xml.pl -i $output_file | perl blast_parser.pl 40 > $parsed_output_file";
				print "\tParsing ublast reference results\n";

				#print "$parse_call\n";
				#exit;
				system($parse_call);
			  }
			foreach my $current_maxaccept (@maxaccept_params)
			  {
				foreach my $current_maxreject (@maxrejects_params)
				  {
					my $maxtargets = $current_maxaccept + $current_maxreject;
					foreach my $alignmentType (@alignment_types)
					  {
						my $parameterCombination;
						my $blast_call =
"$usearch --quiet --query $humanFile1000  --db $currSpeciesFile  --userfields query+target+bits+ql+tl+qlo+qhi+tlo+thi+qs+ts --evalue 10 --maxlen 100000 --minlen 4 --maxaccepts $current_maxaccept  --maxrejects $current_maxreject";
						if ( $alignmentType ne '' )
						  {
							$parameterCombination =
							    "PARAM_"
							  . $current_maxaccept . "_"
							  . $current_maxreject . "_"
							  . $maxtargets . "_"
							  . $alignmentType . ".out";
							$blast_call .= " --$alignmentType";
						  }
						else
						  {
							$parameterCombination =
							    "PARAM_"
							  . $current_maxaccept . "_"
							  . $current_maxreject . "_"
							  . $maxtargets . ".out";
						  }
						my $output_file =
						    "$predDir/"
						  . basename($humanFile1000) . "-"
						  . basename($currSpeciesFile)
						  . "$parameterCombination";
						$blast_call .= " --userout $output_file";
						if ( !-e $output_file && !-s $output_file )
						  {
							print "\tsearching with $current_maxaccept - $current_maxreject\n";

#my $blast_call = "$usearch --quiet --query $humanFile1000  --db $currSpeciesFile  --userfields query+target+evalue+id+bits+ql+tl+qlo+qhi+tlo+thi+qs+ts --userout $output_file --evalue 0.01 --maxlen 100000 --minlen 4 --maxaccepts $current_maxaccept  --maxrejects $current_maxreject";
							print $blast_call. "\n";

							#next;
							#exit;
							my $start_run = time();
							system($blast_call);
							my $end_run  = time();
							my $run_time = $end_run - $start_run;
							print "Job took $run_time seconds\n";

							#system($tmp_inparanoid_call);
							#system($tmpUblastCmd);
							open my $outputFileFH, '>>', $timeFile
							  or croak "Couldn't open '$timeFile': ";
							print {$outputFileFH}
							  "Homo_sapiens - $currSpecies : \"$parameterCombination\" $run_time\n"
							  or croak "Couldn't write '$timeFile': ";
							close $outputFileFH
							  or croak "Couldn't close '$timeFile': ";
						  }
						my $parsed_output_file =
						    "$predDir/"
						  . basename($humanFile1000) . "-"
						  . basename($currSpeciesFile)
						  . "$parameterCombination.parsed";
						print "do the parsing? $parsed_output_file\n";
						if (    !-e $parsed_output_file
							 || !-s $parsed_output_file )
						  {
							my $parse_call =
							  "perl blast2xml.pl -i $output_file | perl blast_parser.pl 40 > $parsed_output_file";
							print "\tParsing ublast  results\n";
							system($parse_call);
						  }
					  }
				  }
			  }
		  }
	  }
	if ( $type eq 'evaluation' )
	  {

		#print "\t\tEvaluating\n";
		# read reference dataset
		my %statistics_hash;
		my %missingHits;
		my %countHitsHash;
		my %different_values_hash;
		my %missingHitsDistribution;

		# read times
		my %method2time_hash;
		readTimes( $timeFile, \%method2time_hash );

		#print Dumper %method2time_hash;
		#exit;
		my $output_file = "$predDir/Homo-Pan.blastReference";

		#&readUBlastOutput($output_file,\%{$statistics_hash{"reference"}});
		#print "\tReading Blast reference\n";
		&readParsedOutput( $output_file, \%{ $statistics_hash{"BlastReference"} } );
		$output_file = "$predDir/Homo-Pan.parsed";

		#print "\tReading UBlast reference\n";
		&readParsedOutput( $output_file, \%{ $statistics_hash{"Reference"} } );

		#my @values = ("evalue","percentage","bit_score");
		my @values = ("bit_score");

		# read parameter datasets
		foreach my $currSpecies (@otherSpecies)
		  {
			$currSpeciesFile = "$fastaDir/$currSpecies.fa";
			foreach my $current_maxaccept (@maxaccept_params)
			  {
				foreach my $current_maxreject (@maxrejects_params)
				  {
					my $maxtargets = $current_maxaccept + $current_maxreject;
					my $parameterCombination =
					  "PARAM_" . $current_maxaccept . "_" . $current_maxreject . "_" . $maxtargets . ".out";
					my $output_file =
					    "$predDir/"
					  . basename($humanFile1000) . "-"
					  . basename($currSpeciesFile)
					  . "$parameterCombination.parsed";

					#print "reading file $output_file\n";
					#&readUBlastOutput($output_file,\%{$statistics_hash{$parameterCombination}});
					&readParsedOutput( $output_file, \%{ $statistics_hash{$parameterCombination} } );

					# last;
				  }

				# last;
			  }
		  }

		#print "\tread everything...\n";
		my $reference = "BlastReference";
		foreach my $query ( keys( %{ $statistics_hash{$reference} } ) )
		  {

			#print "checking $query...\n";
			foreach my $hit ( keys( %{ $statistics_hash{$reference}{$query} } ) )
			  {

				#print "\t\t and $hit\n";
				$countHitsHash{'total_hit_pairs'}++;
				my $noWrongmethods = 0;
				foreach my $method ( sort keys(%statistics_hash) )
				  {
					next if $method eq $reference;

					#other method has same query-hit pair
					if ( exists $statistics_hash{$method}{$query}{$hit} )
					  {
						foreach my $current_value (@values)
						  {

							#evalue
							if ( $statistics_hash{$reference}{$query}{$hit}{$current_value} !=
								 $statistics_hash{$method}{$query}{$hit}{$current_value} )
							  {

#print "\t\tevalue: reference (".$statistics_hash{$reference}{$query}{$hit}{$current_value}.") - method (".$statistics_hash{$method}{$query}{$hit}{$current_value}.")\n";
								$different_values_hash{$method}{$current_value}{"sum_of_differences"} +=
								  abs( $statistics_hash{$reference}{$query}{$hit}{$current_value} -
									   $statistics_hash{$method}{$query}{$hit}{$current_value} );
								$different_values_hash{$method}{$current_value}{"no_of_different_values"}++;

								#print Dumper %different_values_hash;
								#exit;
							  }
						  }

						#print "\t\t\tcount as hit\n";
						$countHitsHash{$method}{'found_hits'} += 1;
					  }
					else
					  {
						foreach my $current_value (@values)
						  {

							#print "\tno prediction. Count/add as missing\n";
							$missingHits{$method}{$current_value}{"sum_of_missing"} +=
							  $statistics_hash{$reference}{$query}{$hit}{$current_value};
							$missingHits{$method}{$current_value}{"no_of_missing_pairs"}++;
							$missingHitsDistribution{$method}{$current_value}
							  { $statistics_hash{$reference}{$query}{$hit}{$current_value} } += 1;

		#$missingHits{$method}{$hit}{$current_value} = $statistics_hash{"reference"}{$query}{$hit}{$current_value};
						  }

						#print "\t\t\tcount as missing\n";
						$countHitsHash{$method}{'missing_hits'} += 1;
					  }

		  #if(!$noWrongmethods){
		  #        print "evalue: reference (".$statistics_hash{"reference"}{$query}{$hit}{"evalue"}.")\n";
		  #        print "percentage: reference (".$statistics_hash{"reference"}{$query}{$hit}{"percentage"}.")\n";
		  #        print "bit_score: reference (".$statistics_hash{"reference"}{$query}{$hit}{"bit_score"}.")\n";
		  #}
		  #last;
				  }

				#last;
			  }

			#print Dumper %missingHits;
			#exit;
		  }

		#print Dumper %countHitsHash;
		print "method\ttime\thits\%\tmissing\%\taverage_diff\taverage_missing\thighest_bit_score\n";
		foreach my $method ( sort keys(%statistics_hash) )
		  {
			next if $method eq $reference;
			my $largest_bitscore;

			#print $method." is method\n";
			my $perc_found_hits =
			  int( ( $countHitsHash{$method}{'found_hits'} * 100 ) / $countHitsHash{'total_hit_pairs'} );
			my $perc_missing_hits =
			  int( ( $countHitsHash{$method}{'missing_hits'} * 100 ) / $countHitsHash{'total_hit_pairs'} );
			my $averagediffValue;
			my $averageMissingValue;
			foreach my $current_value (@values)
			  {
				my $averageDifference =
				  int( $different_values_hash{$method}{$current_value}{"sum_of_differences"} /
					   $different_values_hash{$method}{$current_value}{"no_of_different_values"} );

#print "$method - $current_value: ".$different_values_hash{$method}{$current_value}{"sum_of_differences"}." / ".$different_values_hash{$method}{$current_value}{"no_of_different_values"}."\n";
				$averagediffValue = $averageDifference;
				my $averageMissing =
				  int( $missingHits{$method}{$current_value}{"sum_of_missing"} /
					   $missingHits{$method}{$current_value}{"no_of_missing_pairs"} );
				$averageMissingValue = $averageMissing;

				#my @evalues = sort { $a <=> $b } keys(%{$missingHitsDistribution{$current_value}});
				#my $largest_evalue = $evalues[0];
				#my @percentage = sort { $a <=> $b } keys(%{$missingHitsDistribution{$current_value}});
				#my $largest_percentage = $percentage[0];
				#print Dumper %missingHitsDistribution;
				my @bitscores = sort { $a <=> $b }
				  keys( %{ $missingHitsDistribution{$method}{$current_value} } );

				#print join(",",@bitscores);
				$largest_bitscore = $bitscores[-1];
			  }
			print "$method\t"
			  . $method2time_hash{$method} . "\t"
			  . $perc_found_hits . "\t"
			  . $perc_missing_hits
			  . "\t$averagediffValue\t$averageMissingValue\t$largest_bitscore\n";
		  }
		exit;
	  }

#  my %parameter_combinations = (
#                 # parameter => filename
#                 "--maxlen 100000 --minlen 4 --maxaccepts 0 --maxrejects 0" => "max10000_min4_accpt0_rej0",
#                 "--maxlen 100000 --minlen 4 --maxaccepts 0 --maxrejects 100" => "max10000_min4_accpt0_rej100",
#                 "--maxlen 100000 --minlen 4 --maxaccepts 40" => "max10000_min4_accpt40",
#                 "--maxlen 100000 --minlen 4 --maxaccepts 40 --maxrejects 0" => "max10000_min4_accpt40_rej0",
#                 "--maxlen 100000 --minlen 4 --maxaccepts 40 --maxrejects 100" => "max10000_min4_accpt40_rej100",
#
# # No sorting
#                 #"--maxlen 100000 --minlen 4 --nousort" => "max10000_min4_nousort",
# # U sort
#                # "--maxlen 100000 --minlen 4 --usort" => "max10000_min4_usort",
# #                "--maxlen 100000 --minlen 4 --usort" => "max10000_min4_usort",
#                 "--maxlen 100000 --minlen 4 --usort --maxaccepts 40 --maxrejects 100" => "max10000_min4_accpt40_rej100_usort",
#                 );
#                 #../../bins/usearch5.0.150_i86linux32 --query ../../testdata/inparanoid_data/refGenomes/Homo_sapiens.fa --db ../../testdata/inparanoid_data/refGenomes/Pan_troglodytes.fa --evalue 0.01 --userfields query+target+bits+ql+tl+qlo+qhi+tlo+thi+qs+ts --maxlen 100000 --minlen 4 --maxaccepts 10 --maxrejects 10 --nousort -userout HomoPan.out
#
#                 foreach my $currSpecies (@otherSpecies){
#                         $currSpeciesFile = "$fastaDir/$currSpecies.fa";
#                         my $inparanoid_call = "perl inparanoid.pl -q $humanFile -d $currSpeciesFile ";
#                         #my $USearchCommand = "bins/usearch4.0.38_i86darwin32 --quiet --query $humanFile --db $currSpeciesFile --evalue 0.01 --userfields query+target+bits+ql+tl+qlo+qhi+tlo+thi+qs+ts ";
#                         ParameterCombi:
#                         foreach my $parameterCombination(keys(%parameter_combinations)){
#                                 my $tmp_inparanoid_call = $inparanoid_call;
#
#                                 #my $t0 = [gettimeofday];
#                                 my $outputFile = "$predDir/Homo_sapiens"."-".$currSpecies."_out_".$parameter_combinations{$parameterCombination};
#                                 if(-e $outputFile && -s $outputFile){
#                                         print "\tOutput file $outputFile already exists. Skipping computation\n";
#                                         next ParameterCombi;
#                                 }
#                                 $tmp_inparanoid_call .= "-p \"$parameterCombination\" -o $outputFile";
#                                 print $tmp_inparanoid_call."\n";
#                                 #exit;
#                                 my $blast_outputAB = $predDir."/Homo_sapiens-$currSpecies";
#                                 my $blast_outputBA = $predDir."/$currSpecies-Homo_sapiens";
#                                 my $blast_outputAA = $predDir."/Homo_sapiens-Homo_sapiens";
#                                 my $blast_outputBB = $predDir."/".$currSpecies . "-" . $currSpecies;
#                                 print "\tdeleting $blast_outputAB\t$blast_outputBA\n$blast_outputAA\n$blast_outputBB\n";
#                                 unlink($blast_outputAB);
#                                 unlink($blast_outputBA);
#                                 unlink($blast_outputAA);
#                                 unlink($blast_outputBB);
#                                 #exit;
#                                 #my $tmpUblastCmd = "time $USearchCommand $parameterCombination --userout $outputFile";
#                                 #print "$tmpUblastCmd\n";
#                                 my $start_run = time();
#                                 system($tmp_inparanoid_call);
#                                 my $end_run = time();
#                                 my $run_time = $end_run - $start_run;
#                                 print "Job took $run_time seconds\n";
#                                 system($tmp_inparanoid_call);
#
#                                 #system($tmpUblastCmd);
#                                 open my $outputFileFH, '>>', $timeFile or croak "Couldn't open '$timeFile': ";
#                                 print {$outputFileFH} "Homo_sapiens - $currSpecies : \"$parameterCombination\" $run_time\n" or croak "Couldn't write '$timeFile': ";
#                                 close $outputFileFH or croak "Couldn't close '$timeFile': ";
#                                 #next;
#
#                                 #exit;
#                                 # measuring time
#
#                                 #my $elapsed = tv_interval ( $t0, [gettimeofday]);
#                                 #print "\tTime elapsed: ".$elapsed."\n";
#                         }
#                 }
#         exit;
#
  }

sub readTimes
  {
	my ( $timeFile, $hashReference ) = (@_);
	if ( !-e $timeFile || !-s $timeFile )
	  {
		print "time file $timeFile does not exist. Skipping\n";
		exit;
	  }
	open my $time_FH, '<', $timeFile or croak "Couldn't open '$timeFile': ";
	while (<$time_FH>)
	  {

		#/Homo_sapiens - Pan_troglodytes : "Reference" 789/;
		/\w+ - \w+ : \"(.*)\" (\d+)/;
		my ( $method, $time ) = ( $1, $2 );
		$hashReference->{$method} = $time;
	  }
	if ( !keys( %{$hashReference} ) )
	  {
		print "\tCould not get time data\n";
		exit;
	  }
  }

sub readParsedOutput
  {
	my ( $ublastfile, $hashReference ) = (@_);

	#
	# hash{'query'}{'hit'}{'evalue'} = $evalue
	# hash{'query'}{'hit'}{'bitscore'} = $bitscore
	#
	if ( !-e $ublastfile || !-s $ublastfile )
	  {
		print "ublast file $ublastfile does not exist. Skipping\n";
		exit;
	  }
	open my $blast_output_FH, '<', $ublastfile
	  or croak "Couldn't open '$ublastfile': ";
	while (<$blast_output_FH>)
	  {

		#while (my $line = <$blast_output_FH>) {
		chomp();
		my $line = $_;
		next if /^\s*#/;
		next if !(/\w+/);
		next if /^query/;

		#print $line."\n";
		my @splitted_line = split;
		my (
			 $query_name, $target_name,       $bit_score, $queryLength,
			 $hitLength,  $query_start,       $query_end, $hit_start,
			 $hit_end,    $query_coordinates, $hit_coordinates
		) = (@splitted_line);
		$hashReference->{$query_name}{$target_name}{"queryLength"}       = $queryLength;
		$hashReference->{$query_name}{$target_name}{"hitLength"}         = $hitLength;
		$hashReference->{$query_name}{$target_name}{"bit_score"}         = $bit_score;
		$hashReference->{$query_name}{$target_name}{"query_coordinates"} = $query_coordinates;
		$hashReference->{$query_name}{$target_name}{"hit_coordinates"}   = $hit_coordinates;

		#print Dumper $hashReference;
		#exit;
	  }
	close $blast_output_FH or die "Could not close $ublastfile\n";
  }

sub readUBlastOutput
  {
	my ( $ublastfile, $hashReference ) = (@_);

	#
	# hash{'query'}{'hit'}{'evalue'} = $evalue
	# hash{'query'}{'hit'}{'bitscore'} = $bitscore
	#
	if ( !-e $ublastfile || !-s $ublastfile )
	  {
		print "ublast file $ublastfile does not exist. Skipping\n";
		exit;
	  }
	open my $blast_output_FH, '<', $ublastfile
	  or croak "Couldn't open '$ublastfile': ";
	while (<$blast_output_FH>)
	  {

		#while (my $line = <$blast_output_FH>) {
		chomp();
		my $line = $_;
		next if /^\s*#/;
		next if !(/\w+/);
		next if /^query/;
		my @splitted_line = split;
		my (
			 $query_name,   $target_name,      $evalue,      $percentage, $bit_score,
			 $query_length, $hit_length,       $query_start, $query_end,  $hit_start,
			 $hit_end,      $query_ali_length, $hit_ali_length
		) = (@splitted_line);

		if ( exists $hashReference->{$query_name}{$target_name}
			 && $bit_score < $hashReference->{$query_name}{$target_name}{"bit_score"} )
		  {
			;
		  }
		else
		  {
			$hashReference->{$query_name}{$target_name}{"evalue"}     = $evalue;
			$hashReference->{$query_name}{$target_name}{"percentage"} = $percentage;
			$hashReference->{$query_name}{$target_name}{"bit_score"}  = $bit_score;
		  }

		#print Dumper $hashReference;
		#exit;
	  }
	close $blast_output_FH or die "Could not close $ublastfile\n";
  }
##### Inparanoid Overlap
sub inparanoid_testing
  {
	my ($arg_ref) = @_;

	# Read Reference predictions : Inparanoid predictions
	# $orthologous_to{'geneA'} = geneB
	my $inparanoid_predictions_directory =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/inparanoid_data/refGenomesPredictions/";
	my $inparanoid_small_predictions =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/inparanoid_data/10species_200genes/predictions/10species.pred";
	my $hieranoid_prediction_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/archived_results/results_inparanoid_consensus/nodes/Eukaryota/Eukaryota.expandedGroups.txt";
	my $hieranoid_fasta_seq_directory =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/inparanoid_data/refGenomes/";
	my $InparanoidCompRawDataPairwise =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/pictures/R_sources/InparanoidCompRawDataPairwise.csv";

	#  Get all IDS used in Hieranoid predictions
	my %allHieranoidID_hash;
	my @allHieranoidIDs = `grep ">" $hieranoid_fasta_seq_directory/*.fa`;
	print "\tgrepping IDs from $hieranoid_fasta_seq_directory/\n";
	foreach (@allHieranoidIDs)
	  {
		/>(.*)/;
		if ( exists $allHieranoidID_hash{$1} )
		  {
			print "\tWarn. Duplicate ID $1 found\n";
			exit;
		  }
		$allHieranoidID_hash{$1} = 1;

		#print "saved $1\n";
		#exit;
	  }
	print "\tFinished reading IDs used in Hieranoid\n";
	if ( !keys(%allHieranoidID_hash) )
	  {
		print "\tCould not fetch all IDs used by Hieranoid ( grep \">\" $hieranoid_fasta_seq_directory/*)\n";
		exit;
	  }

	#exit;
	my %hieranoid_prediction;
	my %ID2Group_hash;
	getPredictionMethodPredictions( $hieranoid_prediction_file, \%hieranoid_prediction, \%ID2Group_hash );
	if ( !keys(%hieranoid_prediction) )
	  {
		print "\tCould not read predictions from $hieranoid_prediction_file (hieranoid_prediction)\n";
		exit;
	  }
	if ( !keys(%ID2Group_hash) )
	  {
		print "\tCould not read predictions from $hieranoid_prediction_file (ID2Group_hash)\n";
		exit;
	  }
	my @inparanoid_predictions = (
								   "sqltable.9598_pan_troglodytes-9606_homo_sapiens",
								   "sqltable.10090_mus_musculus-9606_homo_sapiens",
								   "sqltable.8634_xenopus_tropicalis-9606_homo_sapiens",
								   "sqltable.7165_anopheles_gambiae-9606_homo_sapiens",
								   "sqltable.7227_drosophila_melanogaster-9606_homo_sapiens",
								   "sqltable.6239_caenorhabditis_elegans-9606_homo_sapiens",
								   "sqltable.45351_nematostella_vectensis-9606_homo_sapiens",
								   "sqltable.4932_saccharomyces_cerevisiae-9606_homo_sapiens",
								   "sqltable.3702_arabidopsis_thaliana-9606_homo_sapiens",
	);
	my %statistics_hash;
	## {comparison}{tp} = 4
	##             {fp} = 4
	foreach my $currentPredictionFile (@inparanoid_predictions)
	  {
		my %id_to_cluster_hash;
		my %cluster_to_id_hash;
		my %group_hash;
		read_pairwise_results(
							   {
								 profile_file => $inparanoid_predictions_directory . "/" . $currentPredictionFile,
								 id_to_cluster_hashref => \%id_to_cluster_hash,
								 cluster_to_id_hashref => \%cluster_to_id_hash,
								 group_href            => \%group_hash,
								 no_of_groups          => 100000
							   }
		);
		if ( !keys(%cluster_to_id_hash) )
		  {
			print "\tMissing cluster 2 ID information\n";
			exit;
		  }
		if ( !keys(%cluster_to_id_hash) )
		  {
			print "\tMissing cluster 2 ID information\n";
			exit;
		  }
		print "\tchecking $currentPredictionFile\n";

		#exit;
		my $no_predictions_to_look_at = 0;
		$statistics_hash{$currentPredictionFile}{'totalGroups'} = keys(%cluster_to_id_hash);
	  GROUP:
		foreach my $inparanoidPrediction (
										   sort { $a <=> $b }
										   keys(%cluster_to_id_hash)
		  )
		  {
			my %groupHash;
			$no_predictions_to_look_at++;
			if ( $no_predictions_to_look_at++ > 100000 )
			  {
				print Dumper %statistics_hash;
				exit;
			  }
			print "\t\tchecking Inparanoid group $inparanoidPrediction: "
			  . join( ",", @{ $cluster_to_id_hash{$inparanoidPrediction} } ) . "\n";
			foreach my $inparanoid_gene ( @{ $cluster_to_id_hash{$inparanoidPrediction} } )
			  {
				print "\t\t\tchecking gene $inparanoid_gene\n";

				# Gene was not Hieranoid dataset
				if ( !exists $allHieranoidID_hash{$inparanoid_gene} )
				  {
					print "\t\t\t\t ID $inparanoid_gene not used for hieranoid predictions\n";
					exit;
					next GROUP;
				  }

				#else{
				#print "test this one\n";
				#}
				# Gene was not predicted by Hieranoid --> count as FN
				if ( !exists $ID2Group_hash{$inparanoid_gene} )
				  {
					print "\t\t\t\tno hieranoid prediction for $inparanoid_gene\n";
					next;
				  }
				else
				  {
					print "\t\t\thas: " . $ID2Group_hash{$inparanoid_gene} . "\n";
					$groupHash{ $ID2Group_hash{$inparanoid_gene} }++;
				  }
			  }

			#print Dumper %groupHash;
			#exit;
			if ( !keys(%groupHash) )
			  {
				print "\t\t\t\tFN: " . join( ",", @{ $cluster_to_id_hash{$inparanoidPrediction} } ) . "\n";
				$statistics_hash{$currentPredictionFile}{'FN'}++;
				next;
			  }
			if ( keys(%groupHash) > 1 )
			  {
				print "\t\t\t\tFP: " . join( ",", keys(%groupHash) ) . "\n";
				$statistics_hash{$currentPredictionFile}{'FP'}++;
				next;
			  }
			if ( keys(%groupHash) == 1 )
			  {
				print "\t\t\t\tTP: " . join( ",", keys(%groupHash) ) . "\n";
				$statistics_hash{$currentPredictionFile}{'TP'}++;
				next;
			  }
		  }

		#$statistics_hash{$currentPredictionFile}{'allInparanoidGroups'} = keys(%cluster_to_id_hash);
		#print Dumper %statistics_hash;
		#exit;
	  }
	my $line2print = "Comparison\tFN\tFP\tTP\tAllGroups\n";
	foreach my $comp ( (@inparanoid_predictions) )
	  {

		#print $comp." - \n";
		$line2print .= "$comp";
		foreach my $value ( sort keys( %{ $statistics_hash{$comp} } ) )
		  {

			#	print "\t\t".$value.": ".$statistics_hash{$comp}{$value}."\n";
			$line2print .= "\t" . $statistics_hash{$comp}{$value};
		  }
		$line2print .= "\n";
	  }
	write_to_file( { file_name => $InparanoidCompRawDataPairwise, text => $line2print } );
	return 1;
  }

sub inparanoid_testing_pairwise
  {
	my ($arg_ref) = @_;

	# Read Reference predictions : Inparanoid predictions
	# $orthologous_to{'geneA'} = geneB
	my $inparanoid_predictions_directory =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/inparanoid_data/refGenomesPredictions/";
	my $inparanoid_small_predictions =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/inparanoid_data/10species_200genes/predictions/10species.pred";
	my $hieranoid_prediction_file =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/archived_results/results_inparanoid_consensus/nodes/Eukaryota/Eukaryota.expandedGroups.txt";
	my $hieranoid_fasta_seq_directory =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/inparanoid_data/refGenomes/";
	my $InparanoidCompRawDataPairwise =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/pictures/R_sources/InparanoidCompRawDataPairwise.csv";
	my $debugprint = 1;

	#  Get all IDS used in Hieranoid predictions
	my %allHieranoidID_hash;
	my @allHieranoidIDs = `grep ">" $hieranoid_fasta_seq_directory/*.fa`;
	print "\tgrepping IDs from $hieranoid_fasta_seq_directory/\n";
	foreach (@allHieranoidIDs)
	  {
		/>(.*)/;
		if ( exists $allHieranoidID_hash{$1} )
		  {
			print "\tWarn. Duplicate ID $1 found\n";
			exit;
		  }
		$allHieranoidID_hash{$1} = 1;

		#print "saved $1\n";
		#exit;
	  }
	print "\tFinished reading IDs used in Hieranoid\n";
	if ( !keys(%allHieranoidID_hash) )
	  {
		print "\tCould not fetch all IDs used by Hieranoid ( grep \">\" $hieranoid_fasta_seq_directory/*)\n";
		exit;
	  }

	#exit;
	my %hieranoid_prediction;
	my %ID2Group_hash;
	getPredictionMethodPredictions( $hieranoid_prediction_file, \%hieranoid_prediction, \%ID2Group_hash );
	if ( !keys(%hieranoid_prediction) )
	  {
		print "\tCould not read predictions from $hieranoid_prediction_file (hieranoid_prediction)\n";
		exit;
	  }
	if ( !keys(%ID2Group_hash) )
	  {
		print "\tCould not read predictions from $hieranoid_prediction_file (ID2Group_hash)\n";
		exit;
	  }
	my @inparanoid_predictions = (
								   "sqltable.9598_pan_troglodytes-9606_homo_sapiens",
								   "sqltable.10090_mus_musculus-9606_homo_sapiens",
								   "sqltable.8634_xenopus_tropicalis-9606_homo_sapiens",
								   "sqltable.7165_anopheles_gambiae-9606_homo_sapiens",
								   "sqltable.7227_drosophila_melanogaster-9606_homo_sapiens",
								   "sqltable.6239_caenorhabditis_elegans-9606_homo_sapiens",
								   "sqltable.45351_nematostella_vectensis-9606_homo_sapiens",
								   "sqltable.4932_saccharomyces_cerevisiae-9606_homo_sapiens",
								   "sqltable.3702_arabidopsis_thaliana-9606_homo_sapiens",
	);
	my %statistics_hash;
	## {comparison}{tp} = 4
	##             {fp} = 4
	my %countMissingSpeciesIDs;
	foreach my $currentPredictionFile (@inparanoid_predictions)
	  {
		my %id_to_cluster_hash;
		my %cluster_to_id_hash;
		my %group_hash;
		read_pairwise_results(
							   {
								 profile_file => $inparanoid_predictions_directory . "/" . $currentPredictionFile,
								 id_to_cluster_hashref => \%id_to_cluster_hash,
								 cluster_to_id_hashref => \%cluster_to_id_hash,
								 group_href            => \%group_hash,
								 no_of_groups          => 100000
							   }
		);
		if ( !keys(%cluster_to_id_hash) )
		  {
			print "\tMissing cluster 2 ID information\n";
			exit;
		  }
		if ( !keys(%cluster_to_id_hash) )
		  {
			print "\tMissing cluster 2 ID information\n";
			exit;
		  }
		print "\tchecking $currentPredictionFile\n";

		# print Dumper %cluster_to_id_hash;
		#exit;
		my $no_predictions_to_look_at = 0;
		$statistics_hash{$currentPredictionFile}{'totalGroups'} = keys(%cluster_to_id_hash);
	  GROUP:
		foreach my $inparanoidPredictionOG (
											 sort { $a <=> $b }
											 keys(%cluster_to_id_hash)
		  )
		  {

			#foreach my $inparanoidPredictionOG(sort keys(%refOG2Members)){
			#next if $currentRefOG ne 'RefOG001';
			print "checking $inparanoidPredictionOG \n" if $debugprint;
			my $noOfMissingGenes = 0;
			my $noOfGenePairs    = 0;
			my $groupMismatches  = 0;
			my $noMatchingPairs  = 0;
			my $wrongIDs         = 0;
			my @GenePairs        = combine( 2, @{ $cluster_to_id_hash{$inparanoidPredictionOG} } );
			my $noOfPairs        = scalar(@GenePairs);

			#print join(",",@GenePairs);
			$statistics_hash{$currentPredictionFile}{"RefOGPairs"} += $noOfPairs;

			#$statistics_hash{$currentPredictionFile}{"RefOGAllPairs"} = ;
			#	$statistics_hash{$currentPredictionFile}{$currentRefOG}{"RefOGPairs"} = $noOfPairs;
			my %countedMissingGenes;
			for ( my $i = 0 ; $i < scalar(@GenePairs) ; $i++ )
			  {
				my $pairDestiny = "";
				my ( $IDA, $IDB ) = ( $GenePairs[$i][0], $GenePairs[$i][1] );
				print "\t\t$IDA - $IDB: " if $debugprint;
				my $group4IDA = $ID2Group_hash{$IDA};
				my $group4IDB = $ID2Group_hash{$IDB};

				# Check
				if ( !defined($group4IDA) )
				  {
					$pairDestiny .= "Missing gene (no prediction for IDA $IDA)";
					my $Species4ID = $IDA;
					$Species4ID =~ s/\d+//g;
					if (    !exists $countedMissingGenes{$IDA}
						 && !exists $countedMissingGenes{$IDB} )
					  {
						$countMissingSpeciesIDs{$currentPredictionFile}{$Species4ID}++;
						$noOfMissingGenes++;
						$wrongIDs++;
					  }
					else
					  {
						$pairDestiny .= " (not counted)";
					  }
				  }
				elsif ( !defined($group4IDB) )
				  {
					$pairDestiny .= "Missing gene (no prediction for IDB $IDB)";
					if (    !exists $countedMissingGenes{$IDA}
						 && !exists $countedMissingGenes{$IDB} )
					  {
						$noOfMissingGenes++;
						$wrongIDs++;
					  }
					else
					  {
						$pairDestiny .= " (not counted)";
					  }
				  }
				elsif ( $group4IDA eq $group4IDB )
				  {
					$pairDestiny .= "TP (both map to same gene: $group4IDA";
					$noMatchingPairs++;
				  }
				elsif ( $group4IDA ne $group4IDB )
				  {
					$pairDestiny .= "Missing gene (diff group $group4IDA - $group4IDB)";
					if (    !exists $countedMissingGenes{$IDA}
						 && !exists $countedMissingGenes{$IDB} )
					  {
						$noOfMissingGenes++;
						$groupMismatches++;
					  }
					else
					  {
						$pairDestiny .= " (not counted)";
					  }
				  }
				$countedMissingGenes{$IDA} = 1;
				$countedMissingGenes{$IDB} = 1;
				if ($debugprint)
				  {
					print "-> $pairDestiny\n" if $debugprint;
				  }
			  }

#$statistics_hash{$currentPredictionFile}{$inparanoidPredictionOG}{"RefOGMissingPairs_GroupMismatches"} = $groupMismatches;
#$statistics_hash{$currentPredictionFile}{$inparanoidPredictionOG}{"RefOGMissingPairs_WrongIDs"} = $wrongIDs;
#$statistics_hash{$currentPredictionFile}{$inparanoidPredictionOG}{"RefOGMatchingPairs"} = $noMatchingPairs;
#$statistics_hash{$currentPredictionFile}{$inparanoidPredictionOG}{"RefOGMissingPairs"} = $noOfMissingGenes;
#$OGStatistics{$predictionMethod}{$currentRefOG}{"RefOGPairs"} = $noOfGenePairs;
#$statistics_hash{$currentPredictionFile}{$inparanoidPredictionOG}{"FNQuote"} = int(($noOfMissingGenes*100)/$noOfPairs);
			$statistics_hash{$currentPredictionFile}{"RefOGMatchingPairs"} += $noMatchingPairs;
			$statistics_hash{$currentPredictionFile}{"RefOGMissingPairs"}  += $noOfMissingGenes;

			#$overallStatistics{$predictionMethod}{"RefOGPairs"} += $noOfGenePairs;
			$statistics_hash{$currentPredictionFile}{"RefOGMissingPairs_GroupMismatches"} += $groupMismatches;
			$statistics_hash{$currentPredictionFile}{"RefOGMissingPairs_WrongIDs"}        += $wrongIDs;
		  }

		#last;
	  }

	#print Dumper %statistics_hash;
	#exit;
	my $line2print = "Comparison\tFN\tFP\tTP\tAllGroups\n";
	foreach my $comp ( (@inparanoid_predictions) )
	  {

		#print $comp." - \n";
		$line2print .= "$comp";
		$line2print .= "\t" . $statistics_hash{$comp}{"RefOGMissingPairs"};
		$line2print .= "\t1";
		$line2print .= "\t" . $statistics_hash{$comp}{"RefOGMatchingPairs"};
		$line2print .= "\t" . $statistics_hash{$comp}{"RefOGPairs"};

		#foreach my $value(sort keys(%{$statistics_hash{$comp}})){
		#	print "\t\t".$value.": ".$statistics_hash{$comp}{$value}."\n";
		#	$line2print .= "\t".$statistics_hash{$comp}{$value};
		#}
		$line2print .= "\n";
	  }
	write_to_file( { file_name => $InparanoidCompRawDataPairwise, text => $line2print } );
	return 1;
  }

sub getspecies4sequence
  {
	#### PARAMETER VARIABLES
	my ($arg_ref)    = @_;
	my $folder       = $arg_ref->{folder};
	my $hashref      = $arg_ref->{hashref};
	my @species2read = glob("$folder/*.fa");
	foreach my $species_sequence_file (@species2read)
	  {
		my $species_name = basename($species_sequence_file);
		$species_name =~ s/\.fa//;
		print "sequence2species: $species_name\n";
		if ( !-e $species_sequence_file && !-s $species_sequence_file )
		  {
			next SPECIES_TO_ADD;
		  }
		open my $SEQUENCE_FILE, '<', $species_sequence_file
		  or croak "Couldn't open '$species_sequence_file': $OS_ERROR";
		### READING FROM FILE
		while ( my $line = <$SEQUENCE_FILE> )
		  {
			chomp($line);
			if ( $line =~ /^>/ )
			  {
				my @words = split( /\s/, $line );
				my $seq_id = $words[0];
				$seq_id =~ s/>//;
				$hashref->{$seq_id} = $species_name;

				#print "\$hashref->{$seq_id} = $species_name;\n";
				#exit;
			  }
		  }
		### CLOSING FILE
		close $SEQUENCE_FILE
		  or croak "Couldn't close '$species_sequence_file': $OS_ERROR";

		#        print $sequence_href->{"LCT_PANTR"};
		print("Could not read sequence 2 species information\n")
		  if !( keys(%$hashref) );
	  }
  }

sub species_overlap
  {
	my $directory_to_test =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/nodes/Eukaryota/alignments/";
	my $reference_predictions =
	  "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/treefam/new_mapping.txt";

	# Using here two hashes, to capture both mappings: "geneA--> geneB", "geneB--> geneA"
	my %reference_prediction1 = ();
	my %reference_prediction2 = ();
	print "Reading reference orthology predictions\n";
## Read references
	### OPENING FILE
	open my $reference_predictions_FH, '<', $reference_predictions
	  or croak "Couldn't open '$reference_predictions': $OS_ERROR";
	### READING FROM FILE
	while ( my $line = <$reference_predictions_FH> )
	  {
		next if $line =~ /^$/;
		my ( $geneA, $geneB ) = split( /\s/, $line );
		if ( !defined $geneA || !defined $geneB )
		  {
			print "\tCould not read line $_  ($geneA,$geneB)\n";
			next;
		  }
		$reference_prediction1{$geneA}{$geneB} = 1;
		$reference_prediction2{$geneB}{$geneA} = 1;
	  }

	#	    $reference_prediction1{$geneA} = $geneB;
	#	    $reference_prediction2{$geneB} = $geneA;
	close $reference_predictions_FH
	  or croak "Couldn't close '$reference_predictions': $OS_ERROR";
	print "\tRead " . keys(%reference_prediction1) . " and " . keys(%reference_prediction2) . " keys\n";
	my $trees_to_check = 10;
	my ( $true_positive_total, $false_positive_total ) = (0);
  TREE:
	foreach my $tree2test ( glob("$directory_to_test/*.tre") )
	  {
		print "\ttesting $tree2test\n";
		my ( $true_positive, $false_positive ) = (0);
		my %orthology_predictions = ();
		my @orthologs_array       = ();

		#	Hieranoid_module::species_overlap({input_tree => $tree2test,
		#					   orthology_mapping_href => \%orthology_predictions,
		#					   orthologs_aref => \@orthologs_array});
		### USE ETE for species overlap
		# Call python method
		# read results in hash
		my $ete_arg    = "python ete_species_overlap.py $tree2test";
		my @ete_result = `$ete_arg`;
		if ( !@ete_result )
		  {
			print "Could not execue ETE ($ete_arg)\n";
			exit;
		  }
		foreach (@ete_result)
		  {
			next if !/-/;
			my ( $gene1, $gene2 ) = split( /-/, $_ );
			$gene1 =~ s/\s*//g;
			$gene2 =~ s/\s*//g;
			$orthology_predictions{$gene1}{$gene2} = 1;
			$orthology_predictions{$gene2}{$gene1} = 1;
		  }
		if ( !keys(%orthology_predictions) )
		  {
			print "\tno predictions for $tree2test ($ete_arg).\n";
			exit;
			next TREE;
		  }

		# print Dumper %orthology_predictions;
		#	exit;
		my %already_checked = ();
		print "Checking " . keys(%orthology_predictions) . " predictions\n";
		foreach my $geneA ( keys(%orthology_predictions) )
		  {
			foreach my $geneB ( keys( %{ $orthology_predictions{$geneA} } ) )
			  {
				my $is_ortholog = 0;
				next if exists $already_checked{$geneA}{$geneB};
				print "\tchecking predicted orthology $geneA - $geneB\n";
				$is_ortholog = 1
				  if exists $reference_prediction1{$geneA}{$geneB};
				$is_ortholog = 1
				  if exists $reference_prediction1{$geneB}{$geneA};
				$is_ortholog = 1
				  if exists $reference_prediction2{$geneA}{$geneB};
				$is_ortholog = 1
				  if exists $reference_prediction2{$geneB}{$geneA};
				$already_checked{$geneB}{$geneA} = 1;
				$true_positive++        if $is_ortholog;
				$true_positive_total++  if $is_ortholog;
				$false_positive++       if !$is_ortholog;
				$false_positive_total++ if !$is_ortholog;
			  }
		  }
		print "have $true_positive TP ($true_positive_total) and $false_positive FP ($false_positive_total)\n";
		exit if $trees_to_check-- == 0;
	  }
	exit;

	#	orthology_mapping_href
  }

sub BrigitteTrees2species
  {
	my @brigitteTreesFastaFiles =
	  glob("/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/BrigitteTrees/paper/*.fa");
	my $ID2speciesMapping =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/BrigitteTrees/paper/mappingTable.txt";
	my $outputDirectory = "/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/BrigitteTrees/";
	my %ID2speciesHash  = ();
	my $species;
	open my $ID2speciesMapping_FH, '<', $ID2speciesMapping
	  or croak "Couldn't open '$ID2speciesMapping': ";
	while (<$ID2speciesMapping_FH>)
	  {
		next if /^$/;
		my ( $ID, $species1, $species2, $dbID ) = split(/\s+/);
		if (    !defined $ID
			 || !defined $species1
			 || !defined $species2
			 || !defined $dbID )
		  {
			print "Parsing-problem in line $_\n";
			exit;
		  }
		$species = $species1 . "_" . $species2;
		print "save $ID = $species\n";
		$ID2speciesHash{$ID} = $species;
	  }
	close $ID2speciesMapping_FH
	  or croak "Couldn't close '$ID2speciesMapping': ";
	if ( !keys(%ID2speciesHash) )
	  {
		print "\tCould not read id2species mappings\n";
		exit;
	  }
	print Dumper %ID2speciesHash;

	#exit;
	my %ID2seq_hash = ();
	foreach my $currentFile (@brigitteTreesFastaFiles)
	  {
		my $seq_in = Bio::SeqIO->new( -format => 'fasta',
									  -file   => $currentFile, );
		while ( my $seq = $seq_in->next_seq() )
		  {
			my $id = $seq->id;
			if ( !exists $ID2speciesHash{$id} )
			  {
				print "could not find $id\n";

				#                                    exit;
				$id =~ /.*_(\w+)/;
				my $species2map = $1;
				if ( !exists $convert_hash{$species2map} )
				  {
					print "could not find $species2map\n";
					exit;
				  }
				$species = $convert_hash{$species2map};
			  }
			else
			  {
				$species = $ID2speciesHash{$id};
			  }
			print "saving $id in " . $species . "\n";

			#                            $ID2seq_hash{$species}{$id} = $seq->seq;
			$ID2seq_hash{$species} .= ">$id\n" . $seq->seq . "\n";

			#exit;
		  }
	  }
	foreach my $species ( keys(%ID2seq_hash) )
	  {
		my $speciesFile  = "$outputDirectory/$species.fa";
		my $string2print = "";
		open my $speciesFile_FH, '>', $speciesFile
		  or croak "Couldn't open '$speciesFile': ";
		print {$speciesFile_FH} $ID2seq_hash{$species}
		  or croak "Couldn't write '$speciesFile': ";
		close $speciesFile_FH or croak "Couldn't close '$speciesFile': ";
	  }
	exit;
	return 1;
  }

sub test_score_symmetry
  {
	my $AB =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/results/nodes/Eutheria/Euarchontoglires-Laurasiatheria";
	my $BA =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/results/nodes/Eutheria/Laurasiatheria-Euarchontoglires";
	my %stat_hash = ();
	open my $ab_file, '<', $AB or croak "Couldn't open '$AB': $!\n";
	while ( my $line = <$ab_file> )
	  {
		my @splitted_line = split( /\s+/, $line );
		my (
			 $query_name, $target_name, $score,      $queryLength, $hitLength, $qLongSeq,
			 $hLongSeq,   $qTotLenght,  $hTotLenght, $query_coord, $target_coord
		  )
		  = (
			  $splitted_line[0], $splitted_line[1], $splitted_line[2], $splitted_line[3],
			  $splitted_line[4], $splitted_line[5], $splitted_line[6], $splitted_line[7],
			  $splitted_line[8], $splitted_line[9], $splitted_line[10]
		  );
		$stat_hash{$query_name}{$target_name}{"score"}        = $score;
		$stat_hash{$query_name}{$target_name}{"queryLength"}  = $queryLength;
		$stat_hash{$query_name}{$target_name}{"hitLength"}    = $hitLength;
		$stat_hash{$query_name}{$target_name}{"qLongSeq"}     = $qLongSeq;
		$stat_hash{$query_name}{$target_name}{"hLongSeq"}     = $hLongSeq;
		$stat_hash{$query_name}{$target_name}{"qTotLenght"}   = $qTotLenght;
		$stat_hash{$query_name}{$target_name}{"hTotLenght"}   = $hTotLenght;
		$stat_hash{$query_name}{$target_name}{"query_coord"}  = $query_coord;
		$stat_hash{$query_name}{$target_name}{"target_coord"} = $target_coord;
	  }
	close $ab_file or croak "Couldn't close '$AB': $!\n";
	open my $ba_file, '<', $BA or croak "Couldn't open '$BA': $!\n";
	while ( my $line = <$ba_file> )
	  {
		my @splitted_line = split( /\s+/, $line );
		my (
			 $query_name, $target_name, $score,      $queryLength, $hitLength, $qLongSeq,
			 $hLongSeq,   $qTotLenght,  $hTotLenght, $query_coord, $target_coord
		  )
		  = (
			  $splitted_line[0], $splitted_line[1], $splitted_line[2], $splitted_line[3],
			  $splitted_line[4], $splitted_line[5], $splitted_line[6], $splitted_line[7],
			  $splitted_line[8], $splitted_line[9], $splitted_line[10]
		  );
		$stat_hash{$query_name}{$target_name}{"score"}        = $score;
		$stat_hash{$query_name}{$target_name}{"queryLength"}  = $queryLength;
		$stat_hash{$query_name}{$target_name}{"hitLength"}    = $hitLength;
		$stat_hash{$query_name}{$target_name}{"qLongSeq"}     = $qLongSeq;
		$stat_hash{$query_name}{$target_name}{"hLongSeq"}     = $hLongSeq;
		$stat_hash{$query_name}{$target_name}{"qTotLenght"}   = $qTotLenght;
		$stat_hash{$query_name}{$target_name}{"hTotLenght"}   = $hTotLenght;
		$stat_hash{$query_name}{$target_name}{"query_coord"}  = $query_coord;
		$stat_hash{$query_name}{$target_name}{"target_coord"} = $target_coord;
	  }
	close $ba_file or croak "Couldn't close '$BA': $!\n";
	foreach my $query ( keys(%stat_hash) )
	  {
		foreach my $hit ( keys( %{ $stat_hash{$query} } ) )
		  {
			if (    !exists $stat_hash{$query}{$hit}
				 || !exists $stat_hash{$query}{$hit} )
			  {
				print "Just one-sided match\n";
				next;
			  }
			print "score not defined $query - $hit\n"
			  if ( !defined $stat_hash{$query}{$hit}{"score"} );
			print "score not defined $hit - $query\n"
			  if ( !defined $stat_hash{$hit}{$query}{"score"} );
			if ( $stat_hash{$query}{$hit}{"score"} != $stat_hash{$hit}{$query}{"score"} )
			  {
				print "score: "
				  . $stat_hash{$query}{$hit}{"score"} . " != "
				  . $stat_hash{$hit}{$query}{"score"} . "\n";
			  }

#  if($stat_hash{$query}{$hit}{"queryLength"} != $stat_hash{$hit}{$query}{"hitLength"}){print "queryLength: ".$stat_hash{$query}{$hit}{"queryLength"}." != ".$stat_hash{$hit}{$query}{"hitLength"}."\n";}
#  if($stat_hash{$query}{$hit}{"hitLength"} != $stat_hash{$hit}{$query}{"queryLength"}){print "hitLength: ".$stat_hash{$query}{$hit}{"hitLength"}." != ".$stat_hash{$hit}{$query}{"queryLength"}."\n";}
#  if($stat_hash{$query}{$hit}{"qLongSeq"} != $stat_hash{$hit}{$query}{"hLongSeq"}){print "qLongSeq: ".$stat_hash{$query}{$hit}{"qLongSeq"}." != ".$stat_hash{$hit}{$query}{"hLongSeq"}."\n";}
#  if($stat_hash{$query}{$hit}{"hLongSeq"} != $stat_hash{$hit}{$query}{"qLongSeq"}){print "hLongSeq: ".$stat_hash{$query}{$hit}{"hLongSeq"}." != ".$stat_hash{$hit}{$query}{"qLongSeq"}."\n";}
#  if($stat_hash{$query}{$hit}{"qTotLenght"} != $stat_hash{$hit}{$query}{"hTotLenght"}){print "qTotLenght: ".$stat_hash{$query}{$hit}{"qTotLenght"}." != ".$stat_hash{$hit}{$query}{"hTotLenght"}."\n";}
#  if($stat_hash{$query}{$hit}{"hTotLenght"} != $stat_hash{$hit}{$query}{"qTotLenght"}){print "hTotLenght: ".$stat_hash{$query}{$hit}{"hTotLenght"}." != ".$stat_hash{$hit}{$query}{"qTotLenght"}."\n";}
#            print Dumper    $stat_hash{$query}{$hit};
#             print Dumper    $stat_hash{$hit}{$query};
		  }

		#      exit;
	  }
	print "\tNow comparing the two\n";
  }

sub build_small_extended_orthobench
  {
	my $orthobench_groups =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/tmp_stuff/orthobenchDownloads/IDs_OrthoBench2.txt";
	my $fastaInputDir =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/ensemblv60_12species/longestTranscript/";
	my $fastaOutputDir =
"/Users/fab/Documents/documents/work/current_projects/hieranoid/testdata/ensemblv60_12species/longestTranscript_smallSubset";
	mkdir($fastaOutputDir) if !-e $fastaOutputDir;
	my ( %ID2refOG, %ID2refOG, %refOG_hash, %refOGGenes_hash, %refOG2Members, %pair2RefOGGene_hash,
		 %ID2Gene_hash );
	## OrthoBench groups
	print "\tReading OrthoBench data...";
	getOrthoBenchGroupsPairwise(
								 {
								   orthobench_groups_file => $orthobench_groups,
								   ID2refOG_href          => \%ID2refOG,
								   refOG_href             => \%refOG_hash,
								   refOGGenes_href        => \%refOGGenes_hash,
								   refOG2Members_href     => \%refOG2Members,
								   pair2RefOGGene_href    => \%pair2RefOGGene_hash,
								   ID2Gene_href           => \%ID2Gene_hash
								 }
	);
	print "done\n";

	#my @allFastaFiles = glob("$fastaOutputDir/*.fa");
	#die "No input files\n" if !@allFastaFiles;
	foreach my $fastaFile ( glob("$fastaInputDir/*.fa") )
	  {
		my $newFastaFile = "$fastaOutputDir/" . basename($fastaFile);
		next if -e $newFastaFile && -s $newFastaFile;
		my %id2sequences_hash;
		print "reading sequences from $fastaFile into memory...\n";
		readSequenceFileIntoHash(
								  {
									sequence_file => $fastaFile,
									sequence_href => \%id2sequences_hash
								  }
		);
		print "\t\t\tdone\n";
		my $sequenceString;
		foreach my $seqID ( keys(%id2sequences_hash) )
		  {

			if ( exists $ID2Gene_hash{$seqID} )
			  {
				$sequenceString .= ">$seqID\n" . $id2sequences_hash{$seqID} . "\n";
			  }
		  }
		write_to_file( { file_name => $newFastaFile, text => $sequenceString } );
	  }

	# now add close hits
	foreach my $newFastaFile ( glob("$fastaOutputDir/*.fa") )
	  {
		print "" . basename($newFastaFile) . "\n";

		# read all already existing ids -> to avoid duplications
		my @all_headers = `cat $fastaOutputDir/*.fa | grep ">" `;
		my %all_headers_hash;
		foreach (@all_headers)
		  {
			/>(.*)/;
			$all_headers_hash{$1} = 1;
		  }
		print "\talready having " . keys(%all_headers_hash) . " headers\n";

		#exit;
		my %seqs2keep;
		foreach my $fastaFileAsDB ( glob("$fastaInputDir/*.fa") )
		  {
			my %seqs2keep;

			#next if basename($fastaFileAsDB) eq basename($newFastaFile);
			print "\tadding sequences from " . basename($fastaFileAsDB) . "\n";
			my $FastaFile2Add2 = "$fastaOutputDir/" . basename($fastaFileAsDB);
			print "\t\treading sequences from $fastaFileAsDB into memory...";
			my %id2sequences_hash;
			readSequenceFileIntoHash(
									  {
										sequence_file => $fastaFileAsDB,
										sequence_href => \%id2sequences_hash
									  }
			);
			print "\t\t\tdone\n";
			my $usearch_call =
"time usearch5.1.221_i86osx32 -quiet -query $newFastaFile -db $fastaFileAsDB -evalue 0.01 --userfields query+target+bits+ql+tl+qlo+qhi+tlo+thi+qs+ts --userout bla --maxlen 100000 --minlen 4 --evalue 0.01 --maxrejects 1000 --maxaccepts 100";
			system($usearch_call);
			die "\tno hits for $usearch_call?" if !-e 'bla';
			my @hits = `cat bla`;

			foreach (@hits)
			  {
				chomp;
				next if /query/;
				my @splitted_line = split;
				my (
					 $query_name, $target_name,      $bit_score, $query_length,
					 $hit_length, $query_start,      $query_end, $hit_start,
					 $hit_end,    $query_ali_length, $hit_ali_length
				) = (@splitted_line);
				if ( $bit_score < 200
					 || exists $all_headers_hash{$target_name} )
				  {

					#print "\tskipped $query_name, $target_name with $bit_score\n";
					next;
				  }
				$seqs2keep{$target_name} = $id2sequences_hash{$target_name};
			  }
			print "\tadding " . keys(%seqs2keep) . " sequence from " . basename($fastaFileAsDB) . "\n";
			my $string2print;
			foreach my $id ( keys(%seqs2keep) )
			  {
				$string2print .= ">$id\n" . $seqs2keep{$id} . "\n";
			  }
			attach_to_file( { file_name => $FastaFile2Add2, text => $string2print } );

			#last;
			#exit;
		  }
		exit;
	  }
  }

sub write_to_file
  {
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

sub attach_to_file
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
