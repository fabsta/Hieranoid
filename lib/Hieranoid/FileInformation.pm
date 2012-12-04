package FileInformation;
use Moose;
use Carp;
use English;
use Hieranoid::Tree::Daughter;
use Hieranoid::Tree::GroupFile;

#extends 'Daughter';
## ATTRIBUTES
has 'orthologGroups', is => 'rw', isa => 'Object';
has 'sourceFile', is => 'rw', isa => 'Str';
has 'orthologPrediction4Node', is => 'rw', isa => 'Str';
has 'sequenceInputFormat', is => 'rw', isa => 'Str';
has 'speciesFilesDirectory', is => 'rw', isa => 'Str';
has 'sequenceSearchInputFile', is => 'rw', isa => 'Str',required => 1;
has 'profileDirectory', is => 'rw', isa => 'Str';
has 'profileFile', is => 'rw', isa => 'Str';
has 'hmmFile', is => 'rw', isa => 'Str';
has 'csdbFile', is => 'rw', isa => 'Str';
has 'consensusDirectory', is => 'rw', isa => 'Str';
has 'consensusFile', is => 'rw', isa => 'Str';
has 'alignmentFile', is => 'rw', isa => 'Str';
has 'profileSearchInputFile', is => 'rw', isa => 'Str';
has 'outputDirectory', is => 'rw', isa => 'Str';


1;