package ResultsSimilaritySearch;
use Moose;
use Carp;
use English;
use Log::Log4perl qw(:easy);
use Data::Dumper;
#use Daughter;


## ATTRIBUTES
has 'nodeObject', is => 'rw', isa => 'Object';
has 'configuration', is => 'rw', isa => 'Object';

has 'simAA', is => 'rw', isa => 'Str', required => 1;
has 'simAB', is => 'rw', isa => 'Str', required => 1;
has 'simBA', is => 'rw', isa => 'Str', required => 1;
has 'simBB', is => 'rw', isa => 'Str', required => 1;
has 'simAC', is => 'rw', isa => 'Str';
has 'simBC', is => 'rw', isa => 'Str';




sub BUILD {
      my $self = shift;
      DEBUG("\t\tResults from similarity search saved\n");
      }      
1;