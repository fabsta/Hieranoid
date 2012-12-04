package InparanoidPrediction;
=head1 InparanoidPrediction

hieranoid::Comparison::OrthologySearch::OrthologyPrediction::InparanoidPrediction - Container of taxon objects

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
use Carp;
use English;
use Log::Log4perl qw(:easy);
use Data::Dumper;
#use Daughter;
use Hieranoid::Comparison::ResultsSimilaritySearch;

## ATTRIBUTES
has 'nodeObject', is => 'rw', isa => 'Object', required => 1;
has 'configuration', is => 'rw', isa => 'Object';
has 'resultsObject', is => 'rw', isa => 'Object', required => 1;

sub BUILD {
;
}
=item start()

to be written...

 Title   : start
 Usage   : ...
 Function: ...
 Returns : 1 on success
 Args: -




=cut      
sub start{
        my $self = shift;
 #        print "\tChecking input files\n";
        INFO("\t\tStarting Orthology prediction\n");
        my $inparanoid_call = $self->configuration->perl." inparanoid.pl -q ".$self->nodeObject->rightDaughter->fileInformation->sequenceSearchInputFile." -d ".$self->nodeObject->leftDaughter->fileInformation->sequenceSearchInputFile;
        if($self->nodeObject->outgroupDaughter){
                $inparanoid_call .= " -g ".$self->nodeObject->outgroupDaughter->fileInformation->sequenceSearchInputFile."";
        }
        INFO("\tStarting Inparanoid..");
        DEBUG("\tInparanoid call: $inparanoid_call\n");
        `$inparanoid_call`;
        INFO("\t\tOrthology prediction finished\n");
}
1;