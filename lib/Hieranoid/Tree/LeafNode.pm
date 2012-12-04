package LeafNode;
=head1 LeafNode

hieranoid::Tree::LeafNode - Container of taxon objects

=head1 SYNOPSIS

 Abstract class

=head1 DESCRIPTION

Abstract class to handle how orthologous group are summarized

=head1 METHODS

=cut
use Moose;
use Log::Log4perl qw(:easy);

with 'Daughter';

1;