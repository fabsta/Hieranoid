package Daughter;
=head1 Daughter

hieranoid::Tree::Daughter - Container of taxon objects

=head1 SYNOPSIS

 Abstract class

=head1 DESCRIPTION

Abstract class to handle how orthologous group are summarized

=head1 METHODS

=cut
use Moose::Role;
use Log::Log4perl qw(:easy);



## ATTRIBUTES
has 'nodeObject', is => 'rw', isa => 'Object', required => 1;
has 'name', is => 'rw', isa => 'Str', required => 1;
has 'fileInformation', is => 'rw', isa => 'Object';


1;