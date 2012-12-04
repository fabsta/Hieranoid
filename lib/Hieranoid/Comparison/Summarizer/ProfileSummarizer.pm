package ProfileSummarizer;
=head1 ProfileSummarizer

hieranoid::Comparison::Summarizer::ProfileSummarizer - Container of taxon objects

=head1 SYNOPSIS

 Abstract class

=head1 DESCRIPTION

Abstract class to handle how orthologous group are summarized

=head1 METHODS

=cut
use Moose;
use Hieranoid::Tree::InnerNode;
use Hieranoid::FileInformation;
with 'Summarizer';

1;