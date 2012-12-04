Hieranoid
Hierarchical orthology inference

This software is freely available at https://github.com/fabsta/Hieranoid 

Please read the COPYING file before using this software.  

If you use Hieranoid please cite:
Fabian Schreiber and Erik L.L. Sonnhammer
Hieranoid: Hierarchical orthology inference 
Journal of Molecular Biology (submitted)


Installation
Checking system requirements
-	Perl
-	BioPerl
-	HHSearch
-	Muscle
-	Kalign

Perl modules
Hieranoid uses Perl modules that are installed in the 'lib' sub-directory of the Hieranoid installation.
Please add that directory to your PERL5LIB enviromental variable. In bash you do that the following way:
	export PERL5LIB=$PERL5LIB:hieranoid_directory/lib

Configuring Hieranoid
Once you have downloaded Hieranoid, you have to tell Hieranoid where to find the required programs. 
This is done by adapting the settings in the Configuration file (Configurations/Configuration.pm).
All options should be explained in greater detail in the Configuration file


Starting hieranoid
Once, the configuration file has been adapted, the analysis can be started by typing:
	perl hieranoid.pl

Note that depending on the input size, the analysis can take some time


Command line arguments
Additional options:
	-c configurationFile
	-t guide tree file


Test dataset:
To test the installation, simply type “perl hieranoid.pl”.
This will run a quick analysis on the sequences in data/sequences
using species tree in data/tree.
Program output:
All program output is stored in a directory ‘project1’. For each inner node of the provided species tree, a folder will be created in ‘project1/nodes’.
Those folders contain several files: The input query files (e.g. Homo_sapiens.fa), Blast results (e.g. Homo_sapiens-Pan_troglodytes), the InParanoid output (sqltable.Homininae) and the orthology predictions (e.g. Homininae.OGTree.txt).

Where to put my files:
Hieranoid requires a set of sequences as well as a guide tree. By default the Hieranoid searches for the sequences in data/sequences and for the guide tree
in data/tree. These options can be changed in the Configuration file.


Output:
The following is an example of the Hieranoid output format:
HieranoidOG-Eukaryota3:((At1g50030.1:1.000,((ENSPTRP00000000280:1.000,ENSP00000354558:1.000)Homininae1:1.000,ENSMUSP00000099510:1.000)Euarchontoglires1:1.000)Eukaryota3:1.000)root;
The first is the name of the Hieranoid group followed by a newick string with inparalog scores
 




