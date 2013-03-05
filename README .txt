
COGIC RELEASE NOTES
====================

COGIC VERSION 1.0



INSTALLATION
============

	$ tar -xzf cogic_pkg.tar.gz

	$ cd cogic_pkg/


EXAMPLE USAGE
=============


If your Python version is updated to the release 2.7 you can use the program providing arguments from the command line (to check Python version try in your terminal "python --version" ). 

For older Pyhon versions you have two options:

1 - you can try to install yourself the built-in Python module "argparse" (see http://pypi.python.org/pypi/argparse);

2 - you can use the configuration file "configure.txt" (self explanatory) and all the arguments provided in the comand line are neglected 



Calculate the Information Content from SwissProt or GoA:
-------------------------------------------------------- 

	$ python calc_IC.py -obo gene_ontology_test.obo -db uniprot_sprot_test.dat

	$ python calc_IC.py -obo gene_ontology_test.obo -db gene_association_test.goa_uniprot


Where "-obo" indicates the Gene Ontology OBO file and "-db" the annotation database: SwissProt in the first example ("uniprot_sprot_test.dat") or GoA in the second ("gene_association_test.goa_uniprot").

The program without other arguments generates: "ic.csv", "annotation.short" and "outCOGIC.log". 
The first file contains the Information Content for each go term and is necessary to calculate the Cogic score. The second file is useful to recalculate the IC (i.e. with a different Gene Ontology or with different parameters), and the third is a log containing some statistics. 



Calculate the Cogic score:
--------------------------

	$ python calc_COGIC.py -obo gene_ontology_test.obo -inIC ic.csv -ref reference_test.txt -pred predictions_test.txt
	
Where "-obo" indicates the Gene Ontology OBO file, "-inIC" the Information Content file previously generated, "-ref" the reference set of annotation and "-pred" the prediction set.

The program by default generates two files: "outCOGIC.txt" with the Cogic score for each target and "outCOGIC.log" with some statistics related to the calculation.

	 




