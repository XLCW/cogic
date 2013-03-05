#! /usr/bin/env python
# -*- coding: latin-1 -*-

import sys
from classFunc import *

# Get parameters
argumentsObject=argumentConfiguration(program="IC")
args=argumentsObject.arguments
varStr=argumentsObject.argumentsString
if not args.quiet:
	print varStr

# Parse the Gene Ontology given alternativeId, obsolete, partOf flags
dataGo=GO(args.GeneOntologyOBO,args.linksDag)
statGo=dataGo.statGo()
if not args.quiet:
	print statGo

# Information Content calculation
dbData=DB(args.evidenceCodeList,args.ontologyList,args.alternativeId,args.obsolete)
# Parse the annotation input file (SwissProt, GOA or "short")
dbData.parser(args.inputDataBase,dataGo,args.outputShortFile)
statIc=dbData.statIC()
if not args.quiet:
	print statIc


# Write the IC calculation
dbData.writeIC(args.outputInformationContent)


# Write the log file
stat=varStr
stat+=statGo
stat+=statIc
args.logFile.write(stat)



