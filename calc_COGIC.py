#! /usr/bin/env python
# -*- coding: latin-1 -*-

import sys
from classFunc import *

##################################################
	
# Get parameters
argumentsObject=argumentConfiguration(program="COGIC")
args=argumentsObject.arguments
varStr=argumentsObject.argumentsString
if not args.quiet:
	print varStr


# Parse the Gene Ontology given alternativeId, obsolete, partOf flags
dataGo=GO(args.GeneOntologyOBO,args.linksDag)
statGo=dataGo.statGo()
if not args.quiet:
	print statGo


# Parse the Information Content file
header,goIC=readIC(args.inputInformationContent)
statIc=statIc(goIC,header)
if not args.quiet:
	print statIc


###################################################

# Initialize score class
dataScores=calcScores()
msg='#'*60+'\nGenereating the COGIC scores...\n'
if not args.quiet:
	print msg


# Parse annotations

# dataScores.parseAnnotation() OPTIONS:
# separator			[None,string]		# Separator character as Python split() function do
# posAcc			[None,integer]		# Target name column (first column = 0)
# posGO				[None,integer]		# Go term column (first column = 0)
# posScore			[None,integer]		# Score column (None if is reference)
# labelSet			[None,set()]		# set of labels to filter annotations (i.e. testable/non-testable)
# posLabel			[None,integer]		# label column
# targetList=None	[None,set()]		# restrict the target set
referenceFileFormat=[None,0,1,None,None,None,None] 
referenceSet=dataScores.parseAnnotation(args.referenceSet,*referenceFileFormat)

predictionFileFormat=[None,0,1,2,None,None,None]
predictedSet=dataScores.parseAnnotation(args.predictedSet,*predictionFileFormat)



# Expand ancestors

# dataScores.parseAnnotation() OPTIONS:
# removeRoot			[True,False]		# Separator character as Python split() function do
# goFilter				[None,set()]		# Filter go terms
# expandAncestors		[True,False]		# Expand ancestors
referenceSetExpanded=dataScores.expandAncestorsList(dataGo,referenceSet,removeRoot=True,depthThreshold=args.rootDistance,goFilter=None,expandAncestors=True)

predictedSetExpanded=dataScores.expandAncestorsList(dataGo,predictedSet,removeRoot=False,depthThreshold=args.rootDistance,goFilter=None,expandAncestors=True)

dataScores.referenceSet=referenceSetExpanded
dataScores.predictedSet=predictedSetExpanded




####################################################

# Calculate the COGIC score

# dataScores.calculateCOGIC() OPTIONS:
# filterMissing 		[True,False]	# filter reference targets with only go terms without IC
# outputDetail			[True,False]	# return cogic subvalues (not implemented yet)
# alternativeMethod		[True,False]	# non standrd way to consider TP go terms
# goPredictable			[None,set()]	# set of predictable terms for the predictor ...
#										# ... to count reference targets that contain only NOT ...
#										# ... predictable go terms
dataScores.calculateCOGIC(goIC,dataGo,filterMissing=False,outputDetail=False,alternativeMethod=False,goPredictable=None)


# Write the score output file
scoreString='\n'.join([ target+'\t'+ontology+','+str(dataScores.cogic[target][ontology]) for target in dataScores.cogic for ontology in dataScores.cogic[target].keys()])+'\n'
args.outputScore.write(scoreString)
args.outputScore.close()



# Write output log & statistics
stat=varStr+"\n"
stat+=statGo+"\n"
stat+=statIc+"\n"
stat+=dataScores.statScore()+"\n"		
args.logFile.write(stat)
args.logFile.close()


