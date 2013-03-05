# -*- coding: latin-1 -*-

import sys
from math import log
import copy

# The three roots
#rootCodes={'GO:0008150':'P','GO:0005575':'C','GO:0003674':'F'}	
#codesRoot={'P':'GO:0008150','C':'GO:0005575','F':'GO:0003674'}
ontDict={'biological_process':'P','molecular_function':'F','cellular_component':'C'}
	
class fileArgs():
	def __init__(self,program,configFileName):
		f=open(configFileName,"r")
		for line in f.readlines():
			if line.strip():
				if line[0]!="#":
					exec(line.strip())
		f.close()
		
		self.GeneOntologyOBO=open(oboFileName,"rt")
		self.linksDag=links
		self.obsolete=obsoleteFlag
		self.alternativeId=alternativeFlag
		self.quiet=quietFlag
		self.logFile=outputLogFileName
		
		if program=="IC":
			self.inputDataBase=open(dbFileName,"rt")
			self.ontologyList=ontologyList
			self.evidenceList=evidenceList
			self.outputShortFile=outputDBFileName
			self.outputInformationContent=outputICFileName
		elif program=="COGIC":
			self.inputInformationContent=open(inputICFileName,"rt")
			self.referenceSet=open(inputReference,"rt")
			self.predictedSet=open(inputPredicted,"rt")
			self.rootDistance=rootDistanceFlag
			self.outputScore=outputCogicScore
		# Unknown program
		else:
			print "Missing program:",program
			sys.exit(1)
		

class argumentConfiguration:
	def __init__(self,configFileName=None,program=None):
		
		self.argumentsString=None
		self.arguments=None
		
		
		# Command line argument parsing require the Python built-in "argparse" module

		# Argparse module check
		try:
			
			import argparse
			args=self.argparse_Parser(program)

		except ImportError:
			
			print '\nWARNING: Argparse not available, your Python version is not up to date !!\n\nEdit "configure.txt" to define input parameters ...\n'
			
			args=fileArgs(program,"configure.txt")

		self.arguments=self.processArguments(args)
		self.argumentsString=self.makeArgumentStr(args)

	def argparse_Parser(self,program):
		"""
		Parse the argument list with Python Argparse module
		"""
		import argparse
		
		# Parser implementation
		parserCommon=argparse.ArgumentParser(add_help=False)

		# Common parser (GO parser)
		groupOBO = parserCommon.add_argument_group('Gene Ontology parsing')
		groupOBO.add_argument('-obo','--GeneOntologyOBO',required=True,help='Gene Ontology OBO file',type=argparse.FileType('rt'))
		groupOBO.add_argument('-links','--linksDag',help='Type of links used when building the DAG',choices=['is_a','part_of','both'],default='is_a')
		groupOBO.add_argument('-obs','--obsolete',help='Obsolete flag. Add the obsolete terms remapping ancestors from the corresponding terms of "consider" and "replaced_by" fields in the OBO file',  default=True)
		groupOBO.add_argument('-altId','--alternativeId',help='Alternative ID flag. Add terms from the "alt_id" field in the OBO file', default=True)
		parserCommon.add_argument('-q','--quiet',help='Remove output verbosity', action='store_true')
		parserCommon.add_argument('--version', action='version', version='%(prog)s 1.0')
		parserCommon.add_argument('-log','--logFile',help='Print some statistics in a file',default='outCOGIC.log')

		
		# Subparser
		parser = argparse.ArgumentParser(parents=[parserCommon])
		#subparsers = parser.add_subparsers(help='Programs')


		# DB parser
		if program=="IC":
		
			groupSP = parser.add_argument_group('SwissProt parsing')
			groupSP.add_argument('-db','--inputDataBase',required=True,help='Database of annotation: SwissProt flat file (.dat), GoA database or the "short" file (the "short" file is automatically generated after the first run)',type=argparse.FileType('rt'))
			groupSP.add_argument('-ont','--ontologyList',help='Ontology categories to be considered. F = Molecular Function, P = Biological Process, C = Cellular Component',default='FPC')
			
			group = groupSP.add_mutually_exclusive_group()
			group.add_argument('-evL','--evidenceList',help='Evidence codes considered (i.e. TAS,IC,IDA)')
			group.add_argument('-evX','--evidenceExclude',help='Exclude these evidence codes')
			group.add_argument('-evSet','--evidenceSet',help='Preset of evidence codes', choices=['experimental','computational','author','curator','automatic','cafa'])
			
			groupSPout = parser.add_argument_group('Output files')
			groupSPout.add_argument('-outShort',dest='outputShortFile',help='"short" version output file. Can be used as input instead of the annotation database file to speed up calculations')
			groupSPout.add_argument('-outIC','--outputInformationContent',help='Information Content output file',default='ic.csv')

		# COGIC parser
		elif program=="COGIC":
			
			groupCogic = parser.add_argument_group('COGIC parameters')
			
			groupCogic.add_argument('-inIC','--inputInformationContent',required=True,help='Infomation Content file. Use IC program to generate it',type=argparse.FileType('rt'))
			groupCogic.add_argument('-ref','--referenceSet',required=True,help='Reference annotation file',type=argparse.FileType('rt'))
			groupCogic.add_argument('-pred','--predictedSet',required=True,help='Predicted annotation file',type=argparse.FileType('rt'))
			groupCogic.add_argument('-d','--rootDistance',help='Filter GO term based on the distance from root',default=None)
			
			groupCOGICout = cogic_parser.add_argument_group('Output files')
			groupCOGICout.add_argument('-outS','--outputScore',help='Score output file',default='outCOGIC.txt')

		# Unknown program
		else:
			print "Missing program:",program
			sys.exit(1)

		# Initialize the parser object
		args=parser.parse_args()

		return args
		

	def argparse_MergedParser(self):
		"""
		Parse the argument list with Python Argparse module
		"""
		import argparse
		
		# Parser implementation
		parserCommon=argparse.ArgumentParser(add_help=False)

		# Common parser (GO parser)
		groupOBO = parserCommon.add_argument_group('Gene Ontology parsing')
		groupOBO.add_argument('-obo','--GeneOntologyOBO',required=True,help='Gene Ontology OBO file',type=argparse.FileType('rt'))
		groupOBO.add_argument('-links','--linksDag',help='Type of links used when building the DAG',choices=['is_a','part_of','both'],default='is_a')
		groupOBO.add_argument('-obs','--obsolete',help='Obsolete flag. Add the obsolete terms remapping ancestors from the corresponding terms of "consider" and "replaced_by" fields in the OBO file',  default=True)
		groupOBO.add_argument('-altId','--alternativeId',help='Alternative ID flag. Add terms from the "alt_id" field in the OBO file', default=True)
		parserCommon.add_argument('-q','--quiet',help='Remove output verbosity', action='store_true')
		parserCommon.add_argument('--version', action='version', version='%(prog)s 1.0')
		parserCommon.add_argument('-log','--logFile',help='Print some statistics in a file',default='outCOGIC.log')

		
		# Subparser
		parser = argparse.ArgumentParser()
		subparsers = parser.add_subparsers(help='Programs')

		# DB parser
		db_parser = subparsers.add_parser('IC', help='Calculate the Information Content',parents=[parserCommon],formatter_class=argparse.ArgumentDefaultsHelpFormatter)
		groupSP = db_parser.add_argument_group('SwissProt parsing')
		groupSP.add_argument('-db','--inputDataBase',required=True,help='Database of annotation: SwissProt flat file (.dat), GoA database or the "short" file (the "short" file is automatically generated after the first run)',type=argparse.FileType('rt'))
		groupSP.add_argument('-ont','--ontologyList',help='Ontology categories to be considered. F = Molecular Function, P = Biological Process, C = Cellular Component',default='FPC')
		group = groupSP.add_mutually_exclusive_group()
		group.add_argument('-evL','--evidenceList',help='Evidence codes considered (i.e. TAS,IC,IDA)')
		group.add_argument('-evX','--evidenceExclude',help='Exclude these evidence codes')
		group.add_argument('-evSet','--evidenceSet',help='Preset of evidence codes', choices=['experimental','computational','author','curator','automatic','cafa'])
		groupSPout = db_parser.add_argument_group('Output files')
		groupSPout.add_argument('-outShort',dest='outputShortFile',help='"short" version output file. Can be used as input instead of the annotation database file to speed up calculations')
		groupSPout.add_argument('-outIC','--outputInformationContent',help='Information Content output file',default='ic.csv')

		# COGIC parser
		cogic_parser = subparsers.add_parser('COGIC', help='Calculate the COGIC score',parents=[parserCommon],formatter_class=argparse.ArgumentDefaultsHelpFormatter,description='The program calculates the COGIC score between two set of GO terms (predicted/reference)')
		cogic_parser.add_argument('-inIC','--inputInformationContent',required=True,help='Infomation Content file. Use IC program to generate it',type=argparse.FileType('rt'))
		cogic_parser.add_argument('-ref','--referenceSet',required=True,help='Reference annotation file',type=argparse.FileType('rt'))
		cogic_parser.add_argument('-pred','--predictedSet',required=True,help='Predicted annotation file',type=argparse.FileType('rt'))
		cogic_parser.add_argument('-d','--rootDistance',help='Filter GO term based on the distance from root',default=False)
		groupCOGICout = cogic_parser.add_argument_group('Output files')
		groupCOGICout.add_argument('-outS','--outputScore',help='Score output file',default='outCOGIC.txt')


		# Initialize the parser object
		args=parser.parse_args()

		return args


	def processArguments(self,args):
		# Preset of Gene Ontology evidence codes
		evidenceDict={'cafa':set(['EXP','TAS','IC','IDA','IMP','IGI','IEP']),\
					'experimental':set(['EXP','IDA','IPI','IMP','IGI','IEP']),\
					'computational':set(['ISS','ISO','ISA','ISM','IGC','IBA',\
					'IBD','IKR','IRD','RCA']),\
					'author':set(['TAS','NAS']),\
					'curator':set(['IC','ND']),\
					'automatic':set(['IEA']),\
					'all':set(['IC', 'IKR', 'EXP', 'TAS', 'ISS', 'IBA', 'IPI', 'IBD', 'ND', 'NAS', 'IGC', 'IGI', 'ISM', 'IEP', 'IEA', 'ISO', 'ISA', 'IRD', 'RCA', 'IMP', 'IDA'])}
		allEvidenceCodes=set([cod for ele in evidenceDict for cod in evidenceDict[ele]])

		# Setting ontology list
		if hasattr(args, 'ontologyList'):
			args.ontologyList=set([c.upper() for c in args.ontologyList])
			if args.evidenceList:
				evidenceCodeList=set(args.evidenceList.split(','))
			elif args.evidenceSet:
				evidenceCodeList=evidenceDict[args.evidenceSet]
			elif args.evidenceExclude:
				evidenceCodeList=allEvidenceCodes-set(args.evidenceExclude.split(','))
			else:
				evidenceCodeList=allEvidenceCodes
			args.evidenceCodeList=evidenceCodeList

		# DEBUG: to avoid unwanted files overwriting, open file in writing mode only when corresponding functions are called
		# Open files
		args.logFile=open(args.logFile,'wt')
		if hasattr(args, 'outputInformationContent'):
			args.outputInformationContent=open(args.outputInformationContent,'wt')
			if not args.outputShortFile:
				args.outputShortFile=args.inputDataBase.name.split('/')[-1]+'.short'
		if hasattr(args, 'outputScore'):
			args.outputScore=open(args.outputScore,'wt')
			
		return args
		


	def makeArgumentStr(self,args):
		# List all argument in a formatted string
		varStr='\nCOGIC 1.0\n\n'+'#'*60+'\nVariable assignments:\n'+'\n'.join(['{0:>24} = {1}'.format(arg,str(args.__dict__[arg])) for arg in args.__dict__ if not args.__dict__[arg]==None and not isinstance(args.__dict__[arg],file)])+'\n'+'\n'.join(['{0:>24} = {1}'.format(arg,str(args.__dict__[arg].name)) for arg in args.__dict__ if not args.__dict__[arg]==None and isinstance(args.__dict__[arg],file)])+'\n'

		return varStr



#######################################################################################

class GO:
	"""
	Parse the Gene Ontology OBO file and create an object
	"""	
	def __init__(self,oboFile,linksDag):

		self.oboFile=oboFile

		self.linksDag=linksDag	

		# { go : [ ontology, set(is_a parents), set(part_of parents), originalId, depth, set(ancestors) ] }
		self.goData={}
		# { go : [ ontology, set("replaced by" terms), set("consider" terms) ] }
		self.obsolete={}		
		self.alternativeId={}				# { synonym_ID : original_ID}		
		self.roots={}

		self.nTerms={'F':0,'P':0,'C':0}		# Number of terms per ontology
		self.missingRoot={'F':set(),'P':set(),'C':set()}	# Terms without path to root

		self._parseGoOBO() 	# Parse the OBO file and expand the ancestors

	def _parseGoOBO(self):
		"""
		Parse the Gene Ontology OBO file reading the followinf fields: id, namespace,
		alt_id, is_obsolete, consider, replaced_by, is_a, relashionship.
		Update:
		- self.alternativeId
		- self.obsolete
		- self.goData
		- self.leaves
		- rootCodes (gloabal variable)
		"""
		accession=None	
		parentsIs=set()
		parentsPart=set()
		synonym=set()
		obsoleteReplaced=set()
		obsoleteConsider=set()
		ontology=None
		obsolete=False
		start=False
		
		f=self.oboFile
		line=f.readline() 
		while line:
			if len(line.strip())==0 and start==True:
				if accession:
					if not obsolete:
						self.goData[accession]={'ontology':ontDict[ontology],\
								"parentsIs":parentsIs,"parentsPart":parentsPart,"originalId":None,\
								"ancestors":set(),"isRoot":False,"name":name,"minDepth":0,\
								"childrenIs":set(),"childrenPart":set(),"isLeaf":False}
						if not parentsIs and not parentsPart:
							self.goData[accession]["isRoot"]=True
							self.roots[accession]=ontDict[ontology]

						if synonym:
							for synonymAccession in synonym:
								self.alternativeId[synonymAccession]=accession
								
					else:
						self.obsolete[accession]={"ontology":ontDict[ontology],"replaced":obsoleteReplaced,"consider":obsoleteConsider}

				accession=None	
				name=None
				parentsIs=set()
				parentsPart=set()
				synonym=set()
				obsoleteReplaced=set()
				obsoleteConsider=set()
				ontology=None
				obsolete=False
				start=False

			elif line.startswith('[Term]'):
				start=True
			else:
				if line.startswith('id:'):
					accession=line.strip().split()[-1]
				elif line.startswith('name:'):
					name=line.strip()[6:]
				elif line.startswith('namespace:'):
					ontology=line.strip().split()[-1]
				elif line.startswith('alt_id:'):
					synonym.add(line.strip().split()[-1])
				elif line.startswith('is_obsolete: true'):
					obsolete=True
				elif line.startswith('consider:'):
					obsoleteConsider.add(line.strip().split()[1])
				elif line.startswith('replaced_by:'):
					obsoleteReplaced.add(line.strip().split()[1])
				elif line.startswith('is_a'):
					parentsIs.add(line.split()[1])
				elif line.startswith('relationship:'): # or line.startswith('intersection_of:'):
					parentsPart.add(line.split()[2])
			line=f.readline()
		f.close()

		## Count the number of terms for each GO branch
		for go in self.goData:
			self.nTerms[self.goData[go]["ontology"]]+=1

		## Add ancestors (append minimun root distance and ancestors to goData)
		self._buildDag()

		## DEBUG: consider only is_a
		## Find leaf terms
		allParents=set([parent for go in self.goData for parent in self.goData[go]["parentsIs"]])
		for go in self.goData:
			if go not in allParents and not self.goData[go]["isRoot"]:
				self.goData[go]["isLeaf"]=True
		
	
	def _buildDag(self):
		"""
		The function retrieves the ancestors of a term given the goData structure. 
		The function simply iteratively looks for parents of each term untill 
		the root node is reached (the node without any parent)
		Update goData dictionary appending two elements for each go term:
		- the minimum distance from the root
		- the set of ancestors
		"""
		for go in self.goData:
			#print go,self.goData[go]["parentsPart"],self.goData[go]["parentsIs"]	

			parents=set()
			if self.linksDag=='is_a':
				parents=copy.deepcopy(self.goData[go]["parentsIs"])
			elif self.linksDag=='part_of':
				parents=copy.deepcopy(self.goData[go]["parentsPart"])
			elif self.linksDag=='both':
				parents=copy.deepcopy(self.goData[go]["parentsIs"].union(self.goData[go]["parentsPart"]))
			ancestors=copy.deepcopy(parents)
			
			# Build the paths to the root and measure the minimum root distance
			minDepth=0 # Minumum distance from the root
			depth=[]
			c=0
			while parents:
				c+=1
				# Evaluate if the root is reached (to calculate root distance)
				if [goR for goR in parents if self.goData[goR]["isRoot"]]:
					depth.append(c)

				newParents=set()
				for goP in parents:
					if self.linksDag=='is_a':
						newParents=copy.deepcopy(self.goData[goP]["parentsIs"])
					elif self.linksDag=='part_of':
						newParents=copy.deepcopy(self.goData[goP]["parentsPart"])
					elif self.linksDag=='both':
						newParents=copy.deepcopy(self.goData[goP]["parentsIs"].union(self.goData[goP]["parentsPart"]))
					
				parents=copy.deepcopy(newParents)
				ancestors.update(parents)

			if depth:
				minDepth=min(depth)
			
			# Add the root when missing (when using only "part_of" links)
			if not [goR for goR in ancestors if self.goData[goR]["isRoot"]]:
				if go not in self.roots: 
					ontology=self.goData[go]["ontology"]
					
					root=[root for root in self.roots if self.roots[root]==ontology][0]
					ancestors.add(root)
					self.missingRoot[ontology].add(go)

			self.goData[go]["minDepth"]=minDepth
			self.goData[go]["ancestors"]=ancestors
			
			


	def statGo(self):
		sortedList=['Terms','Leaves','Max term depth','Average term depth',\
		'Average leaves depth','Min leaves depth','Obsolete terms',\
		'Obsolete, only "consider" field','Obsolete, only "replaced_by" field',\
		'Obsolete, "consider"+"replaced_by"','Obsolete not mapped','Terms without root path']
		stat={}
		stat['Terms']=self.nTerms
		stat['Max term depth']={'P':0,'C':0,'F':0}	
		stat['Average term depth']={'P':0,'C':0,'F':0}
		stat['Leaves']={'P':0,'C':0,'F':0}
		stat['Average leaves depth']={'P':0,'C':0,'F':0}
		stat['Min leaves depth']={'P':0,'C':0,'F':0}
		stat['Obsolete terms']={'P':0,'C':0,'F':0}		
		stat['Obsolete, only "consider" field']={'P':0,'C':0,'F':0}	
		stat['Obsolete, only "replaced_by" field']={'P':0,'C':0,'F':0}	
		stat['Obsolete, "consider"+"replaced_by"']={'P':0,'C':0,'F':0}	
		stat['Obsolete not mapped']={'P':0,'C':0,'F':0}		
		stat['Terms without root path']={'P':0,'C':0,'F':0}		
		
		for accession in self.obsolete:
			ontology=self.obsolete[accession]['ontology']
			
			replaced=self.obsolete[accession]['replaced']
			consider=self.obsolete[accession]['consider']
			
			stat['Obsolete terms'][ontology]+=1
			if not replaced and consider:
				stat['Obsolete, only "consider" field'][ontology]+=1
			elif replaced and not consider:
				stat['Obsolete, only "replaced_by" field'][ontology]+=1
			elif replaced and consider:
				stat['Obsolete, "consider"+"replaced_by"'][ontology]+=1
			else:
				stat['Obsolete not mapped'][ontology]+=1

		for ontology in self.missingRoot:
			stat['Terms without root path'][ontology]=len(self.missingRoot[ontology])

		depth={'P':[],'C':[],'F':[]}
		depthAll={'P':[],'C':[],'F':[]}
		for go in self.goData:
			
			ontology=self.goData[go]['ontology']
			minDepth=self.goData[go]["minDepth"]
			
			if self.goData[go]["isLeaf"]:
				stat['Leaves'][ontology]+=1
				depth[ontology].append(minDepth)
			depthAll[ontology].append(minDepth)
		
		for ontology in stat['Max term depth']:
			stat['Min leaves depth'][ontology]=min(depth[ontology])
			stat['Max term depth'][ontology]=max(depth[ontology])
			stat['Average leaves depth'][ontology]=round(float(sum(depth[ontology]))/len(depth[ontology]),2)


		for ontology in stat['Average term depth']:
			stat['Average term depth'][ontology]=round(float(sum(depthAll[ontology]))/len(depthAll[ontology]),2)
		
		result='#'*60+'\nGene Ontology OBO statistics:\n'+'\n'.join(['{0:>34} = {1}'.format(k,' '.join(['{0:>5} {1}'.format(str(stat[k][ont]),ont) for ont in stat[k]])) for k in sortedList])+'\n'

		return result

###################################################################################


###################################################################################
class DB:
	"""
	Calculate the Information Content of GO terms 
	from a database of annotation (SwissProt/Trembl or GOA)
	"""
	def __init__(self,evidenceCodeList,ontologyList,alternativeId=True,obsolete=True):
		
		self.evidenceCodeList=evidenceCodeList
		self.ontologyList=ontologyList
		self.alternativeFlag=alternativeId
		self.obsoleteFlag=obsolete
		self.fileSource=None	
		
		
				
		self.nR={'F':0,'P':0,'C':0}			# Number of root occurences
		self.nN={'F':{},'P':{},'C':{}}		# Number of term occurences AFTER expantion
		self.goIC={'F':{},'P':{},'C':{}}	# Information Content	
		self.annotation={}							# Annotation data	

		self.ontologyDist={'F':0,'P':0,'C':0}	# Number of terms per ontology
		self.goMissing=set()				# Terms not present in the Gene Ontology
		self.evidenceCodeDist={}			# Numebr of terms per evidence code
		self.nAccession=0					# Number of annotated entries


	def parser(self,dbFile,dataGo,outputShortFile=None,dataFlag=None):
		"""
		Select the correct parser for the annotation database and call 
		the function to calculate the Information Content
		"""
		self.fileSource=dbFile.name
		sprotFlag=self._dbType(dbFile)

		# Collect GO annotations
		if sprotFlag=='sprot':
			self.parseSwissProt(dbFile,dataGo,outputShortFile,dataFlag) 
		elif sprotFlag=='goa':
			self.parseGoa(dbFile,dataGo,outputShortFile,dataFlag)
		elif sprotFlag=='short':
			self.readShortDb(dbFile,dataGo,outputShortFile,dataFlag)
		else:
			'\nERROR: Database input file in a bad format! Use SwissProt/GoA/"Short"\n'
			sys.exit(1)

		# Calculate the IC 
		self._calculate_IC(dataGo)
		

	def _dbType(self,sprotFile):
		"""
		Determine the input database format from the first two lines
		- GOA
		- SwissProt/Trembl : text ".dat" file
		- "Short" : see format specification
		"""
		fileFormat=None
		line=sprotFile.readline()
		if line[0:2]=='ID' and line[-2]=='.':
			line=sprotFile.readline()
			if line[0:2]=='AC': 
				fileFormat='sprot'
		elif line[0:2]=='##':
			line=sprotFile.readline()
			line=line.strip().split('\t')
			if len(line[0].split(',')[0])==6 and line[1].split()[0].split(':')[0]=='GO':
				fileFormat='short'
		elif line[0:4]=='!gaf':
			line=sprotFile.readline()
			line=line.strip().split('\t')
			if line[0]=='UniProtKB' and len(line[1])==6:
				fileFormat='goa'
		sprotFile.seek(0)

		return fileFormat


	def parseSwissProt(self,sprotFile,dataGo,shortFile,dataFlag):
		"""
		Parse the SwissProt flat file (".dot" format) 
		Calculate the Information Content
		Write annotation in a file ("short" format)
		The GO annotation is extracted from the "DR" Database cross-references fields
		when the "RESOURCE_ABBREVIATION" is equal to "GO;" 
		(See "http://web.expasy.org/docs/userman.html#DR_line")	
		"""	

		## Write the header of the short version output file
		if shortFile:
			f=open(shortFile,'wt')
			f.write('## Source:'+self.fileSource+'\tEvidenceCodes:'+','.join(list(self.evidenceCodeList))+'\tOntology:'+','.join(list(self.ontologyList))+'\n')


		## Parse the SwissProt/Trembl file
		## Update self.nN and self.nR for the IC calculation
		goSet=[]
		synonym=set()
		start=True
		accession=None
		synonymMapping={}	# { accession : set( synonym accessions ) }	

		line=sprotFile.readline() 
		while line:
			if line.startswith('//'):
				if goSet:	
					self.nAccession+=1
					# Update self.nR and self.nN for the IC calculation
					self._expandAncestorCount(goSet,dataGo)
					# Update the synonym dictionary
					if synonym:
						synonymMapping[accession]=synonym
					# Write the "short" version file
					if shortFile:
						f.write(accession+'\t'+' '.join([','.join(list(ele)) for ele in goSet])+'\n')
					if dataFlag:
						self.annotation[accession]=dict([(goTmp,1.0) for goTmp,ont,evcode in goSet])	
				goSet=[]
				synonym=set()
				start=True
			else:
				if line.startswith('AC   '):
					if start==True:
						accession=line.split()[1][0:-1]
						synonym.update(set([ele[0:-1] for ele in line.split()[2:]]))
						start=False
					else:
						synonym.update(set([ele[0:-1] for ele in line.split()[1:]]))
				elif line.startswith('DR   GO;'):
					line=line.split(';')
					evidenceCode=line[-1].split(':')[0].strip()
					ontology=line[2].split(':')[0].strip()
					# filter annotations based on evidenceCodeList and ontologyList
					if evidenceCode in self.evidenceCodeList and ontology in self.ontologyList:
						goSet.append((line[1].strip(),ontology,evidenceCode))
			line=sprotFile.readline()
		sprotFile.close()
		if shortFile:
			f.close()
		
		
	def parseGoa(self,sprotFile,dataGo,shortFile,dataFlag):
		"""
		Parse the GOA flat file and calculate the Information Content 
		and write the annotation in a file ("short" format).
		("ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz")	
		"""
		## Write the header in the short version output file
		if shortFile:
			f=open(shortFile,'wt')
			f.write('## Source:'+self.fileSource+'\tEvidenceCodes:'+','.join(list(self.evidenceCodeList))+'\tOntology:'+','.join(list(self.ontologyList))+'\n')

		oldAccession=None
		goSet=set()

		line=sprotFile.readline()		
		while line:
			if line[0:9]=='UniProtKB':
				line=line.strip().split('\t')
				accession=line[1]
				go=line[4]
				evidenceCode=line[6]
				ontology=line[8]
				if evidenceCode in self.evidenceCodeList and ontology in self.ontologyList:
					if accession!=oldAccession:					
						if oldAccession:
							if goSet:
								self.nAccession+=1
	
								# Write the "short" version file
								if shortFile:
									f.write(oldAccession+'\t'+' '.join([','.join(list(ele)) for ele in goSet])+'\n')
								
								if dataFlag:
									self.annotation[accession]=dict([(goTmp,1.0) for goTmp,ont,evcode in goSet])	
								# Update self.nR and self.nN for the IC calculation
								self._expandAncestorCount(goSet,dataGo)			

							goSet=set()
					goSet.add((go,ontology,evidenceCode))
				oldAccession=accession
			line=sprotFile.readline()
		sprotFile.close()
		if goSet:
			self.nAccession+=1
			# Write the "short" version file (last record)
			if shortFile:
				f.write(oldAccession+'\t'+' '.join([','.join(list(ele)) for ele in goSet])+'\n')
				f.close()
			if dataFlag:
				self.annotation[accession]=dict([(goTmp,1.0) for goTmp,ont,evcode in goSet])	
			# Update self.nR and self.nN for the IC calculation (last record)
			self._expandAncestorCount(goSet,dataGo)	

	def readShortDb(self,shortFileInput,dataGo,shortFile,dataFlag):
		"""
		Read the "short" version of the annotation database. Generated when parsing the
		original databases (GOA or SwissProt/Trembl). See format specification.
		Update: 
		- self.nAccession
		- self.filesource
		"""
		

		line=shortFileInput.readline().strip().split('\t')
		
		self.fileSource=line[0].split(':')[1]
		evidenceCodeList=set(line[1].split(':')[1].split(','))

		if not self.evidenceCodeList.issubset(evidenceCodeList):
			print '\nERROR: Annotations missing in the "short" file (missing evidence codes)! Parse again the original SwissProt/GoA database:\nPresent:'+','.join(list(evidenceCodeList))+'\nRequested:'+','.join(list(self.evidenceCodeList))+'\nMissing:'+','.join([ele for ele in (self.evidenceCodeList-evidenceCodeList)])+'\n'
			sys.exit(1)


		## Write the header of the short version output file
		if shortFile:
			f=open(shortFile,'wt')
			f.write('## Source:'+self.fileSource+'\tEvidenceCodes:'+','.join(list(self.evidenceCodeList))+'\tOntology:'+','.join(list(self.ontologyList))+'\n')

		line=shortFileInput.readline()
		
		while line:
			self.nAccession+=1
			goSet=set()
			accession=line.split('\t')[0]
			for ele in line.split('\t')[1].split():
				go,ontology,evidenceCode=ele.split(',')
				if ontology in self.ontologyList and evidenceCode in self.evidenceCodeList:
					goSet.add((go,ontology,evidenceCode))
			
			# Write the "short" version file
			if shortFile and goSet:
				f.write(accession+'\t'+' '.join([','.join(list(ele)) for ele in goSet])+'\n')
			if dataFlag:
				self.annotation[accession]=dict([(goTmp,1.0) for goTmp,ont,evcode in goSet])

			self._expandAncestorCount(goSet,dataGo)
			line=shortFileInput.readline()

		shortFileInput.close()


		
	def _expandAncestorCount(self,goSet,dataGo):
		"""
		Retrieve ancestors for each GO term of goSet and store them in a list 
		KEEPING duplicated ancestors. 
		Update: 
		- self.nN
		- self.nR
		- self.evidenceCodeDist
		- self.ontologyDist
		- self.goMissing
		"""
		# update altId in goSet
		if self.alternativeFlag:
			goSet=self._updateAlternative(goSet,dataGo)

		# update obsolete in goSet
		if self.obsoleteFlag:
			goSet=self._updateObsolete(goSet,dataGo)

		ancestorList=[]
		# Generate the ancestor list for each GO term
		for ele in goSet:
			go,ontology,evidence=ele

			ancestorList.append(ele)
			ancestorList.extend([(ancestor,ontology,evidence) for ancestor in dataGo.goData[go]["ancestors"]])

			# Statistic x evidence code
			if evidence not in self.evidenceCodeDist:
				self.evidenceCodeDist[evidence]=0
			self.evidenceCodeDist[evidence]+=1
			self.ontologyDist[ontology]+=1

		# Count the terms (self.nN) and the roots (self.nR)
		# for the IC calculation
		for ele in ancestorList:
			go,ontology,evidence=ele
			if dataGo.goData[go]["isRoot"]:
				self.nR[ontology]+=1
			
			if go not in self.nN[ontology]:
				self.nN[ontology][go]=0
			self.nN[ontology][go]+=1


	def _updateAlternative(self,goSet,dataGo):
		# Check if it is an alt_id and replace it
		goSetNew=set()
		for ele in goSet:
			go,ontology,evidence=ele
			if go in dataGo.alternativeId:
				ele=dataGo.alternativeId[go],ontology,evidence
			goSetNew.add(ele)
		return goSetNew
	
	def _updateObsolete(self,goSet,dataGo):
		goSetNew=set()
		for ele in goSet:
			go,ontology,evidence=ele
			if go not in dataGo.goData:
				if go in dataGo.obsolete: 
					for goNew in dataGo.obsolete[go]["consider"].union(dataGo.obsolete[go]["replaced"]):
						if goNew in dataGo.goData:
							ele=goNew,ontology,evidence	
							goSetNew.add(ele)
				else:
					if go not in dataGo.alternativeId:
						self.goMissing.add(go)
					# DEBUG: Add go term not mapped???
					#goSetNew.add(ele)
			else:
				goSetNew.add(ele)
		return goSetNew
					


	def _calculate_IC(self,dataGo):
		"""
		Calculate the IC content for each GO term.
		The root terms are removed from the final dictionary.
		"""

		# Calculate IC
		for ontology in self.nN:
			for go in self.nN[ontology]:
				self.goIC[ontology][go]=-log(float(self.nN[ontology][go])/self.nR[ontology])
							
		# Set root IC to 0
		for ontology in self.goIC:
			for go in self.goIC[ontology]:
				if dataGo.goData[go]["isRoot"]:
					self.goIC[ontology][go]=0.0

			
	


	def writeIC(self,fileIC):
		"""
		Write the IC CSV file (with first line header starting with "##")
		"""
		fileIC.write('## Source:'+self.fileSource+'\tEvidenceCodes:'+','.join(list(self.evidenceCodeList))+'\tOntology:'+','.join(list(self.ontologyList))+'\n')
		for ontology in self.goIC:
			for go in self.goIC[ontology]:
				fileIC.write(go+','+ontology+','+str(self.goIC[ontology][go])+'\n')
		fileIC.close()


	def statIC(self):
		stat={}
		stat['Go terms']=self.ontologyDist
		
		result='#'*60+'\nAnnotation DB statistics:\n'+'{0:>21} = {1}'.format('Accessions',self.nAccession)+'\n'+'\n'.join(['{0:>21} = {1}'.format(k,' '.join(['{0:>6} {1}'.format(str(stat[k][ont]),ont) for ont in stat[k]])) for k in stat])+'\n'
		keys=[m[1] for m in sorted([(self.evidenceCodeDist[k],k)  for k in self.evidenceCodeDist], reverse=True)]
		result+='{0:>21} = {1} {2}'.format("Terms x Evidence Code",keys[0],self.evidenceCodeDist[keys[0]])+'\n'

		
		result+='\n'.join(['{0:>27} {1}'.format(k,self.evidenceCodeDist[k]) for k in keys[1:] ])

		return result


###################################################################################

def readIC(fileIC):
	"""
	Read the Information Content CSV file
	"""
	header={'fileSource':None,'evidenceCodeList':None,'ontologyList':None}	
	goIC={}		
	for line in fileIC.readlines():
		if line[0:2]=='##':
			line=line.strip().split('\t')
			header['fileSource']=line[0].split(':')[1].strip()
			header['evidenceCodeList']=set(line[1].split(':')[1].strip().split(','))
			header['ontologyList']=set(line[2].split(':')[1].strip().split(','))
		else:
			go,ontology,ic=line.strip().split(',')
			goIC[go]=(ontology,float(ic))
	fileIC.close()
	return header,goIC


def statIc(goIC,header):
	stat={}
	stat['Terms']=str(len(goIC))
	stat['Ontologies']=','.join(list(header['ontologyList']))	
	stat['Evidence Code']=','.join(list(header['evidenceCodeList']))	
	
	result='#'*60+'\nInformation Content file:\n'+'\n'.join(['{0:>20} = {1}'.format(k,stat[k]) for k in sorted(stat.keys())])+'\n'

	return result


###################################################################################
class calcScores:
	"""
	Calculate scores.
	"""

	def __init__(self):
		self.referenceSet={}	# { accession : { ontology : { go : score } } } 
		self.predictedSet={}	# { accession : { ontology : { go : score } } }

		self.cogic={}			# { target : { ontology : { go : ( cogic_score, data ) } } }

		self.predictedRef={'P':0,'C':0,'F':0}				# Reference target predicted

		self.missingGoIc={'pred':{'P':set(),'C':set(),'F':set()},'ref':{'P':set(),'C':set(),'F':set()}}	# Missing IC value of reference GO terms
		self.targetMissing={'P':set(),'C':set(),'F':set()}	# Sets of not predicted target
		self.targetNotPredictable={'P':set(),'C':set(),'F':set()}	
		
		self.totalGo={'pred':{'P':set(),'C':set(),'F':set()},'ref':{'P':set(),'C':set(),'F':set()}}
		self.goMapping={'alternative':{'P':set(),'C':set(),'F':set()},'obsolete':{'P':set(),'C':set(),'F':set()},'missing':set()}

		self.precRecall={}
		
		
	def parseAnnotation(self,inputFile,separator=None,posAcc=0,posGo=1,\
										posScore=None,posLabel=None,labelSet=None,targetList=None):
		"""
		Parse an annotation file (prediction/reference). The considered
		target may be filtered passing a targetList.
		"""
		data={}
		line=inputFile.readline()
		while line:
			if line[0]!='#':
				line=line.split(separator)
				if len(line)<2:
					print '\nERROR: Can\'t split the input lines. Possible causes are: a bad separator, a bad file format or a program problem parsing the separator string\n'
					sys.exit(1)
					
				accession=line[posAcc]
				go=line[posGo]
				
				flag=False
				if targetList:
					if accession in targetList:
						flag=True
				else:
					flag=True
					
				if flag:
					if accession not in data:
						data[accession]={}
					data[accession][go]={'score':None,'label':None}
					
					if posScore:
						if len(line)>=posScore+1:
							data[accession][go]['score']=float(line[posScore])
						else:
							print '\nERROR: Score position too big\n'
							print "ERROR LINE:\n",line,posScore
							sys.exit(1)
					else:
						data[accession][go]['score']=1.0
					
					if posLabel:
						if len(line)<posLabel+1:
							print '\nERROR: Label position too big\n'
							print "ERROR LINE:\n",line,posLabel
							sys.exit(1)
						else:
							data[accession][go]['label']=line[posLabel]
				
				
			line=inputFile.readline()
		inputFile.close()
		
		if labelSet:
			for accession in data.keys():
				for go in data[accession].keys():
					if data[accession][go]['label'] not in labelSet:
						del data[accession][go]
				if not data[accession]:
					del data[accession]
		
		return data




	def _updateAlternative(self,dataGo,goSet):
		# Check if it is an alt_id and replace it
		for go in goSet.keys():
			if go in dataGo.alternativeId: 
				goNew=dataGo.alternativeId[go]
				goSet[goNew]=copy.deepcopy(goSet[go])
				self.goMapping['alternative'][dataGo.goData[goNew]['ontology']].add((go,goNew))
				del goSet[go]
				
	def _updateObsolete(self,dataGo,goSet):
		for go in goSet.keys():
			if go not in dataGo.goData:
				if go in dataGo.obsolete: 
					
					#for goNew in dataGo.obsolete[go]["consider"].union(dataGo.obsolete[go]["replaced"]):
					#	if goNew in dataGo.goData:
					#		goSet[goNew]=copy.deepcopy(goSet[go])
					#		self.goMapping['obsolete'][dataGo.goData[goNew]['ontology']].add((go,goNew))

					goNewList=list(dataGo.obsolete[go]["replaced"])
					if len(goNewList)==1:
						goNew=goNewList[0]
						if goNew in dataGo.goData:
							goSet[goNew]=copy.deepcopy(goSet[go])
							self.goMapping['obsolete'][dataGo.goData[goNew]['ontology']].add((go,goNew))

					del goSet[go]
					

	def _filterByDistance(self,dataGo,goSet,depthThreshold):
		for go in goSet.keys():
			depth=dataGo.goData[go]["depth"]
			isLeaf=dataGo.goData[go]["isLeaf"]
			if depthThreshold=='leaves':
				if not isLeaf:
					del goSet[go]
			elif depth<int(depthThreshold):
				del goSet[go]
		
			

	def expandAncestorsList(self,dataGo,targetsGo,removeRoot=False,depthThreshold=None,goFilter=None,expandAncestors=True):
		"""
		Expand the ancestors of each GO term associated to a target. If the
		obsoleteFlag is true it remaps each obsolete term with the new ones 
		listed in the "replaced_by" and "consider" fields of the OBO file.
 
		The ancestor scores are inherited from the predicted one, if a target 
		is endowed with two different go terms with a common ancestor it is 
		taken only once and its score is inherited from the children highest score

		If depthThreshold is set all terms are filtered based on root distance 
		and only kept terms are expanded.
		"""	
		
		targetGoExpanded={}
		for accession in targetsGo:
			goSet=targetsGo[accession]
			
			# Update alternative id
			if goSet:
				self._updateAlternative(dataGo,goSet)
			
			# Update obsolete id
			if goSet:
				self._updateObsolete(dataGo,goSet)
			
			# Filter based on depth
			if goSet and depthThreshold:
				self._filterByDistance(dataGo,goSet,depthThreshold)

			# Remove missing terms
			if goSet:
				for go in goSet.keys():
					if go not in dataGo.goData:
						self.goMapping['missing'].add(go)
						del goSet[go]
					

			if goSet:
				ancestorList=[]
				for go in goSet:

					data=targetsGo[accession][go]
					ontology=dataGo.goData[go]["ontology"]
										
					ancestorList.append((ontology,go,data))
					
					if expandAncestors:
						ancestors=dataGo.goData[go]["ancestors"]
						ancestorList.extend([(ontology,goAncestor,data) for goAncestor in ancestors])
	
				# filter by goFilter and removeRoot
				if ancestorList:
					
					if goFilter:
						ancestorList=[(ont,go,data) for ont,go,data in ancestorList if go in goFilter]
					if removeRoot:
						ancestorList=[(ont,go,data) for ont,go,data in ancestorList if not dataGo.goData[go]['isRoot']]
							
				# Filter common ancestors retaining the best score inherited
				if ancestorList:
					ancestorDict={}
					for ont,go,data in ancestorList:
						if ont not in ancestorDict:
							ancestorDict[ont]={}
						if go not in ancestorDict[ont]:
							ancestorDict[ont][go]=data
						else:
							if data['score']>ancestorDict[ont][go]['score']:
								ancestorDict[ont][go]=data
					if ancestorDict:
						targetGoExpanded[accession]=ancestorDict
			
		return targetGoExpanded
	
	
	def _calculateSimGIC(self,goIC,goSetA,ontology,target): 
		"""
		Calculate the simGIC score
		- goSetA = prediction
		- goSetB = referenece
		"""
		goSetB=set(self.referenceSet[target][ontology].keys())
		
		sumIntersectionIC=0.0
		sumUnionIC=0.0
		
		intersection=goSetA.intersection(goSetB)
		union=goSetA.union(goSetB)

		for go in intersection:
			if go in goIC:
				sumIntersectionIC+=goIC[go][1]
		for go in union:
			if go in goIC:
				sumUnionIC+=goIC[go][1]

		# union != 0.0 filter references with only root terms
		if sumUnionIC!=0.0:
			# filter when the only common term is the root
			if len(intersection)>1: 
				simGIC=sumIntersectionIC/sumUnionIC
				return {'simGIC':simGIC,'lenI':len(intersection),'lenU':len(union),'icI':round(sumIntersectionIC,3),'icU':round(sumUnionIC,3)}
	
	def _calculateSimGIC_alternative(self,goIC,dataGo,goSetA,ontology,target): 
		"""
		Calculate the simGIC score
		- goSetA = prediction
		- goSetB = referenece
		"""
		goSetB=set(self.referenceSet[target][ontology].keys())
		
		
		
		"""
		C_Set=set()
											
		# Check if the predicted go is among the reference ancestors
		for go in self.referenceSet[accession][ontology]:
			for ancestor in dataGo.goData[go]['ancestors']:
				if ancestor in bins[th]:
					C_Set.add(go)
					
		# Check if the reference go are among the predicted ancestors
		for go in bins[th]:
			for ancestor in dataGo.goData[go]['ancestors']:
				if ancestor in self.referenceSet[accession][ontology]:
					C_Set.add(ancestor)
		"""
		
		
		# Calculate intersection in an alternative way
		# A term is correct if it is an ancestor of the reference terms
		# or it is a children of at least one reference term
		
		union=goSetA.union(goSetB)
		intersection=set()
		
		for go in goSetB:
			for ancestor in dataGo.goData[go]['ancestors']:
				if ancestor in goSetA:
					intersection.add(go)
				if ancestor in union:
					union.remove(ancestor)	
					
		for go in goSetA:
			for ancestor in dataGo.goData[go]['ancestors']:
				if ancestor in goSetB:
					intersection.add(ancestor)
					
		intersection.update(goSetB.intersection(goSetA))
		union.update(intersection)
		
		

		sumIntersectionIC=0.0
		sumUnionIC=0.0
		
		for go in intersection:
			if go in goIC:
				sumIntersectionIC+=goIC[go][1]
		for go in union:
			if go in goIC:
				sumUnionIC+=goIC[go][1]

		# union != 0.0 filter references with only root terms
		if sumUnionIC!=0.0:
			# filter when the only common term is the root
			if len(intersection)>1: 
				simGIC=sumIntersectionIC/sumUnionIC
				return {'simGIC':simGIC,'lenI':len(intersection),'lenU':len(union),'icI':round(sumIntersectionIC,3),'icU':round(sumUnionIC,3)}
		
			
			
	def _getConfidenceSubset(self,dataGo,ontology,target,confidenceThreshold):
		"""
		Extract a subset of go terms based on a score threshold
		"""
		subset={}
		targetPredictedAnnotation=self.predictedSet[target][ontology]
		for go in targetPredictedAnnotation:
			if targetPredictedAnnotation[go]['score']>=confidenceThreshold:
				subset[go]=targetPredictedAnnotation[go]
				for ancestor in dataGo.goData[go]['ancestors']:
					subset[ancestor]=copy.deepcopy(subset[go])		
		
		return subset
	
			
	
	
	def calculateCOGIC(self,goIC,dataGo,filterMissing=True,outputDetail=False,alternativeMethod=False,goPredictable=None):
		"""
		Calculate the COGIC score
		"""
		roots=set([go for go in dataGo.goData if dataGo.goData[go]['isRoot']])
		
		thresholdList=[0.75,0.5,0.25,0.0]
		if not self.referenceSet or not self.predictedSet:
			print '\nERROR: Not annotation found!\n'
			sys.exit(1)

		for target in self.referenceSet:
			if target not in self.cogic:
				self.cogic[target]={}
			for ontology in self.referenceSet[target]:
				self.totalGo['ref'][ontology].update(set(self.referenceSet[target][ontology].keys()))
				if target in self.predictedSet:
					if ontology in self.predictedSet[target]:
						self.predictedRef[ontology]+=1
						self.totalGo['pred'][ontology].update(set(self.predictedSet[target][ontology].keys()))
						
						simGicList=None
						for threshold in thresholdList:
							predictedSubset=self._getConfidenceSubset(dataGo,ontology,target,threshold)

							if alternativeMethod:
								simGIC=self._calculateSimGIC_alternative(goIC,dataGo,set(predictedSubset.keys()),ontology,target)
							else:
								simGIC=self._calculateSimGIC(goIC,set(predictedSubset.keys()),ontology,target)
							# If not, missing IC for all reference GO terms 
							if simGIC:
								if not simGicList:
									simGicList=[]
								simGIC['Th']=threshold
								simGicList.append(simGIC)

						if simGicList:
							# Calculate the COGIC score
							
							cogicScore=sum([i['simGIC'] for i in simGicList])/4
							
							# Store data for visualization
							#self.cogic[target][ontology]=cogicScore,'|'.join(['Sc:'+str(round(ele['simGIC'],3))+','+','.join([k+':'+str(ele[k]) for k in ele]) for ele in simGicList])
							self.cogic[target][ontology]=cogicScore

						else:
							# All go terms (reference + predicted) do not have the IC value except the root term
							if filterMissing:
								self.cogic[target][ontology]=None
								#self.cogic[target][ontology]=0.0
							else:
								self.cogic[target][ontology]=0.0
							# Check if all reference terms are not predictable 	
							if goPredictable:
								if not (set(self.referenceSet[target][ontology].keys())-roots).intersection(goPredictable):
									self.targetNotPredictable[ontology].add(target)
								
					else:
						# Not predicted terms for that target and that ontology
						#self.cogic[target][ontology]=None
								
						# Missing target in predicted set 
						self.targetMissing[ontology].add(target)
				else:
					# Target not predicted
					#self.cogic[target][ontology]=None
					
					# Missing target in predicted set 
					self.targetMissing[ontology].add(target)
							
		for ontology in self.totalGo['ref']:			
			self.missingGoIc['ref'][ontology].update(self.totalGo['ref'][ontology]-set(goIC.keys()))
			self.missingGoIc['pred'][ontology].update(self.totalGo['pred'][ontology]-set(goIC.keys()))
					

					
	def statScore(self):
		stat={}
		stat['1 - Reference targets read']={'P':0,'C':0,'F':0}
		stat['2 - Predicted targets read']={'P':0,'C':0,'F':0}
		stat['3 - Predicted reference targets']=self.predictedRef
		stat['4 - Cogic zero score']={'P':0,'C':0,'F':0}
		stat['5 - Cogic N.A.']={'P':0,'C':0,'F':0}
		stat['6 - Not predictable']={'P':0,'C':0,'F':0}
		
		# Targets in the reference
		for target in self.referenceSet:
			for ontology in self.referenceSet[target]:
				stat['1 - Reference targets read'][ontology]+=1
		
		# Number of target after uploading alternative/obsolete
		for target in self.predictedSet:
			for ontology in self.predictedSet[target]:
				stat['2 - Predicted targets read'][ontology]+=1
		
		# Evaluated targets
		for target in self.cogic:
			for ontology in self.cogic[target]:
				if self.cogic[target][ontology]==0.0:
					stat['4 - Cogic zero score'][ontology]+=1
				elif self.cogic[target][ontology]==None:
					stat['5 - Cogic N.A.'][ontology]+=1
		
		# Target annotated with not predictable terms
		for ontology in self.targetNotPredictable:
			for target in self.targetNotPredictable[ontology]:
				stat['6 - Not predictable'][ontology]+=1
				

		for s in stat.keys():
			stat[s]=','.join(['{0:7} {1}'.format(stat[s][k],k) for k in stat[s]])
		
		result='#'*60+'\nScore statistics TARGETS:\n'+'\n'.join(['{0[0]}{0[1]:>32} = {1}'.format(k.split(' - '),stat[k]) for k in sorted(stat.keys())])+'\n'
		
		
		# Terms statistics
		
		# Terms missing the IC value
		stat={}
		stat['1 - Terms in reference set']=dict([(ontology,len(self.totalGo['ref'][ontology])) for ontology in self.totalGo['ref']])
		stat['2 - Terms in predicted set']=dict([(ontology,len(self.totalGo['pred'][ontology])) for ontology in self.totalGo['pred']])
		stat['3 - Missing reference GO IC']=dict([(ont,len(self.missingGoIc['ref'][ont])) for ont in self.missingGoIc['ref']])
		stat['3 - Missing predicted GO IC']=dict([(ont,len(self.missingGoIc['pred'][ont])) for ont in self.missingGoIc['pred']])
		stat['4 - Obsolete mapped terms']=dict([(ont,len(set([ele[0] for ele in self.goMapping['obsolete'][ont]]))) for ont in self.goMapping['obsolete']])
		
		for s in stat.keys():
			stat[s]=','.join(['{0:7} {1}'.format(stat[s][k],k) for k in stat[s]])
		stat['5 - Not mapped GO terms']=len(self.goMapping['missing'])
		result+='\nScore statistics GO TERMS:\n'+'\n'.join(['{0[0]}{0[1]:>26} = {1}'.format(k.split(' - '),stat[k]) for k in sorted(stat.keys())])+'\n'
		
		# GO without the IC value
		statStr=""	
		for ontology in self.missingGoIc['ref']:
			if self.missingGoIc['ref'][ontology]:
				statStr+=ontology+': '+'\n'.join([' '.join(list(self.missingGoIc['ref'][ontology])[i:i+10]) for i in range(0,len(list(self.missingGoIc['ref'][ontology])),10)])+'\n'
		if statStr:
			result+='\n\n\n'+'#'*60+'\nMissing IC value for the following GO terms in the reference set:\n'+statStr
		statStr=""	
		for ontology in self.missingGoIc['pred']:
			if self.missingGoIc['pred'][ontology]:
				statStr+=ontology+': '+'\n'.join([' '.join(list(self.missingGoIc['pred'][ontology])[i:i+10]) for i in range(0,len(list(self.missingGoIc['pred'][ontology])),10)])+'\n'
		if statStr:
			result+='\n'+'#'*60+'\nMissing IC value for the following GO terms in the predicted set:\n'+statStr

		
		# Targets not predicted cogic=0.0
		statStr=""	
		for ontology in self.targetMissing:
			if self.targetMissing[ontology]:
				statStr+=ontology+': '+'\n'.join([' '.join(list(self.targetMissing[ontology])[i:i+10]) for i in range(0,len(list(self.targetMissing[ontology])),10)])+'\n'
		if statStr:
			result+='\n'+'#'*60+'\nTarget not evaluated (not predicted):\n'+statStr


		# Obsolete terms mapping
		statStr=""	
		for ontology in self.goMapping['obsolete']:
			if self.goMapping['obsolete'][ontology]:
				statStr+=ontology+':\n'+'\n'.join([' '.join(list(ele)) for ele in sorted(list(self.goMapping['obsolete'][ontology]))])+'\n'
		if statStr:
			result+='\n'+'\n#'*60+'\nObsolete mapping ( OLD-ID / NEW-ID ):\n'+statStr

		"""
		# Not mapped terms
		if self.goMapping['missing']:
			result+='#'*60+'\nNot mapped GO terms (obsolete or missing in the OBO file):\n'+'\n'.join([' '.join(list(self.goMapping['missing'])[i:i+10]) for i in range(0,len(list(self.goMapping['missing'])),10)])+'\n'
		"""
		
		return result


	def calculatePrecisionRecallThresholdBins(self,dataGo,filterMissing=False,realCoverage=True):
		
		score={'F':{},'P':{},'C':{}}
		nTerm={'C':0,'P':0,'F':0}	# Number of reference terms
		nR={'C':0,'P':0,'F':0}		# Number of reference targets
		nP={'C':0,'P':0,'F':0}		# Number of predicted targets

		roots=set([go for go in dataGo.goData if dataGo.goData[go]['isRoot']])
		
		# Initialize score
		for ontology in score:
			for threshold in range(1,1001,1):
				score[ontology][threshold*0.001]={'C':0,'P':0,'T':0,'TP':0,'FP':0,'FN':0,'pr':0.0,'rc':0.0,'f':0.0,'cov':0.0,'Ntarget':0}
				
		# Count prediction matches for each target
		for accession in self.referenceSet:
			for ontology in self.referenceSet[accession]:
				refSet=set(self.referenceSet[accession][ontology].keys())-roots
				if refSet:
					# To calculate Rc including not predicted target
					if realCoverage:
						nR[ontology]+=1
						nTerm[ontology]+=len(self.referenceSet[accession][ontology])

					if accession in self.predictedSet:
						if ontology in self.predictedSet[accession]:
							predSet=set(self.predictedSet[accession][ontology].keys())-roots
							if predSet:
								
								flag=True
								if filterMissing:
									if refSet.intersection(predSet):
										flag=True
									else:
										flag=False
								if flag:
									# To calculate Rc excluding not predicted target
									if not realCoverage:
										nR[ontology]+=1
										nTerm[ontology]+=len(self.referenceSet[accession][ontology])
									
									nP[ontology]+=1
									# Calculate Prediction(P), CorrectPrediction(C) for each score threshold
									bins={}
									for th in score[ontology]:
										bins[th]=set()
										
									for go in predSet:
										scoreVal=self.predictedSet[accession][ontology][go]['score']
										for th in bins:
											if float(scoreVal)>=th:
												bins[th].add(go)
												
									for th in bins:
										if bins[th]:
											P=len(bins[th])
											T=len(refSet)
											C=len(bins[th].intersection(refSet))
											score[ontology][th]['C']+=C
											score[ontology][th]['P']+=P
											score[ontology][th]['T']+=T
											score[ontology][th]['Ntarget']+=1
											score[ontology][th]['TP']+=C
											score[ontology][th]['FP']+=P-C
											score[ontology][th]['FN']+=T-C
										
									
		# Calculate precision,recall,F-measure,coverage
		for ontology in score:
			for th in sorted(score[ontology].keys()):
				if nTerm[ontology]!=0: 
					if score[ontology][th]['P']!=0:
						pr=float(score[ontology][th]['C'])/score[ontology][th]['P']
					else:
						pr=0.0
					rc=float(score[ontology][th]['C'])/nTerm[ontology]
					f=2.0*score[ontology][th]['C']/float(score[ontology][th]['P']+nTerm[ontology])
					coverage=float(nP[ontology])/nR[ontology]					
	
					score[ontology][th]['pr']=pr
					score[ontology][th]['rc']=rc
					score[ontology][th]['f']=f
					score[ontology][th]['cov']=coverage
					
		self.precRecall=score


	def calculatePrecisionRecallThresholdBins_ALTERNATIVE(self,dataGo,realCoverage=True):
		
		score={'F':{},'P':{},'C':{}}
		nTerm={'C':0,'P':0,'F':0}	# Number of reference terms
		nR={'C':0,'P':0,'F':0}		# Number of reference targets
		nP={'C':0,'P':0,'F':0}		# Number of predicted targets
		
		roots=set([go for go in dataGo.goData if dataGo.goData[go]['isRoot']])

	
			
		# Count prediction matches for each target
		for accession in self.referenceSet:
			for ontology in self.referenceSet[accession]:
				# remove reference with only the root terms
				refSet=set(self.referenceSet[accession][ontology].keys())-roots
				if refSet:
					
					if realCoverage:
						nR[ontology]+=1
						
					if accession in self.predictedSet:
						if ontology in self.predictedSet[accession]:
							# To exclude terms with only the root terms predicted
							predSet=set(self.predictedSet[accession][ontology].keys())-roots
							if predSet:
								
								if not realCoverage:
									nR[ontology]+=1
									
								nP[ontology]+=1	
								
								# Calculate Prediction(P), CorrectPrediction(C) for each score threshold
								bins={}
								
								for go in predSet:
									scoreVal=self.predictedSet[accession][ontology][go]['score']
									if scoreVal not in bins:
										bins[scoreVal]=set()
									bins[scoreVal].add(go)	
									
							
											
								for th in bins:
									if th not in score[ontology]:
										score[ontology][th]={'C':0,'P':0,'T':0,'TP':0,'FP':0,'FN':0,'pr':0.0,'rc':0.0,'f':0.0,'cov':0.0,'Ntarget':0}
		
									# Calculate TP in an alternative way
									# A term is correct if it is an ancestor of the reference terms
									# or it is a children of at least one reference term
									C_Set=set()
									P_Set=copy.deepcopy(bins[th])
									# Check if the predicted go is among the reference ancestors
									for go in refSet:
										for ancestor in dataGo.goData[go]['ancestors']:
											if ancestor in bins[th]:
												C_Set.add(go)
											if ancestor in P_Set:
												P_Set.remove(ancestor)
												
									# Check if the reference go are among the predicted terms ancestors
									for go in bins[th]:
										for ancestor in dataGo.goData[go]['ancestors']:
											if ancestor in self.referenceSet[accession][ontology]:
												C_Set.add(ancestor)
									P_Set.update(C_Set)
									
									
									
									C=len(C_Set)
									P=len(P_Set)
									T=len(refSet)
									
									
									#if C!=0:
									
									score[ontology][th]['pr']+=float(C)/P
									score[ontology][th]['rc']+=float(C)/T
									score[ontology][th]['Ntarget']+=1
								
									score[ontology][th]['TP']+=C
									score[ontology][th]['FP']+=P-C
									score[ontology][th]['FN']+=T-C
								
									score[ontology][th]['C']+=C
									score[ontology][th]['P']+=P
									score[ontology][th]['T']+=T


		# Porpagating for lower thresholds
		for ontology in score:
			thList=sorted(score[ontology].keys(), reverse=True)
			
			for i in range(1,len(thList)):
				th=thList[i]
				for ele in score[ontology][th].keys():
					score[ontology][th][ele]+=score[ontology][thList[i-1]][ele]
				

	
		# Calculate precision,recall,F-measure,coverage
		for ontology in score:
			for th in sorted(score[ontology].keys()):
				if score[ontology][th]['pr']!=0.0 or score[ontology][th]['rc']!=0.0: 
					score[ontology][th]['pr']=score[ontology][th]['pr']/score[ontology][th]['Ntarget']
					score[ontology][th]['rc']=score[ontology][th]['rc']/nR[ontology]
					
					#print th,score[ontology][th]['pr'],score[ontology][th]['rc'],score[ontology][th]['Ntarget']
					Fmeasure=2.0*score[ontology][th]['pr']*score[ontology][th]['rc']/(score[ontology][th]['pr']+score[ontology][th]['rc'])
					if Fmeasure>score[ontology][th]['f']:
						score[ontology][th]['f']=Fmeasure
					coverage=float(score[ontology][th]['Ntarget'])/nR[ontology]
					if coverage>score[ontology][th]['cov']:
						score[ontology][th]['cov']=coverage
					
		self.precRecall=score

	def calculatePrecisionRecallThresholdBins_CAFA(self,dataGo,filterMissing=True):
		
		score={'F':{},'P':{},'C':{}}
		nTerm={'C':0,'P':0,'F':0}	# Number of reference terms
		nR={'C':0,'P':0,'F':0}		# Number of reference targets
		nP={'C':0,'P':0,'F':0}		# Number of predicted targets

		
		roots=set([go for go in dataGo.goData if dataGo.goData[go]['isRoot']])

		# Initialize score
		for ontology in score:
			for threshold in range(1,1001,1):
				score[ontology][threshold*0.001]={'C':0,'P':0,'T':0,'TP':0,'FP':0,'FN':0,'pr':0.0,'rc':0.0,'f':0.0,'cov':0.0,'Ntarget':0}
				
		# Count prediction matches for each target
		for accession in self.referenceSet:
			for ontology in self.referenceSet[accession]:
				refSet=set(self.referenceSet[accession][ontology].keys())-roots
				# remove reference with only the root terms
				if refSet:
					# To calculate Rc including not predicted target
					nR[ontology]+=1
					if accession in self.predictedSet:
						if ontology in self.predictedSet[accession]:
							predSet=set(self.predictedSet[accession][ontology].keys())-roots
							# To exclude terms with only the root terms predicted
							if predSet:
								flag=True
								if filterMissing:
									if predSet.intersection(refSet):
										flag=True
									else:
										flag=False
								if flag:
									nP[ontology]+=1	
								
									# Calculate Prediction(P), CorrectPrediction(C) for each score threshold
									bins={}
									for th in score[ontology]:
										bins[th]=set()
									for go in predSet:
										scoreVal=self.predictedSet[accession][ontology][go]['score']
										
										for th in bins:
											if float(scoreVal)>=th:
												bins[th].add(go)
												
									for th in bins:
										if bins[th]:
											
											P=len(bins[th])
											T=len(refSet)
											C=len(bins[th].intersection(refSet))
											
											score[ontology][th]['Ntarget']+=1
											score[ontology][th]['pr']+=float(C)/P
											score[ontology][th]['rc']+=float(C)/T
											
											score[ontology][th]['TP']+=C
											score[ontology][th]['FP']+=P-C
											score[ontology][th]['FN']+=T-C
											

							
								
		# Calculate precision,recall,F-measure,coverage
		for ontology in score:
			for th in sorted(score[ontology].keys()):
				if score[ontology][th]['pr']!=0.0 or score[ontology][th]['rc']!=0.0: 
				
					score[ontology][th]['pr']=score[ontology][th]['pr']/score[ontology][th]['Ntarget']
					score[ontology][th]['rc']=score[ontology][th]['rc']/nR[ontology]
					
					Fmeasure=2.0*score[ontology][th]['pr']*score[ontology][th]['rc']/(score[ontology][th]['pr']+score[ontology][th]['rc'])
					if Fmeasure>score[ontology][th]['f']:
						score[ontology][th]['f']=Fmeasure
						
					coverage=float(score[ontology][th]['Ntarget'])/nR[ontology]
					if coverage>score[ontology][th]['cov']:
						score[ontology][th]['cov']=coverage
					
		self.precRecall=score

	def calculatePrecisionRecallThresholdBins_CAFA_new(self,dataGo,filterMissing=True):
		
		score={'F':{},'P':{},'C':{}}
		nTerm={'C':0,'P':0,'F':0}	# Number of reference terms
		nR={'C':0,'P':0,'F':0}		# Number of reference targets
		nP={'C':0,'P':0,'F':0}		# Number of predicted targets

		
		roots=set([go for go in dataGo.goData if dataGo.goData[go]['isRoot']])
		"""
		# Initialize score
		for ontology in score:
			for threshold in range(1,1001,1):
				score[ontology][threshold*0.001]={'C':0,'P':0,'T':0,'TP':0,'FP':0,'FN':0,'pr':0.0,'rc':0.0,'f':0.0,'cov':0.0,'Ntarget':0}
		"""		
		# Count prediction matches for each target
		for accession in self.referenceSet:
			for ontology in self.referenceSet[accession]:
				refSet=set(self.referenceSet[accession][ontology].keys())-roots
				# remove reference with only the root terms
				if refSet:
					# To calculate Rc including not predicted target
					nR[ontology]+=1
					if accession in self.predictedSet:
						if ontology in self.predictedSet[accession]:
							predSet=set(self.predictedSet[accession][ontology].keys())-roots
							# To exclude terms with only the root terms predicted
							if predSet:
								flag=True
								if filterMissing:
									if predSet.intersection(refSet):
										flag=True
									else:
										flag=False
								if flag:
									nP[ontology]+=1	
								
									# Calculate Prediction(P), CorrectPrediction(C) for each score threshold
									bins={}
									
									for go in predSet:
										scoreVal=self.predictedSet[accession][ontology][go]['score']
										if scoreVal not in bins:
											bins[scoreVal]=set()
										bins[scoreVal].add(go)
												
									for th in bins:
										if th not in score[ontology]:
											score[ontology][th]={'C':0,'P':0,'T':0,'TP':0,'FP':0,'FN':0,'pr':0.0,'rc':0.0,'f':0.0,'cov':0.0,'Ntarget':0}
		
										P=len(bins[th])
										T=len(refSet)
										C=len(bins[th].intersection(refSet))
										
										score[ontology][th]['Ntarget']+=1
										score[ontology][th]['pr']+=float(C)/P
										score[ontology][th]['rc']+=float(C)/T
										
										score[ontology][th]['TP']+=C
										score[ontology][th]['FP']+=P-C
										score[ontology][th]['FN']+=T-C
											

		# Porpagating for lower thresholds
		for ontology in score:
			thList=sorted(score[ontology].keys(), reverse=True)
			
			for i in range(1,len(thList)):
				th=thList[i]
				for ele in score[ontology][th].keys():
					score[ontology][th][ele]+=score[ontology][thList[i-1]][ele]
				
				
								
		# Calculate precision,recall,F-measure,coverage
		for ontology in score:
			for th in sorted(score[ontology].keys()):
				if score[ontology][th]['pr']!=0.0 or score[ontology][th]['rc']!=0.0: 
				
					score[ontology][th]['pr']=score[ontology][th]['pr']/score[ontology][th]['Ntarget']
					score[ontology][th]['rc']=score[ontology][th]['rc']/nR[ontology]
					
					Fmeasure=2.0*score[ontology][th]['pr']*score[ontology][th]['rc']/(score[ontology][th]['pr']+score[ontology][th]['rc'])
					if Fmeasure>score[ontology][th]['f']:
						score[ontology][th]['f']=Fmeasure
						
					coverage=float(score[ontology][th]['Ntarget'])/nR[ontology]
					if coverage>score[ontology][th]['cov']:
						score[ontology][th]['cov']=coverage
					
		self.precRecall=score
		
		
	def calculatePrecisionRecallThresholdBins_SIFTER(self,dataGo):
		
		score={'F':{},'P':{},'C':{}}
		
		rT={'C':set(),'P':set(),'F':set()}	# Reference Terms
		pT={'C':{},'P':{},'F':{}}

		roots=set([go for go in dataGo.goData if dataGo.goData[go]['isRoot']])

		# Initialize score
		for ontology in score:
			for threshold in range(1,1001,1):
				score[ontology][threshold*0.001]={'C':0,'P':0,'T':0,'TP':0,'FP':0,'FN':0,'pr':0.0,'rc':0.0,'f':0.0,'cov':0.0,'Ntarget':0}
				pT[ontology][threshold*0.001]=set()
				
		# Count prediction matches for each target
		for accession in self.referenceSet:
			for ontology in self.referenceSet[accession]:
				# remove reference with only the root terms
				if not set(self.referenceSet[accession][ontology].keys()).issubset(roots):
					# To calculate Rc including not predicted target
					#nR[ontology]+=1
					#nTerm[ontology]+=len(self.referenceSet[accession][ontology])
					if accession in self.predictedSet:
						if ontology in self.predictedSet[accession]:
							# To exclude terms with only the root terms predicted
							if not set(self.predictedSet[accession][ontology].keys()).issubset(roots):
								
								# Reference terms
								rT[ontology].update(self.referenceSet[accession][ontology])
								
								# Predicted terms
								for go in self.predictedSet[accession][ontology]:
									scoreVal=self.predictedSet[accession][ontology][go]['score']
									for th in pT[ontology]:
										if float(scoreVal)>=th:
											pT[ontology][th].add(go)
											
								
							
								
		# Calculate precision,recall,F-measure,coverage
		for ontology in pT:
			for th in sorted(pT[ontology].keys()):
				if pT[ontology][th]: 
					C=len(rT[ontology].intersection(pT[ontology][th]))
					P=len(pT[ontology][th])
					T=len(rT[ontology])
					score[ontology][th]['pr']=float(C)/P
					score[ontology][th]['rc']=float(C)/T
					score[ontology][th]['TP']=C
					score[ontology][th]['FP']=P-C
					score[ontology][th]['FN']=T-C
					
					Fmeasure=2.0*C/(P+T)
					if Fmeasure>score[ontology][th]['f']:
						score[ontology][th]['f']=Fmeasure
					
					
		self.precRecall=score

	def calculatePrecisionRecall_DEBUG(self):
		
		score={'F':{},'P':{},'C':{}}
		nTerm={'C':0,'P':0,'F':0}	# Number of reference terms
		nR={'C':0,'P':0,'F':0}		# Number of reference targets
		nP={'C':0,'P':0,'F':0}		# Number of predicted targets

					
		# Count prediction matches for each target
		for accession in self.referenceSet:
			for ontology in self.referenceSet[accession]:				
				# To calculate Rc including not predicted target
				#nR[ontology]+=1
				#nTerm[ontology]+=len(self.referenceSet[accession][ontology])
				if accession in self.predictedSet:
					if ontology in self.predictedSet[accession]:
						nP[ontology]+=1	
						# To calculate Rc excluding not predicted target
						nR[ontology]+=1
						nTerm[ontology]+=len(self.referenceSet[accession][ontology])
						
						
						# Calculate Prediction(P), CorrectPrediction(C) for each score threshold
						bins=dict([(th,set()) for th in score[ontology]])
						
						for go in self.predictedSet[accession][ontology]:
							scoreVal=self.predictedSet[accession][ontology][go]['score']
								
							if scoreVal not in bins:	
								bins[scoreVal]=set()
								score[ontology][scoreVal]={'C':0,'P':0,'pr':0.0,'rc':0.0,'f':0.0,'cov':0.0}
							
							for th in bins:
								if scoreVal>=th:
									bins[th].add(go)
							
							
						for th in score[ontology]:
							P=len(bins[th])
							T=len(self.referenceSet[accession][ontology])
							C=len(bins[th].intersection(set(self.referenceSet[accession][ontology].keys())))
							score[ontology][th]['C']+=C
							score[ontology][th]['P']+=P

							
		# Calculate precision,recall,F-measure,coverage
		for ontology in score:
			for th in sorted(score[ontology].keys()):
				if nTerm[ontology]!=0: 
					if score[ontology][th]['P']!=0:
						pr=float(score[ontology][th]['C'])/score[ontology][th]['P']
					else:
						pr=0.0
					rc=float(score[ontology][th]['C'])/nTerm[ontology]
					f=2.0*score[ontology][th]['C']/float(score[ontology][th]['P']+nTerm[ontology])
					coverage=float(nP[ontology])/nR[ontology]					
	
					score[ontology][th]['pr']=pr
					score[ontology][th]['rc']=rc
					score[ontology][th]['f']=f
					score[ontology][th]['cov']=coverage
					
		self.precRecall=score

def readCOGIC(fileCogic):
	dataCogic={}
	for line in fileCogic.readlines():
		if line[0]!="#":
			accession,data=line.strip().split("\t")
			ontology,score=data.split(",")
			if accession not in dataCogic:
				dataCogic[accession]={}
			try:
				dataCogic[accession][ontology]=float(score)
			except:
				dataCogic[accession][ontology]=None
	return dataCogic
				
			
		
		
