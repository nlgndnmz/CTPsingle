
import sys
import math

class TumorSample:
	'''Contains data about a particular tumor sample'''		
	def __init__(self):
		self.purity = 0.0
		self.nsubclones = 0
		self.tree = 0
		self.clusters = {}		
		self.mutations = {}
		
	def setUp(self, num_subclones, tumor_purity, picked_tree):
		self.purity = tumor_purity
		self.nsubclones = num_subclones
		self.tree = picked_tree
		
	def addCluster(self, cluster_name, cancer_cell_fraction, num_mutations):
		self.clusters[cluster_name] = (cancer_cell_fraction, num_mutations)		
		
	def addMutation(self, chromosome, position, cluster_name):
		self.mutations[(chromosome, position)] = cluster_name
		
		
def readVCF(file_prefix, truth):	
	with open(file_prefix + '.snv.vcf', 'r') as vcf:
		for line in vcf:
			columns = line.split()
			if columns[6] == 'PASS':
				truth.addMutation(columns[0], columns[1], columns[7])			
	
	
def readGroundTruth(file_prefix):
	groundTruth = open(file_prefix + '.ground.txt', 'r')	
	ln = groundTruth.readline()		# just the header
	header = ln.split()
	ln = groundTruth.readline()
	values = ln.split()
	truth = TumorSample()
	truth.setUp(int(values[0]), float(values[1]), int(values[2]))
	j = 1
	for i in range(truth.nsubclones+1):		
		j += 2
		truth.addCluster(header[j], float(values[j]), int(values[j+1]))			
	groundTruth.close()
	readVCF(file_prefix, truth)
	return truth
	
	
def getPurity(logfilename):	
	purity = 0.0
	with open(logfilename, 'r') as logfile:
		for line in logfile:
			columns = line.split()
			if len(columns) > 0 and columns[0] == 'Estimated':
				purity = float(columns[2])
	return purity
	

def readAssignments(file_prefix, prediction):
	frequencies = {}
	nmutations = {}
	nsubclones = 0
	with open(file_prefix + '_cluster_assignments.txt', 'r') as res:
		for line in res:
			columns = line.split()
			if columns[0] == 'Chromosome':
				continue			
			prediction.addMutation(columns[0], columns[1], columns[4])
			if columns[4] in frequencies:
				nmutations[columns[4]] += 1
			else:				
				frequencies[columns[4]] = float(columns[5])
				nmutations[columns[4]] = 0				
	for key in frequencies:		
		prediction.addCluster(key, frequencies[key], nmutations[key])
		if int(key) > 0:
			nsubclones += 1
	return nsubclones-1
			
			
def readPrediction(file_prefix, logfilename):
	prediction = TumorSample()
	prediction.purity = getPurity(logfilename)
	prediction.nsubclones = readAssignments(file_prefix, prediction)
	return prediction
	
	
def computeRMSD(difflist):
	rmsd = 0.0
	for i in range(len(difflist)):
		rmsd += difflist[i]**2
	rmsd = rmsd / len(difflist)
	return math.sqrt(rmsd)

	
def comparePrediction2Truth(truth, prediction):
	'''Compares the predicted values to the ground truth'''
	freqDiff = []
	for key in prediction.mutations:
		if not key in truth.mutations:		# could be a false positive mutation
			continue
		trueCluster = truth.mutations[key]
		predictedCluster = prediction.mutations[key]
		trueFreq = truth.clusters[trueCluster][0]
		predictedFreq = prediction.clusters[predictedCluster][0]
		freqDiff.append(trueFreq-predictedFreq)	
	return computeRMSD(freqDiff)
	

if __name__ == "__main__":

	if len(sys.argv) < 8:
		print 'Usage: python', sys.argv[0], '<method> <input_folder> <results_folder> <run_folder> <sim_prefix> <num_simulations> <num_samples>'
		sys.exit(2)
		
	method_name = sys.argv[1]
	input_folder = sys.argv[2]
	results_folder = sys.argv[3]
	run_folder = sys.argv[4]	
	sim_name = sys.argv[5]	
	num_simulations = int(sys.argv[6])
	num_samples = int(sys.argv[7])
	
	print 'method simulation trueNsubclones truePurity trueTree predNsubclones predPurity predTree freqRMSD'		
	for i in range(num_simulations):							
		for j in range(num_samples):
			simulation = sim_name + str(i) + '_N_' + str(j)
			truth = readGroundTruth(input_folder + '/' + simulation)
			prediction = readPrediction(results_folder + '/' + simulation, run_folder + '/run_' + simulation)
			difference = comparePrediction2Truth(truth, prediction)			
			print method_name, simulation, truth.nsubclones, truth.purity, truth.tree, prediction.nsubclones, prediction.purity, prediction.tree, difference
		
	
	