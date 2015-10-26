
import sys
import numpy

def drawDirichlet(num_items, num_samples, size=10):
	pre_alphas = []
	for i in range(num_items):
		pre_alphas.append(numpy.random.random_integers(size))
	return numpy.random.dirichlet(pre_alphas, num_samples)

	
def generateCNVs(num_samples, chromosomes, file_prefix, cnv_rate):	
	for n in range(num_samples):
		cnvFile = open(file_prefix + '_N_' + str(n) + '.cnv.vcf', 'w')
		cnvFile.write('#CRHOM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n')
		for i in range(len(chromosomes)):
			if cnv_rate > 0.0:
				endpos = int(cnv_rate * chromosomes[i])
				cnvFile.write(str(i+1) + '\t1\t.\t.\t.\t.\t.\tSVTYPE=CNV;END=' + str(endpos) + '\tGT:TCN:MCN\t./.:2:1\t./.:3:1\n')
				cnvFile.write(str(i+1) + '\t' + str((endpos+1)) + '\t.\t.\t.\t.\t.\tSVTYPE=CNV;END=' + str(chromosomes[i]) + '\tGT:TCN:MCN\t./.:2:1\t./.:2:1\n')
			else:
				cnvFile.write(str(i+1) + '\t1\t.\t.\t.\t.\t.\tSVTYPE=CNV;END=' + str(chromosomes[i]) + '\tGT:TCN:MCN\t./.:2:1\t./.:2:1\n')
		cnvFile.close()
	
	
def generateSNVs(num_samples, num_subclones, num_mutations, betas, tumor_purity, coverage_depth, chromosomes, file_prefix, false_positive_rate):	
	snvFiles = []
	for n in range(num_samples):
		snvFiles.append(open(file_prefix + '_N_' + str(n) + '.snv.vcf', 'w'))
		snvFiles[n].write('#CRHOM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n')
	# now generate the mutations
	mutations = set()
	events_fraction = drawDirichlet(num_subclones, 1)
	mean_coverage = 0.0
	muts = []
	for i in range(num_subclones):
		events_count = int(events_fraction[0][i] * num_mutations)
		muts.append(events_count)
		for j in range(events_count):
			chr = numpy.random.random_integers(1, len(chromosomes))
			pos = numpy.random.random_integers(1, chromosomes[chr-1])
			while pos in mutations:
				pos = numpy.random.random_integers(1, chromosomes[chr-1])
			mutations.add(pos)	
			for n in range(num_samples):
				num_reads = numpy.random.binomial(coverage_depth, numpy.random.beta(2, 5))
				filter = 'PASS'
				if num_reads < 5:
					num_reads = 4
					filter = 'REJECT'
					
				cell_fraction = 0.5*tumor_purity[n]*betas[n][i]				
				# decide if this mutation will be a false positive				
				if false_positive_rate > 0:
					prob = numpy.random.uniform()
					if prob < false_positive_rate:		
						cell_fraction = numpy.random.uniform()
						filter = 'FALSEPASS'						# this will be treated as PASS by CTPsingle but will be discarded by the evaluation script
										
				alt_count = numpy.random.binomial(num_reads, cell_fraction)
				ref_count = num_reads - alt_count			
				freq = (alt_count)/(num_reads + 0.0)				
				mean_coverage += num_reads
				tumor_counts = '0|1:' + str(ref_count) + ':' + str(alt_count) + ':0:0:0:0:0:0:' + str(freq)
				normal_counts = '0|0:' + str(coverage_depth) + ':0:0:0:0:0:0:0:0.0e+00'		
				
				snvFiles[n].write(str(chr) + '\t' + str(pos) + '\t.\tA\tC\t.\t' + filter + '\tcluster_' + str(i) + '\tGT:FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ:PM\t' + normal_counts + '\t' + tumor_counts + '\n')
	mean_coverage = mean_coverage/(num_samples*sum(muts))
	clonal = muts[0]
	subclonal = sum(muts)-muts[0]
	for n in range(num_samples):
		snvFiles[n].close()		
		print 'none', file_prefix + '_N_' + str(n), num_subclones-1, tumor_purity[n], mean_coverage, clonal, subclonal, clonal/(clonal+subclonal+0.0)
	return muts
		
		
def writeGroundTruth(num_samples, num_subclones, tumor_purity, picked_tree, betas, muts):
	for n in range(num_samples):
		groundTruth = open(file_prefix + '_N_' + str(n) + '.ground.txt', 'w')
		groundTruth.write('num_subclones\ttumor_purity\tpicked_tree')
		for i in range(num_subclones):
			groundTruth.write('\tcluster_' + str(i) +'\tmut_' + str(i))
		# header is written, now write the rest
		groundTruth.write('\n' + str(num_subclones-1) + '\t' + str(tumor_purity[n]) + '\t' + str(picked_tree))
		for i in range(num_subclones):
			groundTruth.write('\t' + str(betas[n][i]) + '\t' + str(muts[i]))
		groundTruth.write('\n')
		groundTruth.close()		
	
			
def generateSubclones(num_samples, num_subclones, betas):
	if(num_subclones < 2):	# i.e. a clonal tumour
		for i in range(num_samples):
			betas[i][0] = 1.0
		return 0
	gammaFile = open('../GammaAdjMatrices/GammaMatrix' + str(num_subclones) + '.txt', 'r')
	description = gammaFile.readline()
	num_trees = int(description.split()[1])
	picked_tree = numpy.random.random_integers(0, num_trees-1)
	# generate alpha and beta frequencies			
	alphas = drawDirichlet(num_subclones, num_samples)	
	for i in range(num_trees):		
		for j in range(num_subclones):
			ln = gammaFile.readline()
			if i == picked_tree:
				arr = ln.split()
				for k in range(len(arr)):
					for n in range(num_samples):
						betas[n][j] += int(arr[k]) * alphas[n][k]
		gammaFile.readline()	# empty line		
	gammaFile.close()	
	return picked_tree 

	
def writeSimulation(coverage_depth, num_subclones, chromosomes, num_mutations, tumor_purity, num_samples, file_prefix, false_positive_rate, cnv_rate):	
	betas = numpy.zeros((num_samples, num_subclones))
	picked_tree = generateSubclones(num_samples, num_subclones, betas)			
	generateCNVs(num_samples, chromosomes, file_prefix, cnv_rate)
	muts = generateSNVs(num_samples, num_subclones, num_mutations, betas, tumor_purity, coverage_depth, chromosomes, file_prefix, false_positive_rate)	
	writeGroundTruth(num_samples, num_subclones, tumor_purity, picked_tree, betas, muts)
	
	
if __name__ == "__main__":
	
	if len(sys.argv) < 8:
		print 'Usage: python', sys.argv[0], '<output_folder> <sim_prefix> <num_simulations> <num_samples> <false_positive_rate> <cnv_rate> <use_deep_coverage>'
		sys.exit(2)
	
	chromosomes = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566]				
	print 'cancer_type samplename num_subclones purity ploidy num_clonal num_subclonal frac_clonal'	
	
	output_folder = sys.argv[1]
	sim_prefix = sys.argv[2]
	num_simulations = int(sys.argv[3])
	num_samples = int(sys.argv[4])
	false_positive_rate = float(sys.argv[5])
	cnv_rate = float(sys.argv[6])
	use_deep_coverage = int(sys.argv[7])	
	
	if num_simulations < 1 or num_samples < 1 or false_positive_rate > 0.99 or cnv_rate > 0.99:
		print 'Invalid parameters, please check your command!'
		print 'Usage: python', sys.argv[0], '<output_folder> <sim_prefix> <num_simulations> <num_samples> <false_positive_rate> <use_deep_coverage>'
		sys.exit(2)
	
	# default simulation parameters
	min_mutations = 500
	max_mutations = 5000
	min_subclones = 1
	max_subclones = 4
	min_purity = 0.3
	max_purity = 1.0
	min_coverage = 100
	max_coverage = 120		
	base_coverage = 10000
	
	for j in range(num_simulations):				
		if not use_deep_coverage:
			base_coverage = numpy.random.random_integers(min_coverage, max_coverage)				
		num_mutations = numpy.random.random_integers(min_mutations, max_mutations)
		num_subclones = numpy.random.random_integers(min_subclones, max_subclones)
		tumor_purity = numpy.random.uniform(min_purity, max_purity, num_samples)			
		file_prefix = output_folder + '/' + sim_prefix + '_sim_' + str(j)
		writeSimulation(base_coverage, num_subclones, chromosomes, num_mutations, tumor_purity, num_samples, file_prefix, false_positive_rate, cnv_rate)
		
	