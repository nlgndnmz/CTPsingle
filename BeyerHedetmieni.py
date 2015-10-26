
import sys

# G is the gamma matrix
# par is the parent array
# n is the number of nodes
def writeGammaMatrix(gammaFile, G, par, n):
	for i in range(n):
		for j in range(n):
			G[i][j] = 0

	for i in range(n):
		G[i][i] = 1
		j = par[i]-1
		while j > -1:
			G[j][i] = 1
			j = par[j]-1

	for i in range(n):
		for j in range(n):
			gammaFile.write(str(G[i][j]) + ' ')
		gammaFile.write('\n')

	gammaFile.write('\n')


# par is the parent array
# n is the number of nodes
def writeAdjMatrix(adjFile, par, n):
	adjFile.write(str(n-1)) # number of edges
	for i in range(1, n):
		adjFile.write('  ' + str(par[i]-1) + '  ' + str(i) + '  ')
	adjFile.write('\n');


# writes a dot file to be processed with Graphviz
def writeDotFile(par, n, num):
	dotFile = open('./GammaAdjMatrices/dotFile' + str(n) + '_' + str(num-1) + '.dot', 'w')
	dotFile.write('digraph G { \n')
	for i in range(1, n):
		dotFile.write(str(par[i]-1) + ' -> ' + str(i) + ';\n')
	dotFile.write('}\n')
	dotFile.close()


# n is the number of nodes
# k is the max number of children allowed per node
def getRootedTrees(adjFile, gammaFile, n, k, writeDots):
	num = 0
	L = []
	par = []
	levels = []
	children = []
	G = []

	p = n-1
	q = 0

	for i in range(n):
		L.append(i)
		par.append(i)
		levels.append(i+1)
		children.append(0)
		G.append([])
		for j in range(n):
			G[i].append(0)

	while (p > 0):
		for i in range(n):
			children[i] = 0
		for i in range(n):
			children[par[i]] += 1

		if max(children) <= k:
			num += 1
			writeAdjMatrix(adjFile, par, n)
			writeGammaMatrix(gammaFile, G, par, n)
			if writeDots:
				writeDotFile(par, n, num)

		p = 0
		for i in range(n-1, -1, -1):
			if L[i] > 1:
				p = i
				break

		if p == 0:
			break

		for i in range(p-1, -1, -1):
			if L[i] == L[p] - 1:
				q = i
				break

		for i in range(p, n):
			L[i] = L[i-p+q]

		for i in range(1, n):
			x = L[i]
			par[i] = levels[x-1]
			levels[x] = i+1


# n is the number of nodes
# k is the max number of children allowed per node
def getNumTrees(n, k):
	num = 0
	L = []
	par = []
	levels = []
	children = []

	p = n-1
	q = 0

	for i in range(n):
		L.append(i)
		par.append(i)
		levels.append(i+1)
		children.append(0)

	while (p > 0):
		for i in range(n):
			children[i] = 0
		for i in range(n):
			children[par[i]] += 1

		if max(children) <= k:
			num += 1

		p = 0
		for i in range(n-1, -1, -1):
			if L[i] > 1:
				p = i
				break

		if p == 0:
			break

		for i in range(p-1, -1, -1):
			if L[i] == L[p] - 1:
				q = i
				break

		for i in range(p, n):
			L[i] = L[i-p+q]

		for i in range(1, n):
			x = L[i]
			par[i] = levels[x-1]
			levels[x] = i+1

	return num


if __name__ == "__main__":
	
	if len(sys.argv) < 2:
		print 'Usage: python', sys.argv[0], '<max_num_nodes> [<max_branching_factor>]'
		sys.exit(0)
		
	maxNumNodes = 1 + int(sys.argv[1])
	k = maxNumNodes

	if len(sys.argv) == 3:
		k = int(sys.argv[2])

	for i in range(2, maxNumNodes):
		x = getNumTrees(i, k)
		print 'Number of trees with ', i, ' nodes: ', x
		if x > 100000:
			print 'Sorry, that is too many trees to write on file. Aborting now, already written files are ok to use.'
			sys.exit(3)
			
		adjFile = open('./GammaAdjMatrices/AdjacencyMatrix' + str(i) + '.txt', 'w')
		gammaFile = open('./GammaAdjMatrices/GammaMatrix' + str(i) + '.txt', 'w')
		
		adjFile.write(str(i) + ' ' + str(x) + '\n\n')
		gammaFile.write(str(i) + ' ' + str(x) + '\n')
		getRootedTrees(adjFile, gammaFile, i, k, False)
		
		adjFile.close()
		gammaFile.close()

