import os
import random as rd
import heapq as hp

cwd = os.getcwd()
force_dat_path = cwd + '/force.dat'
para50_path = cwd + '/para50.txt'

n_set = 50
remain = 20
mutation_rate = 0.01
p_min = [0  ,  0.4, 0, -10, 0.2, 0, -10, -0.5, 0, -45, 0];
p_max = [0.1,  0.6, 1,   0, 0.6, 1,   0,  0.5, 0,  45, 45];

def new_gen(foils):
    #get best DNAs
    pool = selection(foils) #list of DNA
    #empty children list to replace DNAs
    children = []
    for i in range(n_set):
        #choose two parents
        par1, par2 = rd.choice(pool), rd.choice(pool)
        #create new DNA
        child = crossover(par1, par2)
        #do mutation stuff
        mutation(child)
        #put child into list, while sort by fitness
        hp.heappush(children, child)
        
    foils = children.copy()
    return foils
    
def selection(foils): 
    #選出最好的 X個 DNA序列
    # for i in range(len(DNAs)):
    #     DNAs[i].calc_fitness()
    best = hp.nlargest(remain, DNAs)
    # pool = []
    # for DNA in best:
    #     for i in range(DNA.fitness):
    #         pool.append(DNA)
    return best

                    
def crossover(par1,par2): #two parents, make child
    seq1, seq2 = par1.genes, par2.genes
    
    for i in range(target_str_len): #基因隨機配對
        if rd.random()<0.5:
            seq1[i], seq2[i] = seq1[i], seq2[i]
            
    mid = rd.randint(1 , target_str_len) - 1 #聯會
    for i in range(target_str_len):
        if i < mid:
            child_genes[i] = seq1[i]
        else:
            child_genes[i] = seq2[i]

    return child_genes
def mutation(foil):
	for i in range(len(foil.genes)):
		if rd.random()<mutation_rate:
			foil.genes[i] = p_min + (p_max-p_min)*
def get_best(foils):
    bestDNA = hp.nlargest(1,foils)[0]
    # s = ''.join(bestDNA.genes)
    # bestDNA.calc_fitness()
    return bestDNA

class DNA():
    def __init__(self, genes, fitness):
        self.genes = genes
        self.fitness = fitness

    # def calc_fitness(self):  #calculate fitness (how simular to target?)
    #                          #how many letters same as target?
    #     self.fitness = 0
    #     for i in range(target_str_len):
    #         if self.genes[i] == target_str[i]:
    #             self.fitness += 1
                
    #to make DNA() possible to compare
    def __lt__(self, other):
        return self.fitness < other.fitness


# In para file:
# [generation]
# [11 parameters] *50
# 
# In force file, each line:
# [lift] [drag] 
def get_fitness(force):
	return force[0]/force[0]

def readFoils():
	para50_file = open(para50_path, 'r')
	force_file = open(force_dat_path,'r')
	foil_list = []

	gen = int(para50_file.readline().strip('\n'))
	for i in range(n_set):
		p_line = para50_file.readline().strip('\n')
		f_line = force_file.readline().strip('\n')

		parameters = p_line.split(' ')
		for i in range(len(parameters)):
			parameters[i] = float(parameters[i])

		tmp = f_line.split(' ')
		for i in range(len(tmp)):
			tmp[i] = float(tmp[i])
		fitness = get_fitness(tmp)

		foil_list.append(DNA(parameters, fitness))

	return gen, foil_list

def write_output():
# Read PARSEC_para50.txt
gen,foils = readFoils() # List, contains foil obj
rd.seed(gen)
hp.heapify(self.foils)
new_files = new_gen(foils)
# 