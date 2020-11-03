def rand():
    # function to calculate random number
    global rseed 
    rseed = (rseed * 1103515245 + 12345) & RAND_MAX
    return rseed

def get_count_profile(dna, motifs, k, t):
    # get count matrix
    profile = [[1 for i in range(k)] for j in range(4)]
    idx = {'A':0, 'T':1, 'G':2, 'C':3}
    for i,j in enumerate(motifs):
        for x in range(k):
            profile[idx[dna[i][j+x]]][x] +=1
    return profile


def profile_sub(dna, profile, k, exclude_i, motifi, t):
    # subtract the neucleotide counts from the matrix corresponding to excluded motif 
    idx = {'A':0, 'T':1, 'G':2, 'C':3}
    for i in range(k):
        profile[idx[dna[exclude_i][motifi+i]]][i] -= 1
    return profile

def profile_add(dna, profile, k, exclude_i, k_mer, t):
    # Add the neucleotide counts from the matrix corresponding to random generated motif
    idx = {'A':0, 'T':1, 'G':2, 'C':3}
    for i in range(k):
        profile[idx[dna[exclude_i][k_mer+i]]][i] += 1
    return profile

def inverse_sampling(prob, l):
    summ = sum(prob)
    distribution = []
    for p in range(l):
        distribution += [p]*(int((prob[p]/summ)*100))
    sum_per = len(distribution)
            
    for i in range(sum_per, 100):
        distribution.append(distribution[rand() % sum_per])
    return distribution[rand() % 100]
    
def profile_rand_generate_k_mer(dna, profile, exclude_i, t, n, k):
    profile = [[profile[j][i]/t for i in range(k)] for j in range(4)]
    prob = [0 for i in range(n)]
    idx = {'A':0, 'T':1, 'G':2, 'C':3}
    # storing probability at each ith k-mer in excluced motif
    for i in range(n):
        pr_i = 1
        for j in range(k):
            pr_i *= (profile[idx[dna[exclude_i][j+i]]][j])
        prob[i] = pr_i
    
    return inverse_sampling(prob, n)

def score(dna, k, t, profile):
    # calculate motif score by subtracting the max count value from total count
    score = t*k - sum([max(x) for x in zip(*profile)])
    return score

def GIBBSSAMPLER(dna, k, t, N):
    n = len(dna[0])-k+1 #define seq length
    # Random generated motif
    motifs = [rand()%(n) for x in range(t)] 
    # best motif
    best_motif = copy.copy(motifs) 
    count_profile = get_count_profile(dna, motifs, k, t)
    score_bestmotifs = score(dna, k, t+4, count_profile)
    # 9 Random starts (with 11 random starts the average run time for t = 20 and seq length = 200 is 20.23 sec)
    for r in range(20):
        motifs = [rand()%(n) for x in range(t)]
        count_profile = get_count_profile(dna, motifs, k, t)
          
        for j in range(N):
            exclude_i = rand()%(t) ##### generate random number
            count_profile = profile_sub(dna, count_profile, k, exclude_i, motifs[exclude_i], t-1)
            motifs.pop(exclude_i)
            k_mer = profile_rand_generate_k_mer(dna, count_profile, exclude_i, t+3, n, k)
            count_profile = profile_add(dna, count_profile, k, exclude_i, k_mer, t-1)
            motifs.insert(exclude_i,k_mer)
            score_motifs = score(dna, k, t+4, count_profile)
            if score_motifs<score_bestmotifs:
                best_motif = copy.copy(motifs)
                score_bestmotifs = score_motifs
    
    return best_motif, score_bestmotifs


import copy
from sys import stdin, stdout  
ktn = list(map(int, input().rstrip().split()))
k = ktn[0]
t = ktn[1]
N = ktn[2]
rseed = 0
RAND_MAX = (1 << 31) - 1
DNA = []
for i in range(t):
    dna = str(stdin.readline())
    DNA.append(dna.split('\n')[0])

best_motif, score_bestmotifs = GIBBSSAMPLER(DNA, k, t, N)
for i in range(t):
    stdout.write(str(DNA[i][best_motif[i]:best_motif[i]+k] + '\n'))