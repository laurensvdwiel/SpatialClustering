import numpy as np
from joblib import Parallel, delayed, cpu_count

### Spatial clustering computation ###

def spatial_clustering_parallel_lelieveld(k, d_ij, n, n_permutations, random_state, parallel):
    """Generates random arrays with same length of k and computes the pair-
    wise distance of that array. Then compares it with the pair-wise
    distance of k and computes a p-value
    see Lelieveld et al. 2017 PMID: 28867141"""
    # use the settings as in Lelieveld et al.
    pseudo_count = 1
    normalize = True
    
    # compute the actual cluster d_k
    d_k = compute_geometric_mean_distance(d_ij, k, pseudo_count, normalize)
    
    # create a list of permutations based on the number of variants
    # compute for each permutation the distance metric    
    d_R = []
    if parallel:
        d_R = Parallel(n_jobs=CalculateNumberOfActiveThreads(n_permutations))(delayed(spatial_clustering_lelieveld_single_permutation)(k, d_ij, n, random_state, pseudo_count, normalize) for _ in range(n_permutations))
    else:
        d_R = [spatial_clustering_lelieveld_single_permutation(k, d_ij, n, random_state, pseudo_count, normalize) for _ in range(n_permutations)]

    # compute how often the observed distance is more extreme than the permutation (with an added pseudocount to limit for extremes)
    P_value = (np.sum([d_k > d_r for d_r in d_R])+1)/(n_permutations+1)
    
    return d_k, P_value

def spatial_clustering_lelieveld_single_permutation(k, d_ij, n, random_state, pseudo_count, normalize):
    """Computation of a single permutation"""
    r = random_state.choice(n, size=len(k), replace=True)
    d_r = compute_geometric_mean_distance(d_ij, r, pseudo_count, normalize)
    return d_r

def compute_geometric_mean_distance(d_ij, M, pseudo_count = None, normalize = False):
    """Computes the geometric-mean distance positions in M, base on
    the euclidean distance matrix d_ij. 
    For the settings used in Lelieveld et al. 2017 PMID: 28867141
    use a pseudo_count of 1 and set normalize to True"""
    # specify the distances between the positions in M
    M_d = []
    for i in M:
        for j in M:
            if i < j:
                if np.isnan(d_ij[i][j]):
                    print("WARNING: d_ij for i='"+str(i)+"' and j='"+str(j)+"' is NaN")
                else:
                    M_d.append(d_ij[i][j])
    
    # Convert the M_d array to the numpy format
    M_d = np.array(M_d)
    
    if not pseudo_count is None:
        # add pseudo_count to the array
        M_d += pseudo_count
    
    if normalize:
        # retrieve the length for normalizing
        n = len(d_ij)
        
        if not pseudo_count is None:
            # add pseudo count to the normalization
            n += pseudo_count
        
        # Normalize M_d
        M_d /= n
    
    # set k as number of variants    
    k = len(M)
    
    # compute geometric mean distance as d_k
    d_k = M_d.prod(dtype='float128')**(1.0/k)
    
    # return geometric mean distance
    return d_k

### Parallel Helper Functions ###
def MaxNumberOfThreads():
    return cpu_count()

def CalculateNumberOfActiveThreads(numberOfTasks):
    if(cpu_count() == 2):
        return cpu_count()
    elif numberOfTasks < cpu_count():
        return numberOfTasks
    else:
        return cpu_count()

### Main function  ###

def main(gene_name, variant_cDNA_locations, cDNA_length, n_permutations, parallel, random_seed, correction):
    # Report the input
    print("Computing spatial clustering for gene: "+gene_name)
    print("cDNA_length: "+str(cDNA_length))
    print("variant_cDNA_locations: "+str(variant_cDNA_locations))
    print("Settings: random_seed: "+str(random_seed) +", parallel: "+str(parallel)+", n_permutations: "+str(n_permutations))
    print()

    # Specify the settings
    randomstate = np.random.RandomState(random_seed)

    ## Mutation array for the current gene
    M = np.array([int(x) for x in list(variant_cDNA_locations)])

    # cDNA length of the CDS
    n = cDNA_length

    # Correct the mutation locations for numbering starting at zero instead of one
    M_corrected = M - 1

    # Compute the distance matrix
    n_arr = range(n)
    d_ij = np.zeros(shape=(n,n))
    for i in n_arr:
        for j in n_arr:
            d_ij[i][j] = np.abs(i-j)

    # compute the clustering
    d_k, P_value = spatial_clustering_parallel_lelieveld(M_corrected, d_ij, n, n_permutations, randomstate, parallel)


    # Report the results
    print("Results:")
    print("gene: "+gene_name+', with n variants: '+str(len(M)))
    print("geometric_mean: "+str(d_k))
    print("corrected p-value: "+str(np.minimum(P_value*correction, 1.0)) + ' (Bonferroni correction = '+str(correction)+')')
    print()

### command line argumentation parsing ###

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Compute the spatial clustering of variant locations over cDNA')

    # required arguments
    parser.add_argument('--gene_name', type=str, required=True,
                        help='(Required) Name of the gene of interest, example usage: --gene_name=BRCA1')
    parser.add_argument('--variant_cDNA_locations', type=str, required=True,
                        help='(Required) cDNA based variant locations, example usage: --variant_cDNA_locations=10,50,50,123')
    parser.add_argument('--cDNA_length', type=int, required=True,
                        help='(Required) total cDNA length of the gene (including stop codon), example usage: --cDNA_length=1337')
    
    # optional arguments
    parser.add_argument('--n_permutations', type=int, required=False, default=100000000,
                        help='(Optional) total nunber of permutations, default=100000000 (1.00E+08), example usage: --n_permutations=100')
    parser.add_argument('--parallel', type=bool, required=False, default=True,
                        help='(Optional) should the algorithm make use of parallel computation?, default=True, example usage: --parallel=True')
    parser.add_argument('--random_seed', type=int, required=False, default=1,
                        help='(Optional) The seed used for initialization of the random permutations, default=1, example usage: --random_seed=1')
    parser.add_argument('--correction', type=int, required=False, default=1,
                        help='(Optional) The number of genes the p-value must be corrected for in a Bonferonni manner, default=1, example usage: --correction=1')

    args = parser.parse_args()
    main(gene_name=args.gene_name, variant_cDNA_locations=args.variant_cDNA_locations.split(','), cDNA_length=args.cDNA_length, n_permutations=args.n_permutations, parallel=args.parallel, random_seed=args.random_seed, correction=args.correction)