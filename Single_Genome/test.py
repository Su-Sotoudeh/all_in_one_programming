from multiprocessing import Pool

def nucleotid_composition(seq, organism_name, chr_id):
    A = seq.count('A')
    T = seq.count('T')
    C = seq.count('C')
    G = seq.count('G')
    return organism_name, chr_id, A, T, C, G

if __name__ == "__main__":
    chromosomes_pool = ['ATCG','AAAA','G','CCCCC','ACACA','AAAAAAAAAAAA','A','TT']
    chr_id = [1,2,3,4,5,6,7,8]
    organism = ["Human","Human","Human","Human","Human","Human","Human","Human"]

    # prepare input as list of tuples
    tasks = [(seq, organ,id ) for seq,organ,id in zip(chromosomes_pool,organism,chr_id)]

    with Pool(processes=4) as pool:
        results = pool.starmap(nucleotid_composition, tasks)

    print(results)
