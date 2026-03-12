from Bio.Align import substitution_matrices

from blast_index import build_index

subst_mat = substitution_matrices.load("BLOSUM62")


def find_seeds(query, db, k, t):
    """
    Return all seeds whose BLOSUM62 score against a query k-mer meets
    threshold t, using a prebuilt k-mer index for lookups.

    >>> db = {'KTM': [(1, 0)], 'MKT': [(1, 2)], 'MRT': [(0, 0)], 'RTA': [(0, 1)], 'TAY': [(0, 2)], 'TMK': [(1, 1)]}
    >>> find_seeds("MKT", db, 3, 11)
    [{'db_iseq': 1, 'db_pos': 2, 'q_pos': 0, 'q_kmer': 'MKT', 'db_kmer': 'MKT', 'score': 15.0}, {'db_iseq': 0, 'db_pos': 0, 'q_pos': 0, 'q_kmer': 'MKT', 'db_kmer': 'MRT', 'score': 12.0}]
    """
    query_kmers = [query[i:i+k] for i in range(len(query)-k+1)]
    
    seeds = []
    for db_word, db_index in db.items():

        for q_words,q_pos in query:
            
            a=[q_words[i] for i in range(k)]
            score = sum(subst_mat[0])
            if score >= t:
                for db_iseq, db_pos in db_index:
                    seeds.append(
                        {
                            "db_iseq": db_iseq,
                            "db_pos": db_pos,
                            "q_pos": q_pos,
                            "q_kmer": q_words,  # original query k-mer
                            "db_kmer": db_word,  # matched db k-mer (may differ)
                            "score": score,
                        }
                    )
    return seeds
    

if __name__ == "__main__":
    k = 3
    t = 11
    db = {
        "KTM": [(1, 0)],
        "MKT": [(1, 2)],
        "MRT": [(0, 0)],
        "RTA": [(0, 1)],
        "TAY": [(0, 2)],
        "TMK": [(1, 1)],
    }
    print(find_seeds("MKT", db, k, t))