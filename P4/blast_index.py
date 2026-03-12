def build_index(db_seqs, k):
    """
    Preprocess the database into a k-mer index.
    Keys are sorted alphabetically.

    Returns a dict:
        { kmer: [(indx, pos), ...], ... }

    >>> build_index(["MRTAY", "KTMKT"], 3)
    {'KTM': [(1, 0)], 'MKT': [(1, 2)], 'MRT': [(0, 0)], 'RTA': [(0, 1)], 'TAY': [(0, 2)], 'TMK': [(1, 1)]}
    >>> build_index(["AAA"], 2)
    {'AA': [(0, 0), (0, 1)]}
    """
    db={}
    for i, seq, in enumerate(db_seqs):
        for j in range(len(seq)-k+1):
            kmer=seq[j:j+k]
            if kmer not in db:
                db[kmer]=[]
            db[kmer].append((i,j))
    return dict(sorted(db.items()))


if __name__ == "__main__":
    print(build_index(["MRTAY", "KTMKT"], 3))