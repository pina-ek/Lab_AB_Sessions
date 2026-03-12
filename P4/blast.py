# E. Pina 

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

        for q_words, q_pos in query_kmers:
            a=[q_words[i] for i in range(k)]
            score = sum(subst_mat[a])
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

def extend_seeds(seeds, query, db_seqs, k, X=15):
    """Extend seeds in both directions using the BLOSUM62 matrix. Return a list of hits.

    >>> seeds = [{'db_iseq': 1, 'db_pos': 2, 'q_pos': 0, 'q_kmer': 'MKT', 'db_kmer': 'MKT', 'score': 15.0}, {'db_iseq': 0, 'db_pos': 0, 'q_pos': 0, 'q_kmer': 'MKT', 'db_kmer': 'MRT', 'score': 12.0}]
    >>> db_seqs = ["MRTAY", "KTLKT"]
    >>> extend_seeds(seeds, "MKTAY", db_seqs, 3, X=15)
    [{'db_iseq': 1, 'query_start': 0, 'query_end': 2, 'db_start': 2, 'db_end': 4, 'query_seq': 'MKT', 'db_seq': 'LKT', 'seed': 'MKT'}, {'db_iseq': 0, 'query_start': 0, 'query_end': 4, 'db_start': 0, 'db_end': 4, 'query_seq': 'MKTAY', 'db_seq': 'MRTAY', 'seed': 'MKT'}]
    """
    hits = []  # will contain actual hits
    seen = (
        set()
    )  # to avoid duplicates, stores tuples (db_iseq, q_start, q_end, db_start, db_end)

    for s in seeds:
        q_pos = s["q_pos"]
        db_iseq = s["db_iseq"]
        db_pos = s["db_pos"]
        db_seq = db_seqs[db_iseq]
        # --- Extend right ---
        score, best_score, best_right = 0, 0, 0
        right = 0
        while (
            q_pos + k + right < len(query)
            and db_pos + k + right < len(db_seq)
            and score >= best_score - X
        ):
            score += subst_mat[query[q_pos + k + right], db_seq[db_pos + k + right]]
            if score > best_score:
                best_score = score
                best_right = right + 1
            right += 1

        # --- Extend left ---
        score, best_score, best_left = 0, 0, 0
        left = 0
        while (
            q_pos - left - 1 >= 0 and db_pos - left - 1 >= 0 and score >= best_score - X
        ):
            score += subst_mat[query[q_pos - left - 1], db_seq[db_pos - left - 1]]
            if score > best_score:
                best_score = score
                best_left = left + 1
            left += 1

        q_start = q_pos - best_left
        q_end = q_pos + k + best_right - 1
        db_start = db_pos - best_left
        db_end = db_pos + k + best_right - 1

        key = (db_iseq, q_start, q_end, db_start, db_end)
        if key not in seen:
            seen.add(key)
            hits.append(
                {
                    "db_iseq": db_iseq,
                    "query_start": q_start,
                    "query_end": q_end,
                    "db_start": db_start,
                    "db_end": db_end,
                    "query_seq": query[q_start : q_end + 1],
                    "db_seq": db_seq[db_start : db_end + 1],
                    "seed": s["q_kmer"],
                }
            )

    return hits


def merge_overlapping(hits, query, db_seqs):
    """
    Merge overlapping hits per (database sequence, diagonal).

    Group by diagonal (db_start - query_start) so that only hits on the
    same alignment path are merged. Hits on different diagonals represent
    different alignment positions and must NOT be merged.

    >>> h1 = {"db_iseq": 0, "query_start": 0, "query_end": 4, "db_start": 0, "db_end": 4, "query_seq": "MKTAY", "db_seq": "MKTAY", "seed": "MKT"}
    >>> h2 = {"db_iseq": 0, "query_start": 2, "query_end": 6, "db_start": 2, "db_end": 6, "query_seq": "TAYIA", "db_seq": "TAYIA", "seed": "TAY"}
    >>> merged = merge_overlapping([h1, h2], "MKTAYIAK", ["MKTAYIAK"])
    >>> merged[0]['query_seq']
    'MKTAYIA'
    >>> h1 = {"db_iseq": 1, "query_start": 0, "query_end": 2, "db_start": 2, "db_end": 4, "query_seq": "MKT", "db_seq": "LKT", "seed": "MKT"}
    >>> h2 = {"db_iseq": 0, "query_start": 0, "query_end": 4, "db_start": 0, "db_end": 4, "query_seq": "MKTAY", "db_seq": "MRTAY", "seed": "MKT"}
    >>> merged = merge_overlapping([h1, h2], "MKTAY", ["MKTAY", "KTMKT"])
    >>> len(merged)
    2
    """
    merged_hits = []

    groups = {}
    for h in hits:
        diagonal = h["db_start"] - h["query_start"]
        key = (h["db_iseq"], diagonal)
        if key not in groups:
            groups[key] = []
        groups[key].append(h)

    for (db_iseq,_diagonal,), group in groups.items():  # iterate over each db sequence and diagonal
        group.sort(
            key=lambda x: x["query_start"]
        )  # sort hits by query_start to ensure they are in the correct order for merging
        current = dict(
            group[0]
        )  # start with the first hit in the group as the current hit to merge with subsequent hits

        for h in group[1:]:
            if h["query_start"] <= current["query_end"] + 1:
                current["query_end"] = max(current["query_end"], h["query_end"])
                current["db_end"] = max(current["db_end"], h["db_end"])
                current["query_seq"] = query[
                    current["query_start"] : current["query_end"] + 1
                ]
                current["db_seq"] = db_seqs[db_iseq][current["db_start"]: current["db_end"] + 1
                ]
            else:
                merged_hits.append(current) # store the hit
                current = dict(h) # start a new current hit with the next hit in the group

        merged_hits.append(current)

    return merged_hits