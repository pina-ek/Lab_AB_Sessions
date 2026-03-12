def score_seqs(seq1, seq2, match, mismatch, gap):
    """
    Calculate the score of two sequences based on a matching score,
    a mismatching score, and a gap penalty.

    >>> score_seqs("THEFASTCAT", "THEFATCAT-", 1, -1, -2)
    -1
    >>> score_seqs("THEFASTCAT", "THEFATCA-T", 1, -1, -2)
    1
    >>> score_seqs("THEFA-TCAT", "THEFASTCAT", 1, -1, -2)
    7
    >>> score_seqs("THEFASTCAT", "THE", 1, -1, -2)
    >>> score_seqs("THE-FASTCAT", "THE-FASTCAT", 1, -1, -2)
    10
    """