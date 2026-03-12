## Estibaliz Pina - BBI - AB

def seq_identity(seq1, seq2):
    """
    Return the percentage of sequence identity. Exclude positions with gaps on any sequence
    >>> seq_identity("FASTCAT", "FATCAT")
    
    >>> seq_identity("FASTCAT", "FASTCAT")
    100.0
    >>> seq_identity("FASTCAT", "FASTRAT")
    85.7
    >>> seq_identity("-FASTCAT", "-FASTRAT")
    85.7
    >>> seq_identity("FASTCAT", "FA-TCAT")
    100.0
    >>> seq_identity("FASTCAT", "FAT-CAT")
    83.3
    >>> seq_identity("AFASTCAT", "-FASTRAT")
    85.7
    >>> seq_identity("FASTCAT", "AAAAAAA")
    28.6
    >>> seq_identity("FASTCAT", "AFAAAFA")
    0.0
    """
    m = 0
    total = 0
    if len(seq1)==len(seq2):
        for a, b in zip(seq1, seq2):
            if a == "-" or b == "-":
                continue  
            total += 1
            if a == b:
                m += 1
        return round(m / total * 100, 1)



if __name__ == "__main__":
	import doctest
	doctest.testmod(verbose=True)