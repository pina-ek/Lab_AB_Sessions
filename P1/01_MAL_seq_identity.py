
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
    c=0
    no=0
    l=len(seq1)
    if len(seq1)==len(seq2):
        if len(seq1)==1:
            if seq1[0]==seq2[0]:
                c+=1
            elif seq1[0]!=seq2[0]:
                if seq1[0]=="-" or seq2[0]=="-":
                    c+=1
                else:
                    no+=1
            iden=c/l*100
            print("%.1f" % iden)
        elif len(seq1)>1:
            if seq1[0]==seq2[0]:
                c+=1
            elif seq1[0]!=seq2[0]:
                if seq1[0]=="-" or seq2[0]=="-":
                    c+=1
                else:
                    no+=1
            seq_identity(seq1[1::],seq2[1::])
    



if __name__ == "__main__":
	import doctest
	doctest.testmod(verbose=True)