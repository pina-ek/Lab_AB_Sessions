from ini_score_seqs import score_seqs


def move_seq2(seq1, seq2):
    """
    >>> move_seq2("THEFASTCAT", "THEFATCAT")
    ['THEFASTCAT---------', 'THEFATCAT----------', '-THEFATCAT---------', '--THEFATCAT--------', '---THEFATCAT-------', '----THEFATCAT------', '-----THEFATCAT-----', '------THEFATCAT----', '-------THEFATCAT---', '--------THEFATCAT--', '---------THEFATCAT-', '----------THEFATCAT']
    >>> move_seq2("THEFASTCAT", "AFASTCAT")
    ['THEFASTCAT--------', 'AFASTCAT----------', '-AFASTCAT---------', '--AFASTCAT--------', '---AFASTCAT-------', '----AFASTCAT------', '-----AFASTCAT-----', '------AFASTCAT----', '-------AFASTCAT---', '--------AFASTCAT--', '---------AFASTCAT-', '----------AFASTCAT']
    >>> move_seq2("THEFASTCAT", "THECAT")
    ['THEFASTCAT------', 'THECAT----------', '-THECAT---------', '--THECAT--------', '---THECAT-------', '----THECAT------', '-----THECAT-----', '------THECAT----', '-------THECAT---', '--------THECAT--', '---------THECAT-', '----------THECAT']
    """
    # YOUR CODE HERE


def print_scores(seq1, seq2, match, mismatch, gap):
    """
    Score and print the alignments of two sequences with match, mismatch and gap.
    >>> print_scores('THEFASTCAT', 'AFASTCAT', 1, -1, -2)
    THEFASTCAT--------
    AFASTCAT---------- -12
    -AFASTCAT--------- -12
    --AFASTCAT-------- 2
    ---AFASTCAT------- -15
    ----AFASTCAT------ -16
    -----AFASTCAT----- -19
    ------AFASTCAT---- -22
    -------AFASTCAT--- -27
    --------AFASTCAT-- -28
    ---------AFASTCAT- -33
    ----------AFASTCAT -36
    """