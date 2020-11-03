# Gibbs-Sampler
Gibbs sampling technique to find motifs in a set of sequences

Profile-randomly generated k-mer in a string Text - for each k-mer Pattern in Text, compute the probability Pr(Pattern | Profile), resulting in n = |Text| - k + 1 probabilities (p1, …, pn). These probabilities do not necessarily sum to 1, but we can still form the random number generator Random(p1, …, pn) based on them. GIBBSSAMPLER uses this random number generator to select a Profile-randomly generated k-mer at each step: if the die rolls the number i, then we define the Profile-randomly generated k-mer as the i-th k-mer in Text.

GIBBSSAMPLER(Dna, k, t, N)
    randomly select k-mers Motifs = (Motif1, …, Motift) in each string
        from Dna
    BestMotifs ← Motifs
    r_start ← 40
    for r ← 1 to r_start
        Motifs ← randomly select k-mers Motifs = (Motif1, …, Motift) in each string
        from Dna
        for j ← 1 to N
            i ← Random(t)
            Profile ← profile matrix constructed from all strings in Motifs
                       except for Motifi
            Motifi ← Profile-randomly generated k-mer in the i-th sequence
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
    return BestMotifs
