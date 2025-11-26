enzymes = {
    "EcoRI":   {"seq": "GAATTC", "cut_pos": 1},
    "BamHI":   {"seq": "GGATCC", "cut_pos": 1},
    "HindIII": {"seq": "AAGCTT", "cut_pos": 1},
    "TaqI":    {"seq": "TCGA",   "cut_pos": 1},
    "HaeIII":  {"seq": "GGCC",   "cut_pos": 2}
}

dna = (
    "AGGAGTATTCAACGTGATGAAGTCGCAGGGTTAACGTGGGAATGGTGCTTCTGTCCTAACAGGTTAGGGTATAATGCTGGAACCGTCCCCCAAGCGTTCAG"
    "GGTGGGCTTTGCTACGACTTCCGAGTCCAAAGACTCCCTGTTTTCGAAATTTGCGCTCAAGGGCGAGTATTGGACCTGGCTTACGCCTTAGTACGTAGCAAGGTGACACAAGCACAGTAGATCCTGCCCGCGTTTCCT"
    "ATGTATCAAGTTAGAACTTATGGAATATAATAACATGTGGATGGCCAGTGGTCGGTTGTTACACGCCTGCCGCAACGTTGAAAGACCCGGATTAGACTGGCAAGATCTATGGCGTGAGACCCGTTATGC"
    "TCCATTACGGTCAGTGGGTCACAGCTAGTTGTGGATTGGATTGCCATTCTCCGAGTGTATTACCGTGACGGCCGCACGGGTCCCATATAATGCAATCGTAGTCTACCTGACTGTACTTAGAAATGTG"
    "GCTTCGCCTTTGCCCACGCACCTGATCGCTCCTCGTTTGCTTTTAAGGACCGGACGAACCACAGAGCATTAGAAGAATCTCTAGCTGCTTTACAAAGTGCTGGTTCCTTTTCCAGCGGGATGTTTT"
    "ATCTAAACGCAATGAGAGAGGTATTCCTCAGGCCACATCGCTTCCTAGTTCCGCTGGGATCCATCGTTGGCGGCCGAAGCCGCCATTCCATAGTGAGTTCTTCGTCTGAGTCATTCTGTGCCAGATC"
    "GACTGACAGATAGCGGATCCAGTTTATCCCTCGAAACTATAGACGTACAGGTCGAAATCTTAAGTCAAATCGCGCGTCTAGACTCAGCTCTATTTTAGTGGTCATGGGTTCTGGTCCCCCCGAGC"
    "GGCGCAACCGATTAGGACCATGTAGAACATTACTTATAAGTCATCTTTTAAACACAATCTTCCTGCTCAGTGGTACATGGTTTTCGCTATTGCTAGCCAGCCTCATAAGTAACACCACTACTGCGAC"
)  

def find_cleavages(seq, recog, cut_pos):
    positions = []
    start = 0
    while True:
        idx = seq.find(recog, start)
        if idx == -1:
            break
        positions.append(idx + cut_pos)
        start = idx + 1
    return positions

def fragment_lengths(seq, cleavage_sites):
    if not cleavage_sites:
        return [len(seq)]

    cleavage_sites = sorted(cleavage_sites)
    fragments = []

    prev = 0
    for cut in cleavage_sites:
        fragments.append(cut - prev)
        prev = cut

    fragments.append(len(seq) - prev)
    return fragments


def gel_simulation(all_fragments):
    print("\n===== ELECTROPHORESIS GEL (SIMULATION) =====\n")
    print("Large fragments = top, Small fragments = bottom\n")

    max_len = max(all_fragments)
    scale = 50 / max_len

    sorted_frags = sorted(all_fragments, reverse=True)

    for f in sorted_frags:
        line = "|" + "#" * int(f * scale)
        print(f"{f:5d} bp  {line}")


all_fragments_for_gel = []

for name, data in enzymes.items():
    recog = data["seq"]
    cut_pos = data["cut_pos"]

    print("\n========================================")
    print(f"Enzyme: {name}")
    print("Recognition site:", recog)

    cuts = find_cleavages(dna, recog, cut_pos)
    print("Cleavage positions:", cuts)
    print("Number of cleavages:", len(cuts))

    frags = fragment_lengths(dna, cuts)
    print("Fragments:", frags)
    print("Total fragments:", len(frags))

    all_fragments_for_gel.extend(frags)

gel_simulation(all_fragments_for_gel)

