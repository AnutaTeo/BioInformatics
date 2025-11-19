#Ex1 Make an artificial DNA sequence of 200-400b in length, in which to simulate 3-4 transposable elements.

import random

def random_dna(length):
    return ''.join(random.choice("ATGC") for _ in range(length))

transposons = [
    "ATGCGTACGA",
    "TTACGGAATA",
    "CGTATACGGT",
    "GGAATTCCTA"
]

genome = random_dna(250)

positions = []
for trans in transposons:
    pos = random.randint(0, len(genome))
    positions.append((trans, pos))
    genome = genome[:pos] + trans + genome[pos:]

print("Generated DNA sequence (length={}):".format(len(genome)))
print(genome)
print("\nInserted transposon positions (transposon, start_index):")
for t, p in positions:
    print(t, p)
