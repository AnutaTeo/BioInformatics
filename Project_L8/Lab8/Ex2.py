#Ex2  Implement a software application to detect the positions of these transposable elements (start, end) within the created DNA sequence.

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

def find_transposons(genome, transposons):
    results = []
    for trans in transposons:
        start = 0
        while True:
            idx = genome.find(trans, start)
            if idx == -1:
                break
            results.append((trans, idx, idx + len(trans)))
            start = idx + 1
    return results

detections = find_transposons(genome, transposons)

print("\nDetected transposons (sequence, start, end):")
for item in detections:
    print(item)
