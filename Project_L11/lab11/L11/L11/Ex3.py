seq1 = "AAATTAAA"
seq2 = "AAACCCAAA"

if len(seq1) != len(seq2):
    print("ERROR: Sequences are not aligned (different lengths)")
    print("Length seq1 =", len(seq1))
    print("Length seq2 =", len(seq2))
    exit()

L = len(seq1)

matches = 0
mismatches = 0
match_line = ""

for a, b in zip(seq1, seq2):
    if a == b:
        match_line += "I"
        matches += 1
    else:
        match_line += " "
        mismatches += 1

identity_score = (matches / L) * 100

mismatch_penalty_score = (matches - mismatches) / L

weighted_score = 0
for a, b in zip(seq1, seq2):
    weighted_score += 1 if a == b else -1
weighted_score /= L

print(seq1)
print(match_line)
print(seq2)

print("\nSimilarity scores:")
print(f"1) Identity (%)              = {identity_score:.2f}")
print(f"2) Mismatch penalty score    = {mismatch_penalty_score:.2f}")
print(f"3) Weighted similarity score = {weighted_score:.2f}")
