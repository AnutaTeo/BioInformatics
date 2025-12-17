import math

motifs = [
    "GTCATTACTA",
    "ACACAATAGA",
    "GCGAGGGGTG",
    "GGGGGGGGGG",
    "TTTTTTTTTT",
    "AATCCAAAGA",
    "AAGAACATAA",
    "AGGGTTCAGG",
    "CTATTGTCTT",
]

S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

bases = ["A", "C", "G", "T"]

L = len(motifs[0])
N = len(motifs)
background = 0.25
pseudocount = 1

# 1) Count matrix
count = {b: [0]*L for b in bases}

for m in motifs:
    for i, ch in enumerate(m):
        count[ch][i] += 1


# 2) Weight matrix
pwm = {b: [0]*L for b in bases}
for i in range(L):
    col_total = sum(count[b][i] for b in bases) + 4*pseudocount  # N + 4
    for b in bases:
        pwm[b][i] = (count[b][i] + pseudocount) / col_total


# 3) Relative frequencies matrix

freq = {b: [count[b][i]/N for i in range(L)] for b in bases}


# 4) Log-likelihood matrix: ln(pwm / 0.25)

ll = {b: [math.log(pwm[b][i] / background) for i in range(L)] for b in bases}

def print_matrix(title, mat, decimals=3):
    print("\n" + title)
    header = "pos  " + "  ".join(str(i+1).rjust(5) for i in range(L))
    print(header)
    for b in bases:
        row = b + "   " + "  ".join(f"{mat[b][i]:>{5}.{decimals}f}" for i in range(L))
        print(row)

print_matrix("1) COUNT matrix", count, decimals=0)
print_matrix("2) WEIGHT matrix (PWM probabilities, with +1 pseudocount)", pwm, decimals=3)
print_matrix("3) RELATIVE FREQUENCIES matrix (count/N)", freq, decimals=3)
print_matrix("4) LOG-LIKELIHOODS matrix ln(pwm/0.25)", ll, decimals=3)


# 5) Scan S with sliding windows of length 10

print("\n5) Sliding window scores (window length = 10):")
scores = []
for start in range(len(S) - L + 1):
    window = S[start:start+L]
    score = 0.0
    for i, ch in enumerate(window):
        score += ll[ch][i]
    scores.append((start+1, window, score))

# print all scores
for pos, w, sc in scores:
    print(f"pos {pos:2d}: {w}   score = {sc:.3f}")

# best hit
best = max(scores, key=lambda x: x[2])
print("\nBest window:")
print(f"pos {best[0]}: {best[1]}   score = {best[2]:.3f}")

# show windows that look "motif-like" 
print("\nWindows with positive score:")
for pos, w, sc in scores:
    if sc > 0:
        print(f"pos {pos:2d}: {w}   score = {sc:.3f}")
