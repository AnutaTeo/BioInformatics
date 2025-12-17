import tkinter as tk
from tkinter import filedialog, ttk
import matplotlib.pyplot as plt
import math

def read_fasta_first_sequence(path):
    seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line.upper())
    return "".join(seq)

def smith_waterman_score(s1, s2, match=2, mismatch=-1, gap=-2):
    n = len(s1)
    m = len(s2)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    max_score = 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):

            if s1[i-1] == s2[j-1]:
                diag = H[i-1][j-1] + match
            else:
                diag = H[i-1][j-1] + mismatch

            delete = H[i-1][j] + gap
            insert = H[i][j-1] + gap

            H[i][j] = max(0, diag, delete, insert)

            if H[i][j] > max_score:
                max_score = H[i][j]

    return max_score


def smith_waterman_alignment(s1, s2, match=2, mismatch=-1, gap=-2):
    n = len(s1)
    m = len(s2)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):

            if s1[i-1] == s2[j-1]:
                diag = H[i-1][j-1] + match
            else:
                diag = H[i-1][j-1] + mismatch

            delete = H[i-1][j] + gap
            insert = H[i][j-1] + gap

            H[i][j] = max(0, diag, delete, insert)

            if H[i][j] > max_score:
                max_score = H[i][j]
                max_pos = (i, j)

    aligned1 = []
    aligned2 = []
    i, j = max_pos

    while i > 0 and j > 0 and H[i][j] != 0:
        current = H[i][j]
        diag = H[i-1][j-1]
        up = H[i-1][j]
        left = H[i][j-1]

        if current == diag + (match if s1[i-1] == s2[j-1] else mismatch):
            aligned1.append(s1[i-1])
            aligned2.append(s2[j-1])
            i -= 1
            j -= 1
        elif current == up + gap:
            aligned1.append(s1[i-1])
            aligned2.append('-')
            i -= 1
        elif current == left + gap:
            aligned1.append('-')
            aligned2.append(s2[j-1])
            j -= 1
        else:
            break

    aligned1.reverse()
    aligned2.reverse()
    return "".join(aligned1), "".join(aligned2), max_score

def split_into_windows(seq, window_size, step):
    windows = []
    for start in range(0, len(seq) - window_size + 1, step):
        windows.append((start, seq[start:start + window_size]))
    return windows

def choose_files_window():
    root = tk.Tk()
    root.title("Select Influenza & COVID-19 FASTA Files")
    root.geometry("600x220")

    influenza_var = tk.StringVar()
    covid_var = tk.StringVar()

    def choose_influenza():
        path = filedialog.askopenfilename(
            title="Select Influenza Genome FASTA",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffn *.faa"), ("All files", "*.*")]
        )
        if path:
            influenza_var.set(path)

    def choose_covid():
        path = filedialog.askopenfilename(
            title="Select COVID-19 Genome FASTA",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffn *.faa"), ("All files", "*.*")]
        )
        if path:
            covid_var.set(path)

    frame = ttk.Frame(root, padding=15)
    frame.pack(fill="both", expand=True)

    ttk.Label(frame, text="Influenza genome FASTA:").pack(anchor="w")
    ttk.Entry(frame, textvariable=influenza_var, width=70).pack()
    ttk.Button(frame, text="Browse...", command=choose_influenza).pack()

    ttk.Label(frame, text="\nCOVID-19 genome FASTA:").pack(anchor="w")
    ttk.Entry(frame, textvariable=covid_var, width=70).pack()
    ttk.Button(frame, text="Browse...", command=choose_covid).pack()

    def done():
        if influenza_var.get() and covid_var.get():
            root.destroy()

    ttk.Button(frame, text="Confirm", command=done).pack(pady=12)

    root.mainloop()
    return influenza_var.get(), covid_var.get()


def main():

    influenza_file, covid_file = choose_files_window()

    print("\nReading genomes...\n")

    influenza = read_fasta_first_sequence(influenza_file)
    covid = read_fasta_first_sequence(covid_file)

    print("Influenza length:", len(influenza))
    print("COVID-19 length:", len(covid))

    window_size = 200
    step = 200

    influenza_windows = split_into_windows(influenza, window_size, step)
    covid_windows = split_into_windows(covid, window_size, step)

    n1 = len(influenza_windows)
    n2 = len(covid_windows)

    print("\nWindows (Influenza):", n1)
    print("Windows (COVID-19):", n2)

    similarity = [[0.0 for _ in range(n2)] for _ in range(n1)]

    best_score = -math.inf
    best_i = best_j = 0

    print("\nComputing similarity matrix...\n")
    match = 2
    mismatch = -1
    gap = -2

    for i, (_, w1) in enumerate(influenza_windows):
        for j, (_, w2) in enumerate(covid_windows):

            score = smith_waterman_score(w1, w2, match, mismatch, gap)
            max_possible = match * min(len(w1), len(w2))
            similarity[i][j] = score / max_possible if max_possible else 0

            if score > best_score:
                best_score = score
                best_i = i
                best_j = j

    print("Best local score:", best_score)

    plt.figure(figsize=(8, 6))
    plt.imshow(similarity, cmap="magma", aspect="auto", origin="lower")
    plt.colorbar(label="Normalized Local Alignment Score")
    plt.xlabel("COVID-19 windows")
    plt.ylabel("Influenza windows")
    plt.title("Influenza vs COVID-19 – Local Alignment Similarity (Windowed)")
    plt.tight_layout()
    plt.show()

    start1, w1 = influenza_windows[best_i]
    start2, w2 = covid_windows[best_j]

    aligned1, aligned2, true_score = smith_waterman_alignment(w1, w2)

    match_line = ""
    matches = 0
    for a, b in zip(aligned1, aligned2):
        if a == b and a != "-":
            match_line += "|"
            matches += 1
        else:
            match_line += " "

    identity = matches / len(aligned1) * 100 if aligned1 else 0

    print("\nBest Local Alignment Between Windows:")
    print("Influenza start:", start1)
    print("COVID start    :", start2)
    print("Score          :", true_score)
    print("\n" + aligned1)
    print(match_line)
    print(aligned2)
    print(f"\nIdentity: {identity:.2f}%\n")


if __name__ == "__main__":
    main()
