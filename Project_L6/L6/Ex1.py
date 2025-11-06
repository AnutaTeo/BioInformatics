#Gel electrophoresis is an analysis method implemented in all disciplines of life sciences. 
# The results of gel electrophoresis indicate the relative sizes of fragments, which is useful for restriction mapping and analyzing PCR fragments. 
#1. Take an arbitrary DNA sequence from the NCBI (National Center for Biotechnology), between 1000 and 3000 nucleotides (letters).
#2. Take 10 random samples from this sequence, between 100-3000 bases.
#3. Store these samples in an array.
#4. Simulate the migration of these DNA segments on the electrophoresis gel, based on their molecular weights - however, their length should be sufficient for this exercise (show a visual representation).


import random
import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt

def read_fasta_via_dialog():
    root = tk.Tk()
    root.withdraw()

    filepath = filedialog.askopenfilename(
        title="Select FASTA file",
        filetypes=(("FASTA files", "*.fasta"), ("All files", "*.*"))
    )
    if not filepath:
        print("No file selected.")
        return ""

    try:
        with open(filepath, "r") as f:
            lines = f.readlines()
        seq = "".join(line.strip() for line in lines if not line.startswith(">")).upper()
        return seq
    except Exception as e:
        messagebox.showerror("Error", f"Failed to read file:\n{e}")
        return ""

def simulate_and_plot(seq, n_fragments=10):
    n = len(seq)
    if n < 100:
        messagebox.showerror("Error", "Sequence too short (need at least 100 bp).")
        return

    fragments = []
    lengths = []
    for _ in range(n_fragments):
        frag_len = random.randint(100, min(3000, n))
        start = random.randint(0, n - frag_len)
        frag = seq[start:start + frag_len]
        fragments.append(frag)
        lengths.append(len(frag))

    print("Fragment lengths (bp):", lengths)

    K = 1000.0
    distances = [K / L for L in lengths]

    plt.figure(figsize=(6, 8))
    lane_x = 1.0

    for L, d in zip(lengths, distances):
        plt.plot([lane_x - 0.4, lane_x + 0.4], [d, d], linewidth=6)
        plt.text(lane_x + 0.5, d, f"{L} bp", va="center", fontsize=9)

    plt.gca().invert_yaxis()
    plt.title("Simulated Gel Electrophoresis")
    plt.xlabel("Lane")
    plt.ylabel("Migration distance (a.u.)")
    plt.xlim(0, 2)
    plt.ylim(0, max(distances) + 10)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    sequence = read_fasta_via_dialog()
    if sequence:
        simulate_and_plot(sequence, n_fragments=10)
