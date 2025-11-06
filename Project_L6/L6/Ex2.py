#download 10 influenza virus genomes and apply the electrophosis gel simulation on each of them. 
# make a comparison between the 10 electroph gel simulation and show which of the influenza genomes show the most DNA segments. 
# you can plot them in the same graph, but also separately, because the lines may overlap. as the main restriction enzyme phase use ECOR1

import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt

ECO_RI_SITE = "GAATTC"

def read_fasta(path):
    with open(path, "r") as f:
        lines = f.readlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">")).upper()
    return seq

def digest_sequence(seq, site=ECO_RI_SITE):
    fragments = seq.split(site)
    fragment_lengths = [len(f) for f in fragments if len(f) > 0]
    return fragment_lengths

def simulate_gel(lengths, lane_number, label):
    K = 1000.0
    distances = [K / L for L in lengths]
    for d, L in zip(distances, lengths):
        plt.plot([lane_number - 0.3, lane_number + 0.3], [d, d], linewidth=5)
        plt.text(lane_number + 0.4, d, f"{L} bp", fontsize=7, va="center")
    plt.text(lane_number, max(distances) + 20, label, ha="center", fontsize=8)

def main():
    root = tk.Tk()
    root.withdraw()

    num_genomes = 10
    all_fragment_data = []

    print(f"Select {num_genomes} FASTA files (one by one):")
    for i in range(num_genomes):
        filepath = filedialog.askopenfilename(
            title=f"Select genome {i+1}",
            filetypes=(("FASTA files", "*.fasta"), ("All files", "*.*"))
        )
        if not filepath:
            print("File selection cancelled.")
            return

        seq = read_fasta(filepath)
        fragments = digest_sequence(seq)
        all_fragment_data.append((filepath.split("/")[-1], fragments))

    most_fragments = max(all_fragment_data, key=lambda x: len(x[1]))
    print("\nGenome with most EcoRI fragments:")
    print(f"→ {most_fragments[0]} ({len(most_fragments[1])} fragments)\n")

    plt.figure(figsize=(10, 8))
    plt.title("EcoRI Restriction Digest – All Genomes")
    plt.xlabel("Genome Lane")
    plt.ylabel("Migration Distance (a.u.)")

    for lane, (name, fragments) in enumerate(all_fragment_data, start=1):
        simulate_gel(fragments, lane, f"G{lane}")

    plt.gca().invert_yaxis()
    plt.xlim(0, num_genomes + 1)
    plt.tight_layout()
    plt.show()

    for lane, (name, fragments) in enumerate(all_fragment_data, start=1):
        plt.figure(figsize=(5, 7))
        plt.title(f"EcoRI Digest – Genome {lane}\n{name}")
        plt.xlabel("Lane")
        plt.ylabel("Migration Distance (a.u.)")
        simulate_gel(fragments, 1, f"G{lane}")
        plt.gca().invert_yaxis()
        plt.xlim(0, 2)
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    main()
