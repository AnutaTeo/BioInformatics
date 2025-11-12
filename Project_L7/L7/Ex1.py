import os
from collections import Counter
import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt


def read_fasta_file(path):
    seq = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)

def sanitize(seq):
    return "".join(ch for ch in seq if ch in "ACGTN")

def find_repeats(sequence, min_len=6, max_len=10):
    seq = sanitize(sequence)
    counts = Counter()
    n = len(seq)

    for L in range(min_len, max_len + 1):
        if n < L:
            continue
        for i in range(n - L + 1):
            kmer = seq[i:i+L]
            if "N" in kmer:
                continue
            counts[kmer] += 1

    repeats = Counter({k: v for k, v in counts.items() if v > 1})
    return repeats

def plot_repeats(freqs, title, top=15, save_png=False, out_dir="plots"):
    if not freqs:
        messagebox.showinfo("No repeats", f"No 6-10 bp repeats found in {title}.")
        return

    most_common = freqs.most_common(top)
    labels, values = zip(*most_common)

    plt.figure(figsize=(10, 5))
    plt.bar(labels, values)
    plt.xticks(rotation=90)
    plt.title(f"Top {min(top, len(labels))} repeats: {title}")
    plt.xlabel("Repeat (6-10 bases)")
    plt.ylabel("Frequency")
    plt.tight_layout()

    if save_png:
        os.makedirs(out_dir, exist_ok=True)
        safe = os.path.splitext(os.path.basename(title))[0]
        out_path = os.path.join(out_dir, f"{safe}_repeats.png")
        plt.savefig(out_path, dpi=150)
        print(f"[SAVED] {out_path}")

    plt.show()

def plot_all_combined(all_freqs, file_names, top=10):
    n = len(all_freqs)
    if n == 0:
        messagebox.showinfo("No data", "No repeats found in selected files.")
        return

    cols = 2
    rows = (n + 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(14, rows * 4))
    axes = axes.flatten()

    for i, (freqs, name) in enumerate(zip(all_freqs, file_names)):
        if not freqs:
            axes[i].set_title(f"{name}\n(no repeats)")
            axes[i].axis("off")
            continue

        top_items = freqs.most_common(top)
        labels, values = zip(*top_items)
        axes[i].bar(labels, values)
        axes[i].set_title(name)
        axes[i].tick_params(axis='x', rotation=90)
        axes[i].set_ylabel("Freq")

    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.suptitle("Top Repeats Across All Influenza Genomes", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

def analyze_single_sequence():
    path = filedialog.askopenfilename(
        title="Select ONE FASTA (1000-3000 nt)",
        filetypes=(("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*"))
    )
    if not path:
        print("No file selected.")
        return

    seq = read_fasta_file(path)
    length = len(seq)
    print(f"[INFO] Loaded {os.path.basename(path)} length={length} bp")

    if length < 1000 or length > 3000:
        messagebox.showwarning(
            "Length check",
            f"Sequence is {length} bp (expected 1000-3000). "
            "We'll still analyze it."
        )

    reps = find_repeats(seq, 6, 10)
    print(f"[INFO] Found {len(reps)} repeated motifs (count > 1).")
    plot_repeats(reps, title=os.path.basename(path), top=15, save_png=True)

def analyze_multiple_influenza():
    paths = filedialog.askopenfilenames(
        title="Select influenza genome FASTA files",
        filetypes=(("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*"))
    )
    if not paths:
        print("No files selected.")
        return

    all_freqs = []
    names = []

    print(f"[INFO] Selected {len(paths)} files.")
    for idx, path in enumerate(paths, start=1):
        name = os.path.basename(path)
        try:
            seq = read_fasta_file(path)
            print(f"\n[{idx}/{len(paths)}] {name}: length={len(seq)} bp")
            reps = find_repeats(seq, 6, 10)
            print(f"    repeats found: {len(reps)} motifs with count > 1")

            all_freqs.append(reps)
            names.append(name)

            plot_repeats(reps, title=name, top=15, save_png=True)

        except Exception as e:
            messagebox.showerror("Error", f"Failed on {name}:\n{e}")

    plot_all_combined(all_freqs, names, top=10)

def main():
    root = tk.Tk()
    root.title("DNA Repeats (6-10 bp) - Simple Analyzer")

    tk.Label(root, text="Pick what you want to do:", font=("Segoe UI", 11, "bold")).pack(padx=12, pady=10)

    tk.Button(root, text="1) Analyze ONE sequence",
              command=analyze_single_sequence, width=40).pack(padx=12, pady=6)

    tk.Button(root, text="2) Analyze MULTIPLE genomes",
              command=analyze_multiple_influenza, width=40).pack(padx=12, pady=6)

    tk.Button(root, text="Quit", command=root.destroy, width=20).pack(padx=12, pady=10)

    root.mainloop()

if __name__ == "__main__":
    main()
