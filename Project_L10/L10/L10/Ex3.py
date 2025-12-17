import os
import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt

WINDOW = 30
OUTPUT_FOLDER = "ODS"

def read_fasta_multi(path):
    seqs = {}
    name = None
    cur = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(cur).upper()
                name = line[1:].strip()
                cur = []
            else:
                cur.append(line)
        if name:
            seqs[name] = "".join(cur).upper()

    return seqs

def cg_percent(win):
    return 100 * sum(b in "CG" for b in win) / len(win)

def kappa_ic(win):
    n = len(win)
    total = 0

    for shift in range(1, n):
        matches = 0
        L = n - shift
        for i in range(L):
            if win[i] == win[i + shift]:
                matches += 1
        total += (matches / L) * 100

    return total / (n - 1)

def compute_stain(seq):
    cg_vals, ic_vals = [], []

    for start in range(0, len(seq) - WINDOW + 1):
        w = seq[start:start + WINDOW]
        cg_vals.append(cg_percent(w))
        ic_vals.append(kappa_ic(w))

    return cg_vals, ic_vals

def generate_ods():
    fasta_path = fasta_entry.get()

    if not fasta_path or not os.path.exists(fasta_path):
        messagebox.showerror("Error", "Please select a valid FASTA file.")
        return

    promoters = read_fasta_multi(fasta_path)

    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    centers = []

    for name, seq in promoters.items():
        print("Processing:", name)

        cg_vals, ic_vals = compute_stain(seq)

        ods_path = os.path.join(OUTPUT_FOLDER, f"{name}_ODS.txt")
        with open(ods_path, "w") as f:
            f.write("CG%\tIC\n")
            for x, y in zip(cg_vals, ic_vals):
                f.write(f"{x:.3f}\t{y:.3f}\n")

        plt.scatter(cg_vals, ic_vals, s=5)
        plt.xlabel("C+G %")
        plt.ylabel("Kappa IC")
        plt.title(f"ODS — {name}")
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_FOLDER, f"{name}_plot.png"))
        plt.close()

        cx = sum(cg_vals) / len(cg_vals)
        cy = sum(ic_vals) / len(ic_vals)
        centers.append((name, cx, cy))

    plt.figure()
    for name, cx, cy in centers:
        plt.scatter(cx, cy)
        plt.text(cx + 0.3, cy + 0.3, name, fontsize=6)

    plt.xlabel("C+G % (center)")
    plt.ylabel("Kappa IC (center)")
    plt.title("Centers of Objective Digital Stains")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_FOLDER, "centers_plot.png"))
    plt.show()

    messagebox.showinfo("Done", "ODS files and plots saved in the ODS folder!")

def browse_fasta():
    path = filedialog.askopenfilename(title="Select Promoter FASTA File",
                                      filetypes=[("FASTA files", "*.fasta *.fa *.txt")])
    if path:
        fasta_entry.delete(0, tk.END)
        fasta_entry.insert(0, path)

root = tk.Tk()
root.title("Promoter ODS Generator")
root.geometry("600x200")

label = tk.Label(root, text="Select promoter FASTA file:", font=("Arial", 12))
label.pack(pady=10)

fasta_entry = tk.Entry(root, width=60)
fasta_entry.pack()

browse_btn = tk.Button(root, text="Browse", command=browse_fasta)
browse_btn.pack(pady=5)

run_btn = tk.Button(root, text="Generate ODS", command=generate_ods, width=20)
run_btn.pack(pady=20)

root.mainloop()
