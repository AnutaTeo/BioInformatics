import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt

def cg_percent(window):
    count = sum(1 for b in window if b in "CG")
    return (count / len(window)) * 100

def kappa_ic(window):
    n = len(window)
    total = 0

    for shift in range(1, n):
        matches = 0
        L = n - shift
        for i in range(L):
            if window[i] == window[i+shift]:
                matches += 1
        total += (matches / L) * 100

    return total / (n - 1)

def read_fasta(path):
    seq = ""
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                seq += line.upper()
    return seq

def compute_stain(seq, window=30):
    cg_vals = []
    ic_vals = []
    for start in range(0, len(seq) - window + 1):
        w = seq[start:start + window]
        cg_vals.append(cg_percent(w))
        ic_vals.append(kappa_ic(w))
    return cg_vals, ic_vals

def select_files():
    files = filedialog.askopenfilenames(title="Select up to 10 influenza genomes",
                                        filetypes=[("FASTA files", "*.fasta *.fa *.txt")])

    if not files:
        return

    if len(files) > 10:
        messagebox.showerror("Error", "Please select only 10 genomes.")
        return

    file_list.delete(0, tk.END) 
    for f in files:
        file_list.insert(tk.END, f)

def run_analysis():
    if file_list.size() == 0:
        messagebox.showerror("Error", "Please select genome FASTA files first.")
        return

    centers = []

    for i in range(file_list.size()):
        filepath = file_list.get(i)
        seq = read_fasta(filepath)
        name = filepath.split("/")[-1]

        print("Processing:", name)

        cg_vals, ic_vals = compute_stain(seq)

        plt.scatter(cg_vals, ic_vals, s=5)
        plt.title(f"Digital Stain — {name}")
        plt.xlabel("C+G %")
        plt.ylabel("Kappa IC")
        plt.show()

        cx = sum(cg_vals) / len(cg_vals)
        cy = sum(ic_vals) / len(ic_vals)
        centers.append((name, cx, cy))
    plt.figure()
    for name, cx, cy in centers:
        plt.scatter(cx, cy)
        plt.text(cx+0.3, cy+0.3, name, fontsize=8)
    plt.title("Centers of Digital Stains")
    plt.xlabel("C+G % (center)")
    plt.ylabel("Kappa IC (center)")
    plt.grid(True)
    plt.show()

root = tk.Tk()
root.title("Influenza Digital Stain Generator")
root.geometry("600x400")

label = tk.Label(root, text="Select up to 10 FASTA genomes:", font=("Arial", 12))
label.pack(pady=10)

file_list = tk.Listbox(root, width=80, height=10)
file_list.pack(pady=5)

btn_select = tk.Button(root, text="Choose FASTA Files", command=select_files, width=20)
btn_select.pack(pady=10)

btn_run = tk.Button(root, text="Generate Digital Stains", command=run_analysis, width=25)
btn_run.pack(pady=20)

root.mainloop()
