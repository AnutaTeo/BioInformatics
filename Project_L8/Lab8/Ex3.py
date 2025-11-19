#Ex3 Download from NCBI 3 bacterial genoms of your choosing. 
#try to find in these genoms possible transposons. for this one must detect possible inverted repeats without prior knowledge about their existance in the sequence.
#the inverted repeat should have a min length of 4 letters and a max of 6 letters.

import tkinter as tk
from tkinter import filedialog, messagebox

def reverse_complement(seq):
    comp = str.maketrans("ATGC", "TACG")
    return seq.translate(comp)[::-1]

def load_fasta(path):
    seq = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)

def find_inverted_repeats(genome, min_len=4, max_len=6):
    results = []
    n = len(genome)

    for length in range(min_len, max_len + 1):
        for i in range(n - 2 * length):
            left = genome[i:i+length]
            right = genome[i+length:i+2*length]

            if left == reverse_complement(right):
                results.append({
                    "length": length,
                    "repeat": left,
                    "left_start": i,
                    "left_end": i + length - 1,
                    "right_start": i + length,
                    "right_end": i + 2*length - 1
                })
    return results

def open_file():
    filepath = filedialog.askopenfilename(
        title="Select Genome FASTA File",
        filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")]
    )

    if not filepath:
        return

    messagebox.showinfo("Loading", f"Loading genome:\n{filepath}")

    try:
        genome = load_fasta(filepath)
    except Exception as e:
        messagebox.showerror("Error", f"Could not read file:\n{e}")
        return

    print("Genome loaded. Length:", len(genome), "bp")
    print("Detecting inverted repeats (4–6 bp)...")

    hits = find_inverted_repeats(genome)

    print("Found:", len(hits), "inverted repeats")

    out_file = filepath + "_inverted_repeats.txt"
    with open(out_file, "w") as f:
        for h in hits:
            f.write(
                f"{h['length']}\t{h['repeat']}\t"
                f"{h['left_start']}\t{h['left_end']}\t"
                f"{h['right_start']}\t{h['right_end']}\n"
            )

    messagebox.showinfo(
        "Done",
        f"Found {len(hits)} inverted repeats.\n\n"
        f"Results saved to:\n{out_file}"
    )

root = tk.Tk()
root.title("Transposon Finder – Inverted Repeat Detector")
root.geometry("400x200")

label = tk.Label(root, text="Select a bacterial genome FASTA file", font=("Arial", 12))
label.pack(pady=20)

btn = tk.Button(root, text="Choose File", command=open_file, font=("Arial", 12))
btn.pack(pady=20)

root.mainloop()
