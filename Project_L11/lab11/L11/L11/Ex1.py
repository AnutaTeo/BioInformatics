import tkinter as tk
from tkinter import ttk, scrolledtext

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
    n = len(seq1)
    m = len(seq2)

    score = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        score[i][0] = i * gap
    for j in range(1, m + 1):
        score[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            delete = score[i - 1][j] + gap
            insert = score[i][j - 1] + gap
            score[i][j] = max(diag, delete, insert)

    aligned1 = ""
    aligned2 = ""
    path = []

    i, j = n, m
    while i > 0 or j > 0:
        path.append((i, j))
        current = score[i][j]

        diag = score[i - 1][j - 1] if i > 0 and j > 0 else None
        up = score[i - 1][j] if i > 0 else None
        left = score[i][j - 1] if j > 0 else None

        if i > 0 and j > 0 and current == diag + (match if seq1[i - 1] == seq2[j - 1] else mismatch):
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = seq2[j - 1] + aligned2
            i -= 1
            j -= 1
        elif i > 0 and current == up + gap:
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1
        else:
            aligned1 = "-" + aligned1
            aligned2 = seq2[j - 1] + aligned2
            j -= 1

    path.append((0, 0))
    return score, aligned1, aligned2, path

def run_alignment():
    seq1 = entry_seq1.get().strip()
    seq2 = entry_seq2.get().strip()

    if not seq1 or not seq2:
        output_box.insert("end", "Please enter both sequences.\n")
        return

    try:
        match = int(entry_match.get())
        mismatch = int(entry_mismatch.get())
        gap = int(entry_gap.get())
    except ValueError:
        output_box.insert("end", "Error: scoring values must be integers.\n")
        return

    matrix, a1, a2, path = needleman_wunsch(seq1, seq2, match, mismatch, gap)

    matches = sum(1 for i in range(len(a1)) if a1[i] == a2[i])
    similarity = matches / len(a1) * 100
    match_line = "".join("|" if a1[i] == a2[i] else " " for i in range(len(a1)))

    output_box.delete("1.0", "end")
    output_box.insert("end", "Show Alignment:\n")
    output_box.insert("end", a1 + "\n")
    output_box.insert("end", match_line + "\n")
    output_box.insert("end", a2 + "\n\n")
    output_box.insert("end", f"Matches = {matches}\n")
    output_box.insert("end", f"Length  = {len(a1)}\n")
    output_box.insert("end", f"Similarity = {similarity:.0f} %\n")

    update_matrix_plot(matrix)
    update_traceback_plot(matrix, path)


def update_matrix_plot(matrix):
    fig_matrix.clear()
    ax = fig_matrix.add_subplot(111)

    im = ax.imshow(matrix, cmap="magma")
    ax.set_title("Score matrix (heatmap)")
    ax.set_xlabel("Seq2")
    ax.set_ylabel("Seq1")
    fig_matrix.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    canvas_matrix.draw()


def update_traceback_plot(matrix, path):
    rows = len(matrix)
    cols = len(matrix[0])

    path_mat = [[0] * cols for _ in range(rows)]
    for (i, j) in path:
        path_mat[i][j] = 1

    fig_trace.clear()
    ax = fig_trace.add_subplot(111)

    im = ax.imshow(path_mat, cmap="Reds")
    ax.set_title("Traceback path")

    ax.set_xticks(range(cols))
    ax.set_yticks(range(rows))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.grid(True, which="both", linewidth=0.2, color="black")

    canvas_trace.draw()

root = tk.Tk()
root.title("DNA Alignment - Needleman–Wunsch")


left_frame = ttk.Frame(root, padding=10)
left_frame.grid(row=0, column=0, sticky="nsw")

seq_frame = ttk.LabelFrame(left_frame, text="Sequences", padding=5)
seq_frame.grid(row=0, column=0, sticky="ew", pady=5)

ttk.Label(seq_frame, text="Sq 1:").grid(row=0, column=0, sticky="w")
entry_seq1 = ttk.Entry(seq_frame, width=25)
entry_seq1.grid(row=0, column=1, padx=5, pady=2)
entry_seq1.insert(0, "ACCGTGAAGCCAATAC")

ttk.Label(seq_frame, text="Sq 2:").grid(row=1, column=0, sticky="w")
entry_seq2 = ttk.Entry(seq_frame, width=25)
entry_seq2.grid(row=1, column=1, padx=5, pady=2)
entry_seq2.insert(0, "AGCGTGCAGCCAATAC")

param_frame = ttk.LabelFrame(left_frame, text="Parameters", padding=5)
param_frame.grid(row=1, column=0, sticky="ew", pady=5)

ttk.Label(param_frame, text="Gap =").grid(row=0, column=0, sticky="w")
entry_gap = ttk.Entry(param_frame, width=5)
entry_gap.grid(row=0, column=1, padx=3, pady=2)
entry_gap.insert(0, "-1")

ttk.Label(param_frame, text="Match =").grid(row=1, column=0, sticky="w")
entry_match = ttk.Entry(param_frame, width=5)
entry_match.grid(row=1, column=1, padx=3, pady=2)
entry_match.insert(0, "1")

ttk.Label(param_frame, text="Mismatch =").grid(row=2, column=0, sticky="w")
entry_mismatch = ttk.Entry(param_frame, width=5)
entry_mismatch.grid(row=2, column=1, padx=3, pady=2)
entry_mismatch.insert(0, "-1")

btn_align = ttk.Button(left_frame, text="Align", command=run_alignment)
btn_align.grid(row=2, column=0, pady=10, sticky="ew")

right_top = ttk.Frame(root, padding=10)
right_top.grid(row=0, column=1, sticky="nsew")

root.columnconfigure(1, weight=1)
root.rowconfigure(1, weight=1)

frame_matrix = ttk.LabelFrame(right_top, text="Graphic representation of the alignment matrix")
frame_matrix.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

frame_trace = ttk.LabelFrame(right_top, text="Traceback path deviation from optimal alignment")
frame_trace.grid(row=0, column=1, padx=5, pady=5, sticky="nsew")

right_top.columnconfigure(0, weight=1)
right_top.columnconfigure(1, weight=1)

fig_matrix = Figure(figsize=(4, 4))
canvas_matrix = FigureCanvasTkAgg(fig_matrix, master=frame_matrix)
canvas_matrix.get_tk_widget().pack(fill="both", expand=True)

fig_trace = Figure(figsize=(4, 4))
canvas_trace = FigureCanvasTkAgg(fig_trace, master=frame_trace)
canvas_trace.get_tk_widget().pack(fill="both", expand=True)

bottom_frame = ttk.LabelFrame(root, text="Show Alignment", padding=5)
bottom_frame.grid(row=1, column=0, columnspan=2, padx=10, pady=5, sticky="nsew")

output_box = scrolledtext.ScrolledText(bottom_frame, width=80, height=10)
output_box.pack(fill="both", expand=True)

root.mainloop()
