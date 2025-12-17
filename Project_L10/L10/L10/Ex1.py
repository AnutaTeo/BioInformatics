import matplotlib.pyplot as plt

S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
S = S.upper()

window_size = 30

def cg_percent(win):
    count = 0
    for b in win:
        if b == "C" or b == "G":
            count += 1
    return (count / len(win)) * 100

def kappa_ic(win):
    n = len(win)
    total = 0

    for shift in range(1, n):
        matches = 0
        length = n - shift
        for i in range(length):
            if win[i] == win[i + shift]:
                matches += 1
        total += (matches / length) * 100

    return total / (n - 1)

cg_values = []
ic_values = []

for start in range(0, len(S) - window_size + 1):
    window = S[start : start + window_size]
    cg_values.append(cg_percent(window))
    ic_values.append(kappa_ic(window))


cg_whole = cg_percent(S)
print("CG% (whole sequence) =", round(cg_whole, 2)) 

ic_whole = kappa_ic(S)
print("Kappa IC (whole sequence) =", round(ic_whole, 2))

plt.scatter(cg_values, ic_values)
plt.xlabel("C+G %")
plt.ylabel("Kappa IC")
plt.title("DNA Pattern")
plt.show()


center_x = sum(cg_values) / len(cg_values)
center_y = sum(ic_values) / len(ic_values)

print("Center of weight:", round(center_x, 2), round(center_y, 2))

plt.scatter([center_x], [center_y], color="red")
plt.xlabel("C+G %")
plt.ylabel("Kappa IC")
plt.title("Center of Weight")
plt.show()

promoter = input("Enter a promoter DNA sequence (A/C/G/T only): ").upper()

if promoter.strip() != "":
    cg_vals_prom = []
    ic_vals_prom = []

    for start in range(0, len(promoter) - window_size + 1):
        win = promoter[start : start + window_size]
        cg_vals_prom.append(cg_percent(win))
        ic_vals_prom.append(kappa_ic(win))

    plt.scatter(cg_vals_prom, ic_vals_prom)
    plt.xlabel("C+G %")
    plt.ylabel("Kappa IC")
    plt.title("Pattern of Your Promoter")
    plt.show()

    cx = sum(cg_vals_prom) / len(cg_vals_prom)
    cy = sum(ic_vals_prom) / len(ic_vals_prom)
    print("Center of promoter pattern:", round(cx, 2), round(cy, 2))

    plt.scatter([cx], [cy], color="red")
    plt.xlabel("C+G %")
    plt.ylabel("Kappa IC")
    plt.title("Promoter Center of Weight")
    plt.show()

    print("\nNow open PromKappa and compare your pattern.")
