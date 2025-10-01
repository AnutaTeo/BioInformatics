#use the ai to adapt your current algorithm in order to make an app that takes a FASTA file and read the seq content from it and display the rel. percentages for the
#symbols present in the alphabet of seq. Note: FASTA represents a file format that contains DNA, ARN or proteins seq. Thus, it contains the information for your input

with open(r"C:/Users/amamt/Desktop/BioInf/Lab1/sequence.fasta", "r") as f:
    lines = f.readlines()

s = "".join(line.strip() for line in lines if not line.startswith(">")).upper()

n = len(s)

alphabet = set(s)

for letter in alphabet:
    nr = s.count(letter)
    percentage = (nr / n) * 100
    print(letter, ": ", percentage)

