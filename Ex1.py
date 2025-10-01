#make an application that is able to find the alphabet of a sequence of text. This seq may be an ARN seq or ADN seq or protein seq

s = "ACGGGCATATGCGC"
n = len(s)
alphabet = set(s)

for letter in alphabet:
    print(letter)