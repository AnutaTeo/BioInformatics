#A DNA seq is given: s=acgggcatatgcgc. make an app which is able to show the percentage of the component from the alphabet of the seq s. 
# In other words, the input of the seq s and the output is the alphabet of the seq and the percentage of each letter in the alphabet found in seq s.

s = "ACGGGCATATGCGC"
#aaba
n = len(s)
alphabet = set(s)
    
for letter in alphabet:
    nr = s.count(letter)
    percentage = (nr / n) * 100
    print(letter, ": ", percentage)

