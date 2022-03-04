# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 10:53:44 2022

@author: Dmytro Olyva
"""

# Part 1
def s(dna):
    dna_count={'A':0,'C':0,'G':0,'T':0}
    for i in range (0,len(dna)):
        if dna[i]=='A':
            dna_count['A']+=1
        elif dna[i]=='C':
            dna_count['C']+=1
        elif dna[i]=='G':
            dna_count['G']+=1
        elif dna[i]=='T':
            dna_count['T']+=1
    return dna_count


# Part 2
def dna2rna(dna):
    rna=''
    for i in range (0,len(dna)):
        if dna[i]=='A':
            rna+='A'
        elif dna[i]=='C':
            rna+='C'
        elif dna[i]=='G':
            rna+='G'
        elif dna[i]=='T':
            rna+='U'
    return rna

#Part 3
def reverse_compliment(dna):
    dna=dna[::-1]
    rna=''
    for i in range (0,len(dna)):
        if dna[i]=='A':
            rna+='T'
        elif dna[i]=='C':
            rna+='G'
        elif dna[i]=='G':
            rna+='C'
        elif dna[i]=='T':
            rna+='A'
    return rna

"""
Created on Thurs Mar 3 19:59:27 2022

@author: Zach Bates
"""
# Part 4
def mendals_law(hom, het, rec):
    x = hom
    y = het
    z = rec
    total = x+y+z
    twoRecess = (z/total)*((z-1)/(total-1))
    twoHetero=(y/total)*((y-1)/(total-1))
    heteroRecess=(z/total)*(y/(total-1))+(y/total)*(z/(total-1))
    recessProb = twoRecess + twoHetero*1/4 + heteroRecess*1/2
    print(1-recessProb)
    
# Part 5
def fibonacci_rabbits(n,k):
    f1, f2=1, 1
    for i in range(n-1):
        f2,f1=f1,f1+(f2*k)
    print(f2)
    return f2
fibonacci_rabbits(5,3)

# Part 6
def gc_content(dna_list):
    max_gc = 0.0
    index = -1
    for i in range(len(dna_list)):
        count = 0
        for j in range(len(dna_list[i])):
            if dna_list[i][j] == 'G' or dna_list[i][j] == 'C':
                count += 1
        gc_con = ((count * 1.0) / len(dna_list[i])) * 100.0
        if gc_con > max_gc:
            index = i
            max_gc = gc_con
    return index, round(max_gc, 6)


"""
Created on Mon Feb 28 21:00:54 2022

@author: markv_mv4bosn
"""

#Part 7: RNA to Codons
def translate(rna):
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',

        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',

        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', }
        
    return genetic_code.get(rna.upper(),'Invalid')     

def rna2codon(rna):
    amino=''
    for i in range(0,int(len(rna)/3)):
        amino += translate(rna[3*i:3*i+3])
    return amino 
        

#Part 8: Locate Substring
def locate_substring(dna_snippet,dna):
    location=[]
    for i in range(len(dna)):
        if dna.startswith(dna_snippet,i):
            location.append(i)
            i+=1
    return location

    
#Part 9: Hamming Distance
def hamming_dist(dna1,dna2):
    i=0
    distance=0
    while i in range(len(dna1)):
        if dna1[i] != dna2[i]:
            distance += 1
        i +=1
    return distance

"""
Created on Sat Feb 26 15:24:44 2022

@author: Matt Schmitt
"""
#Part 10
def count_dom_phenotype(genotypes):
    counter = 0
    dom_child_count = 0
    while counter < 6:
        if counter == 0 or counter == 1 or counter == 2:
            dom_child_count += 2 * float(genotypes[counter])
        elif counter == 3:
            dom_child_count += 1.5 * float(genotypes[counter])
        elif counter == 4:
            dom_child_count += 1 * float(genotypes[counter])
        else:
            dom_child_count += 0
        counter += 1
    return dom_child_count

#Part 11
def source_rna(protein):
    comboes = 1
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',

        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',

        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    for letter in protein:
        counter = 0
        for triplet in genetic_code:
            if genetic_code[triplet] == letter:
                counter += 1
        comboes = comboes * counter
    comboes = comboes * 3 #always 3 stop codons
    return comboes


#Part 12
def splice_rna(dna, intron_list):
    for intron in intron_list:
        if intron in dna:
            dna = dna.replace(intron, '')
    rna_string = dna2rna(dna)
    protein = rna2codon(rna_string)
    return protein
