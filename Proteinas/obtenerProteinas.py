# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 16:22:45 2015

@author: eboza
"""

import sys, getopt, os

def obtenerProteinas( mRNA):
    listaProteinas = {
        'UUU' : 'Phe',
        'UUC' : 'Phe',
        'UUA' : 'Leu',
        'UUG' : 'Leu',
        'CUU' : 'Leu',
        'CUC' : 'Leu',
        'CUA' : 'Leu',
        'CUG' : 'Leu',   
        'AUU' : 'Ile',
        'AUC' : 'Ile',
        'AUA' : 'Ile',
        'AUG' : 'Met',       
        'GUU' : 'Val',
        'GUC' : 'Val',
        'GUA' : 'Val',
        'GUG' : 'Val',        

        'UCU' : 'Ser',
        'UCC' : 'Ser',
        'UCA' : 'Ser',
        'UCG' : 'Ser',
        'CCU' : 'Pro',
        'CCC' : 'Pro',
        'CCA' : 'Pro',
        'CCG' : 'Pro',   
        'ACU' : 'Thr',
        'ACC' : 'Thr',
        'ACA' : 'Thr',
        'ACG' : 'Thr',       
        'GCU' : 'Ala',
        'GCC' : 'Ala',
        'GCA' : 'Ala',
        'GCG' : 'Ala',    
       
        'UAU' : 'Tyr',
        'UAC' : 'Tyr',
        'UAA' : 'Stop',
        'UAG' : 'Stop',
        'CAU' : 'His',
        'CAC' : 'His',
        'CAA' : 'Gln',
        'CAG' : 'Gln',   
        'AAU' : 'Asn',
        'AAC' : 'Asn',
        'AAA' : 'Lys',
        'AAG' : 'Lys',       
        'GAU' : 'Asp',
        'GAC' : 'Asp',
        'GAA' : 'Glu',
        'GAG' : 'Glu',       
       
        'UGU' : 'Cys',
        'UGC' : 'Cys',
        'UGA' : 'Stop',
        'UGG' : 'Trp',
        'CGU' : 'Arg',
        'CGC' : 'Arg',
        'CGA' : 'Arg',
        'CGG' : 'Arg',   
        'AGU' : 'Ser',
        'AGC' : 'Ser',
        'AGA' : 'Arg',
        'AGG' : 'Arg',       
        'GGU' : 'Gly',
        'GGC' : 'Gly',
        'GGA' : 'Gly',
        'GGG' : 'Gly',             

    };
    indice = 0
    proteinas = ''
    while (indice < len(mRNA)-2):
        codon = mRNA[indice:indice+3]
        if ( codon in listaProteinas):
            proteinas += listaProteinas[codon] + ' '
            indice += 3
        else:
            indice += 1
    return proteinas.rstrip();
    

def main(argv):
    ### Poner argumentos en variables
    nombreScript = os.path.basename(__file__)
    archivoSecuencia = ''
    try:
        opts, args = getopt.getopt(argv,"hf:",["archivoSecuencia="])
    except getopt.GetoptError:
        print nombreScript +' -f <Archivo con Secuencia Genetica>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print nombreScript +' -f <Archivo con Secuencia Genetica>'
            sys.exit()
        elif opt in ("-f", "--archivoSecuencia"):
            archivoSecuencia = arg
    if archivoSecuencia == '':
        print 'ERROR: Debe especificar el nombre del archivo con la secuencia.'
        print nombreScript +' -f <Archivo con Secuencia Genetica>'
        sys.exit(1)

    ### Primera Parte: Leer la secuencia del archivo y obtener la transcripcion
    ### RNA. Almacenarla en memoria en un script largo.
    with open(archivoSecuencia, 'r') as infile:
        mRNA = ''
        for line in infile:
            ### De acuerdo al formato FASTA, si la línea empieza con '>'
            ### corresponde a un identificador y comentarios descriptivos
            ### En este caso solo le agregamos la cadena "- Proteinas".
            if line.startswith('>'):
                if mRNA != '':
                    print obtenerProteinas(mRNA)
                    mRNA = ''
                print line.rstrip() + ' - Proteinas'
            else:
                """
                ### Para obtener la transcripción en RNA (ARN) se intercambian
                ### los caracteres 'A -> U', 'T -> A', 'C -> G' y 'G -> C' 
                for char in line.rstrip():
                    if char == 'A':
                        char2 = 'U'
                    elif char == 'T':
                        char2 = 'A'
                    elif char == 'C':
                        char2 = 'G'
                    elif char == 'G':
                        char2 = 'C'
                    mRNA = mRNA + char2
                """
                ### Para obtener la transcripción en RNA (ARN) se intercambian
                ### los caracteres 'T -> U' 
                for char in line.rstrip():
                    if char == 'T':
                        char = 'U'
                    mRNA = mRNA + char
    print obtenerProteinas(mRNA)
    return;
    
if __name__ == "__main__":
   main(sys.argv[1:])