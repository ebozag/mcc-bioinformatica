# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 16:22:45 2015

@author: eboza
"""

import sys, getopt, os

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

    ### Primera Parte: Leer la secuencia del archivo y obtener la complementaria
    archivoSecuenciaComplementaria = archivoSecuencia + '-complementaria'
    archivoSalida = open(archivoSecuenciaComplementaria, 'w')
    with open(archivoSecuencia, 'r') as infile:
        for line in infile:
            ### De acuerdo al formato FASTA, si la línea empieza con '>'
            ### corresponde a un identificador y comentarios descriptivos
            ### En este caso solo le agregamos la cadena "- Secuencia
            ### Complementaria".
            if line.startswith('>'):
                line2 = line.rstrip() + ' - Secuencia Complementaria'
            else:
                ### Para obtener la secuencia complementaria de ADN se intercambian
                ### los caracteres 'A' con 'T' y 'C' con 'G'
                line2 = ''
                for char in line.rstrip():
                    if char == 'A':
                        char2 = 'T'
                    elif char == 'T':
                        char2 = 'A'
                    elif char == 'C':
                        char2 = 'G'
                    elif char == 'G':
                        char2 = 'C'
                    line2 = line2 + char2
            archivoSalida.write(line2+'\n')
    archivoSalida.close()

 
if __name__ == "__main__":
   main(sys.argv[1:])