# -*- coding: utf-8 -*-
"""
Created on Sun Nov  22 

@author: eboza
"""

import sys, getopt, os

### La función obtenerRNA recibe una cadena mRNA en la cual buscará los ORF.
### La función también recibe los parámetros "flagEucariota" que indica si la 
### cadena recibida es de una eucariota y de esta manera buscar los codones de 
### inicio adecuados.  También recibe la bandera "flagReversed" que indica si 
### la cadena corresponde a la reversa de la cadena complementaria, en cuyo 
### caso, los números de los marcos de lectura seran negativos.
def obtenerORF(mRNA, flagEucariota, flagReversed):
    listaAminoacidos = {
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

    ### Si es la cadena reversa complementaria, los identificadores de los 
    ### marcos son negarivos.
    signal = ''
    if (flagReversed > 0):
        signal = '-'
    

    #### Desactivo la comparación del flagEucariota y defino por defecto 3 marcos de lectura, 
    #### Con esto puedo hacer el mismo ejercicio que en la página:
    #### http://vlab.amrita.edu/?sub=3&brch=273&sim=1432&cnt=1
    '''
    ### Eucariota 1 marco
    ### Procariota 3 marcos.
    if (flagEucariota < 1):
        marcos = range(3)
    else:
        marcos = range(1)
    '''
    marcos = range(1)

    proteinas = ''
    for marco in marcos:   ### Se repite el proceso por cada marco.
        indice = marco
        proteinas = proteinas + "\nMarco de Lectura # " + signal + str(marco+1)
        orf = ''
        while (indice < len(mRNA)-2):  
            codon = mRNA[indice:indice+3]   ### Se toman los codones 
            if (codon in listaAminoacidos):
                ### Se verifica los codones de inicio específicos para eucariotas y procariotas.
                if ((codon in {'GUG', 'UUG'} and flagEucariota < 1) or (codon in {'AUG'} and flagEucariota > 0)):
                    orf = '\n' + listaAminoacidos[codon] + ' '
                else:
                   if (orf != ''):   ## Si el orf está vacío significa que no ha encontrado un codón de inicio.
                      orf += listaAminoacidos[codon] + ' '  ## Se agrega el aminoacio a la lista.
                      ## Si se encuentra un codón de parada, se agrega la secuencia al listado final.
                      if (listaAminoacidos[codon] == 'Stop'): 
                          if (len(orf) > 10):  ## Para evitar imprimir Met Stop sin aminoacidos dentro.
                              proteinas += orf 
                          orf = ''
                indice += 3
            else:
                indice += 1

    return proteinas.rstrip();
    
### Función principal.
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

    ### Setear flagEcuariota = 0, significa que es Procariota.
    flagEucariota = 1
    char2 = ""
    
    ### Primera Parte: Leer la secuencia del archivo y obtener la transcripcion
    ### RNA. Almacenarla en memoria en un script largo.
    with open(archivoSecuencia, 'r') as infile:
        mRNA = ''
        mRNA_complementary_reversed = ''
        for line in infile:
            ### De acuerdo al formato FASTA, si la línea empieza con '>'
            ### corresponde a un identificador y comentarios descriptivos
            ### En este caso solo le agregamos la cadena "- Proteinas".
            if line.startswith('>'):
                if mRNA != '':
                    ### Obtengo los ORF de la cadena original y de la cadena 
                    ### Complementaria reversa.
                    print obtenerORF(mRNA, flagEucariota, 0)
                    print obtenerORF(mRNA_complementary_reversed, flagEucariota, 1)
                    mRNA = ''
                print line.rstrip() + ' - Proteinas'
            else:
                ### Para obtener la transcripción en RNA (ARN) se intercambian
                ### los caracteres 'T -> U' 
                ### al mismo tiempo que obtengo el mRNA de la cadena original,
                ### también voy obteniendo la cadena complementaria reversa.
                for char in line.rstrip():
                    if char == 'T':
                        mRNA += 'U'
                    else:
                        mRNA += char
                    ### Armo la cadena mRNA de la complementaria de la original.
                    if char == 'A':
                        char2 = 'U'
                    elif char == 'T':
                        char2 = 'A'
                    elif char == 'C':
                        char2 = 'G'
                    elif char == 'G':
                        char2 = 'C'
                    mRNA_complementary_reversed = char2 + mRNA_complementary_reversed
                    
                    ### Obtengo los ORF de la cadena original y de la cadena 
                    ### Complementaria reversa.
    ##print mRNA
    print obtenerORF(mRNA, flagEucariota, 0)
    print obtenerORF(mRNA_complementary_reversed, flagEucariota, 1)
    return;
    
if __name__ == "__main__":
   main(sys.argv[1:])