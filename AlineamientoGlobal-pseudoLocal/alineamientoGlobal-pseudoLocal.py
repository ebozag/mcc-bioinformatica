# -*- coding: utf-8 -*-
"""
Created on Sun Jan  10 2016

@author: eboza
"""

import sys, getopt, os
import numpy as np


def leerSecuenciaArchivo(nombreArchivo):
    ### Se asume formato FASTA y una sola secuencia por cada archivo.
    ### Se asume una línea de identificador por cada secuencia, empieza con ">"
    
    ### Inicializo las variables que serán retornadas
    etiqueta = ""    
    secuencia = ""
    
    ### Abro el archivo, leo el contenido y lo almaceno en la variable correspondiente.
    with open(nombreArchivo, 'r') as infile:
        for line in infile:
            ### De acuerdo al formato FASTA, si la línea empieza con '>'
            ### corresponde a un identificador y comentarios descriptivos
            if line.startswith('>'):
                ### Si la variable etiqueta no está vacía, significa que ya se leyó una secuencia
                ### y se termina la función.
                if etiqueta != "":
                    break
                etiqueta = line.rstrip()
            else:
                ### Concatenar la secuencia obtenida
                secuencia += line.rstrip()
    
    ### Retorna la secuencia obtenida del archivo
    return etiqueta, secuencia;


def main(argv):
    ### Poner argumentos en variables
    nombreScript = os.path.basename(__file__)
    archivoSecuencia1 = ''
    archivoSecuencia2 = ''
    try:
        opts, args = getopt.getopt(argv,"hf:g:",["archivoSecuencia1=", "archivoSecuencia2="])
    except getopt.GetoptError:
        print nombreScript +' -f <Archivo con Secuencia Genetica 1> -g <Archivo con Secuencia Genetica 2>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print nombreScript +' -f <Archivo con Secuencia Genetica 1> -g <Archivo con Secuencia Genetica 2>'
            sys.exit()
        elif opt in ("-f", "--archivoSecuencia1"):
            archivoSecuencia1 = arg
        elif opt in ("-g", "--archivoSecuencia2"):
            archivoSecuencia2 = arg
    if archivoSecuencia1 == '' or archivoSecuencia2 == '' :
        print 'ERROR: Debe especificar el nombre del archivo con la secuencia.'
        print nombreScript +' -f <Archivo con Secuencia Genetica 1> -g <Archivo con Secuencia Genetica 2>'
        sys.exit(1)


    ### Obtengo las secuencias de los 2 archivos
    etiquetaSecuencia1, secuencia1_original = leerSecuenciaArchivo(archivoSecuencia1)
    etiquetaSecuencia2, secuencia2_original = leerSecuenciaArchivo(archivoSecuencia2)

    ### Recorte de la parte inicial de la secuencia más larga
    if len(secuencia1_original) >= len(secuencia2_original):
        offset = secuencia1_original.index(secuencia2_original[0])
        secuencia1 = secuencia1_original[offset:]
        secuencia2 = secuencia2_original
    else:
        offset = secuencia2_original.index(secuencia1_original[0])
        secuencia1 = secuencia1_original
        secuencia2 = secuencia2_original[offset:]

    ### Calculo las dimensiones de las secuencias, se usa como dimensiones de los arreglos
    m = len(secuencia1)
    n = len(secuencia2)
    
    ### Inicializo los valores para los puntajes (premio y penalización)
    premioMatch = 1
    penaMismatch = -1
    penaGap = -1
    
    ### Defino tabla dinámica, con una columna y una fila adicional a las de las 
    ### secuencias, para inicializar valores.
    score = np.empty((m+1, n+1), dtype=object)
    path = np.empty((m+1, n+1), dtype=object)
    
    ### Inicializo los valores de la primera fila y primera columna. Corresponde
    ### a los puntajes para cuando una de las cadenas es vacío, mientras mayor el
    ### largo de la segunda cadena, menor el puntaje total.
    for i in range(0,m+1):
        score[i,0] = penaGap * i 
    for j in range(0,n+1):
        score[0,j] = penaGap * j 

    ### Inicializo variable de puntaje máximo en 0 y una variable para sus
    ### coordenadas
    maxScore = 0
    maxScore_i = 0
    maxScore_j = 0
    
    ### Proceso principal de comparación de las cadenas y asignación de puntajes
    for i in range(1,m+1):
        for j in range(1,n+1):
            if secuencia1[i-1] == secuencia2[j-1]:
                score[i,j] = score[i-1,j-1] + premioMatch
                path[i,j] = 1
            else:
                if score[i-1,j-1] >= score[i-1,j]:
                    if score[i-1,j-1] >= score[i,j-1]:
                        score[i,j] = score[i-1,j-1] + penaMismatch
                        path[i,j] = 1
                    else:
                        score[i,j] = score[i,j-1] + penaGap
                        path[i,j] = 2
                elif score[i-1,j] >= score[i,j-1]:
                        score[i,j] = score[i-1,j] + penaGap
                        path[i,j] = 3
                else:
                    score[i,j] = score[i,j-1] + penaGap
                    path[i,j] = 2
                
            ### Verifico si el nuevo puntaje es el mayor y reemplazo la variable.
            if  score[i,j] >= maxScore:
                maxScore = score[i,j]
                maxScore_i = i
                maxScore_j = j
    
    ### Construcción de Cadenas alineadas
    ### Reconstruir la ruta
    i = maxScore_i
    j = maxScore_j
    secuenciaAlineada1 = ""
    secuenciaAlineada2 = ""
    relacionador = ""
    while ((i > 0) and (j > 0)):
        if path[i][j] == 1:
            secuenciaAlineada1 = secuencia1[i-1] + secuenciaAlineada1
            secuenciaAlineada2 = secuencia2[j-1] + secuenciaAlineada2
            if secuencia1[i-1] == secuencia2[j-1]:
                relacionador = "|" + relacionador
            else:
                relacionador = " " + relacionador
            i -= 1
            j -= 1
        elif path[i][j] == 3:
            secuenciaAlineada1 = secuencia1[i-1] + secuenciaAlineada1
            secuenciaAlineada2 = "-" + secuenciaAlineada2
            relacionador = " " + relacionador
            i -= 1                
        elif path[i][j] == 2:
            secuenciaAlineada1 = "-" + secuenciaAlineada1
            secuenciaAlineada2 = secuencia2[j-1] + secuenciaAlineada2
            relacionador = " " + relacionador
            j -= 1


    ### Cálculo del porcentaje de cobertura
    if len(secuencia1_original) >= len(secuencia2_original):
        cobertura = float(len(secuenciaAlineada2)) / float(len(secuencia2_original)) * 100
    else:
        cobertura = float(len(secuenciaAlineada1)) / float(len(secuencia1_original)) * 100
      
    ### Imprimo resultados
    print "#################################################"
    print "### RESULTADOS DE ALINEAMIENTO DE LAS SECUENCIAS:"
    print secuencia1_original
    print secuencia2_original
    print
    print "Score: " + str(maxScore)
    print
    print "### Secuencias alineadas:"
    print secuenciaAlineada1
    print relacionador
    print secuenciaAlineada2
    print
    print "### Porcentaje de cobertura: " + str(cobertura) + "%."
 
if __name__ == "__main__":
   main(sys.argv[1:])