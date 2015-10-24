#!/bin/bash

filename=$1

echo "Procesando secuencia en el archivo: $filename"

grep -v ">" $filename | tr "ATGC" "TACG" > $filename-comp2

echo "El archivo con la secuencia complementaria resultante es: $filename-comp2"
