#!/bin/bash

filename=$1

echo "Procesando secuencia en el archivo: $filename"

tail -n +2 $filename | tr "ATGC" "TACG" > $filename-comp2

echo "El archivo con la secuencia complementaria resultante es: $filename-comp2"
