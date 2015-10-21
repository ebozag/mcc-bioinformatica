#!/bin/bash

filename=$1

echo "Procesando secuencia en el archivo: $filename"

sed 's/A/B/g' $filename > $filename-comp
sed -i 's/T/A/g' $filename-comp
sed -i 's/B/T/g' $filename-comp
sed -i 's/C/D/g' $filename-comp
sed -i 's/G/C/g' $filename-comp
sed -i 's/D/G/g' $filename-comp 

echo "El archivo con la secuencia complementaria resultante es: $filename-comp"
