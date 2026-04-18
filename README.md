[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/gxKL4zJk)
# TP3 - Aligneur de séquences

## Description

Ce projet implémente un aligneur de séquences génomiques en C++.

Il permet de :

* Lire un génome de référence
* Lire des séquences (reads)
* Trouver les meilleures correspondances dans le génome
* Retourner les résultats d’alignement

## Structure du projet

```
src/
├── AlignmentResult
├── BaseLUT
├── BestResultsHeap
├── GenomicPosition
├── HashTable
├── LinkedList
├── Node
├── main.cpp
```

## Compilation

### Avec Makefile

```
make
```

### Avec CMake

```
mkdir build
cd build
cmake ..
make
```

## Exécution

```
./genome_mapper data/chr21_small.fa data/chr21_small_input_reads_100.fa > output.txt
```

## Auteur

Jean Awé
Université Laval

