import numpy as np
import time
import random
import matplotlib.pyplot as plt
import math
from math import gcd
from Bio.Align.substitution_matrices import Array

### NOTE: Feng-Doolitle
def levenshteinDistance(seqs, dicc):
  lista = list()
  mn = [len(seq) for seq in seqs].index(min([len(seq) for seq in seqs]))
  t = seqs.pop(mn)
  s = seqs[0]

  n = len(s)
  m = len(t)

  previa = [j*dicc['GAP'] for j in range(m+1)]
  actual = [0 for j in range(m+1)]

  for i in range(1, n+1):
    actual[0] = i*dicc['GAP']
    for j in range(1, m+1):
      if s[i-1] == t[j-1]:
        actual[j] = previa[j-1] + dicc['MATCH']
      else:
        diagonal = dicc['MISMATCH'] + previa[j-1]
        izq = dicc['GAP'] + previa[j]
        arriba = dicc['GAP'] + actual[j-1]
        actual[j] = min(izq, arriba, diagonal)
    previa = actual.copy()

  return actual[m]

def inicializacion_NWS(n_rows,n_cols,gap):
  matriz = np.full([n_rows, n_cols], 0)
  for i in range(1, n_rows):
    matriz[i,0] = matriz[i-1, 0] + gap
  for j in range(1, n_cols):
    matriz[0,j] = matriz[0, j-1] + gap
  return matriz

def rellenado_NWS(matriz, A, B, n_rows, n_cols, gap, mismatch):
  for i in range(1, n_rows):
    for j in range(1, n_cols):
      izquierda = matriz[i,j-1] + gap
      arriba = matriz[i-1, j] + gap
      if A[j-1] == B[i-1] or A[j-1] == 'X' or B[i-1] == 'X':
        diagonal = matriz[i-1, j-1]
      else:
        diagonal = matriz[i-1, j-1] + mismatch

      matriz[i,j] = min([arriba, izquierda, diagonal])

  return matriz

def vuelta_atras_NWS(matriz, i, j, A, B, gap):
  alin_A = str()
  alin_B = str()

  while i != 0 or j != 0:

    if matriz[i,j] == matriz[i,j-1] + gap:      #desplazamiento en horizontal
      alin_A = A[j-1] + alin_A
      alin_B = "-" + alin_B
      j = j - 1
    elif matriz[i,j] == matriz[i-1, j] + gap:   #desplazamiento en vertical
      alin_A = "-" + alin_A
      alin_B = B[i-1] + alin_B
      i = i - 1
    else:                                       #desplazamiento en diagonal
      alin_A = A[j-1] + alin_A
      alin_B = B[i-1] + alin_B
      i = i - 1
      j = j - 1

  return alin_A, alin_B

def NWSFengDoolittle(A, B, gap, mismatch):
  n_cols = len(A) + 1
  n_rows = len(B) + 1
  m = inicializacion_NWS(n_rows, n_cols, gap)
  m = rellenado_NWS(m,A,B,n_rows, n_cols, gap, mismatch)
  i = n_rows - 1
  j = n_cols - 1
  dist = m[i,j]
  a, b = vuelta_atras_NWS(m,i,j,A, B, gap)
  return dist,a,b

### NOTE: Fitch&Margoliash algorithm
def min_val(M1):
    filas, columnas = M1.shape
    min_val = math.inf
    minx, miny = 0, 0
    for i in range(0, filas):
        for j in range(i+1, columnas):
            if M1[i, j] < min_val:
                min_val = M1[i, j]
                minx, miny = i, j
    return minx, miny

def fusiona_etiquetas(seqs, i, j):
  if i > j:
    i,j = j,i
  seqs[i] = '(' + seqs[i] + ',' + seqs[j] + ')'
  del seqs[j]

def calcula_distancias(M1, seqs, a, b):
  filas, columnas = M1.shape
  sum_i = 0
  sum_j = 0
  count = 0
  for j in range(columnas):
    if j == b or j == a:
      continue
    sum_i += M1[a, j]
    sum_j += M1[b, j]
    count += 1

  avg_i = sum_i/count
  avg_j = sum_j/count

  d_i = (M1[a,b] + (avg_i - avg_j))/2
  d_j = M1[a,b] - d_i

  seqs[a] += ':'+str(d_i)
  seqs[b] += ':'+str(d_j)

def fusiona_matriz(M1, a, b):
  if a > b:
    a,b = b,a
  M2 = M1.copy()
  M2 = np.delete(M2, b, axis=0)
  M2 = np.delete(M2, b, axis=1)

  r = M2.shape[0]
  for j2 in range(r):
    if j2 >= b:
      j1 = j2 + 1
    else:
      j1 = j2

    if a == j2:
      continue

    M2[a,j2] = (M1[a,j1] + M1[b,j1] - M1[a,b])/2
    M2[j2,a] = M2[a,j2]

  return M2

def fusion_final(seqs, M1):
  d_0 = (M1[0,2]+M1[0,1]-M1[1,2])/2
  d_1 = (M1[1,2]+M1[0,1]-M1[0,2])/2
  d_2 = (M1[0,2]+M1[1,2]-M1[0,1])/2
  seqs[0] += ':'+str(d_0)
  seqs[1] += ':'+str(d_1)
  seqs[2] += ':'+str(d_2)
  return seqs

def FitchMargoliash(M2, seqs):
  matriz = M2

  while len(seqs) > 3:
    x,y = min_val(matriz)
    calcula_distancias(matriz, seqs, x, y)
    fusiona_etiquetas(seqs, x, y)
    matriz = fusiona_matriz(matriz, x, y)
  l = fusion_final(seqs, matriz)
  l = ','.join(seqs)
  l = '('+l+')'
  return l

### NOTE: UPGMA
def min_val(M):
    filas, columnas = M.shape
    min_val = math.inf
    minx, miny = 0, 0
    for i in range(0, filas):
        for j in range(i+1, columnas):
            if M[i, j] < min_val:
                min_val = M[i, j]
                minx, miny = i, j
    return minx, miny

def fusiona_etiquetas(seqs, indices, alturas, i, j, M):
  dist = M[i,j]/2
  alt_i = alturas[i]
  alt_j = alturas[j]
  alturas[i] = dist
  seqs[i] = '(' + seqs[i] +':'+str(dist - alt_i) +',' + seqs[j]+':'+str(dist-alt_j) + ')'
  lista = list()
  for index in indices[i]:
    lista.append(index)
  for index in indices[j]:
    lista.append(index)

  indices[i] = lista
  del seqs[j]
  del alturas[j]
  del indices[j]


def fusiona_matriz(M, M1, a, b, indices):
  M2 = M1.copy()
  M2 = np.delete(M2, b,axis=0)
  M2 = np.delete(M2, b,axis=1)
  filas, columnas = M2.shape
  for j in range(columnas):
    if a == j:
      continue
    cont = 0
    suma = 0
    for ind1 in indices[a]:
      for ind2 in indices[j]:
        cont += 1
        suma += M[ind1, ind2]
    M2[a,j] = suma/cont
    M2[j,a] = M2[a,j]

  return M2

def UPGMA(M2, seqs, indices, alturas):
  matriz = M2
  matriz_original = M2
  while len(seqs) > 1:
    #print(matriz)
    x,y = min_val(matriz)
    fusiona_etiquetas(seqs, indices, alturas, x, y, matriz)
    matriz = fusiona_matriz(matriz_original, matriz, x, y, indices)
  return seqs