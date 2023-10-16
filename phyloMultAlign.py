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

#NOTE: Suma de los pares
def sPares(secuencias:dict, matrix:Array, gap:int):
  ids = [k for k,v in secuencias.items()]
  p = 0
  for i in range(0, len(ids)):
    for j in range(i+1, len(ids)):
      for n in range(len(secuencias[ids[i]])):
        n1 = secuencias[ids[i]][n]
        n2 = secuencias[ids[j]][n]
        if n1 == '-' and n2 == '-':
          p += 0
        elif (n1 == '-' and n2 != n1) or (n2 == '-' and n2 != n1):
          p += gap
        else:
          p += matrix[n1][n2]
  return p

def sumaPares(secuencias:list, matrix:Array):
  #ids = [k for k,v in secuencias.items()]
  p = 0
  for i in range(0, len(secuencias)):
    for j in range(i+1, len(secuencias)):
      for n in range(len(secuencias[secuencias[i]])):
        n1 = secuencias[secuencias[i]][n]
        n2 = secuencias[secuencias[j]][n]
        p += matrix[n1][n2]
  return p


#NOTE: Alineamiento a un perfil
from collections import Counter
up_arrow = "\u2191"
left_arrow = "\u2190"
up_left_arrow = "\u2196"

def makeProfile(secuencias:dict):
  '''
  Entrada: 
    secuencias: Diccionario de secuencias (id:seq)
  Salida: perfil, lista por posición de diccionarios que contienen la frecuencia para cada caracter
  '''
  max_length = max(len(seq) for seq in secuencias.values())
  array = np.empty((1, len(secuencias), max_length), dtype='U1')
  for i, seq in enumerate(secuencias.values()):
    array[0, i, :len(seq)] = list(seq)
  sequence_list = array.tolist()
  position_frequencies = []

  total_sequences = len(secuencias)
  for i in range(max_length):
      position = [seq[i] for seq in sequence_list[0]]
      position_counts = {
          'A': position.count('A') / total_sequences,
          'C': position.count('C') / total_sequences,
          'G': position.count('G') / total_sequences,
          'T': position.count('T') / total_sequences,
          '-': position.count('-') / total_sequences,
      }
      position_frequencies.append(position_counts)

  return position_frequencies

def seqToProfile(p:list, alignment:list, s:str, pen:dict, mat:Array):
  '''
  Entrada:
    p: perfil
    alignment: alineamiento multiple que da lugar al perfil
    s: secuencia a alinear
    pen: diccionario de penalizaciones
    mat: matriz de sustitucion
  Salida:
    new_alignment: lista con el nuevo alineamiento
    alin_B: nueva secuecia alineada
    matrizS: matriz de puntuación del alineamiento
    matrizT: matriz de la vuelta atrás
  '''
  n_rows = len(s) + 1
  n_cols = len(p) + 1
  matrizS, matrizT = initStoP(n_rows,n_cols,pen['GAP'],p)
  matrizS, matrizT = fillitStoP(matrizS, matrizT, s, p, n_rows, n_cols, pen, mat)
  new_alignment, alin_B = traceBackStoP(matrizT, alignment, s, n_rows, n_cols)
  return new_alignment, alin_B, matrizS, matrizT


def initStoP(n_rows,n_cols,gap,p):
  '''
  Entrada:
    n_rows: longitud de la secuencia + 1
    n_cols: longitud del alineamiento + 1
    gap: penalización por hueco
    p: perfil del alineamiento
  Devuelve:
    matrizS: matriz de puntuación del alineamiento inicializadas
    matrizT: matriz de la vuelta atrás inicializadas
  '''
  matrizS = np.full([n_rows, n_cols], 0, dtype=float)
  matrizT = np.full([n_rows, n_cols], '-', dtype=str)
  global up_arrow
  global left_arrow
  for i in range(1, n_rows):
    matrizS[i,0] = matrizS[i-1, 0] + gap
    matrizT[i,0] = up_arrow
  for j in range(1, n_cols):
    # Inicializar la primera fila consiste en restar la proporción de hueco a la proporción de letra 
    # en esa posición y multiplicarlo por el valor del hueco.
    matrizS[0,j] = matrizS[0, j-1] + gap*(1-p[j-1]['-'])
    matrizT[0,j] = left_arrow
  return matrizS, matrizT

def fillitStoP(matrizS, matrizT, s, p, n_rows, n_cols,pen, mat):
  '''
  Entrada:
    matrizS: matriz de puntuaciones incializada
    matrizT: matriz de vuelta atrás inicializada
    s: secuencia a alinear
    p: perfil
    n_rows: longitud de la secuencia + 1
    n_cols: longitud del alineamiento + 1
    pen: penalización por hueco
    p: perfil del alineamiento
  Devuelve:
    matrizS: matriz de puntuación del alineamiento completa
    matrizT: matriz de la vuelta atrás completa
  '''
  global up_arrow
  global left_arrow
  global up_left_arrow
  eq = {0:up_arrow,1:left_arrow,2:up_left_arrow}
  for i in range(1, n_rows):
    for j in range(1, n_cols):
      izquierda = matrizS[i,j-1] + pen['GAP']*(1-p[j-1]['-'])
      arriba = matrizS[i-1, j] + pen['GAP']
      b = s[i-1]
      pd = 0
      for a, f in p[j-1].items():
        pd += mat[b,a]*f
      diagonal = matrizS[i-1,j-1] + pd
      score = min([arriba, izquierda, diagonal])

      matrizS[i,j] = score
      matrizT[i,j] = eq[[arriba, izquierda, diagonal].index(score)]

  return matrizS, matrizT

def traceBackStoP(matrizT, alignment, s, n_rows, n_cols):
  '''
  Entrada:
    matrizT: matriz de vuelta atrás completa
    alignment: alineamiento multiple previo
    s: secuencia a alinear
    n_rows: longitud de la secuencia + 1
    n_cols: longitud del alineamiento + 1
  Devuelve:
    new_al: nuevo alineamiento multiple
    alin_B: nueva secuencia alineada
  '''
  i = n_rows - 1
  j = n_cols - 1
  new_al = [str() for _ in alignment]
  alin_B = str()
  global up_arrow
  global left_arrow
  global up_left_arrow
  eq = {0:up_arrow,1:left_arrow,2:up_left_arrow}
  while i != 0 or j != 0:
    if matrizT[i,j] == left_arrow:      #desplazamiento en horizontal
      for k in  range(len(alignment)):
        new_al[k] = alignment[k][j-1] + new_al[k]
      alin_B = '-' + alin_B
      j = j - 1
    elif matrizT[i,j] == up_arrow:   #desplazamiento en vertical
      for k in  range(len(alignment)):
        new_al[k] = '-' + new_al[k]
      alin_B = s[i-1] + alin_B
      i = i - 1
    else:                                       #desplazamiento en diagonal
      for k in  range(len(alignment)):
        new_al[k] =  alignment[k][j-1] + new_al[k]
      alin_B = s[i-1] + alin_B
      i = i - 1
      j = j - 1

  #new_al.append(alin_B)
  return new_al, alin_B