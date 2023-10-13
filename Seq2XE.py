import numpy as np
import time
import random
import matplotlib.pyplot as plt
import math
from math import gcd
from Bio.Align.substitution_matrices import Array


##Matriz sustitucion nucleotidos
def m_sust_nt(id=0.55, ps=0.25, nt='ACGT'):
  idx = 1-id
  nt_list = list(nt)
  m_sust = Array(nt, dims=2)
  for nt1 in nt_list:
      for nt2 in nt_list:
        if nt1 == nt2:
          m_sust[nt1, nt2] = round(2*math.log2((id/4)/ps**2))
        else:
          m_sust[nt1, nt2] = round(2*math.log2((idx/12)/ps**2))
  a = abs(int(m_sust[0,0]))
  b = abs(int(m_sust[1,0]))
  MCD = gcd(a,b)
  m_sust /= MCD

  return m_sust

def m_sust_DL(nt='ACGT-'):
  nt_list = list(nt)
  m_sust = Array(nt, dims=2)
  for nt1 in nt_list:
      for nt2 in nt_list:
        if nt1 == nt2:
          m_sust[nt1, nt2] = 0
        else:
          m_sust[nt1, nt2] = 1
  return m_sust
#NW Original
def inicializacion_NWO_m(A, B, m, eq):
  matriz = matriz = np.full([len(B),len(A)], 0)
  eq = {'A':0, 'R':1, 'N':2, 'D':3, 'C':4, 'Q':5, 'E':6, 'G':7, 'H':8, 'I':9,'L':10,
          'K':11, 'M':12, 'F':13, 'P':14, 'S':15, 'T':16, 'W':17, 'Y':18, 'V':19}
  for j in range(len(A)):
    for i in range(len(B)):
      a = A[j]
      b = B[i]
      matriz[i,j] = m[eq[a], eq[b]]
  return matriz

def inicializacion_NWO(A, B, match):
  matriz = matriz = np.full([len(B),len(A)], 0)
  for j in range(len(A)):
    for i in range(len(B)):
      if A[j] == B[i]:
        matriz[i,j] = match
      else:
        matriz[i,j] = 0
  return matriz

def rellenado_NWO(matriz, n_cols, n_rows):
  for i in range(n_rows-2,-1,-1):       #las ultimas fila y columna no varian(n_rows - 2,)-> desde, (,-1,)-> hasta el 0, (,-1)-> con paso
    for j in range(n_cols-2,-1,-1):
      max = 0
      for j2 in range(j+1, n_cols):     #la columna siguiente a la actual (j+1,)-> desde, final de la matriz (,n_cols)-> hasta
        if matriz[i+1, j2] > max:
          max = matriz[i+1, j2]

      for i2 in range(i+1, n_rows):
        if matriz[i2, j+1] > max:
          max = matriz[i2, j+1]

      matriz[i,j] = matriz[i,j] + max
  return matriz

def coordenadas_inicio_NWO(n_cols, n_rows, matriz):
  i_max = 0
  j_max = 0

  for i in range(1, n_rows):
    if matriz[i, 0] >= matriz[i_max, j_max]:
      i_max = i
      j_max = 0
  for j in range(1, n_cols):
    if matriz[0, j] >= matriz[i_max, j_max]:
      i_max = 0
      j_max = j

  return i_max, j_max

def vuelta_atras_NWO(matriz, i_max, j_max, A, B, match):
  i3 = i_max
  j3 = j_max
  resultado = list()

  while True:
    i = i3
    j = j3
    resultado.append((i,j))
    if i == len(B) - 1 or j == len(A) - 1:        #condicion de parada, cuando se llega a un extremo de la matriz
      break

    else:
      for i2 in range(i+1, len(B)):       #consulta si en la siguiente fila se encuentra un valor que
        if matriz[i,j] == matriz[i2, j+1] + match and A[j] == B[i]:       #sea igual al anterior más el valor del match
          i3 = i2
          j3 = j + 1
          break

        elif matriz[i,j] == matriz[i2, j+1] and A[j] != B[i]:       #o que sea igual al anterior en un mismatch
          i3 = i2
          j3 = j + 1
          break

      for j2 in range(j+1, len(A)):       #si no se cumple para la siguiente fila, se realiza lo mismo para la siguiente columna
        if matriz[i,j] == matriz[i+1, j2] + match and A[j] == B[i]:
          i3 = i + 1
          j3 = j2
          break

        elif matriz[i,j] == matriz[i+1, j2] and A[j] != B[i]:
          i3 = i + 1
          j3 = j2
          break
  return resultado

def alineamiento_NWO(tuplas, A, B):
  (start_B, start_A) = tuplas[0]
  ali_A = A[start_A]
  ali_B = B[start_B]
  past_A = start_A
  past_B = start_B

  for (i, j) in tuplas[1:]:

     if i - past_B > 1:
       for n in range(i-past_B-1):
         ali_A = ali_A + "-"
         ali_B = ali_B + B[past_B+n+1]

     elif j - past_A > 1:
       for n in range(j-past_A-1):
         ali_B = ali_B + "-"
         ali_A = ali_A + A[past_A+n+1]

     ali_A = ali_A + A[j]
     ali_B = ali_B + B[i]
     past_A = j
     past_B = i

  return ali_A, ali_B

def NWorig(A,B,match):
  m = inicializacion_NWO(A, B, match)
  m = rellenado_NWO(m, len(A), len(B))
  i_max, j_max = coordenadas_inicio_NWO(len(A), len(B), m)
  resultado = vuelta_atras_NWO(m, i_max, j_max, A, B, match)
  a,b = alineamiento_NWO(resultado, A, B)
  return a,b

#NW-Sellers
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
      if A[j-1] == B[i-1]:
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

def NWSellers(A, B, gap, mismatch):
  n_cols = len(A) + 1
  n_rows = len(B) + 1
  m = inicializacion_NWS(n_rows, n_cols, gap)
  m = rellenado_NWS(m,A,B,n_rows, n_cols, gap, mismatch)
  i = n_rows - 1
  j = n_cols - 1
  a, b = vuelta_atras_NWS(m,i,j,A, B, gap)
  return a,b

### GOTOH

def f_gap(gap, ext, x):
    '''
    Parametros de entrada:
    gap: puntuación asociada para abrir un hueco
    ext: puntuación asociada por extender el hueco 1 posición
    x: longitud del hueco

    Devuelve:
    wx: puntuación para un tamaño de hueco x
    '''
    wx = (gap + ext*(x-1))
    return wx

def inicializacion_gotoh(n_rows, n_cols, gap, ext):
    '''
    Inicialización de D, H y V. Incialización de una estructura de 3 matrices, donde:
    D--> matriz[0, i, j] permite movimientos en diagonal desde cualquiera de las 3 matrices
    H--> matriz[1, i, j] permite movimientos en horizontal en H y en diagonal en D
    V--> matriz[2, i, j] permite movimeitnos en vertical en V y en diagonal en D
    '''
    matriz = np.full([3, n_rows, n_cols], float(0))
    #Inicialización del origen de D,H y V
    matriz[0, 0, 0] = 0
    matriz[1, 0, 0] = -math.inf
    matriz[2, 0, 0] = -math.inf

    #Rellenado de las primeras filas y columnas de D, H y V
    for i in range(1, n_rows):
        matriz[0, i, 0] = -math.inf
        matriz[1, i, 0] = matriz[0, 0, 0] + f_gap(gap, ext, i)
        matriz[2, i, 0] = -math.inf

    for j in range(1, n_cols):
        matriz[0, 0, j] = -math.inf
        matriz[1, 0, j] = -math.inf
        matriz[2, 0, j] = matriz[0, 0, 0] + f_gap(gap, ext, j)

    return matriz

def rellenado_gotoh(matriz, A, B, n_cols, n_rows, match, mismatch, gap, ext):
    '''
    Rellenado de la matriz. En sentido de descentende de izquierda a derecha.
    '''

    for i in range(1, n_rows):
        for j in range(1, n_cols):

        #Condiciones para rellenar D. Se calcula el máximo valor de venir en diagonal desde cualquiera de las
        #3 matrices, haya una coincidencia (match) o no (mismatch).
            if A[j - 1] == B[i - 1]:
                val_d = matriz[0, i-1, j-1] + match
                val_h = matriz[2, i-1, j-1] + match
                val_v = matriz[1, i-1, j-1] + match
                matriz[0,i,j] = max([val_d, val_h, val_v])

            elif A[j-1] != B[i-1]:
                val_d = matriz[0, i-1, j-1] + mismatch
                val_h = matriz[2, i-1, j-1] + mismatch
                val_v = matriz[1, i-1, j-1] + mismatch
                matriz[0,i,j] = max([val_d, val_h, val_v])

            #Rellenado de H. Se evalua el máximo de venir en diagonal desde D o de venir en horizontal desde H
            val_d_h = matriz[0, i, j-1] + gap
            val_h = matriz[1, i, j-1] + ext
            matriz[1, i, j] = max([val_d_h, val_h])

            #Rellenado de V. Se evalua el máximo de venir en diagonal desde D o de venir en vertical desde V
            val_d_v = matriz[0, i-1, j] + gap
            val_v = matriz[2, i-1, j] + ext
            matriz[2, i, j] = max([val_d_v, val_v])

    return matriz

def vuelta_atras_gotoh(matriz, A, B, n_cols, n_rows, match, mismatch, gap, ext):
    '''Vuelta atrás. Se realiza en sentido ascendente de derecha a izquierda, comprobando en cada paso del alineamiento
    desde que matriz nos hemos movido.
    '''

    alin_A = str()
    alin_B = str()

    #Búsqueda del valor máximo del alineamiento global. Por definición debe estar en matriz[x, n_rows-1, n_cols-1],
    #por lo que buscamos en qué matriz (D,H o V) está el valor máximo por el que se incia el alineamiento.

    init_list = [matriz[0, n_rows-1, n_cols-1], matriz[1, n_rows-1, n_cols-1], matriz[2, n_rows-1, n_cols-1]]
    m_actual = init_list.index(max(init_list))

    i = n_rows - 1
    j = n_cols - 1

    while i != 0 or j != 0:     #Condición de parada es llegar a la primera fila o columna en alguna de las matrices

        if m_actual == 0:         #Si nos encontramos en la matriz diagonal
            #Puede que vengamos en diagonal desde cualquiera de las 3 matrices y que exista un match o un mismatch
            for x in range(3):
                if matriz[0, i, j] == matriz[x, i-1, j-1] + match and A[j-1] == B[i-1]:
                    m_actual = x        #Se actualiza la variable con la matriz actual para la proxima iteración del while

                elif matriz[0, i, j] == matriz[x, i-1, j-1] + mismatch and A[j-1] != B[i-1]:
                    m_actual = x

            #En cualquiera de los casos, los caracteres de A (superior) se alinean con los B (izquierda)
            alin_A = A[j-1] + alin_A
            alin_B = B[i-1] + alin_B

            #Se actualizan los índices para un movimiento en diagonal (restando 1 a ambos)
            i -= 1
            j -= 1
            continue

        #Si nos encontramos en la matriz horizontal
        elif m_actual == 1:     #Solo existen 2 posibilidades, venir en horizontal desde la propia matriz o desde D
            if matriz[1, i, j] == matriz[0, i, j-1] + gap:    #Si venimos de la diagonal la penalización es de apertura
                m_actual = 0

            elif matriz[1, i, j] == matriz[1, i, j-1] + ext:  #Si venimos de la propia matriz la penalizacion es de extension
                m_actual = 1

            #Los movimientos en horizontal provocan huecos en la secuencia B(izquierda) alineados con caracteres en la sequencia A (superior)
            alin_A = A[j-1] + alin_A
            alin_B = '-' + alin_B

            #Se actualiza solo el índice j, nos movemos una columna a la izquierda, manteniendo la misma fila (índice i)
            j -= 1
            continue

        #Si nos encontramos en la matriz vertical
        elif m_actual == 2:       #Solo existen 2 posibilidades, venir en vertical desde la propia matriz o desde D
            if matriz[2, i, j] == matriz[0, i-1, j] + gap:
                m_actual = 0
            elif matriz[2, i, j] == matriz[2, i-1, j] + ext:
                m_actual = 2

            #Los movimientos en vertical provocan huecos en la secuencia A(superior) alineados con caracteres en la sequencia B (izquierda)
            alin_A = '-' + alin_A
            alin_B = B[i-1] + alin_B

            #Se actualiza solo el índice i, nos movemos una fila hacia arriba, manteniendo la misma columna (índice j)
            i -= 1
    return alin_A, alin_B

def gotoh(A,B,match,mismatch,gap,ext):
    n_cols = len(A) + 1
    n_rows = len(B) + 1
    m = inicializacion_gotoh(n_rows, n_cols, gap, ext)
    m = rellenado_gotoh(m, A, B, n_cols, n_rows, match, mismatch, gap, ext)
    a,b = vuelta_atras_gotoh(m, A, B, n_cols, n_rows, match, mismatch, gap, ext)
    return a,b

# NOTE: NW Generalizado

def inicializacion_gen(n_rows, n_cols, gap, ext):
    matriz = np.full([n_rows, n_cols], 0)       #se genera una matriz de ceros con las dimensiones deseadas

    for i in range(1,n_rows):                   #se iteran las celdas de la primera columna
        matriz[i,0] = matriz[0, 0] + f_gap(gap, ext, i)   #al valor de la matriz en 0,0 se le aplica f_gap() en cada iteración

    for j in range(1,n_cols):
        matriz[0,j] = matriz[0, 0] + f_gap(gap, ext, j)
#np.set_printoptions(formatter={'float': lambda x: "{0:0.0f}".format(x)})
    return matriz

def rellenado_gen(matriz, n_rows, n_cols, A, B, match, mismatch, gap, ext):
    for i in range(1, n_rows):   #se fija un fila
        for j in range(1, n_cols):    #se iteran las columnas

            #diagonal
            if A[j - 1] == B[i - 1]:     #si la letra en la diagonal coincide, sumale el valor del match
                d = matriz[i-1, j-1] + match
            else:                        #si no, el del mismatch
                d = matriz[i-1, j-1] + mismatch

            maximo = d                  #se actualiza el valor del máximo

            #arriba
            for z in range(i):          #para cada fila z desde la actual i
                k = i - z                   #k es el tamaño del hueco en vertical
                valor = matriz[z, j] + f_gap(gap, ext, k)   #se calcula el valor a partir de la celda de partida utilizando f_gap()
                if valor > maximo:          #si el valor es mayor que el máximo se almacena
                    maximo = valor

            #izquierda
            for z in range(j):        #para cada columna z desde la actual j
                l = j - z                 #l es el tamaño del hueco en horizontal
                valor = matriz[i, z] + f_gap(gap, ext, l)
                if valor > maximo:
                    maximo = valor


            #el valor en la matriz[i,j] se calcula como el máximo de los dos valores
            matriz[i,j] = max([maximo, d])

    return matriz

def vuelta_atras_gen(matriz, n_rows, n_cols, A, B, match, mismatch, gap, ext):
    '''
    Vuelta atrás. El objetivo es buscar un alineamiento global, por esto se realiza un recorrido
    desde la esquina inferior derecha de la matriz hasta la esquina superior izquierda.
    Hasta llegar a las coordenadas (0,0) de la matriz se mira en cada paso desde que celda
    se produce el movimiento
    '''
    alin_A = str()
    alin_B = str()

    i = n_rows - 1
    j = n_cols - 1

    while i != 0 or j != 0:
        if matriz[i,j] == matriz[i-1, j-1] + match and A[j-1] == B[i-1]:  #si el valor es el de la diagonal + match y las letras coinciden
            alin_A = A[j-1] + alin_A        #añade la letra correspondiente a ambos alineamientos
            alin_B = B[i-1] + alin_B
            i -= 1                          #actualiza los índices para la fila y columna
            j -= 1

        elif matriz[i,j] == matriz[i-1, j-1] + mismatch and A[j-1] != B[i-1]: #si el valor es el de la diagonal + mismatch y las letras no coinciden
            alin_A = A[j-1] + alin_A        #añade las letras correspondientes en cada secuencia
            alin_B = B[i-1] + alin_B
            i -= 1                          #actualiza los índices para la fila y columna
            j -= 1

        else:     #si no se dan cualquiera de las dos condiciones anteriores
            for z in range(i):  #Movimientos en vertical
                k = i - z
                if matriz[i,j] == matriz[z,j] + f_gap(gap, ext, k):
                    alin_A = k*'-' + alin_A         #se añaden k huecos al alineamiento en la secuencia superior
                    alin_B = B[i-k:i] + alin_B      #se alinean con las k letras desde la posicion i de partida en la secuencia de la izq
                    i = i - k                       #se actualiza el valor de la fila restando k (número de posiciones del desplazamiento en vertical)
                    continue

            for z in range(j):  #Movimientos en horizontal
                l = j - z
                if matriz[i,j] == matriz[i,z] + f_gap(gap, ext, l):
                    alin_A = A[j-l:j] + alin_A      #se alinean las l letras desde la posición j de partida en la secuencia superior
                    alin_B = l*'-' + alin_B         #se añaen l huecos al alineamiento en la secuencia de la izq
                    j = j - l                       #se actualiza el valor de la columna restando l (número de posiciones del desplazamineto en horizontal)
                    continue

    return alin_A, alin_B

def NWGeneralizado(A, B, match, mismatch, gap, ext):
    n_cols = len(A) + 1
    n_rows = len(B) + 1
    m = inicializacion_gen(n_rows, n_cols, gap, ext)
    m = rellenado_gen(m, n_rows, n_cols, A, B, match, mismatch, gap, ext)
    a,b = vuelta_atras_gen(m, n_rows, n_cols, A, B, match, mismatch, gap, ext)
    return a,b

#SW_sublocal aligments

def inicializacion_SW_global(A, B, matrizPuntuacion, T, gap, n_cols, n_rows):
  '''
  Parametros:
  A (string) = cadena de caracteres de la secuencia A (en la cual se hace la busqueda).
  B (string) = cadena de caracteres de la secuencia B.
  matrizPuntuacion (matrix) = matriz de puntuacion utilizada.
  T (int) = valor numerico del umbral.
  gap (int) = valor numerico del hueco.
  n_cols (int) = numero de columnas (longitud de la secuencia A).
  n_rows (int) = numero de filas (longitud de la secuencia B).
  '''
  matriz = np.zeros([n_rows, n_cols])       #se genera una matriz de ceros con las dimensiones deseadas
  #rellenado de la matriz
  for i in range(1, n_cols):   #se fija un col

    for j in range(0, n_rows):    #se iteran las filas
      if j == 0:
        maximo = 0
        for j_ in range(1,n_rows):
          if matriz[j_,i-1]-T > maximo:
            maximo = matriz[j_,i-1] - T
        matriz[0,i] = max(maximo, matriz[0,i-1])
      else:

        diagonal = matriz[j-1, i-1] + matrizPuntuacion[A[i-1], B[j-1]]

        arriba = matriz[j,i-1] + gap
        izquierda = matriz[j-1,i] + gap

        matriz[j,i] = max([diagonal, arriba, izquierda, matriz[0,i]])
  return matriz

def vuelta_atras_SW_global(T, gap, matrizPuntuacion, matriz, n_cols, n_rows, A, B):
  '''
  Parametros:
  T (int) = valor numerico del umbral.
  gap (int) = valor numerico del hueco.
  matrizPuntuacion (matrix) = matriz de puntuación utilizada.
  resultadoMatriz (matrix) = matriz resultante de realizar el rellenado.
  n_cols (int) = numero de columnas (longitud de la secuencia A).
  n_rows (int) = numero de filas (longitud de la secuencia B).
  '''
  for j_ in range(1,n_rows):
    if matriz[0,n_cols-1] < matriz[j_,n_cols-1] - T:
      j = j_
    else:
      j = 0
    i = n_cols - 1
    alin_B = str()

  while True:
    if j == 0 and i == 0:
      break
    if j == 0:
      if matriz[0,i] == matriz[0,i-1]:
        alin_B = '.' + alin_B
        i -= 1
      else:
        for j_ in range(1, n_rows):
          if matriz[j_,i-1] == matriz[0,i] + T:
            maximo = matriz[j_,i-1]
            j = j_
            i -= 1
            alin_B = '.' + alin_B
            break
    else:
      if matriz[j,i] == matriz[j,i-1] + gap:
        alin_B = '-' + alin_B
        i -= 1
      elif matriz[j,i] == matriz[j-1,i-1] + matrizPuntuacion[A[i-1],B[j-1]]:
        alin_B = B[j-1] + alin_B
        i -= 1
        j -= 1
      elif matriz[j,i] == matriz[j-1,i]:
        j = 0
  return alin_B

def SW_global(A, B, matrizPuntuacion, T, gap):
  '''
  Parametros:
  A (string) = cadena de caracteres de la secuencia A (en la cual se hace la busqueda).
  B (string) = cadena de caracteres de la secuencia B.
  matrizPuntuacion (matrix) = matriz de puntuacion utilizada.
  T (int) = valor numerico del umbral.
  gap (int) = valor numerico del hueco.
  n_cols (int) = numero de columnas (longitud de la secuencia A).
  n_rows (int) = numero de filas (longitud de la secuencia B).
  '''
  n_cols = len(A) + 1
  n_rows = len(B) + 1
  matriz = inicializacion_SW_global(A, B, matrizPuntuacion, T, gap, n_cols, n_rows)
  b = vuelta_atras_SW_global(T, gap, matrizPuntuacion, matriz, n_cols, n_rows,A,B)
  print(A)
  print(b)
  return b

###k-band

def inicializacion_kb_global(n_rows, n_cols, k, gap):
    #Incialización de la matriz
    matriz = np.full([n_rows, n_cols], float(0))  #se genera una matriz de ceros
    for j_2 in range(k+1, n_cols):    #rellenado de la primera fila, valor de la celda de la izq + el gap, desde k+1 hasta el final
        matriz[0, j_2] = matriz[0, j_2 - 1] + gap

    for i in range(1, k + 1):
        matriz[i, k-i] = matriz[i-1, k-i+1] + gap

    return matriz


def rellenado_kb_global(A,B,matriz,k, mp, gap):
    for i in range(1, len(B)+1):  #Rango de las filas
        for j_2 in range(max(0, k-i+1), min(2*k+1, len(A)+k+1-i)):  #Rango de las columnas, varía según la fila
            #incialización de los valores de la nueva casilla i, j_2
            diagonal = -math.inf
            izq = -math.inf
            arriba = -math.inf

            if j_2 <= 2*k-1:  #no se puede venir en diagonal desde valores de columna fuera de este rango
                diagonal = matriz[i-1, j_2 + 1] + gap

            if j_2 >= 1:      #no se puede venir desde la izquierda si el valor de columna no es superior o igual 1
                izq = matriz[i, j_2-1] + gap

            b = B[i-1]
            a = A[j_2 + i - k - 1]
            arriba = matriz[i-1, j_2] + mp[a,b]

            matriz[i, j_2] = max([arriba, izq, diagonal]) #el valor de la matriz se calcula como el máximo de los tres valores

    return matriz

def rellenado_kb_local(A,B,matriz,k, mp, gap):
    i_max = 0
    j_max = 0
    for i in range(1, len(B)+1):  #Rango de las filas
        for j_2 in range(max(0, k-i+1), min(2*k+1, len(A)+k+1-i)):  #Rango de las columnas, varía según la fila
            #incialización de los valores de la nueva casilla i, j_2
            diagonal = 0
            izq = 0
            arriba = 0

            if j_2 <= 2*k-1:  #no se puede venir en diagonal desde valores de columna fuera de este rango
                diagonal = matriz[i-1, j_2 + 1] + gap

            if j_2 >= 1:      #no se puede venir desde la izquierda si el valor de columna no es superior a 1
                izq = matriz[i, j_2-1] + gap

            b = B[i-1]
            a = A[j_2 + i - k - 1]

            arriba = matriz[i-1, j_2] + mp[a,b]

            matriz[i, j_2] = max([arriba, izq, diagonal, 0]) #el valor de la matriz se calcula como el máximo de los tres valores
            if matriz[i, j_2] > matriz[i_max, j_max]:
                i_max = i
                j_max = j_2

    return matriz, i_max, j_max

def vuelta_atras_kb_global(A,B,matriz, k, mp, gap):
    i = len(B)
    j_2 = 2*k - (k -(len(B)-len(A)))

    alin_A = str()
    alin_B = str()

    while i != 0 or j_2 != k - i:
        #izquierda(k-band) == izquierda(NW)
        if matriz[i, j_2] == matriz[i, j_2 -1] + gap and j_2 >= max(1,(k-i)+1):
            alin_A = A[j_2 + i - k - 1] + alin_A
            alin_B = '-' + alin_B
            j_2 -= 1

        #diagonal(k-band) == arriba(NW)
        elif matriz[i, j_2] == matriz[i-1, j_2 + 1] + gap and j_2 < 2*k-1:
            alin_A = '-' + alin_A
            alin_B = B[i-1] + alin_B
            i -= 1
            j_2 += 1

        #arriba(k-band match) == diagonal(match NW)
        elif matriz[i, j_2] == matriz[i-1, j_2] + mp[A[j_2 + i -k - 1],B[i-1]]:
            alin_B = B[i-1] + alin_B
            alin_A = A[j_2 + i -k - 1] + alin_A
            i -= 1

    return alin_A, alin_B

def vuelta_atras_kb_local(A,B,matriz, k, mp, gap, i_max, j_max):
    i = i_max
    j_2 = j_max

    alin_A = str()
    alin_B = str()

    while matriz[i, j_2] != 0:
        #izquierda(k-band) == izquierda(NW)
        if matriz[i, j_2] == matriz[i, j_2 -1] + gap and j_2 >= max(1,(k-i)+1):
            alin_A = A[j_2 + i - k - 1] + alin_A
            alin_B = '-' + alin_B
            j_2 -= 1

        #diagonal(k-band) == arriba(NW)
        elif matriz[i, j_2] == matriz[i-1, j_2 + 1] + gap and j_2 < 2*k-1:
            alin_A = '-' + alin_A
            alin_B = B[i-1] + alin_B
            i -= 1
            j_2 += 1

        #arriba(k-band match) == diagonal(match NW)
        elif matriz[i, j_2] == matriz[i-1, j_2] + mp[A[j_2 + i -k - 1],B[i-1]]:
            alin_B = B[i-1] + alin_B
            alin_A = A[j_2 + i -k - 1] + alin_A
            i -= 1

    return alin_A, alin_B
### FASTP
def k_band_fastp_results(alin_A, alin_B, idA, idB, score):
    print('aligment: {0}\t ID: {1}'.format(alin_A, idA))
    print('aligment: {0}\t ID: {1}'.format(alin_B, idB))
    print('Optimized Score: {}'.format(score))

def k_band_fastp(A,idA,B,idB,k,mp,gap):
    if len(B) > len(A):
        B, A = A, B
    n_cols = 2*k + 1
    n_rows = len(B) + 1
    m = np.full([n_rows, n_cols], float(0))
    m, i_max, j_max = rellenado_kb_local(A,B,m,k, mp, gap)
    score = m[i_max, j_max]
    a,b = vuelta_atras_kb_local(A,B,m, k, mp, gap, i_max, j_max)
    k_band_fastp_results(a, b, idA, idB, score)
    return score

def k_band_global(A,B,k,mp,gap):
    if len(B) > len(A):
        B, A = A, B
    n_cols = 2*k + 1
    n_rows = len(B) + 1
    m = inicializacion_kb_global(n_rows, n_cols, k, gap)
    m = rellenado_kb_global(A,B,m,k, mp, gap)
    a,b = vuelta_atras_kb_global(A,B,m, k, mp, gap)
    print(a)
    print(b)

###DISTANCIA DE LEVENSHTEIN
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