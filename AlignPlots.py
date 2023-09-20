import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from IPython.core.display import HTML,display
import pandas as pd

def tabula(matrix, seq1, seq2):
    row_names = [i for i in '-'+seq1]
    col_names = [i for i in '-'+seq2]
    df = pd.DataFrame(matrix,index=row_names,columns=col_names)
    table_html = df.to_html()
    return HTML(table_html)
#### Algoritmo Needleman-Wunsch-Sellers

def opt_path_nws(x,y,alinA, alinB, p_matriz, A, B, ax, arrowprops, dicc):
  ax.grid(True, which='minor')
  if (x,y) == (0,0):
   #print(alinA)
   #print(alinB)
   #print('.'*len(alinA))
   return
  if p_matriz[x,y] == p_matriz[x,y-1] + dicc['GAP']:
    ax.annotate('',xy=(x,y-1), xytext=(x,y), arrowprops=arrowprops)
    opt_path_nws(x, y-1,A[y-1] + alinA, "-" + alinB, p_matriz, A, B, ax, arrowprops, dicc)
  if p_matriz[x,y] == p_matriz[x-1, y] + dicc['GAP']:
    ax.annotate('',xy=(x-1,y),xytext=(x,y), arrowprops=arrowprops)
    opt_path_nws(x-1,y,"-" + alinA,B[x-1] + alinB, p_matriz, A, B, ax, arrowprops, dicc)
  if (p_matriz[x,y] == p_matriz[x-1, y-1] + dicc['MATCH'] and A[y-1] == B[x-1]) or (p_matriz[x,y] == p_matriz[x-1, y-1] + dicc['MISMATCH'] and A[y-1] != B[x-1]):
    ax.annotate('',xy=(x-1,y-1),xytext=(x,y), arrowprops=arrowprops)
    opt_path_nws(x-1,y-1, A[y-1] + alinA, B[x-1] + alinB, p_matriz, A, B, ax, arrowprops, dicc)

def plot_nw(A,B, dicc):
  plt.rcParams["figure.figsize"] = 0.8*(len(B)+1), 0.8*(len(A) + 1)
  param = {"grid.linewidth": 1.6,
           "grid.color":"lightgray",
           "axes.linewidth": 1.6,
          "axes.edgecolor": "lightgray",
          "font.size": 11}
  desp = 0.3
  desp2 = 0.2
  plt.rcParams.update(param)
  arrowprops = dict(facecolor='blue', alpha=0.5, lw=0,
                    shrink=0, width=0.5, headwidth=5, headlength=5)
  fig, ax = plt.subplots()
  n_cols = len(A) + 1
  n_rows = len(B) + 1
  ax.set_xlim(-1.5, n_rows - 0.5)
  ax.set_ylim(-1.5, n_cols - 0.5)

  ax.invert_yaxis()
  p_matriz = np.full([n_rows, n_cols], 0, dtype=int)

  for i, l in enumerate(B):
    ax.text(i + 1, -1, l, ha="center", va="center", fontweight="semibold")
  for i, l in enumerate(A):
    ax.text(-1, i + 1, l, ha="center", va="center", fontweight="semibold")

  ax.tick_params(labelbottom=False)
  ax.tick_params(labelleft=False)
  ax.set_xticks([])
  ax.set_yticks([])
  ax.text(0,0,p_matriz[0,0],ha='center',va='center')

  for i in range(1,n_rows):
    p_matriz[i,0] = p_matriz[i-1, 0] + dicc['GAP']
    ax.text(i,0,p_matriz[i,0],ha='center', va='center')
    ax.annotate('',xy=(i-1+desp,0), xytext=(i-desp,0), arrowprops=arrowprops)

  for j in range(1,n_cols):
    p_matriz[0,j] = p_matriz[0, j-1] + dicc['GAP']
    ax.text(0,j,p_matriz[0,j],ha='center', va='center')
    ax.annotate('',xy=(0,j-1+desp), xytext=(0,j-desp), arrowprops=arrowprops)

  ax.xaxis.set_minor_locator(ticker.FixedLocator(np.arange(-1.5, p_matriz.shape[0] - .5, 1)))
  ax.yaxis.set_minor_locator(ticker.FixedLocator(np.arange(-1.5, p_matriz.shape[1] - .5, 1)))
  plt.tick_params(axis='both', which='both', bottom='off', top='off', left="off", right="off", labelbottom='off', labelleft='off')
  ax.grid(True, which='minor')

  for i in range(1, n_rows):    #se fija una fila
    for j in range(1, n_cols):    #se iteran sus columnas
      izquierda = p_matriz[i,j-1] + dicc['GAP']     #el valor de venir de la casilla de la izquierda + el hueco
      arriba = p_matriz[i-1, j] + dicc['GAP']       #el valor de venir de la casiila de arriba + el hueco

      if A[j-1] == B[i-1]:                #si las letras son iguales
        diagonal = p_matriz[i-1, j-1] + dicc['MATCH']       #valor de la casilla diagonal superior + una coincidencia
      else:                               #si las letras difieren
        diagonal = p_matriz[i-1, j-1] + dicc['MISMATCH']   #valor de la casilla diagonal superior + mismatch

      p = max([arriba, izquierda, diagonal])    #se calcula el máximo de la 3 valores y se actualiza el valor de matriz[i,j]
      p_matriz[i,j] = p
      ax.text(i,j,p_matriz[i,j],ha='center', va='center')
      if p == arriba:
        ax.annotate('', xy=(i-1+desp,j),xytext=(i-desp,j), arrowprops=arrowprops)
      if p == izquierda:
        ax.annotate('', xy=(i,j-1+desp),xytext=(i,j-desp), arrowprops=arrowprops)
      if p == diagonal:
        ax.annotate('', xy=(i-1+desp2,j-1+desp2),xytext=(i-desp2,j-desp2), arrowprops=arrowprops)

  arrowprops = dict(facecolor='orange', alpha=0.2, lw=0,
                    shrink=0, width=10, headwidth=10, headlength=10)
  opt_path_nws(n_rows - 1, n_cols - 1, '', '', p_matriz ,A, B, ax, arrowprops, dicc)
  plt.gca().set_aspect('auto')
  plt.show()

  ### Nussinov
  def nussinov_tabula(matrix, seq):
    row_names = [i for i in seq]
    df = pd.DataFrame(matrix,index=row_names,columns=row_names)
    table_html = df.to_html()
    return HTML(table_html)