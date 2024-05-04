# =============================================================================
# MIT License
# 
# Copyright (c) 2024 Gael Rigaud  <https://www.f08.uni-stuttgart.de/organisation/team/Rigaud/>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# =============================================================================


import numpy as np
import time



def filterdata(g,filtertype='ram-lak'):
    """
    filters CT-data for the FBP algorithm.

    Parameters
    ----------
    g : the data as a numpy.ndarray type matrix. 
    filtertype : 
        'ram-lak' (Default), 'shepp-logan', 'cosine', 'hamming', 'hann'

    Returns
    -------
    a numpy.ndarray type matrix

    """

    Np = g.shape[0]
    Npadded = int(2**np.ceil(np.log2(Np)))   
    x = np.concatenate((np.arange(1, Npadded / 2 + 1, 2, dtype=int),
            np.arange(Npadded / 2 - 1, 0, -2, dtype=int)))
    y = np.zeros(Npadded)
    y[0] = 0.25
    y[1::2] = -1 / (np.pi * x) ** 2
    H = 2 * np.real(np.fft.fft(y))[:, np.newaxis]
    w = 2/Npadded*np.pi*np.concatenate((np.arange(0, Npadded / 2 + 1, 1, dtype=int),
            np.arange(Npadded / 2 - 1, 0, -1, dtype=int)))
    w = w[:,np.newaxis]
    
    
    if filtertype == 'ram-lak': pass
    elif filtertype == 'shepp-logan': H = H * np.sinc(w/2) 
    elif filtertype == 'cosine': H = H * np.cos(w/2)
    elif filtertype == 'hamming': H = H * (0.54+0.46*np.cos(w))
    elif filtertype == 'hann': H = H * (1+np.cos(w))/2
            
    g = np.pad(g, ((0, Npadded - Np), (0, 0)), mode="constant", constant_values=0)
    
    return np.real(np.fft.ifft(np.fft.fft(g, axis=0) * H, axis=0)[:Np, :])



def bilinear_interp(F,x,y):
    """
    Bilinear interpolation to estimate F(x,y).

    Parameters
    ----------
    F : 2D numpy.ndarray
    x : interpolation points (column/abscissa)
    y : interpolation points (row/ordinate)
    Returns
    -------
    numpy.ndarray matrix.

    """
    x1    = [int(i) for i in x]
    y1    = [int(i) for i in y]
    x2    = [i+1 for i in x1]
    y2    = [i+1 for i in y1]        

    return ((y2-y)*(F[y1,x1]*(x2-x) + F[y1,x2]*(x-x1)) + (y-y1)*(F[y2,x1]*(x2-x) + F[y2,x2]*(x-x1))) 
    

def get_intersection_with_grid(x,y,xgrid,ygrid,d1,d2,N):
    """
    computes the intersection of the line with grid.

    Parameters
    ----------
    x,y : original coordinates
    xgrid, ygrid : coordinates of the meshgrid
    d1,d2 : direction vector of the line
    N : grid size

    Returns
    -------
    the coordinates and the pixel index of the intersection.

    """

    t1 = (xgrid-x) /d1 #vertical
    t2 = (ygrid-y) /d2 #horizontal
    x1 = x + t1*d1
    x2 = x + t2*d1
    y1 = y + t1*d2
    y2 = y + t2*d2
    r1 = abs(x1-N/2) + abs(y1-N/2)
    r2 = abs(x2-N/2) + abs(y2-N/2)


    if r2>r1: return x1,y1,N*int(y1)+int(x1)
    else:     return x2,y2,N*int(y2)+int(x2)
    
    
def get_intersection_with_pixel(x,y,xpixel,ypixel,d1,d2,pixel,shiftx,shifty):
    """
    computes the intersection of the line with a given pixel.

    Parameters
    ----------
    x,y : original coordinates of the original point
    xpixel, ypixel : coordinates of the pixel
    d1,d2 : direction vector of the line
    pixel : index of the targeted pixel
    shiftx,shifty : relative indices for the final index.

    Returns
    -------
    distance to the intersection, 
    the coordinates of the intersection
    and the pixel index of the intersection.

    """
    t1 = (xpixel-x) /d1 #vertical
    t2 = (ypixel-y) /d2 #horizontal
    
    if t2>t1:
        return t1,x+t1*d1,y+t1*d2,pixel+shiftx
    else:
        return t2,x+t2*d1,y+t2*d2,pixel+shifty



def progressbar(i,N,strg_loop,bar_lengh=10,progress_char='#'):
    """
    progression bar
    
    """
    percentage = int((i/N)*100)                                                
    progress = int((bar_lengh * i ) / N)                                     
    loadbar = strg_loop + ":" + " [{:<{len}}]{}%".format(progress*progress_char,percentage,len = bar_lengh)
    print("\r" +  loadbar, end="", flush=True)    
    
    
    
########################################
### The following addresses sparsity ###
########################################

def switchsparsity(A):
    """
    Switch sparsity from column to row. Used to compute the transpose for 
    sparse matrices.

    Parameters
    ----------
    A : a sparse matrix (dict)

    Returns
    -------
    B : a sparse matrix (dict)

    """
    
    index_i = list(A.keys())
    B = {}
    
    for i in index_i:
        for j in A[i].keys():
            if j not in B: B[j] = {}
            B[j][i] = A[i][j]
    
    return B



def sparsematrixproduct(A,X,N):
    """
    Compute the matrix product for sparse matrices (dict). Y = A@X

    Parameters
    ----------
    A : a sparse matrix (dict)
    X : a numpy.ndarray
    N : number of rows in A

    Returns
    -------
    Y : a np.ndarray matrix 

    """
    
    index_i = list(A.keys())
    M = X.shape
    
    if len(M) == 1: 
        M1 = M[0]
        M2 = 1
        X  = X[:,np.newaxis]
    else:  M1,M2 = M
    
    Y = np.zeros((N,M2))
    B = np.zeros(M1)    
    
    for j in range(M2):
        for i in index_i:
            B[list(A[i].keys())] = list(A[i].values())
            Y[i,j] = B.dot(X[:,j])
            B    = 0*B
                  
    return Y

def getMatrixfromSparse(A,n,m):
    """
    Converts a sparse matrix (dict) to a numpy.ndarray.

    Parameters
    ----------
    A : a sparse matrix (dict)
    n : number of rows
    m : number of columns
    Returns
    -------
    B : a numpy.ndarray

    """
    
    B = np.zeros((n,m))
    index_i = list(A.keys())
    for i in index_i:
        B[i,list(A[i].keys())] = list(A[i].values())
    
    
    return B

def getColumnfromSparse(A,n,m,k):
    """
    Extract the k column into a sparse matrix A 

    Parameters
    ----------
    A : a sparse matrix (dict)
    n : number of rows
    m : number of columns
    k : column index

    Returns
    -------
    b : the column (numpy.ndarray)

    """
    
    b = np.zeros(n)
    index_i = list(A.keys())
    for i in index_i:
        if k in A[i].keys(): b[i] = A[i][k]   
    
    return b

def getRowfromSparse(A,n,m,k):
    """
    Extract the k row into a sparse matrix A 

    Parameters
    ----------
    A : a sparse matrix (dict)
    n : number of rows
    m : number of columns
    k : row index

    Returns
    -------
    b : the row (numpy.ndarray)

    """
    
    b = np.zeros(m)
    index_i = list(A.keys())
    if k in index_i: b[list(A[k].keys())] = list(A[k].values())
    
    
    return b

def getAtA(A,blocksize,n,m):
    """
    Computes the normal operator AtA for A sparse

    Parameters
    ----------
    A : a sparse matrix (dict)
    blocksize : parameter of the block multiplication
    n : number of rows
    m : number of columns

    Returns
    -------
    B : TYPE
        DESCRIPTION.

    """
    
    A = switchsparsity(A)
    index_i = list(A.keys())
    indb = np.arange(0,n+1,blocksize,dtype=int)
    B = np.zeros((n,n)) 
    index_block = [1]*(len(indb)-1)
    start_time = time.time() 

    for i in range(len(indb)-1):
        index_block[i] = [k for k in range(n) if (k>=indb[i] and k<indb[i+1])] 
    
    for i in range(len(indb)-1):
        progressbar(i,len(indb)-2,'Auxiliary matrix RtR')
        index_i_block = [k for k in index_i if  (k>=indb[i] and k<indb[i+1])]
        subA = dict((k-indb[i],A[k]) for k in index_i_block)
        subM = getMatrixfromSparse(subA, blocksize, m)
        
        for j in range(len(indb)-1):
            index_j_block = [k for k in index_i if  (k>=indb[j] and k<indb[j+1])]
            subA = dict((k-indb[j],A[k]) for k in index_j_block)
            subMj = getMatrixfromSparse(subA, blocksize, m).transpose()
            B[index_block[i][0]:index_block[i][-1]+1,index_block[j][0]:index_block[j][-1]+1] = subM @ subMj
            
    end_time = time.time()
    elapsed_time = int(10*(end_time - start_time))/10
    print('  '+f'Time elapsed: {elapsed_time} seconds',end='\n',flush=True)               
    return B

    

