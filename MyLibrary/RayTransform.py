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
import numpy.matlib
import time 

if __name__ == '__main__':
    import utils
else: 
    import MyLibrary.utils as utils

from scipy import special



class RayTransform(): 
    
    """  
RayTransform
=====

Provides
  1. An object to handle and manipulate X-ray and Radon transforms
  2. Easy to handle your 2D-data in computerized tomography
  3. Different geometries: parallel, fanbeam (ring, plane, any)
  4. Projection matrix, integral representation, FBP, ART, Approximate Inverse, 
      Landweber
  5. and assumes square images.
  
How to create an instance
----------------------------
    >>> parameters = (p,phi,N,center)
    >>> Xobj       = RT.RayTransform(parameters,modality='parallel',datatype ='matrix')

    Parameters
    ----------
    parameters : tuple containing the different parameters of the class
        depending on the type of modality/geometry.
    modality : indicates the type of CT-scan.
        'parallel' (Default): Standard parallel geometry. 
                            Leads to the so-called Radon transform.
        'fanbeam-ring': Detectors on an annulus.
        'fanbeam-plane': Detectors on a plane (line).
        'fanbeam': offers more freedom and requires only the positions of the 
                    detector sets and source all along the acquisition.
    datatype : of the projection matrix
        'matrix' (Default) uses a numpy-array 
        'sparse' uses a dictionary architecture. To prefer for large datasets.
                
    if modality == 'parallel': 
        
        parameters is a tuple of length 4 containing
            parameterRot : numpy.ndarray with one dimension. 
                        Gives the detector parameters.
            parameterDet : numpy.ndarray with one dimension. 
                        Gives the rotation parameters. 
            N            : the width in pixel of the targeted resolution.
            center       : list or array of size 2. Gives the center of the image domain.
        
    if modality == 'fanbeam-ring' or modality == 'fanbeam-plane': 
    
        parameters is a tuple of length 6 containing
            parameterRot : numpy.ndarray with one dimension. 
                        Gives the detector parameters.
            parameterDet : numpy.ndarray with one dimension. 
                        Gives the rotation parameters. 
            N            : the width in pixel of the targeted resolution.
            center       : list or array of size 2. 
                        Gives the center of the image domain.
            distOS       : gives the distance between the center and the source.
            distOD       : gives the distance between the center and the detector set.
    
    if modality == 'fanbeam': 
        
        parameters is a tuple of length 4 containing
            parameterD : numpy.ndarray with two dimension. 
                    Gives the detector positions.
            parameterS : numpy.ndarray with one dimension. 
                    Gives the source positions. 
            N          : the width in pixel of the targeted resolution.
            center     : list or array of size 2. 
                    Gives the center of the image domain.
        
  
Functions
----------------------------
    computeProjectionMatrix : computes the projection matrix R
    ----------
    getData : computes g = R@f with f an image
    ----------
    RadonTransform : computes the Radon transform (parallel geometry) 
            via line integrals.
    ----------
    FBP : computes the standard "inverse" of the Radon transform, 
            a.k.a. the filtered backrpojection (only parallel geometry)  
    ----------
    ApproximateInverse : computes the approximate inverse for the Gaussian 
            mollifier. (all geometries except "fanbeam").
    ----------
    ART : computes the Kaczmarz algorithm to solve R@f = g
    ----------
    Landweber : computes the Landweber iterate to solve R@f = g. 
            Includes the Tikhonov regularization.
    ----------
    
See corresponding documentations

    Example:
    >>> Xobj.computeProjectionMatrix(method='pixel')
    >>> g     = Xobj.getData(f) 
    >>> gRT   = Xobj.RadonTransform(f)
    >>> fFBP  = Xobj.FBP(gRT,filtertype='hann')
    >>> fK    = Xobj.ART(g,sweeps=5)
    
References
----------
    Toft, P.; "The Radon Transform - Theory and Implementation", 
        	Ph.D. thesis, Department of Mathematical Modelling, Technical 
            University of Denmark, June 1996.

    Natterer, F.; "The mathematics of computerized tomography Classics 
            in Mathematics",
            Society for Industrial and Applied Mathematics, New York, 2001.

    Louis, A. K.; "Combining Image Reconstruction and Image Analysis with 
                Application to 2D-tomography", 
                M J. Imaging Sciences 1 pp. 188-208, 2008.
    
    """
    
    
    def __init__(self,parameters,modality='parallel', datatype ='matrix'): 
        """
    Parameters
    ----------
    parameters : tuple containing the different parameters of the class
        depending on the type of modality/geometry.
    modality : indicates the type of CT-scan.
        'parallel' (Default)
        'fanbeam-ring' and 'fanbeam-plane'
        'fanbeam'.

        """
                
        if modality == 'fanbeam': 
            """
            parameters is a tuple of length 4 containing
                parameterD : numpy.ndarray with two dimension. Gives the detector positions.
                parameterS : numpy.ndarray with one dimension. Gives the source positions. 
                N          : the width in pixel of the targeted resolution.
                center     : list or array of size 2. Gives the center of the image domain.
            """
            
            self.parameterD   = parameters[0]
            self.parameterS   = parameters[1]
            self.N            = parameters[2]
            self.center       = parameters[3]
            self.matrix       = 0
            self.geometry     = modality
            self.datatype     = datatype
            self.datasize     = self.parameterD.shape[0:2]
            
        elif modality == 'fanbeam-ring' or modality == 'fanbeam-plane': 
            """
            parameters is a tuple of length 6 containing
                parameterRot : numpy.ndarray with one dimension. Gives the detector parameters.
                parameterDet : numpy.ndarray with one dimension. Gives the rotation parameters. 
                N            : the width in pixel of the targeted resolution.
                center       : list or array of size 2. Gives the center of the image domain.
                distOS       : gives the distance between the center and the source.
                distOD       : gives the distance between the center and the detector set.
            """
            self.parameterDet = parameters[0]
            self.parameterRot = parameters[1]
            self.N            = parameters[2]
            self.center       = parameters[3]
            self.distOS       = parameters[4]
            self.distOD       = parameters[5]
            self.matrix       = 0
            self.geometry     = modality
            self.datatype     = datatype
            self.datasize     = (self.parameterDet.shape[0],self.parameterRot.shape[0])
                    
        elif modality == 'parallel': 
            """
            parameters is a tuple of length 4 containing
                parameterRot : numpy.ndarray with one dimension. Gives the detector parameters.
                parameterDet : numpy.ndarray with one dimension. Gives the rotation parameters. 
                N            : the width in pixel of the targeted resolution.
                center       : list or array of size 2. Gives the center of the image domain.
            """
            self.parameterDet = parameters[0]
            self.parameterRot = parameters[1]
            self.N            = parameters[2]
            self.center       = parameters[3]
            self.matrix       = 0
            self.geometry     = modality
            self.datatype     = datatype
            self.datasize     = (self.parameterDet.shape[0],self.parameterRot.shape[0])
            
        else :
            raise AssertionError('The modality should be either: parallel, fanbeam, fanbeam-plane or fanbeam-ring.')

    def __str__(self):
        return str(self.matrix)
    
    def __repr__(self):
        return repr(self.matrix)
    
    def computeProjectionMatrix(self,gamma=0.5,method='pixel'):
    
        """
        Compute the Projection Matrix and update the attribute self.matrix
        
        Parameters
        ----------
        gamma : stands for the regularization parameter used in the Gaussian mollifier. 
                        Only used with the method 'gaussian'. Default is set to 0.5.
        method : 'stands for the way the projection matrix is computed. 
                        'pixel' -- length of the intersection pixel/line (Default)
                        'gaussian' -- analytical expression of the projection matrix for a Gaussian mollifier
        datatype : decides the type of the projection matrix.
                        'matrix' -- np.ndarray (Default)
                        'sparse' -- dictionary with a sparse architecture
        """
        
        geometry = self.geometry
        N = self.N
        center = self.center
        start_time = time.time()  
        
        if geometry == 'fanbeam': pass
                        
        elif geometry == 'fanbeam-ring' or geometry == 'fanbeam-plane':
            param_rota = self.parameterRot
            P = self.getPfromD() 
            Beta = np.arcsin(P/self.distOS) 

        elif geometry == 'parallel': 
            P   = self.parameterDet
            Beta = np.zeros(P.shape)
            param_rota = self.parameterRot
            
        else: raise Exception('Wrong geometry.')
                    
        if geometry == 'fanbeam': Np, Nphi = self.parameterD.shape[0:2]
        else : Nphi, Np = (param_rota.shape[0], P.shape[0])
        
        N_proj = Np*Nphi
        gridsize = (N,N)
        
        if method == 'pixel':      

            if self.datatype == 'sparse': self.matrix={}
            else:     self.matrix = np.zeros((N_proj,N*N))    
    
            for j in range(Nphi):   
                
                if self.datatype == 'sparse': utils.progressbar(j,Nphi-1,'Sparse Projection Matrix')
                else: utils.progressbar(j,Nphi-1,'Projection Matrix')
                
                for k in range(Np):
                    line    = Nphi*k+j
    
                    if geometry == 'fanbeam':
                        D = self.parameterD[k,j,:]
                        S = self.parameterS[j]
                        distOS = np.sqrt((S-center).dot(S-center))
                        distDS = np.sqrt((D-S).dot(D-S))
                        L      = (S-center).dot(D-S)
                        p      = (distOS**2 - (L/distDS)**2)
                        
                        if p<0: p = 0
                        else : p = np.sqrt(p)
                        
                        if (S - D  + L/distOS**2*(S-center)).dot([S[1],-S[0]]) <0: p = -p
                        
                        d1,d2 = tuple((D-S)/np.sqrt((D-S).dot(D-S)))
                        n1,n2 = (-d2,d1)
                        phi = np.arctan2(n1,n2)
                        
                    else : 
                        p,phi   = (P[k], param_rota[j] + Beta[k])  
                        n1,n2   = (np.cos(phi), np.sin(phi))
                        d1,d2   = (n2, -n1)
                    
                    xmin    = center[1]+0.5 + p*n1 - N*d1 
                    ymin    = center[0]+0.5 + p*n2 - N*d2 
                    xmax    = center[1]+0.5 + p*n1 + N*d1
                    ymax    = center[0]+0.5 + p*n2 + N*d2
                    pixels = []
                    values = []        
                    eps = 10**(-12)
            
                    if abs(phi % (np.pi/2)) < eps: 
                        if abs(d1) < abs(d2):
                            pixel = int(min(xmin,xmax))
                            pixels = [pixel+N*i for i in range(N) if pixel < N and pixel > -1]
                            values = [1]*N                
                        if abs(d1) > abs(d2):
                            pixel = N*int(min(ymin,ymax))
                            pixels = [pixel+i for i in range(N)  if pixel > -1 and  pixel < N**2]
                            values = [1]*N      
                        if len(pixels)>0:
                            if self.datatype == 'sparse': self.matrix[line] = {pixel:value for pixel,value in zip(pixels,values)}
                            else : 
                                self.matrix[line,pixels] = values
                        
                    else: 
                        if d1>0 and d2>0:
                            xgridmin,ygridmin,xgridmax,ygridmax = (eps,eps,N-eps,N-eps)
                            xpixel,ypixel,xshift,yshift         = (1,1,1,N)
                        if d1<0 and d2<0:
                            xgridmin,ygridmin,xgridmax,ygridmax = (N-eps,N-eps,eps,eps)
                            xpixel,ypixel,xshift,yshift         = (0,0,-1,-N)
                        if d1>0 and d2<0:
                            xgridmin,ygridmin,xgridmax,ygridmax = (eps,N-eps,N-eps,eps)
                            xpixel,ypixel,xshift,yshift         = (1,0,1,-N)
                        if d1<0 and d2>0:
                            xgridmin,ygridmin,xgridmax,ygridmax = (N-eps,eps,eps,N-eps)
                            xpixel,ypixel,xshift,yshift         = (0,1,-1,N)
                            
                        x,y,pixel = utils.get_intersection_with_grid(xmin,ymin,xgridmin,ygridmin,d1,d2,N)
                        xtest = int(x/eps)*eps
                        ytest = int(y/eps)*eps
            
                        if (xtest>=0 and xtest<=N and ytest>=0 and ytest<=N ): 
                    
                            xmax,ymax,pixelmax = utils.get_intersection_with_grid(xmax,ymax,xgridmax,ygridmax,-d1,-d2,N)
            
                            while ((x-xmax)**2+(y-ymax)**2>10**(-10)):    
                                
                                pixels.append(pixel)
                                x2 = pixel % N
                                y2 = (pixel - x2)/N
                                t,x,y,pixel_next = utils.get_intersection_with_pixel(x,y,x2+xpixel,y2+ypixel,d1,d2,pixel,xshift,yshift)
                                pixel = pixel_next
                                values.append(t)
                            if self.datatype == 'sparse': self.matrix[line] = {pixel:value for pixel,value in zip(pixels,values)}
                            else : self.matrix[line,pixels] = values
                
        elif method == 'gaussian': 
            x = np.arange(0,gridsize[1],1,dtype='int')
            y = np.arange(0,gridsize[0],1,dtype='int') 
            X,Y = np.meshgrid(y,x)
            X = X.flatten()
            Y = Y.flatten()
            Pixel = list(gridsize[1]*Y+X)
                
            if self.datatype == 'sparse': 
                self.matrix={}
                valuemax = np.exp(-4*(gamma)**2/(2*gamma**2))/(np.sqrt(2*np.pi)*gamma)
            else:     self.matrix = np.zeros((N_proj,gridsize[0]*gridsize[1]))    
            
            for j in range(Nphi):   
                
                if self.datatype == 'sparse': utils.progressbar(j,Nphi-1,'Sparse Projection Matrix')
                else: utils.progressbar(j,Nphi-1,'Projection Matrix')
                
                for k in range(Np):
                    
                    line    = Nphi*k+j
                    if geometry == 'fanbeam':
                        D = self.parameterD[k,j,:]
                        S = self.parameterS[j]
                        distOS = np.sqrt((S-center).dot(S-center))
                        distDS = np.sqrt((D-S).dot(D-S))
                        L      = (S-center).dot(D-S)
                        p      = (distOS**2 - (L/distDS)**2)
                        
                        if p<0: p = 0
                        else : p = np.sqrt(p)
                        
                        if (S - D  + L/distOS**2*(S-center)).dot([S[1],-S[0]]) <0: p = -p
                        
                        d1,d2 = tuple((D-S)/np.sqrt((D-S).dot(D-S)))
                        n1,n2 = (-d2,d1)
                        phi = np.arctan2(n1,n2)
                    else : 
                        p,phi   = (P[k], param_rota[j] + Beta[k])  
                        n1,n2   = (np.cos(phi), np.sin(phi))
      
                    
                    dist_pixel_line  = p - (X-center[1])*n1 - (Y-center[0])*n2 #n1*(center[0]-0.5) + n2*(center[1]-0.5)
                    values  = np.exp(-(dist_pixel_line)**2/(2*gamma**2))/(np.sqrt(2*np.pi)*gamma)
                    
                    if self.datatype == 'sparse':
                        self.matrix[line] = {pixel:value for pixel,value in zip(Pixel,values) if value > valuemax}     
            
                    else: self.matrix[line,Pixel] = values
                
        else: raise Exception('Wrong method. Only length and gaussian are supported.')  
        
        end_time = time.time()
        elapsed_time = int(10*(end_time - start_time))/10
        print('  '+f'Time elapsed: {elapsed_time} seconds',end='\n',flush=True)


    def getData(self,f):
        """
        Generates CT-data by evaluating R@f with R the projection matrix
        Parameters
        ----------
        f : a signal (typically an image) of numpy.ndarray type.

        Returns
        -------
        The data g as a numpy.ndarray type matrix. 

        """
        
        f = f.flatten()
                
        if isinstance(self.matrix,numpy.ndarray): g = np.matmul(self.matrix, f) 
        else: g = utils.sparsematrixproduct(self.matrix, f,self.datasize[0]*self.datasize[1]) 
        
        return  g.reshape(self.datasize)

    
    def RadonTransform(self,f):
        """
        Computes the line integrals (Radon transform) of an image f.

        Parameters
        ----------
        f : a square matrix of numpy.ndarray type.

        Raises
        ------
        works only for parallel geometry

        Returns
        -------
        The data g as a numpy.ndarray type matrix. 

        """
        
        if self.geometry != 'parallel': raise Exception('Geometry is not supported.')
        
        p       = self.parameterDet
        phi     = self.parameterRot
        N       = f.shape[0]
        Np      = p.shape[0]
        Nphi    = phi.shape[0]
        g       = np.zeros((Np,Nphi))
        n1,n2   = (np.cos(phi), np.sin(phi))
        d1,d2   = (n2, -n1)
        center  = self.center
        dq      = 0.5  
        q       = np.arange(-N/np.sqrt(2),N/np.sqrt(2),dq) 
        Nq      = q.shape[0]
        
        start_time = time.time()      
    
        for i in range(Nphi):
            utils.progressbar(i,Nphi-1,'Line integrals')
            
            for j in range(Np):
        
                xline = center[1] + p[j]*n1[i] + q*d1[i] 
                yline = center[0] + p[j]*n2[i] + q*d2[i] 
                index = [k for k in range(Nq) if (xline[k]>=2 and xline[k]<=N-2 and yline[k]>=2 and yline[k]<=N-2)]
    
                if not index: continue
                else:  
                    g[j,i] = sum(utils.bilinear_interp(f, xline[index], yline[index]))

    
        end_time = time.time()
        elapsed_time = int(10*(end_time - start_time))/10
        print('  '+f'Time elapsed: {elapsed_time} seconds',end='\n',flush=True)                       
        return g*dq  
    
    

            
    def FBP(self,g,method='integral',filtertype='ram-lak'):
        """
        Computes the filtered backprojection for parallel geometry.

        Parameters
        ----------
        g : the data as a numpy.ndarray type matrix. 
        method : 
            'integral' (Default): computes the adjoint operator as an integrak operator.
            'matrix': uses the transpose of the projection matrix
        filtertype : 
            'ram-lak' (Default), 'shepp-logan', 'cosine', 'hamming', 'hann'
                    
        Returns
        -------
        a numpy.ndarray type matrix

        """
        
        p       = self.parameterDet
        phi     = self.parameterRot
        N       = self.N
        dphi = (phi[1]-phi[0])
        dp   = (p[1]-p[0])
        
        wtrapz  = np.ones(phi.shape)
        wtrapz[0],wtrapz[-1] = (0.5,0.5) 
        Wtrapz  = np.tile(wtrapz, (p.shape[0],1))
        
        gf   = Wtrapz * utils.filterdata(g,filtertype)
        
        if method == 'integral':
            f  = self.backprojection(gf, method='FBP')
        
        if method == 'matrix':
            ProjMat = self.matrix
            if isinstance(ProjMat, numpy.ndarray): 
                f  = np.matmul(ProjMat.transpose(),gf.flatten()).reshape(N,N) 
                f  = dphi*dp*f 
            else: 
                Rt = utils.switchsparsity(self.matrix)
                f  = dphi*dp*utils.sparsematrixproduct(Rt, gf.flatten(),self.datasize[0]*self.datasize[1]).reshape(N,N)
    
        f[f<0] = 0
     
        return f/(2*dp)        
    
        
    def ApproximateInverseRT(self,g,gamma=1):
        """
        Computes the approximate inverse for all geometries expect 'fanbeam'.

        Parameters
        ----------
        g : the data as a numpy.ndarray type matrix.
        gamma : the standard deviation of the Gaussian distribution. 
                Plays the role of regularization parameter.

        Returns
        -------
        a numpy.ndarray type matrix

        References
        ----------
        Louis, A. K.; "Combining Image Reconstruction and Image Analysis with 
                    Application to 2D-tomography", 
                    M J. Imaging Sciences 1 pp. 188-208, 2008.

        """
        
        if self.geometry == 'fanbeam': raise Exception('Geometry is not supported.')
        
        p       = self.parameterDet
        phi     = self.parameterRot
        N       = self.N
        center  = self.center
        
        if self.geometry == 'parallel':
        
            Np,Nphi = (p.shape[0], phi.shape[0])
            f       = np.zeros((N,N))   
            P       = np.tile(p, (Nphi,1)).transpose()
            wtrapz  = np.ones(phi.shape)
            wtrapz[0],wtrapz[-1] = (0.5,0.5) 
            Phi     = np.tile(phi, (Np,1))
            Wtrapz  = np.tile(wtrapz, (Np,1))
            wtrapz  = Wtrapz.flatten()
            gvec    = g.flatten()
            pvec    = P.flatten()
            phivec  = Phi.flatten()
            n1,n2   = (np.cos(phivec), np.sin(phivec))
            kern    = np.zeros((N**2,Np*Nphi))
            
            start_time = time.time()       
            for x in range(N):
                utils.progressbar(x,N-1,'ApproxInverse')
                for y in range(N):      
                    z = (pvec - (x-center[1])*n1 - (y-center[0])*n2)/(np.sqrt(2)*gamma)
                    kern =  1 - 2*z*special.dawsn(z)            
                    f[y,x] = gvec.dot(wtrapz * kern)
            
            
            f[f<0] = 0
            f = f.reshape((N,N))
                   
            end_time = time.time()
            elapsed_time = int(10*(end_time - start_time))/10
            print('  '+f'Time elapsed: {elapsed_time} seconds',end='\n',flush=True) 
            
            return np.pi/(phi[-1]-phi[0])*2/(2*np.pi*gamma)**2*f*(p[1] - p[0])*(phi[1] - phi[0])
    
        elif self.geometry == 'fanbeam-ring' or self.geometry == 'fanbeam-plane':
            
            Np,Nphi = (p.shape[0], phi.shape[0])
            f = np.zeros((N,N))    
            pvec = np.tile(self.getPfromD(), (Nphi,1)).transpose().flatten()         
            phivec = np.tile(phi, (Np,1)).flatten() + np.arcsin(pvec/self.distOS)     
            gvec   = g.flatten()
            n1,n2  = (np.cos(phivec), np.sin(phivec))
            kern   = np.zeros((N**2,Np*Nphi))
            jacobi = np.tile(self.getDetJacobi(), (Nphi,1)).transpose().flatten()
            start_time = time.time()       
                
            for x in range(N):
                utils.progressbar(x,N-1,'ApproxInverse')
                for y in range(N):
                    if (x-N/2)**2 + (y-N/2)**2 < N**2/4:
                        z = (pvec - (x-N/2)*n1 - (y-N/2)*n2)/(np.sqrt(2)*gamma) 
                        kern =  (1 - 2*z*special.dawsn(z)) * jacobi        
                        f[y,x] = gvec.dot(kern)
            
            f[f<0] = 0
            f = f.reshape((N,N))
                   
            end_time = time.time()
            elapsed_time = int(10*(end_time - start_time))/10
            print('  '+f'Time elapsed: {elapsed_time} seconds',end='\n',flush=True) 
            
            return np.pi/(phi[-1]-phi[0])*2/(2*np.pi*gamma)**2*f*(p[1] - p[0])*(phi[1] - phi[0])
        
        else : raise Exception('This geometry is not supported for the function ApproximateInverseRT()')
        
        
    
    def Landweber(self, g, relax=0.000005, tikhonov_param=0, iterations=100):
        """
        Computes the Landweber iterate (Gradient descent with fixed stepsize).

        Parameters
        ----------
        g : the data as a numpy.ndarray type matrix.
        relax : relaxation parameter.
        tikhonov_param : controls the Tikhonov regularization. The default is 0.
        iterations : Number of iterations. The default is 100.

        Returns
        -------
        a numpy.ndarray type matrix
        
        """

        start_time = time.time() 
        g          = g.flatten()
        N          = self.N
        if self.datatype == 'matrix':
            Rt = self.matrix.transpose()
            Rtg = np.matmul(Rt,g)
            Rtg = Rtg[:,np.newaxis]
            RtR = np.matmul(Rt,self.matrix) + tikhonov_param*np.eye(self.matrix.shape[1])
            
        elif self.datatype == 'sparse':
            N = self.N
            Rt = utils.switchsparsity(self.matrix)
            Rtg = utils.sparsematrixproduct(Rt, g,self.N**2)  
            RtR = utils.getAtA(self.matrix,N,N**2,self.parameterDet.shape[0]*self.parameterRot.shape[0]) +  tikhonov_param*np.eye(N**2)
            
        f = np.zeros((N**2,1))
    
        for i in range(iterations):
            utils.progressbar(i,iterations-1,'Landweber')
            f = f + relax * (Rtg - np.matmul(RtR,f))
            f[f<0] = 0
            
        end_time = time.time()
        elapsed_time = int(10*(end_time - start_time))/10
        print('  '+f'Time elapsed: {elapsed_time} seconds',end='\n',flush=True) 
        return f.reshape(N,N)
    
    
        
    def ART(self,g,sweeps=1,relax=1):
        """
        computes the Kaczmarz algorithm (ART).

        Parameters
        ----------
        g : the data as a numpy.ndarray type matrix.
        sweeps : number of sweeps. The default is 1.
        relax : relaxation parameter. The default is 1.

        Returns
        -------
        a numpy.ndarray type matrix
        
        """
        
        start_time = time.time() 
        g          = g.flatten()
        f          = np.zeros(self.N**2)
        N          = self.N
        if self.geometry == 'fanbeam':  Np, Nphi = self.parameterD.shape[0:2]
        else : Np,Nphi    = (self.parameterDet.shape[0], self.parameterRot.shape[0])
        
        j = 0
        
        
        if self.datatype == 'matrix':
            normR = np.zeros(Np*Nphi)
            for k in range(Np*Nphi): normR[k] = np.dot(self.matrix[k,:],self.matrix[k,:])
        
            for i in range(sweeps):
                for k in np.random.choice(range(Np*Nphi), p=normR/np.sum(normR), size=(Np*Nphi)):
                    j+=1
                    if j%Np == 0: utils.progressbar(j,sweeps*Nphi*Np,'ART')
                    Rvec = self.matrix[k,:]
                    f += relax**i * (g[k]-np.dot(Rvec,f))/normR[k]*Rvec
                    f[f<0] = 0
                    
        elif self.datatype == 'sparse':
            for i in range(sweeps):
                for k in np.random.permutation(Np*Nphi):
                    j+=1
                    if j%Np == 0: utils.progressbar(j,sweeps*Nphi*Np,'ART')
                    Rvec = utils.getRowfromSparse(self.matrix,Np*Nphi,N**2,k)
                    normR = np.dot(Rvec,Rvec)
                    if normR>10**(-5): f += relax**i * (g[k]-np.dot(Rvec,f))/normR*Rvec
                    f[f<0] = 0

        
        end_time = time.time()
        elapsed_time = int(10*(end_time - start_time))/10
        print('  '+f'Time elapsed: {elapsed_time} seconds',end='\n',flush=True) 
        return f.reshape(N,N)
    
    
    
    def switchfanbeam2parallel(self,g,new_parameters):
        """
        transforms the fanbeam (ring or plane) data into Radon data (parallel).

        Parameters
        ----------
        g : the data as a numpy.ndarray type matrix.
        new_parameters : parameters of the parallel geometry.

        Returns
        -------
        newg : a numpy.ndarray type matrix

        """
        
        if self.geometry == 'fanbeam': raise Exception('Geometry is not supported.')
        
        p,phi = new_parameters
        Nphi  = phi.shape[0]
        Np    = p.shape[0]
        newg = np.zeros((Np,Nphi))
        
        param_detector = self.parameterDet
        param_rotation = self.parameterRot
        Nb    = param_rotation.shape[0]
        Na    = param_detector.shape[0]
        Da    = param_detector[1] - param_detector[0]
        amin  = np.min(param_detector)
        Db    = param_rotation[1] - param_rotation[0]
        bmin  = np.min(param_rotation)
        alpha = self.getDfromP(p) 
        beta0 = np.arcsin(p/self.distOS) 
        indexa= (alpha - amin)/Da     
        
        start_time = time.time()
        
        for i in range(Nphi):
            utils.progressbar(i,Nphi-1,'Fanbeam --> Parallel')
        
            beta  = phi[i] - beta0
            indexb = (beta - bmin)/Db
            index = [k for k in range(Na) if (indexa[k]>=1 and indexa[k]<=Na-1 and indexb[k]>=0 and indexb[k]<=Nb)]
    
            if index: 
                newg[index,i] = utils.bilinear_interp(g, indexb[index], indexa[index])
        
        end_time = time.time()
        elapsed_time = int(10*(end_time - start_time))/10
        print('  '+f'Time elapsed: {elapsed_time} seconds',end='\n',flush=True) 
        return newg
    
    
    

#######################################
### The following are inner methods ###
#######################################
    
    def backprojection(self,g,method= 'BP'):
        p       = self.parameterDet
        phi     = self.parameterRot
        N       = self.N
        center = self.center
        Nphi = phi.shape[0]
        f = np.zeros((N,N))
        n1,n2   = (np.cos(phi), np.sin(phi))   
        x = np.arange(0,N,1)
        y = np.arange(0,N,1)
        xx,yy = np.meshgrid(x-center[1],y-center[0])
        
        start_time = time.time()  
    
        for k in range(Nphi):
            utils.progressbar(k,Nphi-1,method)
            f +=  np.interp(n1[k]*xx+n2[k]*yy, p, g[:,k])
                            
        end_time = time.time()
        elapsed_time = int(10*(end_time - start_time))/10
        print('  '+f'Time elapsed: {elapsed_time} seconds',end='\n',flush=True)  
        return (phi[1]-phi[0])*f
    
    
    def getPfromD(self):
        OS,OD = (self.distOS,self.distOD)        
        if self.geometry == 'fanbeam-ring':    
            p = np.divide(OS*OD*np.sin(self.parameterDet),np.sqrt(OS**2+OD**2+2*OS*OD*np.cos(self.parameterDet)))        
        elif self.geometry == 'fanbeam-plane': 
            p = np.divide(OS*self.parameterDet,np.sqrt((OS+OD)**2 + self.parameterDet**2))
        else: raise Exception('Wrong geometry.')
        return p
    
    def getDfromP(self,p):
        OS,OD = (self.distOS,self.distOD)   
        if self.geometry == 'fanbeam-ring': d = np.arcsin(p/OD)+np.arcsin(p/OS)         
        elif self.geometry == 'fanbeam-plane': d = np.divide((OS+OD)*p,np.sqrt(OS**2 - p**2))
        else: raise Exception('Wrong geometry.')
        return d
    
    def getDetJacobi(self):
        OS,OD = (self.distOS,self.distOD)   
        if self.geometry == 'fanbeam-ring': 
            J = (OS+OD*np.cos(self.parameterDet))*(OD+OS*np.cos(self.parameterDet))* \
                                                     pow(OS**2+OD**2+2*OS*OD*np.cos(self.parameterDet),-1.5)       
        elif self.geometry == 'fanbeam-plane': J = OS*(OS+OD)**2*pow((OS+OD)**2+self.parameterDet**2,-1.5)
        else: raise Exception('Wrong geometry.')
        return J
    
    def getLinefromSD(self):
        
        S = self.S 
        D = self.D
        
        Vx = D[0] - S[0]
        Vy = D[1] - S[1]    
        normV = np.sqrt(Vx**2+Vy**2)
        Vx = np.divide(Vx,normV)
        Vy = np.divide(Vy,normV)
        Hx = S[0] - (S[0]*Vx+S[1]*Vy) * Vx
        Hy = S[1] - (S[0]*Vx+S[1]*Vy) * Vy
        p = np.sqrt(Hx**2+Hy**2)
        phi = np.arctan2(-Vy,Vx)%(2*np.pi)-np.pi/2 
        
        return p,phi
    
    
    
    
    
    
       

###############################################################################        

if __name__ == '__main__': pass
        
