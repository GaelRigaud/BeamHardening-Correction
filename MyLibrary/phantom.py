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

def phantom (N):
	"""
	Create a modified Shepp-Logan phantom of size NxN.
	
	"""
	
	ellipses = [
            [   1,   .69,   .92,    0,      0,   0],
	        [-.80, .6624, .8740,    0, -.0184,   0],
	        [-.20, .1100, .3100,  .22,      0, -18],
	        [-.20, .1600, .4100, -.22,      0,  18],
	        [ .10, .2100, .2500,    0,    .35,   0],
	        [ .10, .0460, .0460,    0,     .1,   0],
	        [ .10, .0460, .0460,    0,    -.1,   0],
	        [ .10, .0460, .0230, -.08,  -.605,   0],
	        [ .10, .0230, .0230,    0,  -.606,   0],
	        [ .10, .0230, .0460,  .06,  -.605,   0]]
	
	ph = np.zeros ((N,N))
	ygrid, xgrid = np.mgrid[-1:1:(1j*N), -1:1:(1j*N)]

	for ell in ellipses:
		I   = ell[0]
		a2  = ell[1]**2
		b2  = ell[2]**2
		x0  = ell[3]
		y0  = ell[4]
		phi = ell[5] * np.pi / 180  
		
		cos_p = np.cos (phi) 
		sin_p = np.sin (phi)
		
		ph [((((xgrid - x0) * cos_p + (ygrid - y0) * sin_p)**2) / a2 
              + (((ygrid - y0) * cos_p - (xgrid - x0) * sin_p)**2) / b2) <= 1] += I

	return ph




