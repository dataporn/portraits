#!/usr/bin/env python
# encoding: utf-8

# B_matrix_igraph.py
# Jim Bagrow 
# @DataPornStar (conversion from networkx to igraph)
# Last Modified: 2013-07-03

"""B_matrix.py - Calculates complex network portraits

COPYRIGHT:
    Copyright (C) 2008 Jim Bagrow
    
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
    
    See http://www.gnu.org/licenses/gpl.txt for more details.

ABOUT:
    Plot complex networks portraits, requires python 2.7,
    igraph, and (optionally) matplotlib/pylab for plotting
    
    If this software is used in an article, an acknowledgment would be 
    awesome, as would an email with the article.  Please cite as:
    J. P. Bagrow et al 2008 EPL 81 68004
    
    Dependencies:
    http://www.python.org/
    http://igraph.sourceforge.net/
    http://matplotlib.sourceforge.net/
    
    References:
    doi: 10.1209/0295-5075/81/68004
    http://arxiv.org/abs/cond-mat/0703470v2
    http://people.clarkson.edu/~qd00/B_matrix_site/

USAGE:
    python B_matrix.py input_edgelist.txt output_matrix.txt

Jim Bagrow, 2008-04-21
bagrowjp [at] gmail [dot] com
"""

import sys, os, igraph
from math import log
try:
	import numpy
except:
	import scipy_base as numpy

def elementWiseLog(mat):
	""" Take log of each element+1 in matrix, the +1
	keeps 0 from being a problem.  
	"""
	new_mat = zeros( mat.shape, tCode=float )
	i = 0
	for row in mat:
		j = 0
		for e in row:
			if e !=0:
				new_mat[i,j] = log( e+1 )
			else:
				new_mat[i,j] = 0
			j += 1
		i += 1
	return new_mat


def zeros( shape, tCode=None):
	try:
		return numpy.zeros(shape,dtype=tCode)
	except TypeError:
		return numpy.zeros(shape,typecode='fd') # hardwired to float


def fileMat(fileName, S=None):
	"""Read and write matrices to file at fileName
	if S=None, read S from fileName else write S.
	"""
	if S != None:
		f = open(fileName, 'w')
		for row in S:
			for el in row:
				print >>f, el,
			print >>f
		f.close()
	else:
		f = open(fileName, 'r')
		S = []
		for row in f.readlines():
			S.append( [float(i) for i in row.split(" ")] )
		return numpy.array(S)

def plotMatrix(o_mat, **kwargs):
	""" DOC STRING
	"""
	kwargs['interpolation']='nearest'
	origin = kwargs.get('origin',1); kwargs['origin']='lower'
	showColorBar = kwargs.get('showColorBar',False)
	if kwargs.has_key("showColorBar"): kwargs.pop("showColorBar")
	logColors    = kwargs.get('logColors',False)
	if kwargs.has_key("logColors"):    kwargs.pop("logColors")
	ifShow       = kwargs.get('show',False)
	if kwargs.has_key("show"):         kwargs.pop("show")
	fileName     = kwargs.get('fileName',None)
	if kwargs.has_key("fileName"):     kwargs.pop("fileName")

	mat = o_mat.copy() # don't modify original matrix
	if logColors: mat = elementWiseLog(mat)

	if not kwargs.has_key("vmax"):
		kwargs['vmax'] = float(mat[origin:,origin:].max())

	ax = pylab.axes()#[.05,.05,.9,.9])
	ax.xaxis.tick_top()
	H = pylab.imshow( mat, **kwargs)
	pylab.axis('tight')	
	ax.set_xlim( (origin,mat.shape[1]) )
	ax.set_ylim( (mat.shape[0],origin) )

	if showColorBar: pylab.colorbar()

	if fileName != None:
		pylab.savefig(fileName)

	if ifShow: pylab.show()
	else: pylab.draw_if_interactive()

	return H

def portrait(G):
	""" return matrix where M[i][j] is the number of starting nodes in G
	with j nodes in shell i.
	"""

	from collections import Counter

	dia = G.diameter(directed=False)
	N = G.vcount()
	# B indices are 0...dia x 0...N-1:
	B = zeros( (dia+1,N) ) 
	
	max_path = 1
	for node in G.vs.indices:
		distances = [ d for n,d,p in G.bfsiter(node, mode='ALL', advanced=True) ]

		# build individual distribution:
		distribution = Counter(distances)
		for shell,count in distribution.iteritems():
			B[shell][count] += 1
			
		max_distance = max(distances)
		if max_distance > max_path:
			max_path = max_distance

		# HACK: count starting nodes that have zero nodes in farther shells
		max_shell = dia
		while max_shell > max_distance:
			B[max_shell][0] += 1
			max_shell -= 1
	
	return B[:max_path+1,:]

# based on B_Distance.m
def distance(mat1, mat2):
	""" Distance between B-Matrix mat1 and mat2 """

	ns, ms = mat1.shape
	nl, ml = mat2.shape
	size_to = max([ns,nl])

	# pad smaller matrix with rows of zeros 
	# except column 'zero' is # of nodes:
	if ns < nl:
		mat1 = numpy.vstack((
			mat1,
			numpy.hstack((
				mat1[0,1] * numpy.ones((size_to - ns,1)),
              	zeros((size_to - ns, ms - 1))
			))
		))
	else:
		mat2 = numpy.vstack((
			mat2, 
			numpy.hstack((
				mat2[0,1] * numpy.ones((size_to - nl,1)), 
				zeros((size_to - nl, ml - 1))
			))
		))

	# pad columns with zeros, to align CDFs:
	if ms < ml:
		mat1 = numpy.hstack((mat1, zeros((size_to, ml - ms))))
	else:
		mat2 = numpy.hstack((mat2, zeros((size_to, ms - ml))))

	# Get row-wise test statistic:
	K = zeros((size_to,1))

	C1 = bcdf(mat1)
	C2 = bcdf(mat2)
	for i in range(1,size_to):
		K[i] = max(abs(C1[i,] - C2[i,]))

	# shell weights:
	b = numpy.sum(mat1, axis = 1) + numpy.sum(mat2, axis = 1)
	return numpy.dot(b.conj().transpose(),K) / numpy.sum(b)
		
# based on Bcdf.m
def bcdf(B):
	""" compute the matrix of cumulative distributions of B """
	n,m = numpy.shape(B)
	C = zeros((n,m))

	row_sum = numpy.sum(B, axis = 1)
	for i in range(0,n):
		if row_sum[i].any() > 0:
			C[i,] = B[i,].cumsum() / row_sum[i]
		else: # avoid this by including column 'zero' in B:
			C[i,] = numpy.ones((1,m)); # ones or zeros or 0.5's?

	return C

if __name__ == '__main__':
	G1 = igraph.Graph.Read_Edgelist(sys.argv[1])
	G2 = igraph.Graph.Read_Edgelist(sys.argv[2])
	B1 = portrait(G1)
	B2 = portrait(G2)

	D = distance(B1, B2)
	print D

	sys.exit()

	try: # plot the portrait with pylab, but I prefer matlab:
		import pylab
		plotMatrix(B, origin=1, logColors=True, show=True)
	except ImportError:
		print "pylab failed, no plotting"
	try: 
		print "writing matrix to file...", sys.argv[2]
		fileMat(sys.argv[2], B)
	except:
		pass
