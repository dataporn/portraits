#=
B_matrix.jl - Calculates complex network portraits

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
    
    If this software is used in an article, an acknowledgment would be 
    awesome, as would an email with the article.  Please cite as:
    J. P. Bagrow et al 2008 EPL 81 68004
    
    References:
    doi: 10.1209/0295-5075/81/68004
    http://arxiv.org/abs/cond-mat/0703470v2
    http://people.clarkson.edu/~qd00/B_matrix_site/

Jim Bagrow, 2008-04-21
bagrowjp [at] gmail [dot] com

@DataPornStar, 2014-05 (julia version)

=#

using Graphs
using DataStructures

typealias BMatrix Matrix

function portrait(g::GenericGraph)
	# return matrix where M[i,j] is the number of starting nodes in G
	# with j nodes in shell i.

	# dia = diameter(g)
    dia = 500
	N = num_vertices(g)
	B = zeros(dia + 1, N) 
	
	max_path = 1
	for v in vertices(g)
		distances = gdistances(g, v)

		# build individual distribution:
		distribution = counter(distances)
		for (shell, count) in distribution
            if shell !== -1
			    B[shell+1,count] += 1
            end
        end
			
		max_distance = maximum(distances)
		if max_distance > max_path
			max_path = max_distance
        end

		# HACK: count starting nodes that have zero nodes in farther shells
		max_shell = dia
		while max_shell > max_distance
			B[max_shell,1] += 1
			max_shell -= 1
        end
    end
	return B
end

function distance(mat1::BMatrix, mat2::BMatrix)
    
	# Distance between B-Matrix mat1 and mat2 

	ns, ms = size(mat1)
	nl, ml = size(mat2)
	size_to = maximum([ns, nl])

	# pad smaller matrix with rows of zeros 
	# except column 'zero' is # of nodes:
	if ns < nl
		mat1 = vcat(
			mat1,
			hcat(
				mat1[0,1] * ones(size_to - ns, 1),
              	zeros(size_to - ns, ms - 1)
			)
		)
	else
		mat2 = vcat(
			mat2, 
			hcat(
				mat2[1,1] * ones(size_to - nl,1), 
				zeros(size_to - nl, ml - 1)
			)
		)
    end

	# pad columns with zeros, to align CDFs:
	if ms < ml
		mat1 = hcat(mat1, zeros(size_to, ml - ms))
	else
		mat2 = hcat(mat2, zeros(size_to, ms - ml))
    end
	# Get row-wise test statistic:
	K = zeros(1, size_to)

	C1 = bcdf(mat1)
	C2 = bcdf(mat2)
	for i in 1:size_to
        K[i] = maximum(abs(C1[i:end] - C2[i:end]))
    end
	# shell weights:
	b = cumsum(mat1, 1) + cumsum(mat2, 1)
	return 1 #dot(ctranspose(b),K) / cumsum(b)
end

function bcdf(B::Matrix)
	# compute the matrix of cumulative distributions of B
	n, m = size(B)
	C = zeros(n,m)

	row_sum = cumsum(B, 1)
	for i in 1:n
		if row_sum[i] > 0
            C[i:end] = cumsum(B[i:end]) / row_sum[i]
		else 
            C[i:end] = ones(1,m)
        end
    end

	return C
end

function main()
	G1 = simple_graph(3)
	G2 = simple_graph(3)
    add_edge!(G1, 1, 2)
    add_edge!(G2, 1, 2)
	B1 = portrait(G1)
	B2 = portrait(G2)

	D = distance(B1, B2)
	println(D)
end

main()
