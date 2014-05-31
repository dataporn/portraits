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

#=
function diameter(g::GenericGraph)
    lengths = Array(Int,0)
    weights = zeros(num_edges(g))
    for source in vertices(g)
        for target in vertices(g)
            sp = shortest_path(g, weights, source, target)
            push!(lengths, length(sp))
            println(lengths)
        end
    end
    return maximum(lengths)
end
=#

diameter(g::GenericGraph) = 5

# return matrix where M[i,j] is the number of starting nodes in G
# with j nodes in shell i.
function portrait(g::GenericGraph)

    dia = diameter(g)
	bmatrix = zeros(
        dia + 1, 
        num_vertices(g)
    )
	
    max_path = 1
	for v in vertices(g)
		distances = sort(filter(n -> n !== -1, gdistances(g, v)))

        max_node_distances = maximum(distances)
        curr_max_path = max_node_distances
        if curr_max_path > max_path
            max_path = curr_max_path
        end

		# build individual distribution
		distribution = counter(distances)
		for (shell, count) in distribution
		    bmatrix[shell + 1, count + 1] += 1
        end
			
        max_shell = dia 
        while max_shell > max_node_distances
            bmatrix[max_shell + 1, 1] += 1
            max_shell -= 1
        end
    end
    return bmatrix[1:max_path+1,:]
end

# Distance between B-Matrix mat1 and mat2 
function distance(mat1::BMatrix, mat2::BMatrix)

	ns, ms = size(mat1)
	nl, ml = size(mat2)
	size_to = maximum([ns, nl])

	# pad smaller matrix with rows of zeros 
	# except column 'zero' is # of nodes:
	if ns < nl
		mat1 = vcat(
			mat1,
			hcat(
				mat1[1,2] * ones(size_to - ns, 1),
              	zeros(size_to - ns, ms - 1)
			)
		)
	else
		mat2 = vcat(
			mat2, 
			hcat(
				mat2[1,2] * ones(size_to - nl,1), 
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
	K = zeros(size_to,1)

	C1 = bcdf(mat1)
	C2 = bcdf(mat2)
	for i in 1:size_to
        K[i] = maximum( abs( C1[i,:] - C2[i,:] ) )
    end
	# shell weights:
	b = sum(mat1, 2) + sum(mat2, 2)
    D = b' * K / sum(b)
    return D[1]
end

# compute the matrix of cumulative distributions of B
function bcdf(B::Matrix)
	n, m = size(B)
	C = zeros(n,m)

	row_sum = sum(B, 2)
	for i in 1:n
		if row_sum[i] > 0
            C[i,:] = cumsum(B[i,:]) ./ row_sum[i]
		else 
            C[i,:] = ones(1,m)
        end
    end

	return C
end
