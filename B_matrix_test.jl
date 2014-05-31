# Test of bmatrix generation and distance

using Graphs
using Base.Test
include("B_matrix.jl")

# g1: the example in CLRS (2nd Ed.)
g1 = simple_graph(5)

g1_wedges = [
    (1, 2, 10.),
    (1, 3, 5.),
    (2, 3, 2.),
    (3, 2, 3.),
    (2, 4, 1.),
    (3, 5, 2.),
    (4, 5, 4.),
    (5, 4, 6.),
    (5, 1, 7.),
    (3, 4, 9.) ]

ne = length(g1_wedges)
eweights1 = zeros(ne)
for i = 1 : ne
    we = g1_wedges[i]
    add_edge!(g1, we[1], we[2])
    eweights1[i] = we[3]
end

@assert num_vertices(g1) == 5
@assert num_edges(g1) == 10
@test diameter(g1) == 5

b1 = portrait(g1)
@test b1 == [
    0.  5.  0.  0.  0.
    0.  1.  3.  1.  0.
    0.  3.  2.  0.  0.
    3.  1.  1.  0.  0.
]

# g2: the example in Wikipedia
g2 = simple_graph(6, is_directed=false)

g2_wedges = [
    (5, 6, 9.),
    (5, 4, 6.),
    (6, 3, 2.),
    (4, 3, 11.),
    (6, 1, 14.),
    (3, 1, 9.),
    (3, 2, 10.),
    (4, 2, 15.),
    (1, 2, 7.) ]

ne = length(g2_wedges)
eweights2 = zeros(ne)
for i = 1 : ne
    we = g2_wedges[i]
    add_edge!(g2, we[1], we[2])
    eweights2[i] = we[3]
end
    
@assert num_vertices(g2) == 6
@assert num_edges(g2) == 9
@test diameter(g2) == 5
b2 = portrait(g2)
@test b2 == [
    0.  6.  0.  0.  0.  0.
    0.  0.  1.  4.  1.  0.
    0.  1.  4.  1.  0.  0.
]

@test distance(b1, b2) == 0.32499999999999996

