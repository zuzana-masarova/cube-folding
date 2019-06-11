Code for FOLDING POLYOMINOES WITH HOLES INTO A CUBE

This code is supplementary material for paper ‘Folding Polyominoes with Holes into a Cube’ by Aichholzer, O., Akitaya, H. A., Cheung, K. C., Demaine, E. D., Demaine, M. L., Fekete, S. P., Kleist L., Kostitsyna, I., Löffler, M., Masárová, Z., Mundilova, K., Schmidt, C., CCCG (2019)

WHAT THE CODE DOES:
===================
Given a rectangular grid P in which some squares and slits along grid edges were cut out, the program checks whether such polyomino could fold into a unit cube C (a necessary condition for folding). 
More precisely, when folding into the unit cube, we only allow folds along the grid edges, by 90 or 180 degrees, and we require that all faces of the unit cube get covered by some polyomino square. The code checks whether there exists a mapping between all vertices of squares of P to vertices of C such that every pair of adjacent polyomino squares gets mapped consistently. (This means that, for every adjacent pair of squares s and s’ in P, if vertices of s have been mapped to some 4 vertices of C, then there are two possibilities how to map vertices of s’ onto C: the two vertices shared by s and s’ must be mapped consistently and for the other two vertices of s’ there are two options depending on whether s’ is folded at 90 degree angle to an adjacent face of C, or whether it is folded at 180 degree to the same face of C.)
The code can be used to show non-foldability for small polyominoes: if no consistent mapping exists for a polyomino, then the polyomino cannot be folded onto C. On the other hand, any consistent vertex mapping covering all faces of C obtained by the algorithm that we tried could in practice be turned into a folding. However, we have not been able to prove that this is always the case.

The code
1) Runs a breadth-first-search on the polyomino squares, starting with the leftmost square in the top row of P and continue via adjacent squares. This produces a numbering of polyomino squares in which each but the first square is adjacent to at least one square with smaller number. 
2) Maps vertices of the first square to the bottom face of C. Then extends the mapping one square at a time according to the numbering in up to two consistent ways. All such partial mappings are being tracked.
3) Once all vertices of polyomino’s squares have been consistently mapped to vertices of the unit cube, the code checks whether all faces of the unit cube have been covered and outputs all the consistent mappings that do so.


HOW TO USE THE CODE WITH YOUR OWN POLYOMINOES:
==============================================
Input of the code = polyomino P whose foldability we want to check
Output = all consistent mappings of vertices of P onto the unit cube C and their number

1) How to specify the input polyomino: 
Specify parameter dim_a, dim_b, holes, slits
dim_a = number of rows in polyomino P
dim_b = number of columns in P
Thus P is a rectangular grid (dim_a)x(dim_b) and the following specifies the positions of holes and slits:
holes = list of ordered pairs specifying the coordinates of the cut-out squares from the grid (numbered from 0)
e.g. holes = [[1,1],] means that the square  with coordinates [1,1] in P is a hole
slits = list of square coordinates in the grid whose boundary contains slits, and numbers of square sides where the slit is (down=0,right=1,up=2,left=3)
e.g. slits = [[1,0,0],[2,3,1,2]] means that the grid square with coordinates [1,0] has a slit on side 0 (=down); and square [2,3] has slits on sides 1,2 (right&up)

2) How to interpret the code output:
‘Found valid labelling’ gives the number of found consistent mappings of vertices of P onto the cube and all mappings are listed.

Each consistent mapping consists of (dim_a+1)x(dim_b+1) array A and a dictionary D:
The array A gives the vertex mappings of the original (dim_a)x(dim_b) polyomino. 
e.g. A[2,1] = 100. means that the polyomino vertex with coordinates [2,1] (counting from 0) got mapped to cube vertex with coordinates [1,0,0] in all its incident polyomino squares.
Sometimes it may happen that a vertex v of the polyomino gets mapped to multiple different cube vertices from its multiple incident squares (e.g. if there was a slit through v), in that case v is said to be a special vertex, its label in A is 5. and for its mapped cube vertices see the dictionary D.
All possible labels in A:
0. = coordinates [0,0,0]
1. = [0,0,1]
10. = [0,1,0]
11. = [0,1,1]
100. = [1,0,0]
101. = [1,0,1]
110. = [1,1,0]
111. = [1,1,1]
5. = special vertex, for its labels see the dictionary D
-1 = vertex left unassigned (e.g. there was a hole)

The dictionary D gives the mappings for ’special’ polyomino vertices. 
Keys in D = coordinates of vertex v in array A
Values in D = 2x2 array giving the coordinates of cube vertices to which v was assigned in all its incident polyomino squares
e.g. (2, 2): array([[  -1.,  101.],
       	    [   0.,   -1.]])
means that vertex (2,2) of the polyomino P was left unassigned in polyomino squares (1,1) and (2,2); and it got mapped to cube vertex [1,0,1] in the polyomino square (1,2); and to the cube vertex [0,0,0] in the polyomino square (2,1).