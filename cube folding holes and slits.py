import numpy as np
import copy
import time


#dimensions of the rectangular net
dim_a = 5
dim_b = 4
#coordinates of the square holes and slits
holes = [[1,1],[2,2]] #list of ordered pairs
slits = []#[[1,0,0],[2,3,1,2]]#the first two numbers in each list are the coordinates of a
#square in grid1 and the remaining numbers denote which of its sides have slits (down=0,right=1,up=2,left=3). E.g.
#[2,2,1,2,3] means the square (2,2) in grid1 has slits on its right, upper and left side.
# NOTE: later (when queries to slit_dictionary are made) it is assumed that each item in slits has length >=3 (i.e. when
# a grid1 square is mentioned in slits, it has a slit on a side). Also, it is assumed that there are no slits along the boundary
# of grid1.
#NOTE: algorithm assumes that holes and slits do not disconnect grid1

#coordinates of the square in grid1 to start breadth first search
start_pos = [0,0]

#Create the cube net as an array and an extra (dim_a+1)(dim_b+1) grid2
#uses global variables dim_a, dim_b,holes
def create_grids():
    grid1 = np.ones((dim_a, dim_b))
    for hole in holes:
        grid1[hole[0]][hole[1]] = 0

    grid2 = np.full((3,dim_a+1,dim_b+1),-1) #3 copies of dima+1 x dim_b+1 grid, each copy stores x-,y-,z- labels
    return grid1, grid2

#Create dictionary to store the slits
#keys = tuple of coordinates of squares in grid1, values = list of legth 4, has entry 1 for the sides where the square has a slit (down=list[0],right=list[1],up=list[2],left=list[3])
#uses global variable slits
def create_dictionary_for_slits():
    slit_dictionary = {}
    for item in slits:
        for slit_coor in range(2,len(item)):
            #add 1 in the list in dictionary corresponding to the relevant grid1 square
            if tuple([item[0],item[1]]) not in slit_dictionary.keys():
                slit_dictionary[tuple([item[0],item[1]])] = [0]*4 #list of length 4 full of 0s
            slit_dictionary[tuple([item[0], item[1]])][item[slit_coor]] = 1
            #add one also to the neighbouring square:
            if item[slit_coor] == 0: #then square below has a slit on upper side (2)
                if tuple([item[0]+1, item[1]]) not in slit_dictionary.keys():
                    slit_dictionary[tuple([item[0]+1, item[1]])] = [0] * 4  # list of length 4 full of 0s
                slit_dictionary[tuple([item[0]+1, item[1]])][2] = 1
            if item[slit_coor] == 1:  # then square to the right has a slit on left side (3)
                if tuple([item[0], item[1]+1]) not in slit_dictionary.keys():
                    slit_dictionary[tuple([item[0], item[1]+1])] = [0] * 4  # list of length 4 full of 0s
                slit_dictionary[tuple([item[0], item[1]+1])][3] = 1
            if item[slit_coor] == 2:  # then square above has a slit on lower side (0)
                if tuple([item[0]-1, item[1]]) not in slit_dictionary.keys():
                    slit_dictionary[tuple([item[0]-1, item[1]])] = [0] * 4  # list of length 4 full of 0s
                slit_dictionary[tuple([item[0]-1, item[1]])][0] = 1
            if item[slit_coor] == 3:  # then square left has a slit on right side (1)
                if tuple([item[0], item[1]-1]) not in slit_dictionary.keys():
                    slit_dictionary[tuple([item[0], item[1]-1])] = [0] * 4  # list of length 4 full of 0s
                slit_dictionary[tuple([item[0], item[1]-1])][1] = 1
    return slit_dictionary

#For all vertices, fill in '5' as their x,y,z label in grid2 and create a dictionary of 3x2x2 arrays. Key in dictionary
# = tuple(coordinates of vertex in grid2 (=a list of length 2)), value = 3x2x2 np array to store the 4 possibly different
#labels of the x-,y-,z-coordinate of special vertex (at the beginning all are set -1, the non-existent behind boundary squares of 2x2 grid are set to 100)
#uses global variables dim_a,dim_b,grid2
def create_dictionary_for_vertices():
    vertex_dictionary = {}
    for i in range(dim_a+1):
        for j in range(dim_b+1):
            vertex = [i,j]

            for k in range(3): #i is x-,y-,z-label
                grid2[k][vertex[0]][vertex[1]]=5
            vertex_dictionary[tuple(vertex)] = np.full((3,2,2),-1)
            if i==0: #if vertex [i,j] is on top boundary of grid2, the top two squares of 2x2 grid contain 100
                for l in range(3):
                    vertex_dictionary[tuple(vertex)][l][0][0]=100
                    vertex_dictionary[tuple(vertex)][l][0][1]=100
            if i==dim_a: #if vertex [i,j] is on bottom boundary of grid2, the bottom two squares of 2x2 grid contain 100
                for l in range(3):
                    vertex_dictionary[tuple(vertex)][l][1][0]=100
                    vertex_dictionary[tuple(vertex)][l][1][1]=100
            if j==0: #if vertex [i,j] is on left boundary of grid2, the left two squares of 2x2 grid contain 100
                for l in range(3):
                    vertex_dictionary[tuple(vertex)][l][0][0]=100
                    vertex_dictionary[tuple(vertex)][l][1][0]=100
            if j==dim_b: #if vertex [i,j] is on right boundary of grid2, the right two squares of 2x2 grid contain 100
                for l in range(3):
                    vertex_dictionary[tuple(vertex)][l][0][1]=100
                    vertex_dictionary[tuple(vertex)][l][1][1]=100
    return vertex_dictionary

#Find neighbouring squares of a square in grid1
#Input: square is a list of 2 coordinates
#uses global variable grid1
def find_neighbours(square):
    s=square
    if s[0]<0 or s[0]>dim_a-1 or s[1]<0 or s[1]>dim_b-1:
        print 'Function find_neighbours: GIVEN SQUARE IS OUTSIDE OF GRID1 \n'
    neighbours = [[s[0]+1,s[1]],[s[0],s[1]+1],[s[0]-1,s[1]],[s[0],s[1]-1]] #list of coordinates of squares neighbouring the given square in grid1
    #print 'neighbours at start', neighbours
    for side in range(3,-1,-1): #0=bottom, 1=right, 2=up, 3=left
        #if square s has slit on that side, delete the neighbour from the list
        if tuple(s) in slit_dictionary.keys() and slit_dictionary[tuple(s)][side]==1:
            neighbours.pop(side)
    for square in list(neighbours):
        #print 'square to examine',square
        if square[0]<0 or square[0]>dim_a-1 or square[1]<0 or square[1]>dim_b-1:
            neighbours.pop(neighbours.index(square))
        #print 'neighbours after examining square', neighbours
        if square[0]>-1 and square[0]<dim_a and square[1]>-1 and square[1]<dim_b:
            #print square[0],square[1]
            if grid1[square[0]][square[1]] == 0: #if the neighbouring square is a hole
                #print 'square is a hole'
                neighbours.pop(neighbours.index(square))
    return neighbours # a list of pairs

#breadth first search on grid1
def search_net(start_pos):
    grid = np.full((dim_a,dim_b),-1)
    grid[start_pos[0]][start_pos[1]] = 0
    queue = find_neighbours(start_pos)
    for square in queue:
        grid[square[0]][square[1]] = 1
    step = 1
    while len(queue)>0:
        step = step+1
        from_square = queue.pop(0)
        additional_items_into_queue = []
        for nbh_square in find_neighbours(from_square):
            if grid[nbh_square[0]][nbh_square[1]] == -1:
                additional_items_into_queue.append(nbh_square)
                grid[nbh_square[0]][nbh_square[1]] = step
        queue = queue + additional_items_into_queue
    return grid

#vertex coordinates of a given square in grid2
#Input: [a,b] coordinates of a square in grid1
#Output: [down_left_xy,down_right_xy,up_right_xy,up_left_xy] each element is a list of length 2
def coord_for_sq(square):
    return [[square[0]+1,square[1]],[square[0]+1,square[1]+1],[square[0],square[1]+1],[square[0],square[1]]]

#Given a position of a vertex in a grid1 square (0 = bottom left, 1 = bottom right, 2 = top right, 3 = top left), find which
#position of a 2x2 array (in vertex_dictionary) it corresponds to
#Input: number 0,1,2,or 3
def from_vertex_pos_in_square_to_array_pos(vertex_pos_in_square):
    if vertex_pos_in_square == 0:
        dict_array = [0, 1]
    if vertex_pos_in_square == 1:
        dict_array = [0, 0]
    if vertex_pos_in_square == 2:
        dict_array = [1, 0]
    if vertex_pos_in_square == 3:
        dict_array = [1, 1]
    return dict_array

#Given a position in a 2x2 array (in vertex_dictionary), find which position of a vertex in a grid1 square (0 = bottom left,
# 1 = bottom right, 2 = top right, 3 = top left) it corresponds to
#Input:[0,0],[1,0],[0,1] or [1,1] (positions in 2x2 array)
#Output: number 0,1,2,or 3
def from_array_pos_to_vertex_pos_in_square(array_pos):
    if array_pos == [0,0]:
        vertex_pos_in_square = 1
    if array_pos == [0,1]:
        vertex_pos_in_square = 0
    if array_pos == [1,0]:
        vertex_pos_in_square = 2
    if array_pos == [1,1]:
        vertex_pos_in_square = 3
    return vertex_pos_in_square

#Extract the x,y,z coordinates of a vertex in vertex_dictionary
#Input: vertex_dictionary (storing some partial labelling), vertex=
#[a,b] a=row, b=column of 'grid2', vertex_pos_in_square = which of the 4 square vertices it is = 0 <-> vertex is the bottom left of the square, i.e. in dictionary want up right;
# 1 <-> vertex is the bottom right corner, i.e. in dictionary array look up left, 2 <-> vertex is the top right corner,
# i.e. in dictionary look down left, 3 <-> vertex is the top left corner, i.e. in dictionary look down right
#Output: a list of length 3 (x,y,z label of the vertex)
def label_of_vertex_in_vertex_dictionary(vertex_dictionary,vertex,vertex_pos_in_square): #vertex=[a,b] a=row, b=column of 'grid2'
    dict_array = from_vertex_pos_in_square_to_array_pos(vertex_pos_in_square)
    label = [vertex_dictionary[tuple(vertex)][0][dict_array[0]][dict_array[1]],
             vertex_dictionary[tuple(vertex)][1][dict_array[0]][dict_array[1]],
             vertex_dictionary[tuple(vertex)][2][dict_array[0]][dict_array[1]]]
    return label

#Given a square in grid1 (as a pair of its coordinates), and one of square's vertices (as a pair of its coordinates in grid2),
#it determinates which square vertex the given vertex is.
#Output: 0 = bottom left, 1 = bottom right, 2 = top right, 3 = top left
def find_vertex_position_in_square(square,vertex):
    square_vertices = coord_for_sq(square)
    #print 'square_vertices',square_vertices
    return square_vertices.index(vertex)

#Having a labelled square abcd (a=bottom left, then continue ccw, b=bottom right, etc.) and assuming that the square
#dcfe (i.e. the square 'above' abcd) should be at 90 degrees to abcd, the labels of dcfe can be determined as follows:
#compute label separately for each coordinate (x-,y-,z-axis). [x_a,x_b] are the respective coordinate values of vertices
# a and b, [x_c,x_d] are coordinates of c and d.
#(Example: a = 000, b=100, c=110, d=010, and computing the x-coordinate of dcfe. Clearly, x_d=0, x_c=1. Remaining two
# x-coordinates [x_f,x_e] are determined from P1=[x_a,x_b]=[0,1] and P2=[x_c,x_d]=[1,0] by conditions below to be
# [1,0], i.e. x_f=1 and x_e=0. Repeat to obtain y- and z-coordinates of vertices e and f.)
#Input: P1=[x_a,x_b],P2=[x_c,x_d].
#Output pair3=[x_3,x_4] (top right, top left of the new square)
def deduce_coordinates(p1,p2):#p1=[x_a,x_b], p2=[x_c,x_d].
    p_new=list()
    if p1 == [0,0] and p2 == [1,1]:
        p_new = [1,1]
    if p1 == [1,1] and p2 == [0,0]:
        p_new = [0,0]
    if p1 == [1,1] and p2 == [1,1]:
        p_new = [0,0]
    if p1 == [0,0] and p2 == [0,0]:
        p_new = [1,1]
    if p1 == [0,1] and p2 == [1,0]:
        p_new = [1,0]
    if p1 == [1,0] and p2 == [0,1]:
        p_new = [0,1]
    return p_new


#Label neighbour that is up, given a labelled square below
#Input: list of length 4, each element (vertex labels) is a list of x-, y-, z-axes labels (order of vertices: down left, down right, up right, up left)
#Output: 2 possible labellings (one assuming the new square is at 180 degrees fold to the old square, one assuming it is
# at 90 degrees)
def label_up(from_sq_coord):
    lbling1 = [from_sq_coord[3],from_sq_coord[2],from_sq_coord[1],from_sq_coord[0]] #possibility 1: 180 turn, same unit cube
    vertex3_coord = []
    vertex4_coord = []
    for i in range(3): #compute x-, y-, z- coordinate of vertex 3 and 4
        new_i_pair = deduce_coordinates([from_sq_coord[0][i],from_sq_coord[1][i]],[from_sq_coord[2][i],from_sq_coord[3][i]]) #x-coordinate of vertex 3 and of vertex 4, then y-coordinate of vertex 3 and vertex 4 (of the square being labelled), etc.
        vertex3_coord.append(new_i_pair[0])
        vertex4_coord.append(new_i_pair[1])
    lbling2 = [from_sq_coord[3],from_sq_coord[2],vertex3_coord,vertex4_coord]
    return [lbling1,lbling2]

#label neighbour that is down
def label_down(from_sq_coord): #input= list of length 4, each element is a list of x-, y-, z-coordinates
    lbling1 = [from_sq_coord[3],from_sq_coord[2],from_sq_coord[1],from_sq_coord[0]] #possibility 1: 180 turn, same unit cube
    vertex0_coord = []
    vertex1_coord = []
    for i in range(3): #compute x-, y-, z- coordinate of vertex 3 and 4
        new_i_pair = deduce_coordinates([from_sq_coord[2][i],from_sq_coord[3][i]],[from_sq_coord[0][i],from_sq_coord[1][i]]) #x-coordinate of vertex 3 and of vertex 4, then y-coordinate of vertex 3 and vertex 4 (of the square being labelled), etc.
        vertex0_coord.append(new_i_pair[0])
        vertex1_coord.append(new_i_pair[1])
    lbling2 = [vertex0_coord, vertex1_coord, from_sq_coord[1],from_sq_coord[0]]
    return [lbling1,lbling2]

def label_right(from_sq_coord): #input= list of length 4, each element is a list of x-, y-, z-coordinates
    lbling1 = [from_sq_coord[1],from_sq_coord[0],from_sq_coord[3],from_sq_coord[2]] #possibility 1: 180 turn, same unit cube
    vertex1_coord = []
    vertex2_coord = []
    for i in range(3): #compute x-, y-, z- coordinate of vertex 3 and 4
        new_i_pair = deduce_coordinates([from_sq_coord[3][i],from_sq_coord[0][i]],[from_sq_coord[1][i],from_sq_coord[2][i]]) #x-coordinate of vertex 3 and of vertex 4, then y-coordinate of vertex 3 and vertex 4 (of the square being labelled), etc.
        vertex1_coord.append(new_i_pair[0])
        vertex2_coord.append(new_i_pair[1])
    lbling2 = [from_sq_coord[1],vertex1_coord,vertex2_coord,from_sq_coord[2]]
    return [lbling1,lbling2]

def label_left(from_sq_coord): #input= list of length 4, each element is a list of x-, y-, z-coordinates
    lbling1 = [from_sq_coord[1],from_sq_coord[0],from_sq_coord[3],from_sq_coord[2]] #possibility 1: 180 turn, same unit cube
    vertex0_coord = []
    vertex3_coord = []
    for i in range(3): #compute x-, y-, z- coordinate of vertex 3 and 4
        new_i_pair = deduce_coordinates([from_sq_coord[1][i],from_sq_coord[2][i]],[from_sq_coord[3][i],from_sq_coord[0][i]]) #x-coordinate of vertex 3 and of vertex 4, then y-coordinate of vertex 3 and vertex 4 (of the square being labelled), etc.
        vertex3_coord.append(new_i_pair[0]) #vertex3 and vertex 0 in different order on purpose!
        vertex0_coord.append(new_i_pair[1])
    lbling2 = [vertex0_coord,from_sq_coord[0],from_sq_coord[3],vertex3_coord]
    return [lbling1,lbling2]

#Given a square of a 2x2 array, it finds its two neighbours
#Input: one of [0,0],[0,1],[1,0],[1,1]
#Output: list of length 2, each item is a list of two
def neighbours_in_2x2_array(position):
    return [[(position[0]+1)%2,position[1]],[position[0],(position[1]+1)%2]]

#Given a vertex in grid2 and its position in a square of grid1, compute the square in grid1
#Input: vertex = list of length 2 (=coordinates in grid2), vertex_pos_in_square = which of the 4 square vertices it is (see function below)
#Output: square in grid1 (list of length 2)
def square_from_vertex_and_its_position_in_square(vertex, vertex_pos_in_square):
    if vertex_pos_in_square == 3:
        square = vertex
    if vertex_pos_in_square == 2:
        square = [vertex[0],vertex[1]-1]
    if vertex_pos_in_square == 1:
        square = [vertex[0]-1,vertex[1]-1]
    if vertex_pos_in_square == 0:
        square = [vertex[0]-1,vertex[1]]
    # if square[0]<0 or square[0]>dim_a-1 or square[1]<0 or square[1]>dim_b-1:
    #     print 'Function square_from_vertex_and_its_position_in_square: THIS SQUARE IS OUT OF BOUNDARY OF GRID1'
    return square


#Insert x,y,z label of a vertex into vertex_dictionary - to the corresponding square of 2x2 array and to all other of its
# squares that share a side with it so have to have the same vertex label
#Input: vertex = list of length 2 (coordinates in grid2), vertex_pos_in_square = which of the 4 square vertices it is = 0 <-> vertex is the bottom left of the square, i.e. in dictionary want up right;
# 1 <-> vertex is the bottom right corner, i.e. in dictionary array look up left, 2 <-> vertex is the top right corner,
# i.e. in dictionary look down left, 3 <-> vertex is the top left corner, i.e. in dictionary look down right
#label=list of length 3 (x-,y-,z-label to insert)
def insert_vertex_label_into_vertex_dictionary(vertex_dictionary,vertex,vertex_pos_in_square,label):
    dict_array = from_vertex_pos_in_square_to_array_pos(vertex_pos_in_square)#vertex_pos_in_square => dict_array = position in 2x2 array
    for i in range(3): #i as x-,y-,z-coordinate
        vertex_dictionary[tuple(vertex)][i][dict_array[0]][dict_array[1]]=label[i]

    #propagate the label also to other squares of 2x2 grid that share a side with this square
    #idea: square, dict_array => neighbour_array => neighbour_pos_in_square   +vertex => neighbour_square. Then test whether square
    # and neighbour_square are indeed neighbours in grid1. If so, insert the label in neighbour_array, set square to be
    # the neighbour_square and repeat the procedure
    square = square_from_vertex_and_its_position_in_square(vertex,vertex_pos_in_square) #vertex + vertex_pos_in_square => square in grid1
    neighbrs_in_dict_array = neighbours_in_2x2_array(dict_array) #finds neighbouring squares in 2x2 array of a given square 'dict_array' in the array
    queue_squares = []
    queue_dict_array = []
    # after repeating at most 3 times, if not blocked by slits/holes, we will reach all squares of the array sharing that vertex.
    for j in range(3):
        for nbr_array in neighbrs_in_dict_array:
            nbr_pos_in_square = from_array_pos_to_vertex_pos_in_square(nbr_array)
            nbh_square = square_from_vertex_and_its_position_in_square(vertex,nbr_pos_in_square)
            if nbh_square in find_neighbours(square) and vertex_dictionary[tuple(vertex)][0][nbr_array[0]][nbr_array[1]] != label[0]:
                for i in range(3):
                    vertex_dictionary[tuple(vertex)][i][nbr_array[0]][nbr_array[1]] = label[i]
                queue_squares.append(nbh_square)
                queue_dict_array.append(nbr_array)
        if len(queue_squares)>0: #If ==0, because no new square was found because we were blocked by slits/holes, just repeat with the current square/array, not changing anything.
            square = queue_squares.pop(0)
            neighbrs_in_dict_array = neighbours_in_2x2_array(queue_dict_array.pop(0))
    return vertex_dictionary

#Inserts the 4 labels of a grid1 square into vertex_dictionary (and spreads out the labels in the corresponding 2x2 arrays
#whenever not blocked by slits or holes)
#Input: square_to_label = 2 coordinates of square in grid1 whose labels we want to add, e.g. square [0,0], labels is a
# list of length 4, each element is a list of length 3 (x,y,z coordinates) for vertices of the square_to_label
# (downleft,downright,upright,upleft)
def insert_square_labels_into_vertex_dictionary(vertex_dictionary,square_to_label,labels):
    coord_in_grid2 = coord_for_sq(square_to_label) # list of length 4 - the 4 cells in grid2 to enter the square's labels (each item is a list of length 2), starting down left, going ccw
    for i in range(4): #ith vertex of the square to label (down left, down right, up right, up left)
        vertex_dictionary = insert_vertex_label_into_vertex_dictionary(vertex_dictionary,coord_in_grid2[i],i,labels[i])
    return vertex_dictionary

#Given a square (a list of length 2) in grid1, find in which directions it has neighbours
def get_neighbours_direction(from_square):
    neighbours = find_neighbours(from_square)
    nbh_directions = []
    for nbh in neighbours:
        if nbh[0]>from_square[0]: #down
            nbh_directions.append(0)
        if nbh[1]>from_square[1]: #right
            nbh_directions.append(1)
        if nbh[0]<from_square[0]: #up
            nbh_directions.append(2)
        if nbh[1]<from_square[1]: #left
            nbh_directions.append(3)
    return nbh_directions,neighbours

#Given a square from grid1 (=a pair [x,y]), extract its vertex labels stored in vertex_dictionary.
# Outputs a list of length 4 (bottom left, bottom right,top right,top left vertex of the square), each item is [x-label,y,z] of the vertex
def get_square_labels_from_vertex_dictionary(vertex_dictionary,square):
    vertices_in_grid2 = coord_for_sq(square)
    #print 'vertex_coord',vertices_in_grid2
    labels = []
    for i in range(4):
        labels.append(label_of_vertex_in_vertex_dictionary(vertex_dictionary,vertices_in_grid2[i],i))
    return labels

#Checks whether all vertices of a square are labelled (NOTE: a square may have all vertices labelled before it was
# itself added to a labelling - as a result of its neighbours being labelled before. However, this will still be a consistent
# labelling)
#Input: sq = coord of square in grid1
#Output: 1 if all vertices have labels, 0 otherwise
def is_square_labelled(vertex_dictionary,sq):
    labelled = 1
    square_labels =  get_square_labels_from_vertex_dictionary(vertex_dictionary,sq)
    if [-1,-1,-1] in square_labels:
        labelled = 0
    return labelled

#Given a vertex, the function checks whether its 4 labellings stored in 2x2 array vertex_dictionary[tuple(vertex)] are the same.
#If so, it also outputs that x-,y-,z-label (=list of length 3, e.g. [1,0,0]).
#If some cells store 100 because the vertex is on the boundary of grid1, these labels are ignored.
#If some cells store -1 because the vertex is on the boundary of a hole but the array has some other proper labels as well, -1s are
#ignored.
#If all array cells have only -1's (and 100s), then the function outputs labels = [-1,-1,-1]
#If there are more than 1 proper (!=100,-1) x-,y-,z-label combination, then same = 0 and labels can contain some partial label so just disregard it
#Input: vertex = coordinates of vertex in grid2
#Output: same =0/1, labels = list of length 3
def does_vertex_have_the_same_labels_in_2x2_array(vertex_dictionary,vertex):
    same = 1
    labels = []
    for i in range(3): #checking x-,y-,z-coordinate
        non_boundary_labels = [] #collects labels that are != -1, 100
        for j in range(2):
            for k in range(2):
                if vertex_dictionary[tuple(vertex)][i][j][k]!= -1 and vertex_dictionary[tuple(vertex)][i][j][k]!= 100 and \
                vertex_dictionary[tuple(vertex)][i][j][k] not in non_boundary_labels:
                    non_boundary_labels.append(vertex_dictionary[tuple(vertex)][i][j][k])
        #print 'len(non-boundary labels)',len(non_boundary_labels)
        if len(non_boundary_labels)>1: #if there are multiple proper labels for the given vertex coordinate
            same = 0
        #if there is always only one/zero proper label for the given coordinate, len(non-boundary_labels) is always 1/zero
        if len(non_boundary_labels)==1: #exactly one proper label
            labels.append(non_boundary_labels[0])
        if len(non_boundary_labels)==0: #if all labels are -1 or 100
            labels.append(-1)
    return same,labels


#Given a vertex (as list of length 2, coordinates in grid2) and a position in the 2x2 array, e.g. [0,0] (not position in
# the square!), it collects the x-,y-,z- label of the vertex stored in vertex_dictionary in that array position and outputs
# the label as an integer (e.g. 10 = 010, etc.)
def collect_label_of_vertex_into_int(vertex_dictionary,vertex,array_pos):
    label = ''
    for coor in range(3): #x-,y-,z-coordinate
        label = label + str(int(vertex_dictionary[tuple(vertex)][coor][array_pos[0]][array_pos[1]]))
    # print 'LABEL',label
    if label == '-1-1-1': #if label yet unassigned
        label = -1
    elif label == '100100100':  # if outside of boundary of grid1
        label = -1
    else:
        label = int(label)
    return label

#Given a vertex label as a list of length 3, the function concatenates it into an integer. If the label is yet unassigned
# (=[-1,-1,-1]) or outside of boundary (=[100,100,100]), the output label is set to -1
def concatenate_label_of_vertex_into_int(list_of_3_coordinates):
    label = ''
    for coor in range(3):
        label = label + str(int(list_of_3_coordinates[coor]))
    # print 'LABEL',label
    if label == '-1-1-1':
        label = -1
    elif label == '100100100':  # if outside of boundary of grid1
        label = -1
    else:
        label = int(label)
    return label

#Given a labelling, transform it into an easier to read format
#Output: grid with dimensions as grid2 and a dictionary vt_dict2.
#Idea: vertices that got in vertex_dictionary 4 same labels will have that label written in the grid. Vertices that
# received multiple labels have '5' in the grid and their labels are in vt_dict2. Its keys = coordinates of vertex in grid2,
# value = one 2x2 array (all possible positions of a vertex in a square).
# The labels (both in the grid and dict) are integers concatenating x-,y-,z-labels (e.g. 0.=000,10=010 etc), -1 = this (side of a) vertex was not labelled
# due to a hole or because it is outside of a grid1 boundary
def translate_labelling(vertex_dictionary):
    #create a new grid and a dictionary
    grid = np.full((dim_a+1,dim_b+1),-1)
    vt_dict2 = {}
    for i in range(dim_a+1):
        for j in range(dim_b+1):
            #for each vertex, if it has the same 4 labels, insert the labels into the grid, else into the dictionary
            vertex = [i,j]
            same, labels = does_vertex_have_the_same_labels_in_2x2_array(vertex_dictionary, vertex)
            if same:
                grid[i][j]= concatenate_label_of_vertex_into_int(labels)
            else:
                grid[i][j] = 5
                vt_dict2[tuple(vertex)]=np.full((2,2),-1)
                for row in range(2):
                    for col in range(2):
                        vt_dict2[tuple(vertex)][row][col] = collect_label_of_vertex_into_int(vertex_dictionary,vertex,[row,col])
    return grid, vt_dict2

#Sums given lists elementwise
def elementwise_list_sum(lists): #assumes given lists are of the same length, but any number of lists can be given (inside []). Each list should contain only integers
    summed = []
    for i in range(len(lists[0])):
        sum=0
        for list in lists:
            sum = sum + list[i]
        summed.append(sum)
    return summed

def are_2_lists_same(list1,list2): #assumes equal length of given two lists, elements of lists assumed to be comparable
    same = 1
    for i in range(len(list1)):
        if list1[i]!=list2[i]:
            same = 0
    return same

#Counts number of different lists in a (nonempty) collection of lists (collection = list). Assumes all lists in the collection have
#the same length and their elements are comparable
def number_of_different_lists(lists):
    number = 1
    numbers = np.zeros(len(lists))
    for i  in range(len(lists)):
        if numbers[i]==0: #if the list is different from all previous lists, find lists that are equal to it and change
            # their number to the same one. Then increase the number.
            for j in range(len(lists)):
                if are_2_lists_same(lists[i],lists[j]):
                    numbers[j]=number
            number = number +1
    return max(numbers) #number of different numbers in numbers = number of different lists

#check whether a labelled grid covers the unit cube
#idea: for each square in grid1, sum elementwise its 4 labels (e.g. from labels 000,100,110,010, get a new list 220) and check
#whether after doing this for all squares, we get six different lists (220,... - the six different elementwise sums that
#the six different squares of a unit cube produce)
def does_lbld_grid_cover_cube(vertex_dictionary):
    # if np.amin(grid2) == -1 or True in [-1 in array for array in dictionary.values()]:
    #     print 'There are some unassigned values (-1s) in grid2!!! Possible if there are holes. Or check the program'
        #To have this warning stop appearing when there are just holes, I could modify the grid2 right after it's been created at the start -
        #place 0's to vertices which will never be labelled because they are inside holes. (Also this would be consistent with grid1
        #notation for holes)...then amend also the translate_labellings function please
    sums = []
    for sq_row in range(dim_a):
        for sq_col in range(dim_b):
            if grid1[sq_row][sq_col] != 0: #if the square is not a hole in grid1
                sq_labels = get_square_labels_from_vertex_dictionary(vertex_dictionary,[sq_row,sq_col])
                #print 'square labels',sq_labels
                sums.append(elementwise_list_sum(sq_labels))
    cover = 0
    if number_of_different_lists(sums)>5:
        cover = 1
    return cover

#Extends partial labellings by adding labels of a new square, produces a new list of consistent labellings
#Input: a list of partial labellings, where each item - vt_dict - is a partially labelled vertex dictionary,
# from_sq = square from which we reached the new_sq, new_sq is its neighbour. new_sq = square we
# want to add to labellings (given by coords in grid1)
#Output: a list of partially labelled vertex_dictionaries [vt_dict,vt_dict,...]
def add_square_lblng(partial_lblngs,from_sq,new_sq):
    vertices_of_new_sq_in_grid2 = coord_for_sq(new_sq) #[[r1,c1],[r2,c2],[r3,c3],[r4,c4]], starting bottom left, going ccw around new square
    #print 'vertices of new_sq in grid2',vertices_of_new_sq_in_grid2
    new_partial_lblngs = []
    for vt_dict in partial_lblngs:
        #compute the proposed labellings for the new_sq based on the given net
        dirs,nbhs = get_neighbours_direction(from_sq)
        direction_new_sq_from_sq = dirs[nbhs.index(new_sq)] #0=down, 1=right, 2=up, 3=left
        new_sq_lglngs = []
        if direction_new_sq_from_sq == 0:
            #print '\n','lbl net:',lbl_net
            #print from_sq
            #print get_square_labels_from_grid2(lbl_net,from_sq)
            new_sq_lglngs = label_down(get_square_labels_from_vertex_dictionary(vt_dict,from_sq))
        if direction_new_sq_from_sq == 1:
            new_sq_lglngs = label_right(get_square_labels_from_vertex_dictionary(vt_dict,from_sq))
        if direction_new_sq_from_sq == 2:
            new_sq_lglngs = label_up(get_square_labels_from_vertex_dictionary(vt_dict,from_sq))
        if direction_new_sq_from_sq == 3:
            new_sq_lglngs = label_left(get_square_labels_from_vertex_dictionary(vt_dict,from_sq))
        #print 'new_sq_lblngs',new_sq_lglngs
        #print len(new_sq_lglngs)
        for new_sq_labels in new_sq_lglngs: #extend labelling if consistent with already labelled vertices of new_square
            consistent = 1
            for i in range(4):#vertex i lbl in grid2 is different from the suggested one in new_sq_labels
                label_in_vt_dict = label_of_vertex_in_vertex_dictionary(vt_dict,vertices_of_new_sq_in_grid2[i],i)
                suggested_label = new_sq_labels[i]  #both label_in_grid2 and suggested label are lists of length 3 e.g. [0,1,0]
                #print i,label_in_grid2, suggested_label
                #print label_in_grid2[0]!= -1
                for coor in range(3):
                    if label_in_vt_dict[coor] !=-1 and label_in_vt_dict[coor] != suggested_label[coor]:
                        consistent = 0
            #print consistent
            if consistent: #create a new partial labelling and add it to the list
                new_partial_vt_dict = copy.deepcopy(vt_dict)
                new_partial_vt_dict = insert_square_labels_into_vertex_dictionary(new_partial_vt_dict,new_sq,new_sq_labels)
                new_partial_lblngs.append(new_partial_vt_dict)
                # new_partial_net_and_dict = copy.deepcopy(lbl_net_and_dict)
                # new_partial_net, new_partial_dict = insert_square_labels_into_grid2(new_partial_net_and_dict[0],new_partial_net_and_dict[1],new_sq,new_sq_labels)
                # new_partial_lblngs.append([new_partial_net, new_partial_dict])
    return new_partial_lblngs

#Outputs a list, each item is a valid labelling in the form [grid,vt_dictionary] (keys in dictionary = coordinates of vertices in
# grid2, values = 2x2 array, each array position contains a label (e.g. '100','10',etc.) of that vertex from that position
# in the square). grid contains labels for those vertices that had 4 same labels in the array, all others vertices
# [i,j] have grid[i][j] = 5 and their labels are stored in vt_dict
def find_grid_labellings(): #grid1 is assumed to contain 0 for every hole square
    partial_lblgs = [] #list of grids+dictionaries to store partial labellings found so far
    step = -1
    start_row, start_column = np.where(grid1_BFS == 0)
    #starting square
    from_sq = [start_row[0],start_column[0]]   #grid1[i[0]][j[0]]
#     print start_row, start_column
#     print from_sq
#     #label the first square
    partial_lblng1 = copy.deepcopy(vertex_dictionary)

    #WLOG, the from_sq gets assigned cube face [0,0,1],[1,0,1],[1,1,1],[0,1,1]
    partial_lblng1 = insert_square_labels_into_vertex_dictionary(partial_lblng1,from_sq,[[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
    partial_lblgs.append(partial_lblng1)
    while step < np.amax(grid1_BFS)+1:
        step = step + 1
        i, j = np.where(grid1_BFS == step) #coordinates of squares found at this step in BFS (up to 4 new squares neighbouring some of previously labelled squares)
        from_sqrs = zip(i,j) #extend partial labellings by adding labels of one of these squares at a time
        for from_sq in from_sqrs: #from_sqrs can be empty
            #print 'from_sq',from_sq
            dirs, nbhs = get_neighbours_direction(from_sq)
            for i in range(len(nbhs)):
                if is_square_labelled(partial_lblgs[0],nbhs[i])==0:
                    partial_lblgs = add_square_lblng(partial_lblgs, from_sq, nbhs[i])
    #check which of the found consistent labellings covers the entire unit cube:
    valid_lblngs = []
    for labelling in partial_lblgs:
        #print 'LABELING',labelling
        if does_lbld_grid_cover_cube(labelling):
            valid_lblngs.append(labelling)
    #translate the found labellings into an easier to read format:
    translated_valid_labellings = []
    for labelling in valid_lblngs:
        grid,dict = translate_labelling(labelling)
        translated_valid_labellings.append([grid,dict])
    return translated_valid_labellings #valid_lblngs #partial_lblgs   #translated_labellings

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
time1=time.time()

grid1,grid2 = create_grids()

print 'grid1 \n',grid1,'\n'
print 'grid2 \n',grid2,'\n'

slit_dictionary = create_dictionary_for_slits()
print 'slit_dictionary \n',slit_dictionary,'\n'

vertex_dictionary = create_dictionary_for_vertices()
print 'vertex_dictionary \n',vertex_dictionary,'\n'

grid1_BFS = search_net(start_pos)
print 'grid_BFS \n',grid1_BFS, '\n'

complete_labellings = find_grid_labellings()

for i in range(len(complete_labellings)):
    print complete_labellings[i][0]
    print complete_labellings[i][1],'\n'
print 'Found valid labellings:', len(complete_labellings),'\n'

time2=time.time()
print 'Time:',time2-time1,'\n'




