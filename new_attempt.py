
import itertools
import copy

def BinTreeGen(n):
    # A function which takes an integer as input and returns a list of all
    # binary rooted trees with n leaves.
    if  n <= 1:
        print('Input must be at least 1')
    else:
        return trees(n,0)

# We proceed recursiely: at the root node we consider the subtrees to the left
# and the right respectively, and cycle over the possible sizes of these
# subtrees.

def trees(n,leafnode):
    alltrees = []
    if n == 1:
        label = leafnode + 1
        return [label]       #labels the nodes acording the the standard order.
    else:
        for i in range(1,n):
            L = trees(i,leafnode)
            R = trees(n-i,leafnode+i)

            for j in range(0,len(L)):
                for k in range(0, len(R)):
                    subtree = [L[j],R[k]]
                    alltrees.append(subtree)
    return alltrees

##############################################################################

def interactions(tree):
    # Takes as input a tree with n leaves
    # Returns a list of all points on a binary rooted tree on which interactions
    # can take place, collectiely called locations.
    if type(tree) == int:
        return [tree]

    # Working recursively we compute the list of interactions as the union of
    # the list of interactions from the left and right branches of oour tree,
    # adding the interaction at the root node in the form of the entire tree.
    else:
        L = tree[0]
        R = tree[1]
        interactionlist = interactions(L) + interactions(R)
        interactionlist.append(tree)


        # Add the B-type interctions aster internal nodes
        if type(L) == list:
                interactionlist = interactionlist + [[L]]
        if type(R) == list:
                interactionlist = interactionlist + [[R]]

        # Orders the interactions in interactionlist: no actual application,
        # just makes it easier to check that the code works.
        leaf = []
        internal_edge = []
        internal_vertex = []

        for location in interactionlist:
            if type(location) == int:
                leaf.append(location)
            elif len(location) == 1:
                internal_edge.append(location)
            elif len(location) == 2:
                internal_vertex.append(location)

        interactionlist = leaf + internal_vertex + internal_edge

        return interactionlist

###############################################################################

def TopEdges(T, edgelist):
    # Given any list of edges on a tree T, returns a list of those edges which
    # maximal with respect to other edges in the list.
    temp = []
    for edge1 in edgelist:
        count = 0
        for edge2 in edgelist:
            if edge1 != edge2 and IsAbove(edge2, edge1) == 1:
                count += 1
        if count == 0:
            temp.append(edge1)
    return(temp)

def PaintEdgepath(a_type, location, edgelist, T):
    edges_below = Edge_paths_below(location, edgelist, T)
    monomial = a_type[0]
    partition = list(prodlist(Mon_Part(monomial, edges_below, T)))
    return partition

def PaintTop(T, pot, edgelist):
    if edgelist == ():
        return()
    numvars = len(pot[0])
    top = TopEdges(T, edgelist)
    a_list = A_interactions(edgelist, pot, T)
    # create lits of patterns, called pathlist. A pattern is a tuple
    # ((a_type, edges),...) with one sub tuple for each edge in top, where
    # a_type is a tuple ofthe form (mon, theta_type, psi_type).
    templist = Empty(len(top))
    for item in a_list:
        temp = []
        partition = PaintEdgepath(item[0], item[1][1], edgelist, T)
        # collect together pairs [a_type, edgelist]
        loc = item[1]
        a = item[0] + [loc]
        temp = [[a, p[0]] for p in partition]
        for i in range(0, len(top)):
            if IsAbove(item[1][1],top[i]) == 1:
                templist[i] += temp
    # collect togehter one pattern for each top edge
    pathlist = list(prodlist(templist))

    # check that there are no intersecting edgepaths
    newpathlist = []
    for path in pathlist:
        if OverlapsThetaPaths(Thetapaths(path), T) == 0:
            if Psi_Test(path, T, pot) == 1:
                newpathlist.append(path)

    return(newpathlist)

def Paint(T, pot):
    numvars = len(pot[0])
    edgelist = Edges(T)
    pathlist = PaintTop(T, pot, edgelist)
    newpathlist = []
    # prepare the output to be fed into FullPaint: collect unused edges and mark
    # the patterns with 1 if there are edges still to paint. Also remove any
    # edges which are used twice
    for pattern in pathlist:
        edge_paths = []
        for path in pattern:
            edge_paths += path[1]
        usededges = [path[0] for path in edge_paths]
        if usededges != Reduce(usededges):
                pathlist.remove(pattern)
        else:
            newedges = list(remove_subset(Edges(T), usededges))
            if newedges != []:
                pattern.append(newedges)
                pattern += [1]

    while len(pathlist) > 0:
        for pattern in pathlist:
            # check for he the label assigned after the first painting
            if pattern[-1] != 1:
                newpathlist.append(pattern)
                pathlist.remove(pattern)

            else:
                pattern.remove(1)
                # identify the free edges and then remove them from the data
                newedges = pattern[-1]
                pattern.remove(newedges)
                paint_unused = PaintTop(T, pot, newedges)
                if paint_unused == [[]]:
                    if newedges != []:
                        # newpathlist.append(pattern)
                        pathlist.remove(pattern)
                    else:
                        pathlist.remove(pattern)
                else:
                    templist = []
                    for p in paint_unused:
                        edge_paths = []
                        for path in p:
                            edge_paths += path[1]
                        usededges = [path[0] for path in edge_paths]
                        if usededges == Reduce(usededges):
                            free_edges = list(remove_subset(newedges, usededges))
                            newpattern = pattern + p
                            if OverlapsThetaPaths(Thetapaths(newpattern), T) == 0:
                                if Thetapaths(newpattern) == Reduce(Thetapaths(newpattern)):
                                    if Psi_Test(newpattern, T, pot) == 1:
                                        if free_edges == []:
                                            newpathlist.append(newpattern)
                                        else:
                                            newpattern.append(free_edges)
                                            newpattern += [1]
                                            templist.append(newpattern)

                pathlist.remove(pattern)
                pathlist += templist

    # templist = []
    # for pattern in newpathlist:
    #     if Psi_Test(pattern, T, pot) == 1:
    #         templist.append(pattern)

    return(newpathlist)


def Edge_paths_below(e, edgelist, T):
    edge = e
    templist = []
    for edge1 in edgelist:
        if IsAbove(edge, edge1) == 1:
            #choose a vertex where the theta term is destroyed
            vertices = Vertices_below(edge1, T)
            paths  = [[edge1, v] for v in vertices]
            for item in paths:
                templist.append(item)
    return(templist)

def Vertices_below(loc, T):
    vertices = Vertices(T)
    templist = []
    for v in vertices:
        if IsAbove(loc, v) == 1:
            templist.append(v)
    return(templist)

def A_interactions(edgelist, pot, T):
    numvars = len(pot[0])
    # gives a list of the possible combinations of location and type of A-type
    # interactions above a given edgelist.
    top = TopEdges(T, edgelist)   # Identifies the top edges in the edgelist
    a_types = [list(item) for item in Feynmann(pot, numvars)]

    a_list = []
    for edge in top:
        locations_above = []
        for location in Edges(T) + Leaves(T):
            if IsAbove(location, edge) == 1 and location != edge:
                for v in Vertices_below(location, T):
                    locations_above.append([location, v])

        loc = list(prodlist([a_types, locations_above]))
        a_list = a_list + [a for a in loc]
    return(a_list)

def Mon_Part(mon, list, T):
    templist = []
    for num in mon:
        if num > 0:
            subs = subsets(list, num)
            newsubs = []
            for sub in subs:
                if OverlapsThetaPaths1(sub, T) == 0:
                    newsubs.append(sub)
            templist.append(newsubs)
    return(templist)

def Thetapaths(config):
    # returns the edge paths in a config
    templist = []
    for path in config:
        templist.append([path[0][2], path[0][3]])
        # for index in range(0,len(path[0][0])):
        #     if index != 0:
        #         j = index
        #         break
        j = 0
        for theta_type in range(0, len(path[0][0])):
            index = path[0][0][theta_type]
            edges = [[theta_type, p] for p in path[1][j:j+index]]
            templist += edges
            j += index
    return(templist)



###############################################################################
# MAIN PROCEDURES
###############################################################################
def OverlapsThetaPaths1(theta_paths, T):
    # special version of the proc below for primitive theta paths.
    overlapcount = 0
    for path1 in theta_paths:
        for path2 in theta_paths:
            if path1 != path2:
                if len(Overlaps(path1, path2, T)) > 0:
                    overlapcount += 1
    return(overlapcount)

def OverlapsThetaPaths(theta_paths, T):
    # For each theta path in t, tests that it has no intersections of
    # theta paths of the same type.
    overlapcount = 0
    for path1 in theta_paths:
        for path2 in theta_paths:
            if path1 != path2:
                if len(Overlaps(path1, path2, T)) > 0 and path1[0][1] == path2[0][1]:
                    overlapcount += 1
    return(overlapcount)

def Psi_Test(config, T, potential):
    # Checks that there are enough psi input terms in T to feed the A and C type
    # interactions in in config.
    numvars = len(potential[0])
    A = A_types(config)
    C = C_types(config)
    count = 0
    for ints in A:
        if Enough_Leaves_AboveA(ints, config, T, numvars) == 0:
            count += 1
    for ints in C:
        if Enough_Leaves_AboveC(ints, config, T, numvars) == 0:
            count += 1
    if count == 0:
        return(1)
    else:
        return(0)

def Enough_Leaves_AboveC(interaction, config, T, numvars):
    loc = Location(interaction)

    if IsInput(loc[0]) == 1:
        numleaves = 1
        leaf = loc[0]
        leafints = Interactions_Above(leaf, config)
        if len(leafints) > 0:
            ints = [interaction] + leafints
        else:
            ints = [interaction]
    else:
        numleaves = len(Leaves_Above([loc[0]], T))
        edgeints = Interactions_Above([loc[0]], config)
        ints = [interaction] + edgeints

    if len(ints) > numvars*numleaves:
        return(0)

    typecount = Zeros(numvars)
    for item in ints:
        for i in range(0, numvars):
            if item[1] == i:
                typecount[i] += 1

    if max(typecount) <= numleaves:
        return(1)
    else:
        return(0)

def Enough_Leaves_AboveA(interaction, config, T, numvars):

    loc = Location(interaction)
    ints = Interactions_Above(Location(interaction), config)
    numleaves = len(Leaves_Above(loc, T))
    if len(ints) > numvars*numleaves:
        return(0)

    typecount = Zeros(numvars)
    for item in ints:
        for i in range(0, numvars):
            if item[1] == i:
                typecount[i] += 1

    if max(typecount) <= numleaves:
        return(1)
    else:
        return(0)

###############################################################################
# OPERATIONS ON CONFIGS
###############################################################################

def C_types(pattern):
    # Extracts a list of all C-type interactins in config, giving a list of
    # tuples (location, psi_type).
    templist = []
    for path in pattern:
        c_int = [path[0][3][1], path[0][1]]
        templist.append(c_int)
    edge_paths = Thetapaths(pattern)
    for path in edge_paths:
        c_int = [path[1][1], path[0]]
        templist.append(c_int)
    return(templist)

def A_types(pattern):
    # Extracts a list of all A-type interactions in config, giving a list of
    # tuples (location, psi_type, mon, theta_type)
    templist = []
    for path in pattern:
        a_int = [path[0][3][0], path[0][2], path[0][0], path[0][1]]
        templist.append(a_int)
    return(templist)

def Location(interaction):
    # Given an intreractions of type A, B, or C, returns its location.
    return(interaction[0])

def InteractionType(interaction):
    # Given an A or C type interaction, returns the phi_type of the interaction.
    return(interaction[1])

def Interactions_Above(loc, config):
    # Given an interaction which is part of a configuration config on a
    # tree T, returns a list of the A or C type interactions which lie above
    # it on T. Note: this includes he given interaction.

    # Collect a list of all A or C type interactions above. Note: this includes
    # our given interaction.
    A_or_C = A_types(config) + C_types(config)
    templist = []
    for item in A_or_C:
        if IsAbove(Location(item), loc) == 1:
            templist.append(item)
    return(templist)

def ThetaPaths(config):
    templist = []

    for pattern in config:
        if type(pattern[1][0][1]) == int:
            templist.append(pattern[1])
        else:
            templist = templist + pattern[1]
        templist.append(pattern[2])
    return(templist)

###############################################################################
# MONOMIAL PREPARATION AND MANIPULATION
###############################################################################

def Feynmann(pot, numvars):
    # Give a potential in with standard decomposition (in the form of a tuple
    # of tuples), returns a list of all A-type interactions, i.e. a list of
    # tuples (mon, theta_type, psi_type)
    potential = []
    for i in range(0,len(pot)):
        potential.append(list(pot[i]))
    templist = []
    for i in range(0, len(potential)):
        mon = potential[i]
        for index in range(0,numvars):
            if mon[index] > 0:
                derived_mon = Derive(index, mon)
                templist.append((derived_mon, index, i))
    return(templist)

def valid_A_types(edge, loc, t, T, potential):
    # Given an edge where a boson of type t is destroyed and a chosen vertex
    # loc at which the boson was created, generates a list of possible A-type
    # interactions which could lie at loc.
    numvars = len(potential[0])
    h_loc = Height(loc, T)
    h_edge = Height(edge, T)
    templist = []
    for type in Feynmann(potential, numvars):
        # Check that type produces a monomial with nonzero x_t term and that
        # the order of the monomial is large enough to survive to reach edge
        # sand small enough to be destrpyed by the edges below loc in the tree.
        if type[0][t] > 0 and sum(type[0]) in range(h_loc - h_edge, h_loc + 1):
            templist.append(type)
    return(templist)

# def colour(edgelist, potential, numvars):
#     k = sum(potential)
#     paths = subsets(edgelist, k)
#     temp = []
#     for item in paths:
#         for i in range(0, numvars):
#             temp.append(subsets(item, potential[i]))

def Partition(L, l):
    # L is a list and l is a list of integers. Returns a list of ways to
    # L into subsets of the sizes in l.
    # # if l is None:
    # #     return()
    if sum(l) > len(L):
        return()

    if len(l) == 2:
        if max(l) == 1:
            templist = []
            for item in L:
                G = []
                for i in L:
                    if i != item:
                        G.append(i)
                temp = [[item, 1]] + [[non, 0] for non in G]
                templist.append(temp)
            return(templist)
        else:
            temp = []
            subs = subsets(L, max(l))
            for sub in subs:
                reduced_L = list(remove_subset(L, sub))
                if min(l) == 0:
                    templist = [[sub], [reduced_L, 0]]
                    temp.append(templist)
                else:
                    k = [min(l), 0]
                    templist = [sub] + [item for item in Partition(reduced_L, k)]
                    temp.append(templist)
            return(temp)


    subs = subsets(L, max(l))
    l.remove(max(l))

    temp = []
    for sub in subs:
        for item in partition(remove_subset(L, sub), l):
            templist = [sub] + item
            temp.append(templist)
    return(temp)
# THIS NEEDS TO BE FIXED, DOES WORK FOR TWO VARIABLES THOUGH

###############################################################################
# TREE OPERATIONS
###############################################################################

# Extracts the edges, vartices, and leaves form the location list of the tree.

def Edges(T):
    templist = []
    for loc in interactions(T):
        if IsEdge(loc) == 1:
            templist.append(loc)
    return(templist)

def Vertices(T):
    templist = []
    for loc in interactions(T):
        if IsVertex(loc) == 1:
            templist.append(loc)
    return(templist)

def Leaves(T):
    templist = []
    for loc in interactions(T):
        if IsInput(loc) == 1:
            templist.append(loc)
    return(templist)

# Three procedures which decide whether a given locaiton is a leaf, an internal
# edge, ot an internal vertex respectively.

def IsInput(x):
    if type(x) == int:
        return(1)
    else:
        return(0)

def IsEdge(x):
    if type(x) == list:
        if len(x) == 1:
            return(1)
        else:
            return(0)
    else:
        return(0)

def IsVertex(x):
    if type(x) == list:
        if len(x) == 2:
            return(1)
        else:
            return(0)
    else:
        return(0)


def IsAbove(x,y):
    # Given two locations x and y on a tree T, decides whether x is above (including
    # equal to) y in T.
    if type(y) != list:
        if x == y:
            return(1)
        else:
            return(0)
    if x == y:
        return(1)
    if type(x) == list:
        if len(x) == 1:
            if x[0] == y:
                return(0)
            else:
                return(IsAbove(x[0],y))
    if len(y) == 1:
        return(IsAbove(x,y[0]))
    if len(y) == 2:
        l = IsAbove(x,y[0])
        r = IsAbove(x,y[1])

        if l + r >= 1:
            return(1)
        edge_aboveright = y[1]
        edge_aboveleft = y[0]

        if edge_aboveleft == x or edge_aboveright == x:
            return(1)
        else:
            return(0)

def count_edges_below(u, T):
    # Gives the number of edges below a location u. Note: Does NOT include u, this
    # wluld cause problems otherwise.
    count = 0
    for edge in Edges(T):
        if IsAbove(u,edge) and u != edge:
            count += 1
    return(count)


def IsBetween(a,x,y):
    # Determines whether a location a is between two other locations x and y,
    # i.e. above x and downstrean of y.

    if IsAbove(a,x) == 1 and IsAbove(y,a) == 1:
        if a != x:
            return(1)
    else:
        return(0)

def Overlaps(x,y,T):
    # Determines whether two paths in the tree overlap and gives the locations
    # in both.
    # takes as imput two lists x and y, each a source and target for a theta
    # in a given tree T
    Locations = interactions(T)
    Intersection = []
    for loc in Locations:
        if IsBetween(loc, x[2], x[1]) == 1 and IsBetween(loc, y[2], y[1]):
            Intersection.append(loc)
    return(Intersection)

def Height(loc, T):
    count = 0
    for edge in Edges(T):
        if IsAbove(loc, edge) == 1:
            count = count + 1
    return(count)

def Leaves_Above(v,T):
    # Returns a list of the leaves above a given location
    templist = []
    if IsVertex(v) == 1:
        if type(v[0]) == list:
            return(Leaves_Above([v[0]],T))
        else:
            return(Leaves_Above(v[0],T))
    if IsInput(v) == 1:
        return([v])
    if IsEdge(v) == 1:
        for leaf in Leaves(T):
            if IsAbove(leaf, v) == 1:
                templist.append(leaf)
        return(templist)

###############################################################################
# MISC PROCEDURES
###############################################################################

def Zeros(n):
    # Generates a list of length n contiaining all zeros. Used for initialising
    # lists of a fixed length.
    templist = []
    for i in range(0,n):
        templist.append(0)
    return(templist)

def Empty(n):
    # Generates a list of length n contiaining all zeros. Used for initialising
    # lists of a fixed length.
    templist = []
    for i in range(0,n):
        templist.append([])
    return(templist)


def prodlist(L):
    # Given a list of lists L = [L1, ..., Ln] returns the product list, i.e. a list
    # L1 x ... x Ln whose entries are lists with first entry form L1, second from L2,
    # etc.
    prod = list(itertools.product(*L))
    templist = []
    # Itertools outputs tuples: we want lists for that the data type is mutable.
    for item in prod:
        item = list(item)
        templist.append(item)
    return(templist)

def ItemPrint(list):
    # prints each element of a given list on a new line
    for item in list:
        print(item)

def Reduce(list):
    templist = []
    for item in list:
        count = 0
        for item2 in templist:
            if item == item2:
                count += 1
        if count == 0:
            templist.append(item)
    return(templist)

def Derive(k, L):
    # treating L as an integer vector representing a monomial, this returns
    # the representative of the monomial under the x_k derivative.
    temp = []
    for item in range(0,k):
        temp.append(L[item])
    temp.append(L[k] - 1)
    for item in range(k+1,len(L)):
        temp.append(L[item])
    return(temp)

def subsets(L, n):
    temp = [list(item) for item in itertools.combinations(L, n)]
    return(temp)

def remove_subset(L, s):
    LL = list(L)
    for item in s:
        LL.remove(item)
    return(tuple(LL))
