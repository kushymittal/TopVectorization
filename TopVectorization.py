import numpy as np
import scipy.misc as spm
import matplotlib.pyplot as plt 
import matplotlib.lines as lines 
import pickle
import heapq

def stroke_disam(im):
    (gradx, grady) = np.gradient(im)
    #build M mask
    norm_threshold = .1*np.max(np.sqrt(np.square(gradx) + np.square(grady))) 
    M = np.sqrt(np.square(gradx) + np.square(grady)) >= norm_threshold
    global sh
    sh = np.shape(im)
    count = np.sum(M)
    rad = np.zeros(sh)
    deltax = np.zeros(sh)
    deltay = np.zeros(sh)
    rdeltax = np.zeros(sh)
    rdeltay = np.zeros(sh)
    #while sum of M mask is greater than min
    it = 0
    dt = - .01
    while np.sum(M) > count * .01:
        it += 1
        #move pixels
        deltax += dt*gradx*M
        deltay += dt*grady*M
        shifted_grads = {} 
        #remove finished ones
        for i in range(sh[0]):
            for j in range(sh[1]):
                si = int(i + deltax[i,j])
                sj = int(j + deltay[i,j])
                if si >= sh[0] or sj >= sh[1] or si < 0 or sj < 0:
                    deltax[i,j] -= dt*gradx[i,j]
                    deltay[i,j] -= dt*grady[i,j]
                    si = int(i + deltax[i,j])
                    sj = int(j + deltay[i,j])
                    M[i,j] = 0
                if si not in shifted_grads:
                    shifted_grads[si] = {}
                if sj not in shifted_grads[si]:
                    shifted_grads[si][sj] = []
                if deltax[i,j] != 0 or deltay[i,j] != 0:
                    shifted_grads[si][sj] += [(gradx[i,j], grady[i,j])]  
        for i in range(1,sh[0]-1):
            for j in range(1,sh[1]-1):
                if M[i,j] == 0:
                    continue
                si = int(i + deltax[i,j])
                sj = int(j + deltay[i,j])

                Nij = [[-1,0], [-1,1], [0,1], [1,1], [1, 0], [1,-1], [0,-1], [-1,-1]]
                for k in range(len(Nij)):
                    x = Nij[k]
                    if si + x[0] >= sh[0] or sj + x[1] >= sh[1] or si + x[0] < 0 or sj + x[1] < 0:
                        continue
                    for test in shifted_grads.get(si+x[0], {}).get(sj+x[1], []): 
                        if np.dot(test,[gradx[i,j], grady[i, j]]) < 0 and  \
                           np.dot(x,[gradx[i,j], grady[i,j]]):
                            M[i,j] = 0

    print(np.max(deltax), np.max(deltay), np.sum(M))
    for i in range(1, sh[0]-1):
        for j in range(1, sh[1]-1):
            N = [[deltax[i+1, j], deltay[i+1, j]],
                 [deltax[i, j+1], deltay[i, j+1]],
                 [deltax[i-1, j], deltay[i-1, j]], 
                 [deltax[i, j], deltay[i, j]], 
                 [deltax[i+1, j+1], deltay[i+1, j+1]],
                 [deltax[i-1, j+1], deltay[i-1, j+1]],
                 [deltax[i-1, j-1], deltay[i-1, j-1]],
                 [deltax[i+1, j-1], deltay[i+1, j-1]],
                 [deltax[i, j-1], deltay[i, j-1]]]
            for test in N: 
                if rad[i, j] < np.linalg.norm(test):
                    rad[i, j] = 1.5*np.linalg.norm(test)

    nodes = {}
    for i in range(sh[0]):
        for j in range(sh[1]):
            if deltax[i,j] != 0 or deltay[i,j] != 0:
                nodes[i+deltax[i,j],j+deltay[i,j]] = rad[i,j]

    print("mean width", np.sum(rad)/np.count_nonzero(rad))
    print("node count", len(nodes)) 
    #plot nodes
    arr = np.array(list(nodes.keys()))
    x = arr[:, 1]
    y = sh[0] - arr[:, 0]
    plt.figure(4)
    plt.plot(x,y, 'b.')
    return nodes

def topo_extraction(nodes):
    #build graph
    edgelist = {} 
    k = list(nodes.keys())
    for n in range(len(k)): 
        if n % 500 == 0:
            print(n)
        node = k[n]
        #for di in range(-int(maxdist), int(maxdist)):
        #    for dj in range(-int(maxdist), int(maxdist)):
        #        other = (i+di, j + dj)
        for m in range(n,len(k)):
            other = k[m]
            if other == node: # or other not in nodes:
                continue
            dist = np.sqrt(np.square(node[0]-other[0]) + np.square(node[1]-other[1]))
            if dist < nodes[node]:
                edgelist[key(node, other)] = dist
    with open("dump/edgelist", "wb") as f:
        pickle.dump(edgelist, f)
    """
    with open("dump/edgelist", "rb") as f:
        edgelist = pickle.load(f)
    """
    global sh
    #plot graph 
    fig = plt.figure(5)
    ax = fig.add_subplot(111)
    ax.set_ylim(0, sh[0])
    ax.set_xlim(0, sh[1])
    for edge in edgelist.keys():
        x = [edge[1], edge[3]]
        y = [sh[0] - edge[0], sh[0] - edge[2]]
        line = lines.Line2D(x,y)
        ax.add_line(line)
    #calculate MST
    print("num edges", len(edgelist))
    MST = getMST(nodes, edgelist)
    #plot MST 
    fig = plt.figure(6)
    ax = fig.add_subplot(111)
    ax.set_ylim(0, sh[0])
    ax.set_xlim(0, sh[1])
    for u in MST.keys():
        for v in MST[u].keys():
            x = [u[1], v[1]]
            y = [sh[0] - u[0], sh[0] - v[0]]
            line = lines.Line2D(x,y)
            ax.add_line(line)
    #prune 'twigs'
    pMST = prune(MST, nodes)
    #plot pMST 
    fig = plt.figure(7)
    ax = fig.add_subplot(111)
    ax.set_ylim(0, sh[0])
    ax.set_xlim(0, sh[1])
    for u in pMST.keys():
        for v in pMST[u].keys():
            x = [u[1], v[1]]
            y = [sh[0] - u[0], sh[0] - v[0]]
            line = lines.Line2D(x,y)
            ax.add_line(line)
    plt.show()
    return edgelist, pMST 
    pass

def getMST(nodes, edgelist):
    MST = {} 
    sets = {}
    s = np.argsort(list(edgelist.values()))
    keys = list(edgelist.keys())
    edges = [keys[x] for x in s]
    print(s)
    for node in nodes.keys():
        sets[node] = set([node])
        MST[node] = {}
    for edge in edges:
        u = (edge[0], edge[1])
        v = (edge[2], edge[3])
        if findset(sets, u) is not findset(sets, v):
            MST[u][v] = edgelist[key(u,v)] 
            MST[v][u] = edgelist[key(u,v)]
            sets = union(sets, u, v) 
    count = 0
    for curr_set in sets.values():
        if len(curr_set) != 0:
            count += 1
            #print(len(curr_set))
    print("conn comp", count)

    

    return MST 

def findset(sets, node):
    for s in sets.values():
        if node in s:
            return s
    assert False 
def union(sets, u, v):
    for k in sets.keys():
        if u in sets[k]:
            uk = k
    for k in sets.keys():
        if v in sets[k]:
            vk = k
    sets[uk] = sets[uk] | sets[vk]
    sets[vk] = set() 
    return sets

def prune(MST, nodes):
    heap = []
    biggest_deleted = {}
    for node in MST.keys(): 
        biggest_deleted[node] = 0
        if len(MST[node]) == 1:
            heap += [node]
    while heap != []:
        node = min(heap, key=lambda x: list(MST[x].values())[0])
        heap.remove(node)

        other_coords = list(MST[node].keys())[0]
        maxdist = nodes[node] 
        if MST[node][other_coords] > maxdist:
            continue 
        biggest_deleted[other_coords] = max(biggest_deleted[other_coords], MST[node][other_coords])
        del MST[other_coords][node]
        if len(MST[other_coords]) == 1:
            MST[other_coords][list(MST[other_coords].keys())[0]] += biggest_deleted[other_coords]
            heap += [other_coords]
        elif len(MST[other_coords]) == 0:
            del MST[other_coords]
            heap.remove(other_coords)

        del MST[node]


    #find all terminals and put in minheap based on current length + edgeweight
    #pop and delete edge, add new terminal to min heap(if there is one)
    print(len(MST))
    return MST


def key(a, b):
    if a[0] < b[0] or (a[0] == b[0] and a[1] < b[1]):
        return (a[0], a[1], b[0], b[1])
    return (b[0], b[1], a[0], a[1])

def redraw(MST, edgelist, dims):
    #choose endpoints
    endpoints = []
    for node in MST.keys():
        if len(MST[node]) == 1 or len(MST[node]) > 2:
            endpoints += [node]
    print("endpoints: ", len(endpoints))
    #trace centerlines
    fig = plt.figure(2)
    ax = fig.add_subplot(111)
    count = 0

    el= {}
    for pair in edgelist.keys():
        u = (pair[0],pair[1])
        v = (pair[2],pair[3])
        if u not in el:
            el[u] = {}
        if v not in el:
            el[v] = {}
        el[u][v] = edgelist[pair]
        el[v][u] = edgelist[pair]

    for node in endpoints:
        count += 1
        if count % 5 == 0:
            print(count)
        shortest_path = dijkstra(node, el)       
        for other in endpoints:
            #other_obj = astar(node, other, edgelist)       
            if other not in shortest_path: 
                continue
            x = []
            y = []
            other_obj = shortest_path[other]
            while other_obj != None:
                x += [dims[0] - other_obj.coords[0]]
                y += [other_obj.coords[1]]
                other_obj = other_obj.prev
            line = lines.Line2D(y,x)
            ax.add_line(line)
        pass
    ax.set_ylim(0, dims[0])
    ax.set_xlim(0, dims[1])
    plt.show()
    #draw the lines I guess


    
    # Second, our reverse drawing procedure
    # utilizes the drawing topology to identify ambigu-
    # ous regions (e.g., junctions and sharp corners) and then corrects the
    # centerline estimates by choosing the most likely centerline config-
    # uration among all possible ones.

    pass
"""
def astar(start_node, dest, edeglist):

    node = start_node

    frontier = [Node(node, 0, 0, 0, None)]

    while len(frontier):
        node = min(frontier, key=lambda x: x.dist+x.h)
        #print(node)
        #print(frontier)
        frontier.remove(node)
        if node.coords == dest:
            return node 
        for n in edgelist[node.coords].keys():
            h = np.sqrt((dest[0] - n[0])**2+(dest[1] - n[1])**2) #line dist 
            frontier.append(Node(n, node.dist+node.weight, edgelist[node.coords][n], h, node))
    return None

"""

def dijkstra(start_node, edgelist):

    shortest_path = {}
    node = start_node
    frontier = [Node(node, 0, None)]
    while len(frontier):
        node = min(frontier, key=lambda x: x.dist)
        #print(node)
        #print(frontier)
        frontier.remove(node)

        if node.coords not in shortest_path:
            shortest_path[node.coords] = node 
        else:
            continue

        for n in edgelist[node.coords].keys():
            frontier.append(Node(n, node.dist + edgelist[node.coords][n], node))

    return shortest_path

        
        
class Node:
    #def __init__(self, coords, dist, weight, h, node):
    def __init__(self, coords, dist, node):
        self.coords = coords
        self.dist = dist
        #self.weight = weight
        #self.h = h
        self.prev = node


if __name__ == '__main__':
    #im = spm.imread('img/Screen Shot 2017-12-06 at 6.04.30 PM.png', flatten=True)
    im = spm.imread('img/small.gif', flatten=True)
    plt.figure(1)
    plt.imshow(im, cmap='gray')
    nodes = stroke_disam(im)
    with open("dump/nodes", "wb") as f:
        pickle.dump(nodes, f) 
    #dis = pickle.load(open("dump/nodes", "rb"))
    edgelist, MST = topo_extraction(nodes)
    dims = np.shape(im)
    redraw(MST, edgelist, dims)
    plt.show(block=True)
