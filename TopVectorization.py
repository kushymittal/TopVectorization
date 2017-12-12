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
    count = np.sum(M)
    rad = np.zeros(np.shape(im))
    deltax = np.zeros(np.shape(im))
    deltay = np.zeros(np.shape(im))
    rdeltax = np.zeros(np.shape(im))
    rdeltay = np.zeros(np.shape(im))
    #while sum of M mask is greater than min
    it = 0
    dt = -.01
    while np.sum(M) > count * .01 and it < 1000:
        it += 1
        #move pixels
        deltax += dt*gradx*M
        deltay += dt*grady*M
        shifted_gradx = np.zeros(np.shape(im))
        shifted_grady = np.zeros(np.shape(im))
        #remove finished ones
        for i in range(np.shape(im)[0]):
            for j in range(np.shape(im)[1]):
                si = int(i + deltax[i,j])
                sj = int(j + deltay[i,j])
                if np.sqrt(np.square(gradx[i,j]) + np.square(grady[i,j])) > \
                        np.sqrt(np.square(shifted_gradx[i,j]) + np.square(shifted_grady[i,j])):
                    shifted_gradx[si,sj] = gradx[i,j]  
                    shifted_grady[si,sj] = grady[i,j]  
        for i in range(1,np.shape(im)[0]-1):
            for j in range(1,np.shape(im)[1]-1):
                if M[i,j] == 0:
                    continue
                si = int(i + deltax[i,j])
                sj = int(j + deltay[i,j])

                N = [[shifted_gradx[si+1, sj], shifted_grady[si+1, sj]],
                     [shifted_gradx[si, sj+1], shifted_grady[si, sj+1]], 
                     [shifted_gradx[si-1, sj], shifted_grady[si-1, sj]], 
                     [shifted_gradx[si, sj-1], shifted_grady[si, sj-1]]]
                Nij = [[-1,0], [0,1], [1, 0], [0,-1]]
                for x in range(len(N)):
                    test = N[x]
                    if np.dot(test,[gradx[i,j], grady[i, j]]) < 0 and  \
                       np.dot(Nij[x],[gradx[i,j], grady[i,j]]):
                        M[i,j] = 0

    print(np.sum(M), count)
    print(it)
    for i in range(1,np.shape(im)[0]-1):
        for j in range(1,np.shape(im)[1]-1):
            N = [[deltax[i+1, j], deltay[i+1, j]],
                 [deltax[i, j+1], deltay[i, j+1]], 
                 [deltax[i-1, j], deltay[i-1, j]], 
                 [deltax[i, j-1], deltay[i, j-1]]]
            for test in N: 
                if rad[i, j] < np.linalg.norm(test):
                    rad[i, j] = np.linalg.norm(test)



    img = np.zeros(np.shape(im))
    srad = np.zeros(np.shape(im))
    print(np.sum(rad)/np.count_nonzero(rad))
    for i in range(np.shape(im)[0]):
        for j in range(np.shape(im)[1]):
            if(deltax[i,j] == 0 and deltay[i,j] == 0):
                continue
            si = int(i + deltax[i,j])
            sj = int(j + deltay[i,j])
            img[int(si),int(sj)] = 1 
            srad[int(si),int(sj)] = rad[i,j]

    plt.figure(2)
    plt.imshow(srad, cmap='gray')

    plt.figure(3)
    plt.imshow(img, cmap='gray')
    return (img, srad)

def topo_extraction(im, m):
    #build graph
    nodes = set()
    edgelist = {} 
    for i in range(np.shape(im)[0]):
        for j in range(np.shape(im)[1]):
            if im[i,j] != 0:
                nodes.add((i,j))
    for node in nodes: 
        maxdist = m[node[0],node[1]]
        #for di in range(-int(maxdist), int(maxdist)):
        #    for dj in range(-int(maxdist), int(maxdist)):
        #        other = (i+di, j + dj)
        for other in nodes:
            if other == node: # or other not in nodes:
                continue
            dist = np.sqrt(np.square(node[0]-other[0]) + np.square(node[1]-other[1]))
            if dist < maxdist:
                edgelist[key(node, other)] = dist
    with open("dump/edgelist", "wb") as f:
        pickle.dump(edgelist, f)
    with open("dump/nodes", "wb") as f:
        pickle.dump(nodes, f)
    """
    with open("dump/edgelist", "rb") as f:
        edgelist = pickle.load(f)
    with open("dump/nodes", "rb") as f:
        nodes = pickle.load(f)
    """
    #calculate MST
    print(len(edgelist))
    MST = getMST(nodes, edgelist)
    #prune 'twigs'
    return edgelist, prune(MST, m)
    pass

def getMST(nodes, edgelist):
    MST = {} 
    sets = {}
    s = np.argsort(list(edgelist.values()))
    keys = list(edgelist.keys())
    edges = [keys[x] for x in s]
    print(s)
    for node in nodes:
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
    print(count)

    

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

def prune(MST, m):
    heap = []
    biggest_deleted = {}
    for node in MST.keys(): 
        biggest_deleted[node] = 0
        if len(MST[node]) == 1:
            heap += [node]
    while heap != []:
        node = min(heap, key=lambda x: list(MST[x].values())[0])
        del heap[node]

        other_coords = list(MST[node].keys())[0]
        maxdist = m[node[0],node[1]]
        if MST[node][other_coords] > 2*maxdist:
            if len(MST[other_coords]) == 0:
                heap.remove(other_coords)
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
    plt.show()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    count = 0
    for node in endpoints:
        count += 1
        if count % 50 == 0:
            print(count)
        shortest_path = dijkstra(node, edgelist)       
        for other in endpoints:
            x = []
            y = []
            if other not in shortest_path:
                continue
            other_obj = shortest_path[other]
            while other_obj.coords != node:
                x += [dims[0] - other_obj.coords[0]]
                y += [other_obj.coords[1]]
                other_obj = other_obj.prev
            x += [dims[0] - node[0]] 
            y += [node[1]] 
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

def dijkstra(start_node, el):

    shortest_path = {}
    node = start_node

    frontier = [Node(node, 0, None)]
    edgelist = {}
    for pair in el.keys():
        u = (pair[0],pair[1])
        v = (pair[2],pair[3])
        if u not in edgelist:
            edgelist[u] = {}
        if v not in edgelist:
            edgelist[v] = {}
        edgelist[u][v] = el[pair]
        edgelist[v][u] = el[pair]

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
    def __init__(self, coords, weight, node):
        self.coords = coords
        self.dist = weight
        self.prev = node


if __name__ == '__main__':
    im = spm.imread('img/Screen Shot 2017-12-06 at 6.04.30 PM.png', flatten=True)
    plt.figure(1)
    plt.imshow(im, cmap='gray')
    (dis, rad) = stroke_disam(im)
    with open("dump/disam", "wb") as f:
        pickle.dump(dis, f) 
    with open("dump/rad", "wb") as f:
        pickle.dump(rad, f)
    #dis = pickle.load(open("dump/disam", "rb"))
    #rad = pickle.load(open("dump/rad", "rb"))
    edgelist, MST = topo_extraction(dis, rad)
    dims = np.shape(im)
    redraw(MST, edgelist, dims)
    plt.show(block=True)
