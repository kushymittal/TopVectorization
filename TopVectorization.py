import numpy as np
from scipy.misc import imread
from matplotlib.pyplot import imshow, show
import matplotlib.pyplot as plt 
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
                si = i + deltax[i,j]
                sj = j + deltay[i,j]

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
                    rad[i, j] = 2*np.linalg.norm(test)



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
    imshow(srad, cmap='gray')

    plt.figure(3)
    imshow(img, cmap='gray')
    return (img, srad)

def topo_extraction(im, m):
    #build graph
    nodes = set()
    edgelist = {} 
    for i in range(np.shape(im)[0]):
        for j in range(np.shape(im)[1]):
            if im[i,j] != 0:
                nodes.add((i,j))
    """
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
    """
    with open("dump/edgelist", "r") as f:
    #    pickle.dump(edgelist, f)
        edgelist = pickle.load(f)

    #calculate MST
    print(len(edgelist))
    MST = getMST(nodes, edgelist)
    #prune 'twigs'
    return prune(MST, m)
    pass

def getMST(nodes, edgelist):
    A = {} 
    sets = {}
    s = np.argsort(edgelist.values())
    keys = edgelist.keys()
    edges = [keys[x] for x in s]
    print(s)
    for node in nodes:
        sets[node] = set([node])
        A[node] = {}
    for edge in edges:
        u = (edge[0], edge[1])
        v = (edge[2], edge[3])
        if findset(sets, u) is not findset(sets, v):
            A[u][v] = edgelist[key(u,v)] 
            A[v][u] = edgelist[key(u,v)]
            sets = union(sets, u, v) 
    count = 0
    for curr_set in sets.values():
        if len(curr_set) != 0:
            count += 1
            print(len(curr_set))
    print(count)

    

    return A

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
        heap = sorted(heap, key=lambda x: MST[x].values()[0])
        node = heap[0]
        heap = heap[1:]

        other_coords = MST[node].keys()[0]
        maxdist = m[node[0],node[1]]
        if MST[node][other_coords] > maxdist:
            continue 
        biggest_deleted[other_coords] = max(biggest_deleted[other_coords], MST[node][other_coords])
        del MST[other_coords][node]
        if len(MST[otbreakher_coords]) == 1:
            MST[other_coords][MST[other_coords].keys()[0]] += biggest_deleted[other_coords]
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

def redraw(MST, dim):
    #trace centerlines
    for node in MST.keys():
        for other in MST[node].keys():
            #draw a straight line
            pass
        pass
    
    # Second, our reverse drawing procedure
    # utilizes the drawing topology to identify ambigu-
    # ous regions (e.g., junctions and sharp corners) and then corrects the
    # centerline estimates by choosing the most likely centerline config-
    # uration among all possible ones.

    pass

if __name__ == '__main__':
    im = imread('img/Screen Shot 2017-12-06 at 6.04.30 PM.png', flatten=True)
    plt.figure(1)
    imshow(im, cmap='gray')
    # (dis, rad) = stroke_disam(im)
    # with open("dump/disam", "w") as f:
    #     pickle.dump(dis, f) 
    # with open("dump/rad", "w") as f:
    #     pickle.dump(rad, f)
    dis = pickle.load(open("dump/disam", "r"))
    rad = pickle.load(open("dump/rad", "r"))
    topo = topo_extraction(dis, rad)
    redraw(topo)
    show(block=True)
