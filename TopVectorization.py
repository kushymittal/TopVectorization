import numpy as np
from scipy.misc import imread
from matplotlib.pyplot import imshow, show
import matplotlib.pyplot as plt 

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
            N = [[shifted_gradx[i+1, j], shifted_grady[i+1, j]],
                 [shifted_gradx[i, j+1], shifted_grady[i, j+1]], 
                 [shifted_gradx[i-1, j], shifted_grady[i-1, j]], 
                 [shifted_gradx[i, j-1], shifted_grady[i, j-1]]]
            for test in N: 
                if rad[i, j] < np.linalg.norm(test):
                    rad[i, j] = np.linalg.norm(test)


    plt.figure(2)
    imshow(M, cmap='gray')
    plt.figure(3)
    imshow(rad, cmap='gray')

    img = np.zeros(np.shape(im))
    srad = np.zeros(np.shape(im))
    for i in range(np.shape(im)[0]):
        for j in range(np.shape(im)[1]):
            if(deltax[i,j] == 0 and deltay[i,j] == 0):
                continue
            si = int(i + deltax[i,j])
            sj = int(j + deltay[i,j])
            img[int(si),int(sj)] = 1 
            srad[int(si),int(sj)] = rad[i,j]

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
                print("hello")

    #calculate MST
    print(len(edgelist))
    MST = getMST(nodes, edgelist)
    #prune 'twigs'
    prune(MST)
    pass

def getMST(nodes, edgelist):
    A = []
    sets = {}
    s = np.argsort(edgelist.values())
    keys = edgelist.keys()
    edges = [keys[x] for x in s]
    print(s)
    for node in nodes:
        sets[node] = set([node])
    for edge in edges:
        u = (edge[0], edge[1])
        v = (edge[2], edge[3])
        if findset(sets, u) is not findset(sets, v):
            A += [edge]
            sets[u].union(v)
    print(len(A))
    return A

def findset(sets, node):
    for s in sets:
        if node in s:
            return s

def prune(MST):
    return MST


def key(a, b):
    if a[0] < b[0] or (a[0] == b[0] and a[1] < b[1]):
        return (a[0], a[1], b[0], b[1])
    return (b[0], b[1], a[0], a[1])

def redraw(im):
    pass

if __name__ == '__main__':
    im = imread('img/Screen Shot 2017-12-06 at 6.04.30 PM.png', flatten=True)
    plt.figure(1)
    imshow(im, cmap='gray')
    (dis, rad) = stroke_disam(im)
    topo = topo_extraction(dis, rad)
    redraw(topo)
    show(block=True)
