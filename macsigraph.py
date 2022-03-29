class Graph ():
    """a simple graph library"""

    def __init__ (self):
        """initializes a graph object"""

        import random

        self.compsizes = {}
        self.compnames = {}

        self.N = int (1e3)
        self.z = 5
        self.p = self.z / self.N

        self.C = 0
        self.m = 0

        self.rule = 0 # 0 Erdos-Renyi, 1 Bohman-Frieze, 2 product

    def dfs (self, source, comp):
        """perform a depth-first search that relabels a connected component"""

        S = [source]

        explored = []

        while (len (S) > 0):
            u = S.pop ()
            if (u not in explored):
                explored.append (u)
                self.compnames[u] = comp # labelling the component
                for v in self.nlist[u]:
                    S.append (v)

    def printnlist (self):
        """prints the graph, or more precisely, its neighbour list"""

        for i in self.nlist:
            if (len (self.nlist[i]) > 0):
                print (i, self.nlist[i])
            else:
                print (i, "{}")

    def erdosrenyi_gnp (self):
        """generates an Erdos-Renyi graph, from the Gnp ensemble"""

        import random

        self.nlist = {}

        for u in range (0, self.N):
            self.nlist[u] = set ()

        for u in range (0, self.N):
            for v in range (0, u):
                if (random.uniform (0, 1) < self.p):
                    self.nlist[u].add (v)
                    self.nlist[v].add (u)

    def nullnodes (self):
        """initialises sampled nodes"""

        return 0, 0, 0, 0

    def randnodes (self):
        """sample nodes uniformly at random"""

        import random

        u = random.randrange (0, self.N);
        v = random.randrange (0, self.N);
        x = random.randrange (0, self.N);
        y = random.randrange (0, self.N);

        return u, v, x, y

    def readcomps (self, u, v):
        """simply return component labels of nodes u, v"""

        return self.compnames[u], self.compnames[v]

    def swapvariables (self, u, v):
        """without loss of generality swap what we're calling u and v"""

        if (self.compsizes[self.compnames[u]] < self.compsizes[self.compnames[v]]):
          return v, u
        else:
          return u, v

    def bfrule (self, u, v, x, y):
        """returns an edge resulting from the Bohman-Frieze rule"""

        ucomp = self.compnames[u] 
        vcomp = self.compnames[v]

        if (self.compsizes[ucomp] > 1 or self.compsizes[vcomp] > 1):
            u = x
            v = y

        return u, v

    def prrule (self, u, v, x, y):
        """returns an edge resulting from the product rule"""

        cuv = set (); cxy = set ()

        ucomp = self.compnames[u]; cuv.add (ucomp) 
        vcomp = self.compnames[v]; cuv.add (vcomp)
        xcomp = self.compnames[x]; cxy.add (xcomp)
        ycomp = self.compnames[y]; cxy.add (ycomp)

        if (len (cuv) == 2 and len (cxy) == 2 and cuv != cxy):
            puv = self.compsizes[ucomp] * self.compsizes[vcomp]
            pxy = self.compsizes[xcomp] * self.compsizes[ycomp]
            if (pxy < puv):
                u = x
                v = y
        elif (len (cxy) < 2):
            u = x
            v = y

        return u, v

    def initnewmanziff (self):
        """initialise variables for Newman-Ziff algorithm"""

        self.C = 1                        # size of lcc
        self.m = 0                        # number of edges

        self.compnames = {}               # maps node index to component label
        self.compsizes = {}               # maps component label to component size

        for i in range (0, self.N):
            self.compsizes[i] = 1
            self.compnames[i] = i

        self.nlist = {}                   # stores graph itself as neighbour list

        for i in range (0, self.N):
            self.nlist[i] = set ()

    def sampleedge (self):
        """returns an edge sampled according to some rule, ER by default"""

        u, v, x, y = self.nullnodes ()                                            

        while (u == v or u in self.nlist[v] or x == y or x in self.nlist[y]):
            u, v, x, y = self.randnodes ()
        
        if (self.rule == 1):
            u, v = self.bfrule (u, v, x, y) # Bohman-Frieze rule

        if (self.rule == 2):
            u, v = self.prrule (u, v, x, y) # product rule

        return u, v

    def newmanziff (self):
        """efficient percolation simulation"""

        self.initnewmanziff ()

        while (self.m / self.N < 1.5):

            u, v = self.sampleedge ()

            ucomp, vcomp = self.readcomps (u, v)
            
            if (ucomp != vcomp):  # if u and v not in same component

                u, v = self.swapvariables (u, v)

                ucomp, vcomp = self.readcomps (u, v)

                if ucomp in self.compsizes:
                    self.compsizes[ucomp] += self.compsizes[vcomp]
                else:
                    self.compsizes[ucomp] = 1 

                del self.compsizes[vcomp] 

                self.dfs (v, ucomp)  # merge v's cluster with u's

            self.nlist[u].add (v)
            self.nlist[v].add (u)

            if (self.compsizes[ucomp] > self.C):
                self.C = self.compsizes[ucomp]

            self.m += 1

            print (self.m / self.N, self.C / self.N)  # r, C
