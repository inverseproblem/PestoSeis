
#------------------------------------------------------------------------
#
#    PestoSeis, a numerical laboratory to learn about seismology, written
#    in the Python language.
#    Copyright (C) 2021  Andrea Zunino 
#
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#------------------------------------------------------------------------

"""
Binary heap trees (min/max).

"""

import numpy as NP

#=================================

def index_parent(i):
    return i//2 

def index_leftchild(i):
    return 2*i ## shift bits???

def index_rightchild(i):
    return 2*i+1 ## shift bits???


# #=============================================================
# #               MAX stuff
# #=============================================================

class BinHeapMax():
    """
      Max heap binary tree.
    """
    #==========================================
    
    def __init__(self,values,Nmax,handles):
        
        assert(values.size==handles.size)
        ## Init the structure 
        ## Nmax,Nh,nodes,handles
        Nh=values.size
        ## Init big arrays
        ##Array{Float64}(0),10,Array{Int64}(0)
        values2  = NP.zeros(Nmax)
        handles2 = NP.zeros(Nmax,dtype=NP.int64)-1
        ## copy to avoid modifying input arrays!!
        values2[0:Nh]  = NP.asarray(values.copy())
        handles2[0:Nh] = NP.asarray(handles.copy())
        self.Nh = Nh
        self.Nmax = Nmax
        self.nodes = values2
        self.handles = handles2
        #h = BinHeapMax(Nmax,Nh,values2,handles2)
        ## heapify the structure
        sta = self.Nh//2
        for i in range(sta-1,-1,-1): #sta:-1:1 ## backwards
            self.__max_heapify(i)
          
    #==========================================    
    def __swap_nodes_heap(self,p,q):
        # temporary copy
        ptmp_node = self.nodes[p]
        ptmp_handle = self.handles[p]
        # swap last and first node
        self.nodes[p] = self.nodes[q]
        self.nodes[q] = ptmp_node
        # swap handles too
        self.handles[p] = self.handles[q]
        self.handles[q] = ptmp_handle

    #==========================================
    
    def topval_heap(self):
        return self.nodes[0],self.handles[0]
                  
    ##=============================================================
    
    def __max_heapify(self,i):
        ## Go DOWN the tree (max heap)...
        # get children indices
        l = index_leftchild(i)
        r = index_rightchild(i)
        ## Introduction to algorithms, p. 154
        if (l <= self.Nh-1) and (self.nodes[l]>self.nodes[i]) :
            # largest on left
            largest = l
        else :
            # largest on i
            largest = i
        if (r <= self.Nh-1) and (self.nodes[r]>self.nodes[largest]) :
            # largest on right
            largest = r
        if largest!=i :
            # swap the nodes
            self.__swap_nodes_heap(i,largest)
            # keep going
            self.__max_heapify(largest)

    ##=============================================================
    def update_node_maxheap(self,val,handle) :
        ## find index of node given handle
        idxh = NP.where(self.handles==handle) #find(self.handles==handle)
        #assert size(idxh,1)<=1
        i = idxh[0] ## from 1D array to scalar
        self.nodes[i] = val
    
        # GO UP THE TREE: compare with parents...
        while (i>0) and (self.nodes[index_parent(i)]<self.nodes[i]) :
            ## swap node with parent
            self.__swap_nodes_heap(i,index_parent(i))
            i = index_parent(i)

        # GO DOWN THE TREE: compare with children...
        if i<self.Nh-1 :
            self.__max_heapify(i)

    ##=============================================================
    
    def insert_maxheap(self,val,handle) :
        ## Go UP the tree (max heap) from the bottom...
        # extend heap
        assert self.Nmax>(self.Nh+1)
        ## IMPORTANT: make sure handle is unique
        ## commented because slows things down
        #   @assert ~(handle in self.handles[1:self.Nh])
        # resize
        self.Nh = self.Nh+1
        # add the new node at the bottom
        self.nodes[self.Nh-1] = val
        self.handles[self.Nh-1] = handle # set handle too
        # start from bottom right leaf
        i = self.Nh-1
        # compare with parents...
        while (i>0) and (self.nodes[index_parent(i)]<self.nodes[i]) :
            ## swap node with parent
            self.__swap_nodes_heap(i,index_parent(i))
            i = index_parent(i)

    ##=============================================================
    
    def pop_maxheap(self) :
        # save handle and values of top to return them
        poppedtophandle = self.handles[0]
        poppedval = self.nodes[0]
        # move last elem to top
        self.nodes[0] = self.nodes[self.Nh-1]
        self.handles[0] = self.handles[self.Nh-1]
        ## set out of Nh nodes to 0
        self.nodes[self.Nh-1] = 0.0
        self.handles[self.Nh-1] = 0
        # shorten heap
        self.Nh = self.Nh-1
        # swap the root [1] with largest child, and so on...
        self.__max_heapify(0)
        return poppedtophandle,poppedval

    ##=============================================================


# #=============================================================
# #               MIN stuff
# #=============================================================


class BinHeapMin():
    """
      Min heap binary tree.
    """
    #====================================================
    
    def __init__(self,values,Nmax,handles):

        assert values.size==handles.size
        ## Init the structure 
        ## Nmax,Nh,nodes,handles
        Nh=values.size
        ## Init big arrays
        ##Array{Float64}(0),10,Array{Int64}(0)
        values2  = NP.zeros(Nmax)
        handles2 = NP.zeros(Nmax,dtype=NP.int64)-1
        ## copy to avoid modifying input arrays!!
        values2[0:Nh]  = values.copy()
        handles2[0:Nh] = handles.copy()
        self.Nh = Nh
        self.Nmax = Nmax
        self.nodes = values2
        self.handles = handles2
        #h = BinHeapMax(Nmax,Nh,values2,handles2)
        ## heapify the structure
        sta = self.Nh//2
        for i in range(sta-1,-1,-1): #sta:-1:1 ## backwards
            self.__min_heapify(i)

    #==================================================    

    def __swap_nodes_heap(self,p,q):
        # temporary copy
        ptmp_node = self.nodes[p]
        ptmp_handle = self.handles[p]
        # swap last and first node
        self.nodes[p] = self.nodes[q]
        self.nodes[q] = ptmp_node
        # swap handles too
        self.handles[p] = self.handles[q]
        self.handles[q] = ptmp_handle

    #==================================================    

    def topval_heap(self):
        return self.nodes[0],self.handles[0]

    ##=============================================================
    
    def __min_heapify(self,i) :
        ## Go DOWN the tree (min heap)...
        # get children indices
        l = index_leftchild(i)
        r = index_rightchild(i)
        ## Introduction to algorithms, p. 154
        if (l <= self.Nh-1) and (self.nodes[l]<self.nodes[i]) :
            # smallest on left
            smallest = l
        else :
            # smallest on i
            smallest = i
        if (r <= self.Nh-1) and (self.nodes[r]<self.nodes[smallest]) :
            # smallest on right
            smallest = r
        if smallest!=i :
            # swap the nodes
            self.__swap_nodes_heap(i,smallest)
            # keep going
            self.__min_heapify(smallest)

    ##=============================================================
    
    def update_node_minheap(self,val,handle) :
        ## find index of node given handle
        #idxh = 0
        # for l in range(0,self.Nh) :
        #     if self.handles[l]==handle :
        #         idxh = l
        idxh = NP.where(self.handles==handle)[0] 
        i = idxh[0] ## from 1D array to scalar
        self.nodes[i] = val
        # GO UP THE TREE: compare with parents...
        while (i>0) and (self.nodes[index_parent(i)]>self.nodes[i]) :
            ## swap node with parent
            self.__swap_nodes_heap(i,index_parent(i))
            i = index_parent(i)
        # GO DOWN THE TREE: compare with children...
        if i<self.Nh-1 :
            self.__min_heapify(i)

    ##=============================================================

    def insert_minheap(self,val,handle) :
        ## Go UP the tree (min heap) from the bottom...
        # extend heap
        assert self.Nmax>(self.Nh+1)
        ## IMPORTANT: make sure handle is unique
        ## commented because slows things down
        #   @assert ~(handle in self.handles[1:self.Nh])
        # resize
        self.Nh = self.Nh+1
        # add the new node at the bottom
        self.nodes[self.Nh-1] = val
        self.handles[self.Nh-1] = handle # set handle too
        # start from bottom right leaf
        i = self.Nh-1
        # compare with parents...
        while (i>0) and (self.nodes[index_parent(i)]>self.nodes[i]) :
            ## swap node with parent
            self.__swap_nodes_heap(i,index_parent(i))
            i = index_parent(i)

    ##=============================================================
    
    def pop_minheap(self):
        # save handle and values of top to return them
        poppedtophandle = self.handles[0]
        poppedval = self.nodes[0]
        # move last elem to top
        self.nodes[0] = self.nodes[self.Nh-1]
        self.handles[0] = self.handles[self.Nh-1]
        ## set out of Nh nodes to 0
        self.nodes[self.Nh-1] = 0.0
        self.handles[self.Nh-1] = 0
        # shorten heap
        self.Nh = self.Nh-1
        # swap the root [1] with largest child, and so on...
        self.__min_heapify(0)
        return poppedtophandle,poppedval


    ##=============================================================

