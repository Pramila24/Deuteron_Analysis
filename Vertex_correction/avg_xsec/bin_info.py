# bin information for root histograms
#
# this newer version skips overflow bins
#
# the data file created includes parameters and can reconstruct the full 1d or 2d histogram
# you need to use the module bin_info2_datafile
#
# not NOT mix version one and two !

import sys
#Set proper paths to import ROOT
sys.path.append('../../../../pyroot/')
sys.path.append('/apps/root/PRO/lib/')
sys.path.append('/apps/root/PRO/')

import ROOT as R
import numpy as np


#------------------------------------------------------------
# header information for the bin_info output file
header = \
"""
# bin_info data
#! ib[i,0]/ ix[i,1]/ iy[i,2]/ xb[f,3]/ yb[f,4]/ cont[f,5]/ dcont[f,6]/
"""
#------------------------------------------------------------


# class to store bin information for root histograms
class bin_info:
    def __init__(self, nx = 0, ny = 0, dx = 1, dy = 1):
        self.i = []
        self.ix = []
        self.iy = []
        self.xb = []
        self.yb = []
        self.cont = []
        self.dcont = []
        self.weight = []
        # number of bins
        self.nx = nx
        self.ny = ny
        # bin width
        self.dx = dx
        self.dy = dy
    def add(self, i, ix, iy, xv, yv, cont, dcont, weight ):
        self.i.append(i)
        self.ix.append(ix)
        self.iy.append(iy)
        self.xb.append(xv)
        self.yb.append(yv)
        self.cont.append(cont)
        self.dcont.append(dcont)
        self.weight.append(weight)
    def copy_from(self,from_b):
        self.i = np.copy( from_b.i)
        self.ix = np.copy( from_b.ix)
        self.iy = np.copy( from_b.iy)
        self.xb = np.copy(from_b.xb)
        self.yb = np.copy(from_b.yb)
        self.cont = np.copy( from_b.cont)
        self.dcont = np.copy(from_b.dcont)
        self.weight = np.copy(from_b.weight)
        self.nx = from_b.nx
        self.ny = from_b.ny
        self.dx = from_b.dx
        self.dy = from_b.dy
        self.xmin = from_b.xmin
        self.ymin = from_b.ymin
    def make_arrays(self, reshape= False):
        self.i = np.array(self.i)
        self.counter = np.ones_like(self.i)
        self.ix = np.array(self.ix)
        self.iy = np.array(self.iy)
        self.xb = np.array(self.xb)
        self.yb = np.array(self.yb)
        self.cont = np.array(self.cont)
        self.dcont = np.array(self.dcont)
        self.weight = np.array(self.weight)
        shape = (self.nx, self.ny)
        if (reshape):
            self.i=self.i.reshape(shape)
            self.counter=self.counter.reshape(shape)
            self.ix=self.ix.reshape(shape)
            self.iy=self.iy.reshape(shape)
            self.xb=self.xb.reshape(shape)
            self.yb=self.yb.reshape(shape)
            self.cont=self.cont.reshape(shape)
            self.dcont=self.dcont.reshape(shape)
            self.weight=self.weight.reshape(shape)
        # that's all
    
    def make_ones_like(self, another):
        self.i = np.ones_like(another.i)
        self.counter = np.ones_like(another.counter)
        self.ix = np.ones_like(another.ix)
        self.iy = np.ones_like(another.iy)
        self.xb = np.ones_like(another.xb)
        self.yb = np.ones_like(another.yb)
        self.cont = np.ones_like(another.cont)
        self.dcont = np.ones_like(another.dcont)
        self.weight = np.ones_like(another.weight)
        self.nx = another.nx
        self.ny = another.ny
    def make_zeros_like(self, another):
        self.i = np.zeros_like(another.i)
        self.counter = np.zeros_like(another.counter)
        self.ix = np.zeros_like(another.ix)
        self.iy = np.zeros_like(another.iy)
        self.xb = np.zeros_like(another.xb)
        self.yb = np.zeros_like(another.yb)
        self.cont = np.zeros_like(another.cont)
        self.dcont = np.zeros_like(another.dcont)
        self.weight = np.zeros_like(another.weight)
        self.nx = another.nx
        self.ny = another.ny
    def save(self, fo, zeros = False):
        """
        save data into a file object named fo as a flat file
        """
        # write parameters
        fo.write('# histogram parameters \n')
        fo.write('#\\ dx = {0:}\n'.format(repr(self.dx)))
        fo.write('#\\ dy = {0:}\n'.format(repr(self.dy)))
        fo.write('#\\ nx = {0:}\n'.format(repr(self.nx)))
        fo.write('#\\ ny = {0:}\n'.format(repr(self.ny)))
        fo.write('#\\ xmin = {0:}\n'.format(repr(self.xmin)))
        fo.write('#\\ ymin = {0:}\n'.format(repr(self.ymin)))
        # write header
        fo.write(header)
        # now the data
        # bin_info data
        # flatten the arrays before saving
        ib = self.i.flatten()
        ix = self.ix.flatten()
        iy = self.iy.flatten()
        xb = self.xb.flatten()
        yb = self.yb.flatten()
        cont = self.cont.flatten()
        dcont = self.dcont.flatten()
        if (not zeros) :
            # check where the content is not zero
            nz = np.where( cont != 0.)[0]
        else:
            nz = np.array(range(len(cont)))
        for i, i_b in enumerate(ib[nz]):
            l = "%d %d %d %r %r %r %r\n"%(\
                i_b, \
                ix[nz][i], \
                iy[nz][i], \
                xb[nz][i], \
                yb[nz][i], \
                cont[nz][i], \
                dcont[nz][i])
            fo.write(l)
        fo.close()
#----------------------------------------------------------------------

#end of class definition
#----------------------------------------------------------------------

# get histrogram and extract information
# histo is a root histogram
def get_histo_data(histo):
    # number of real bins (excl. overflows)
    data=bin_info()    
    data.nx = histo.GetNbinsX()
    data.ny = histo.GetNbinsY()
    # bin limits
    xaxis = histo.GetXaxis()
    yaxis = histo.GetYaxis()
    data.dx = xaxis.GetBinWidth(1)
    data.dy = yaxis.GetBinWidth(1)    
    # first bin low edge
    data.xmin = xaxis.GetBinLowEdge(1)
    data.ymin = yaxis.GetBinLowEdge(1)
    #
    # extract bin information
    # range starts at one to skip overflow bins
    for ix in range(1, data.nx+1):
        for iy in range(1, data.ny+1):
            i = histo.GetBin(ix,iy)
            xv  = xaxis.GetBinCenter(ix)
            yv  = yaxis.GetBinCenter(iy)
            cont = histo.GetBinContent(i)
            dcont = histo.GetBinError(i)
            if (dcont > 0.):
                weight = 1./dcont**2
            else:
                weight = 0.
            data.add(i, ix, iy, xv, yv, cont, dcont,weight)
    return data
#
def get_histo_data_arrays(histo):
    data = get_histo_data(histo)
    data.make_arrays()
    return data
#
def get_index(bi, ix, iy):
    #----------------------------------------------------------------------
    # calculate the flat array index as well as the root hist index
    # 
    # i_flat, i_root = get_index( ix, iy)
    #
    # flat array index calculation
    i_flat = bi.ny*(ix-1) + (iy-1)
    # root bin number
    i_root = iy*(bi.nx + 2) + ix
    return (i_flat, i_root)

#----------------------------------------------------------------------
