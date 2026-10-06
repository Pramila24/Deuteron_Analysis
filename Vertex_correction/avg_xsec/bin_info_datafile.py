# bin information for histograms stored as datafiles
#
# this uses a more compact form. Only channels with a non-zero content are stored
#
# this module recreates 1d and 2d histogram created by the modul bin_info2
#

import numpy as np

# class to store bin information for root histograms
class bin_info:
    def __init__(self, nx = 0, ny = 0, dx = 1., dy = 1., xmin = 0., ymin = 0.):
        self.il = []
        self.ixl = []
        self.iyl = []
        self.xbl = []
        self.ybl = []
        self.contl = []
        self.dcontl = []
        self.weightl = []
        self.weight2l = []
        self.nx = nx
        self.ny = ny
        self.dx = dx
        self.dy = dy
        self.xmin = xmin
        self.ymin = ymin
        self.good_parameters = False
    def add(self, i, ix, iy, xv, yv, cont, dcont, weight ):
        self.il.append(i)
        self.ixl.append(ix)
        self.iyl.append(iy)
        self.xbl.append(xv)
        self.ybl.append(yv)
        self.contl.append(cont)
        self.dcontl.append(dcont)
        self.weightl.append(weight)
        self.weight2l.append(weight**2)
    def copy_from(self,from_b):
        self.i = np.copy( from_b.i)
        self.ix = np.copy( from_b.ix)
        self.iy = np.copy( from_b.iy)
        self.xb = np.copy(from_b.xb)
        self.yb = np.copy(from_b.yb)
        self.cont = np.copy( from_b.cont)
        self.dcont = np.copy(from_b.dcont)
        self.weight = np.copy(from_b.weight)
        self.weight2 = np.copy(from_b.weight2)
        self.nx = from_b.nx
        self.ny = from_b.ny
        self.dx = from_b.dx
        self.dy = from_b.dy
        self.xmin = from_b.xmin
        self.ymin = from_b.ymin
    def make_arrays(self, reshape = False):
        if (not self.good_parameters):
            print  'histogram parameters not properly set !'
            return
        # make 2d arrays according to shape
        a_type = np.zeros((self.nx, self.ny))
        # flatten the arrays to get 1d arrays
        a_type_flat = a_type.flatten()
        # setup the 1d arrays
        self.counter = np.ones_like(a_type_flat)
        self.i = np.zeros_like(a_type_flat)
        self.ix = np.zeros_like(a_type_flat)
        self.iy = np.zeros_like(a_type_flat)
        self.xb = np.zeros_like(a_type_flat)
        self.yb = np.zeros_like(a_type_flat)
        # fill the binning information, ignoring the overflow bins
        # of the root histogram
        for ix in range(1, self.nx+1):
            for jy in range(1, self.ny+1):
                # root bin number
                i_bin = jy*(self.nx + 2) + ix
                # flat array index
                ib = self.ny*(ix-1) + (jy-1)
                self.i[ib] = i_bin
                self.ix[ib] = ix
                self.iy[ib] = jy
                # fill the bin centers
                self.xb[ib] = self.xmin + (ix-0.5)*self.dx
                self.yb[ib] = self.ymin + (jy-0.5)*self.dy 
        self.cont = np.zeros_like(a_type_flat)
        self.dcont = np.zeros_like(a_type_flat)
        self.weight = np.zeros_like(a_type_flat)
        self.weight2 = np.zeros_like(a_type_flat)
        # fill the 1d (flat) arrays
        for i, i_bin in enumerate(self.il):
            # flat array index calculation
            ix = self.ixl[i]
            iy = self.iyl[i]
            ib = self.ny*(ix-1) + (iy-1)
            self.i[ib] = i_bin
            self.ix[ib] = ix
            self.iy[ib] = iy
            self.xb[ib] = self.xbl[i]
            self.yb[ib] = self.ybl[i]
            self.cont[ib] = self.contl[i]
            self.dcont[ib] = self.dcontl[i]
            self.weight[ib] = self.weightl[i]
            self.weight2[ib] = self.weight2l[i]
        # if desired reshape the arrays
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
        self.weight2 = np.ones_like(another.weight2)
        self.nx = another.nx
        self.ny = another.ny
        self.dx = another.dx
        self.dy = another.dy
        self.xmin = another.xmin
        self.ymin = another.ymin
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
        self.weight2 = np.zeros_like(another.weight2)
        self.nx = another.nx
        self.ny = another.ny
        self.dx = another.dx
        self.dy = another.dy
        self.xmin = another.xmin
        self.ymin = another.ymin
#end of class definition
#----------------------------------------------------------------------

# get histrogram information from a datafile
# it is important that all information has been stored

# histo is originally a root histogram
def get_histo_data(histo, bincontent = None, binerror = None, weight = None):
    # histo is a datafile containing the histogram information
    # part of the information (that which does not change) is stored in the parameter section
    # read the relevant data
    data=bin_info()
    # get  parameters
    try:
        # x-axis
        data.nx = histo.par.get_value('nx', int)
        data.dx = histo.par.get_value('dx')
        data.xmin = histo.par.get_value('xmin')
        # y-axis
        data.ny = histo.par.get_value('ny', int)
        data.dy = histo.par.get_value('dy')
        data.ymin = histo.par.get_value('ymin')
        data.good_parameters = True
    except:
        print 'parameter data missing, will not be able create arrays'
        print 'you need to set the binning parameters and set good_parameters to True'
        data.good_parameters = False
    # get the rest of the data
    for d in histo.data:
        i = d['ib']
        ix  = d['ix']
        iy  = d['iy']
        xv = d['xb']
        yv = d['yb']
        if bincontent == None:
            cont = 0.
        else:
            cont = d[bincontent]
        if binerror == None:
            dcont = 0.
        else:
            dcont = d[binerror]
        if weight == None:
            if (dcont > 0.):
                weight_val = 1./dcont**2
            else:
                weight_val = 1.
        else:
            weight_val = d[weight]
        data.add(i, ix, iy, xv, yv, cont, dcont, weight_val)
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