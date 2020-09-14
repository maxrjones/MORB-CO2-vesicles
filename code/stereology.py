from abc import ABC, abstractmethod
import numpy as np
import warnings

class VSD():
    '''
    Vesicle size distribution for 2D images, microCT data, or synthetics
    
    Parameters
    ----------
    vesicles : array, sizes of individual vesicles
    '''
    
    def __init__(self, vesicles):
        self.vesicles = vesicles
        self.bins = Bins(self.vesicles, self.params)
        self.bins.bin_data()
        
    def compute_nv(self, voi_area):
        self.nv = np.zeros_like(self.bin_centers)
        for i in np.arange(self.params['nbins']):
            self.nv[i] = len(self.ind[self.ind == i]) / voi_area
        
    def to_lnn(self):
        divisor = self.bins.bin_widths * 1e-3 * self.bins.norm
        if self.params['length_type'] != 'diameter':
            divisor *= 2
        self.n = np.divide(self.nv * 1e9, divisor)
        self.lnn = np.log(self.n)
        
    def plot_data(self, ax, yvar, color, marker):
        xdata = np.multiply(self.bins.bin_centers,self.bins.norm)
        if self.bins.params['length_type'] != 'diameter':
            xdata *= 2
        ax.plot(xdata, getattr(self,yvar),
                marker=marker,color=color,linestyle='None')

class VSDCorrection(ABC, VSD):
    ''' Base class for vesicles size distribution correction methods '''
    def __init__(self, vesicles, roi_area):
        super().__init__(vesicles)
        self.roi_area = roi_area
        self.nv = np.zeros_like(self.bins.bin_centers)
        self.bins.compute_na(roi_area)
        self.bins.compute_hbar()

    @abstractmethod
    def to_nv(self):
        ''' Convert from number per unit area to number per unit volume '''
        pass

class Bins():
    ''' 
    Bins for characterizing vesicle size distributions and/or applying
    stereological corrections
    
    Parameters
    ----------
    vesicles: array, sizes of individual vesicles
    params : dictionary of parameters used for creating bins
        - 'normalized' : Bool , determines whether normalized sizes are used
        - 'bin_type' : str , determines between 'linear' or 'geometric' bins
        - 'nbins' : int , number of bins to create
        - 'hbar_method' : str , method for determine characteristic size 
                                within each ('bin_center', 'median', 'mean')
                            
    '''
    
    def __init__(self, vesicles, params):
        self.params = params
        self.params['bin_order'] = 'desending'
        if self.params['normalized']:
            self.norm = np.max(vesicles)
            self.vesicles = np.divide(vesicles, self.norm)
        else: 
            self.norm = 1
            self.vesicles = vesicles
        self._create_bins()
        
    def bin_data(self):
        self.ind = np.digitize(self.vesicles, self.bin_edges, right=True) - 1
        self.__check_inrange()
        
    def _create_bins(self):
        ''' Create either linear or geometric bins '''
        if self.params['bin_method'] == 'linear':
            self.bin_edges = np.linspace(np.max(self.vesicles) * 1.0000001,
                                    np.min(self.vesicles) * 0.9999999,
                                    num = self.params['nbins'] + 1)
        elif self.params['bin_method'] == 'geometric':
            self.bin_edges = np.max(self.vesicles) \
                * 10 ** (-0.1 * np.arange(self.params['nbins'] + 1))
        self.bin_widths = abs(self.bin_edges[1:]-self.bin_edges[:-1])
        self.bin_centers = self.bin_edges[:-1] - self.bin_widths / 2
    
    def __check_inrange(self):
        if np.max(self.ind) >= self.params['nbins']:
            warnings.warn(
                '{0} vesicles are smaller than any bin edge'.format(
                    len(self.ind[self.ind >= self.params['nbins']])),
                RuntimeWarning, stacklevel=2)
            
    def compute_na(self, roi_area):
        self.na = np.zeros_like(self.bin_centers)
        for i in np.arange(self.params['nbins']):
            self.na[i] = len(self.ind[self.ind == i]) / roi_area
            
    def compute_hbar(self):
        ''' Determine the characteristic vesicle size for each bin '''
        self.hbar = np.zeros_like(self.bin_centers)
        # Add code to check for accepted method
        if self.params['hbar_method'] == 'bin_center':
            self.hbar = self.bin_centers
        else:
            for i in np.arange(self.params['nbins']):
                in_vesicles = self.vesicles[self.ind == i]
                if not in_vesicles.any():
                    warnings.warn(
                        'Bin {} does not contain vesicles, using bin center'. \
                            format(i), RuntimeWarning, stacklevel=2)
                    self.hbar[i] = self.bin_centers[i]
                else:
                    self.hbar[i] = getattr(np, self.params['hbar_method'])(
                        in_vesicles)
    
    
class ChengLemlich(VSDCorrection):
    ''' 
    Stereological corretion using the methods of Cheng and Lemlich (1983)
    
    Parameters
    ----------
    vesicles: array, sizes of individual vesicles
    roi_area : float, area of the planar section of the sample analyzed
    length_type : str, 'radius' or 'diameter'
    nbins : int, number of bins to use in the VSD correction
    
    '''
    
    def __init__(self, vesicles, roi_area, length_type, nbins):
        if length_type != 'radius':
            warnings.warn('Input length for C&L correction should be radius',
                          RuntimeWarning, stacklevel=2)
        self.params = {'bin_method': 'linear', 'hbar_method': 'mean',
                       'length_type': length_type, 'nbins': nbins,
                       'normalized': False}
        super().__init__(vesicles, roi_area)
        self.to_nv()
        self.to_lnn()
    
    def to_nv(self):
        self.nv = np.divide(self.bins.na, self.bins.hbar)
        
              
class SahagianProussevitch(VSDCorrection):
    ''' 
    Stereological corretion using the methods of Sahagian and 
    Proussevitch (1998)
    
    Parameters
    ----------
    vesicles: array, sizes of individual vesicles
    roi_area : float, area of the planar section of the sample analyzed
    length_type : str, 'radius' or 'diameter'
    nbins : int, number of bins to use in the VSD correction
    
    '''
    
    def __init__(self, vesicles, roi_area, length_type, nbins):
        if length_type != 'diameter':
            warnings.warn('Input length for S&P correction should be diameter',
                          RuntimeWarning, stacklevel=2)
        self.params = {'bin_method': 'geometric', 'hbar_method': 'bin_center',
                       'length_type': length_type, 'nbins': nbins,
                       'normalized': False}
        super().__init__(vesicles, roi_area)
        self.intersection_probabilities()
        self.to_nv()
        self.to_lnn()

        
    def intersection_probabilities(self):
        ''' Intersection probabilities valid for geometric bins only '''
        self.inter_prob = np.zeros_like(self.bins.bin_centers)
        for i in np.arange(self.params['nbins']):
            self.inter_prob[i] = (1 / self.bins.bin_edges[0]) * (
                np.sqrt(self.bins.bin_edges[0] ** 2 - 
                        self.bins.bin_edges[i+1] ** 2) - 
                np.sqrt(self.bins.bin_edges[0] ** 2 - 
                        self.bins.bin_edges[i] ** 2))
            
    def to_nv(self):
        for i in np.arange(self.params['nbins']):
            previous = 0
            for j in np.arange(i):
                previous = previous + (self.inter_prob[j+1] * 
                                       self.bins.hbar[j+1] * 
                                       self.nv[i-j-1])
            self.nv[i] = (1 / (self.inter_prob[0] * self.bins.hbar[i])) \
                         * (self.bins.na[i] - previous)


class Saltikov(VSDCorrection):
    ''' 
    Stereological corretion using the methods of Saltikov (1967)
    
    Parameters
    ----------
    vesicles: array, sizes of individual vesicles
    roi_area : float, area of the planar section of the sample analyzed
    length_type : str, 'radius' or 'diameter'
    nbins : int, number of bins to use in the VSD correction
    
    '''
    
    def __init__(self, vesicles, roi_area, length_type, nbins):
        if length_type != 'diameter':
            warnings.warn('Input length for S correction should be diameter',
                          RuntimeWarning, stacklevel=2)
        self.params = {'bin_method': 'linear', 'hbar_method': 'mean',
                       'length_type': length_type, 'nbins': nbins,
                       'normalized': False}
        super().__init__(vesicles, roi_area)
        self.to_nv()
        self.to_lnn()

    def to_nv(self):
        coeff = [1.6461, -0.4561, -0.1162, -0.0415, -0.0173, -0.0079, 
                -0.0038, -0.0018, -0.0010, -0.0003, -0.0002, -0.0002]
        for i in np.arange(self.params['nbins']):
            na_sum = coeff[0] * self.bins.na[i]
            for k in np.arange(1,i+1):
                na_sum = na_sum + coeff[k] * self.bins.na[i-k]
            self.nv[i] = 1 / self.bins.hbar[i] * na_sum
