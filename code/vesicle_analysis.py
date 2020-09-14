import pandas as pd
import numpy as np
import warnings
from importlib import import_module


class Sample():
    ''' Sample with associated 2D and 3D vesicle data '''
    def __init__(self, sample_name):
        self.sample = sample_name
        self.vsd_num = 0; self.vsd_index = {}; self.vsd = []
        self.image_num = 0; self.image_index = {}; self.images = []
        self.scan_num = 0; self.scan_index = {}; self.scans = []
        self.roi_area = 0
        self.image_data = pd.DataFrame()
    
    def add_image(self, image_id, file, meta_file, units, min_diameter):
        ''' Add image data (ImageJ output and metadata) for the sample '''
        self.images.append(Image(file, meta_file, units, min_diameter))
        self.image_data = pd.concat([self.image_data,
                                     self.images[self.image_num].image_data])
        self.roi_area += self.images[self.image_num].roi_area
        self.image_index.update({image_id: self.image_num})
        self.image_num += 1
     
    def add_scan(self, scan_id, i3d_file, ctan_file, min_diameter):
        self.scans.append(Scan(i3d_file, ctan_file, min_diameter))
        self.scan_index.update({scan_id: self.scan_num})
        self.scan_num += 1
        
    def apply_vsd_corr(self, key, method, length_type, nbins):
        ''' 
        Apply stereological correction to the cumulative image data 
        
        Parameters
        ----------
        key : str, identifying key for vsd correction
        method : str, 'ChengLemlich', 'SahagianProussevitch', or 'Saltikov'
        length_type: str, whether 'radius' or 'diameter' is used
        nbins: int, number of bins for the analysis
        
        '''
        stereology = import_module('stereology')
        self.vsd_index.update({key: self.vsd_num})
        self.vsd.append(getattr(stereology, method)(vesicles=
                                            self.image_data[length_type],
                                            roi_area=self.roi_area,
                                            length_type=length_type,
                                            nbins=nbins))
        self.vsd_num += 1

class Scan():
    '''
    3D data from x-ray microtomography scans of basalts
    
    Parameters
    ----------
    i3d_file : str, path to file with results from 3D object analysis (i3d)
    ctan_file : str, path to  file with results from slice analysis (CTAn)
    min_diameter: float or int, minimum vesicle size to include in analysis
    '''
    
    def __init__(self, i3d_file, ctan_file, min_diameter):
        self.i3d_file = i3d_file
        self.ctan_file = ctan_file
        self.params = {'scaled': False,
                       'thresholded': False,
                       'min_diameter': min_diameter,
                       'units': 'micrometers'}
        self.load_ctan_data()
        self.load_i3d_data()
    
    def load_ctan_data(self):
        ind = 0
        with open(self.ctan_file, 'rt', encoding='latin1') as file:
            for ind, ln in enumerate(file,1):
                if ln.startswith('Number of images inside VOI'):
                    self.params.update({'n_slices': 
                                        float(ln.split(',')[1].strip())})
                if ln.startswith('Summary 2D data'):
                    file.readline()
                    self.params.update({'voi_volume': float(
                        file.readline().strip().split(',')[2])})
                    file.readline().strip().split(',')
                    self.params.update({'vesicularity': float(
                        file.readline().strip().split(',')[2])})
                if ln.startswith('File name,Z position'):
                    start_ind = ind + 6
        self.ctan_data = pd.read_csv(self.ctan_file,
                                     skiprows = start_ind,
                                     header=None,
                                     nrows = self.params['n_slices'],
                                     usecols=[0, 1, 2, 3, 4, 5, 9, 10])
        self.ctan_data.columns = \
            ['image_number','z_position','num_objects','total_ROI_area',
             'object_area','percent_object_area','mean_vesicle_area',
             'mean_vesicle_diameter']  
    
    def load_i3d_data(self):
        self.i3d_data = pd.read_csv(self.i3d_file,
                                    skiprows=[0,1,2,4,5],
                                    usecols=[1,2,7,8,9,11,14])
        self.i3d_data.columns = ['volume','surface','x','y','z','diameter',
                                 'sphericity']
        self.i3d_data['radius'] = self.i3d_data['diameter']/2
        #self.i3d_data.columns = 
        
    def threshold_data(self):
        self.i3d_data = self.i3d_data[self.i3d_data['Volume-equivalent sphere diameter'] \
                          > self.params['min_diameter']]
        self.params['thresholded'] = True
        
    
class Image():
    ''' 
    2D image data for basalt samples 
    
    Parameters
    ----------
    file : str, path to the csv file with results from ImageJ analysis
    meta_file : str, path to text file with roi size and pixel resolution
    units : str, units for metadata (usually micrometers)
    min_diameter: float or int, minimum vesicle size to include in analysis
    '''
    
    def __init__(self, file, meta_file, units, min_diameter):
        self.file = file
        self.meta_file = meta_file
        self.params = {'units': units,
                       'scaled': False,
                       'thresholded': False,
                       'min_diameter': min_diameter}
        self.vsd = []
        self.vsd_index = {}
        self.vsd_num = 0
        self.load_data()
        self.load_metadata()
        self.scale_data()
        self.threshold_data()
        
    def load_data(self):
        ''' Load results from ImageJ analyses of vesicles '''
        self.image_data = pd.read_csv(self.file, header=None)
        self.image_data.columns = ['area', 
                                   'x',
                                   'y',
                                   'perimeter',
                                   'bx',
                                   'by',
                                   'width',
                                   'height',
                                   'major_axis',
                                   'minor_axis',
                                   'angle',
                                   'xstart',
                                   'ystart']
        self.image_data['radius'] = np.sqrt(self.image_data['area'] / np.pi)
        self.image_data['diameter'] = np.multiply(self.image_data['radius'],2)
        
    def scale_data(self):
        ''' ImageJ results are in pixels, apply unit conversion from metadata'''
        if not self.params['scaled']:
            self.image_data['area'] = self.image_data['area'] \
                * self.pixel_size ** 2
            for measurement in ['x','y','perimeter','bx','by','width','height',
                                'major_axis','minor_axis','xstart','ystart',
                                'radius','diameter']:
                self.image_data[measurement] = self.image_data[measurement] * \
                    self.pixel_size
            self.params['scaled'] = True
        else:
            warnings.warm('Image data has already been scaled, scale not applied',
                          RuntimeWarning, stacklevel=2)
            
    def load_metadata(self):
        ''' Load metadata (roi size and image pixel resolution) '''
        with open(self.meta_file) as mf:
            mf.readline()
            self.pixel_size = np.array(
                mf.readline().strip().split(',')[1]).astype(float)
            self.roi_area = np.array(
                mf.readline().strip().split(',')[1]).astype(float)
        
    def threshold_data(self):
        ''' Threshold based on the minimum allowed vesicle size '''
        self.image_data = self.image_data[self.image_data['diameter'] \
                          > self.params['min_diameter']]
        self.params['thresholded'] = True
            
            
    def apply_vsd_corr(self, key, method, length_type, nbins):
        ''' 
        Apply stereological correction to the data 
        
        Parameters
        ----------
        key : str, identifying key for vsd correction
        method : str, 'ChengLemlich', 'SahagianProussevitch', or 'Saltikov'
        length_type: str, whether 'radius' or 'diameter' is used
        nbins: int, number of bins for the analysis
        
        '''
        stereology = import_module('stereology')
        self.vsd_index.update({key: self.vsd_num})
        self.vsd.append(getattr(stereology, method)(vesicles=
                                            self.image_data[length_type],
                                            roi_area=self.roi_area,
                                            length_type=length_type,
                                            nbins=nbins))
        self.vsd_num += 1
