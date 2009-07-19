
def ds9(self,regionfile,wcs):
    
    reg = RegionFile()
    reg.parse(regionfile,wcs)
    
    for p in reg.plot():
        self._ax1.add_patch(p)
    
    self.refresh()

def clear_markers(self):
    nmarks = len(self._ax1.patches)
    ncolls = len(self._ax1.collections)
    for i in xrange(nmarks):
        self._ax1.patches.pop()
    for i in xrange(ncolls):
        self._ax1.collections.pop()

    self.refresh()

import numpy as np
import matplotlib as mpl
import matplotlib.patches as mpp
import coords as pywcs
from ds9_util import *

def ds9type(line):
    return line.replace("# ","")

class RegionFile(object):
    
    ds9Dict = {'box':'',}
    
    """
    This class parses DS9 region files. This is used
    for the main APLpy class.
    """
    
    def __init__(self):
        self.coord_sys=None
        pass
    
    def parse(self, fileName, wcs):
        try:
            f = open(fileName,'r')
            self.file = f.read().splitlines()
            f.close()
        except:
            print "This is not a proper ds9 region file. It's spam!"
            raise NameError('Wrong ds9 regions file type')
            return
        
        self.shapes = []
        
        self._wcs = wcs
        for i in xrange(len(self.file)):
            if self.file[i][0] == "#":
                continue
            if self.file[i].split()[0] == 'global':
                line = self.file[i].split()
                """
                INSERT GLOBALS PARSING HERE
                """
            if self.file[i] in  ['fk5','icrs','galactic','image','fk4']:
                self.coord_sys = self.file[i]
                break
            elif self.file[i] in ['physical','ecliptic']:
                raise NotImplementedError("We don't support 'physical' or 'ecliptic' coordinate systems")
        if self.coord_sys is None:
            raise EOFError("Did not find a coordinate system specification in your regions file.  " + \
                    "Valid systems are: fk4,fk5,icrs,galactic,physical,image,ecliptic")
        
        for line in self.file[4:]:
            
            shape = {}
            
            p1 = line.find("(")
            p2 = line.find(")")
            
            line_prefix = line[:p1]
            line_coords = line[p1+1:p2]
            line_format = line[p2:].replace("# ","")
            line_format = line_format[2:]
            format = line_format.split(" ")
            
            shape['edgecolor'] = "g"
            shape['lw'] = None
            shape['ls'] = 'solid'
            
            for elem in format:
                arr = elem.split("=")
                if arr[0] == 'color':
                    shape['edgecolor'] = arr[1]
                if arr[0] == 'width':
                    shape['lw'] = arr[1]
                if arr[0] == 'dash':
                    if arr[1] == "1": shape['ls'] = "dashed"
                if arr[0] == 'point':
                    shape['point'] = arr[1]
                if arr[0] == 'text':
                    L=line_format.find("{")+1
                    R=line_format.find("}")
                    shape['text'] = line_format[L:R]
            shape['type'] = ds9type(line_prefix)
            shape['coord_sys'] = self.coord_sys
            coords = line_coords.split(",")

            # lines necessary to pick conversion type for pywcs coords package
            # (need to use pywcs coords to parse sexagesimal)
            cosys,epoch,degrees=wcs_util.system(wcs)
            if cosys == 'celestial':
                coordfunc = epoch
            else:
                coordfunc = cosys

            # ALL regions have first two values as coordinates
            # may be sexagesimal or degree
            if self.coord_sys == 'image':
                sexagesimal=False
                x1 = float(coords[0])
                y1 = float(coords[1])
            else:
                if coords[0].find(":") != -1:
                    sexagesimal=True
                    pos = pywcs.Position(coords[0]+" "+coords[1])
                else:
                    sexagesimal=False
                    pos = pywcs.Position((coords[0],coords[1]))

                exec("ra1,dec1 = pos."+coordfunc+"()")
                x1,y1=wcs_util.world2pix(wcs,ra1,dec1)

            shape['x'] = x1
            shape['y'] = y1
            
            if shape['type'] == 'circle':
                shape['radius'] = float(coords[2].replace('"',''))
            elif shape['type'] == 'box':
                shape['width'] = float(coords[2].replace('"',''))
                shape['height'] = float(coords[3].replace('"',''))
                shape['angle'] = float(coords[4])
            elif shape['type'] == 'line':
                if self.coord_sys == 'image':
                    x2 = float(coords[2])
                    y2 = float(coords[3])
                else:
                    if sexagesimal:
                        pos2 = pywcs.Position(coords[2]+" "+coords[3])
                    else:
                        pos2 = pywcs.Position((coords[2],coords[3]))
                    exec("ra2,dec2 = pos2."+coordfunc+"()")
                    x2,y2=wcs_util.world2pix(wcs,ra2,dec2)
                shape['x'] = [x1,x2]
                shape['y'] = [y1,y2]
            elif shape['type'] == 'point':
                pass
            elif shape['type'] == 'vector':
                veclen   = float(coords[2].rstrip('\'"')) / wcs_util.arcperpix(wcs)
                vecangle = (float(coords[3])) / 180.0 * np.pi
                dx,dy=veclen*cos(vecangle),veclen*sin(vecangle)
                shape['dx'] = float(dx)
                shape['dy'] = float(dy)
            elif shape['type'] == 'polygon':
                vertices = np.array(coords,float)
                n = vertices.size / 2
                vertices = vertices.reshape((2,n),order='F')
                shape['v'] = vertices.tolist()
            
            self.shapes.append(shape)
    
    def plot(self,**kwargs):
        
        patches = []
        collections = []
        
        for shape in self.shapes:
            
            if shape['type'] == 'circle':
                #shape = dict2pix(shape,self._wcs)
                patches.append(circle_patch(**shape))
            elif shape['type'] == 'box':
                #shape = dict2pix(shape,self._wcs)
                patches.append(box_patch(**shape))
            elif shape['type'] == 'line':
                #shape = dict2pix(shape,self._wcs)
                patches.append(line_patch(**shape))
            elif shape['type'] == 'vector':
                #shape = dict2pix(shape,self._wcs)
                patches.append(mpl.patches.FancyArrow(shape['x'],shape['y'],shape['dx'],shape['dy'],
                        edgecolor=shape['edgecolor'],facecolor=shape['edgecolor'],
                        head_width=8,**kwargs))
            elif shape['type'] == 'point':
                if shape['point'] == 'circle':
                    #patches.append(mpp.Circle(shape['x'],shape['y'],radius=5,**kwargs))
                    collections.append(matplotlib.pyplot.scatter(shape['x'],shape['y'],marker='o',edgecolor=shape['edgecolor'],facecolor='none',**kwargs))
                if shape['point'] == 'box':
                    #patches.append(mpp.RegularPolygon(shape['x'],shape['y'],4,orientation=0,facecolor='none',**kwargs))
                    collections.append(matplotlib.pyplot.scatter(shape['x'],shape['y'],marker='s',edgecolor=shape['edgecolor'],facecolor='none',**kwargs))
                if shape['point'] == 'diamond':
                    collections.append(matplotlib.pyplot.scatter(shape['x'],shape['y'],marker='d',edgecolor=shape['edgecolor'],facecolor='none',**kwargs))
                    #patches.append(mpp.RegularPolygon(shape['x'],shape['y'],4,orientation=(45.0/180.0*pi),facecolor='none',**kwargs))
                if shape['point'] == 'cross':
                    collections.append(matplotlib.pyplot.scatter(shape['x'],shape['y'],marker='+',edgecolor=shape['edgecolor'],facecolor='none',**kwargs))
                    # agh, this needs to be done with scatter.  
                    """
                    if kwargs.has_key('radius'):
                        r = float(kwargs['radius'])
                    else:
                        r = 5.0
                    x,y = shape['x'],shape['y']
                    cross = mpp.Path([[x-r/2.0,y],[x+r/2.0,y],[x,y-r/2.0],[x,y+r/2.0]],
                            [2,2,1,2],**kwargs)
                    patches.append(cross)
                    """
                if shape['point'] == 'x':
                    collections.append(matplotlib.pyplot.scatter(shape['x'],shape['y'],marker='x',edgecolor=shape['edgecolor'],facecolor='none',**kwargs))
                if shape['point'] == 'arrow':
                    collections.append(matplotlib.pyplot.scatter(shape['x'],shape['y'],marker='>',edgecolor=shape['edgecolor'],facecolor='none',**kwargs))
                if shape['point'] == 'boxcircle':
                    collections.append(matplotlib.pyplot.scatter(shape['x'],shape['y'],marker='s',edgecolor=shape['edgecolor'],facecolor='none',**kwargs))
                    collections.append(matplotlib.pyplot.scatter(shape['x'],shape['y'],marker='o',edgecolor=shape['edgecolor'],facecolor='none',**kwargs))
        
        return patches
