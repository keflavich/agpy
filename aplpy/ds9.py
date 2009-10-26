

import numpy as np
import matplotlib as mpl
import matplotlib.patches as mpp
import coords as pywcs
import re
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

        fontre = re.compile('font="(.*?)"')
        
        for line in self.file[4:]:
            
            shape = {}
            
            p1 = line.find("(")
            p2 = line.find(")")
            
            line_prefix = line[:p1]
            line_coords = line[p1+1:p2]
            line_format = line[p2:].replace("# ","")
            line_format = line_format[2:]
            format = line_format.split(" ")
            fontsform = fontre.search(line)
            
            shape['edgecolor'] = "g"
            shape['lw'] = None
            shape['linestyle'] = 'solid'
            shape['fontsize'] = 'large'
            shape['fontstyle'] = 'normal'
            shape['fontweight'] = 'normal'
            shape['zorder'] = 2  # default shapes to draw after whatever draws first, probably contours
            
            for elem in format:
                arr = elem.split("=")
                if arr[0] == 'color':
                    shape['edgecolor'] = arr[1]
                if arr[0] == 'width':
                    shape['lw'] = float(arr[1])
                if arr[0] == 'dash':
                    if arr[1] == "1": shape['linestyle'] = "dashed"
                if arr[0] == 'point':
                    shape['point'] = arr[1]
                if arr[0] == 'zorder':
                    shape['zorder'] = int(arr[1])
                if arr[0] == 'text':
                    L=line_format.find("{")+1
                    R=line_format.find("}")
                    shape['text'] = line_format[L:R]
                if fontsform is not None:
                    fonts = fontsform.groups()[0].split()
                    shape['fontname'] = fonts[0]
                    shape['fontsize'] = fonts[1]
                    if fonts[2] == 'bold':
                        shape['fontweight'] = 'bold'
                    if fonts[2] == 'italic':
                        shape['fontstyle'] = 'italic'
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
                    pos = pywcs.Position((float(coords[0]),float(coords[1])))

                exec("ra1,dec1 = pos."+coordfunc+"()")
                x1,y1=wcs_util.world2pix(wcs,ra1,dec1)

            shape['x'] = x1
            shape['y'] = y1
            
            if shape['type'] == 'circle':
                shape['radius'] = float(coords[2].replace('"',''))
            elif shape['type'] == 'box':
                width  = float(coords[2].replace('"','')) 
                height = float(coords[3].replace('"','')) 
                angle = float(coords[4])
                cosa = cos(angle/180.0 * np.pi)
                sina = sin(angle/180.0 * np.pi)
                dxm = (width/2.0*cosa-height/2.0*sina)/wcs_util.arcperpix(wcs)
                dym = (width/2.0*sina-height/2.0*cosa)/wcs_util.arcperpix(wcs)
                dxp = (width/2.0*cosa+height/2.0*sina)/wcs_util.arcperpix(wcs)
                dyp = (width/2.0*sina+height/2.0*cosa)/wcs_util.arcperpix(wcs)
                xb1,yb1 = x1 + dxm, y1 + dyp
                xb2,yb2 = x1 + dxp, y1 + dym
                xb3,yb3 = x1 - dxm, y1 - dyp
                xb4,yb4 = x1 - dxp, y1 - dym
                shape['box'] = [ [xb1,yb1] , [xb2,yb2] , [xb3,yb3], [xb4,yb4] ]
                if len(coords) == 7:
                    shape['width2'] = float(coords[4].replace('"',''))
                    shape['height2'] = float(coords[5].replace('"',''))
                    print "box panda not implemented"
                    continue
                #print x1,y1,dx,dy,shape,angle
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
            elif shape['type'] == 'ellipse':
                shape['width']  = float(coords[2].rstrip('\'"')) * 2 / wcs_util.arcperpix(wcs)
                shape['height'] = float(coords[3].rstrip('\'"')) * 2 / wcs_util.arcperpix(wcs)
                shape['angle'] = float(coords[4])
            elif shape['type'] == 'vector':
                veclen   = float(coords[2].rstrip('\'"')) / wcs_util.arcperpix(wcs)
                vecangle = (float(coords[3])) / 180.0 * np.pi
                dx,dy=veclen*cos(vecangle),veclen*sin(vecangle)
                shape['dx'] = float(dx)
                shape['dy'] = float(dy)
            elif shape['type'] == 'polygon':
                print "polygons are not implemented yet."
                continue
                vertices = np.array(coords,float)
                n = vertices.size / 2
                vertices = vertices.reshape((2,n),order='F')
                shape['v'] = vertices.tolist()
            
            self.shapes.append(shape)
    
    def plot(self,ax,**kwargs):
        
        patches = []
        collections = []
        fill = [] # HACK to deal with http://www.nabble.com/Alpha-values-not-preserved-with-Ellipse-and-PatchCollection-td23408210.html#a23408210

        headsize = np.abs(np.diff(ax.get_xlim()))[0] * .01
        
        for shape in self.shapes:
            
            if shape['type'] == 'circle':
                #shape = dict2pix(shape,self._wcs)
                patches.append(circle_patch(**shape))
                fill.append(0)
            elif shape['type'] == 'ellipse':
                shape = dict2pix(shape,self._wcs)
                patches.append(ellipse_patch(**shape))
                fill.append(0)
            elif shape['type'] == 'box':
                #shape = dict2pix(shape,self._wcs)
                patches.append(box_patch2(**shape))
                fill.append(0)
            elif shape['type'] == 'line':
                #shape = dict2pix(shape,self._wcs)
                patches.append(line_patch(**shape))
                fill.append(0)
            elif shape['type'] == 'vector':
                #shape = dict2pix(shape,self._wcs)
                patches.append(mpl.patches.FancyArrow(shape['x'],shape['y'],shape['dx'],shape['dy'],
                        edgecolor=shape['edgecolor'],facecolor=shape['edgecolor'],
                        head_width=headsize,**kwargs))
                fill.append(1.0)
            elif shape['type'] == 'point':
                if shape['point'] == 'circle':
                    collections.append(matplotlib.pyplot.scatter(shape['x'],shape['y'],marker='o',edgecolor=shape['edgecolor'],facecolor='none',zorder=shape['zorder'],**kwargs))
                if shape['point'] == 'box':
                    collections.append(matplotlib.pyplot.scatter(shape['x'],shape['y'],marker='s',edgecolor=shape['edgecolor'],facecolor='none',zorder=shape['zorder'],**kwargs))
                if shape['point'] == 'diamond':
                    collections.append(matplotlib.pyplot.scatter(shape['x'],shape['y'],marker='d',edgecolor=shape['edgecolor'],facecolor='none',zorder=shape['zorder'],**kwargs))
                if shape['point'] == 'cross':
                    collections.append(matplotlib.pyplot.scatter(shape['x'],shape['y'],marker='+',edgecolor=shape['edgecolor'],facecolor='none',zorder=shape['zorder'],**kwargs))
                if shape['point'] == 'x':
                    collections.append(matplotlib.pyplot.scatter(shape['x'],shape['y'],marker='x',edgecolor=shape['edgecolor'],facecolor='none',zorder=shape['zorder'],**kwargs))
                if shape['point'] == 'arrow':
                    collections.append(matplotlib.pyplot.scatter(shape['x'],shape['y'],marker='>',edgecolor=shape['edgecolor'],facecolor='none',zorder=shape['zorder'],**kwargs))
                if shape['point'] == 'boxcircle':
                    collections.append(matplotlib.pyplot.scatter(shape['x'],shape['y'],marker='s',edgecolor=shape['edgecolor'],facecolor='none',zorder=shape['zorder'],**kwargs))
                    collections.append(matplotlib.pyplot.scatter(shape['x'],shape['y'],marker='o',edgecolor=shape['edgecolor'],facecolor='none',zorder=shape['zorder'],**kwargs))
            if shape['type'] == 'text':
                matplotlib.pyplot.annotate(shape['text'],(shape['x'],shape['y']),color=shape['edgecolor'],ha='center',va='center',fontsize=shape['fontsize'],fontweight=shape['fontweight'],fontstyle=shape['fontstyle'],zorder=shape['zorder'],**kwargs)
            elif shape.has_key('text'):
                matplotlib.pyplot.annotate(shape['text'],(shape['x'],shape['y']),color=shape['edgecolor'],ha='center',va='center',fontsize=shape['fontsize'],fontweight=shape['fontweight'],fontstyle=shape['fontstyle'],zorder=shape['zorder'],**kwargs)

        
        # DO NOT USE COLLECTIONS!  Collections do not preserve linestyle or opacity.
        # Hack to prevent facecolors from being drawn
        #PC = matplotlib.collections.PatchCollection(patches,match_original=True)
        #PC.set_facecolors(PC.get_facecolors()*np.array([fill]).T)
        #PC._is_filled=False
        #import pdb; pdb.set_trace()
        return patches #PC
