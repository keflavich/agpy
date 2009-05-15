from numpy import sqrt,repeat,indices,newaxis

def extract_aperture(cube,apcen,apwidth):
	sh = cube.shape
	yind,xind = indices(sh[1:3]) # recall that python indices are backwards
	dis = sqrt((xind-apcen[0])**2+(yind-apcen[1])**2)
	dis3d = repeat(dis[newaxis,:,:],sh[0],axis=0)
	spec = (cube* (dis3d < apwidth)).sum(axis=2).sum(axis=1)
	return spec



