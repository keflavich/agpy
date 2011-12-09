
def drizzle(tstomap,ts,weights,mapshape):
    """
    Drizzle a timestream onto a map
    """
    newmap = numpy.zeros(mapshape)
    wm = numpy.zeros(mapshape)
    ts_to_index = masktozero((ts*weights).ravel())
    weights_to_index = masktozero((weights).ravel())

    if len(tstomap.shape) > 1:
        tstomap = tstomap.ravel()
    tsmapped = numpy.bincount(tstomap,ts_to_index)
    tsweights = numpy.bincount(tstomap,weights_to_index)
    if tsmapped.shape[0] < product(mapshape):
        tsmapped = numpy.concatenate([tsmapped,
            numpy.zeros(product(mapshape)-tsmapped.shape[0])])
        tsweights = numpy.concatenate([tsweights,
            numpy.zeros(product(mapshape)-tsweights.shape[0])])
    xinds,yinds = (a.ravel() for a in numpy.indices(mapshape))
    newmap[xinds,yinds] = tsmapped
    wm[xinds,yinds] = tsweights

    return newmap/wm
