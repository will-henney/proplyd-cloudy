"""
Plots of bivariate distributions

Creates images of two-dimensional joint probability density functions
18 Aug 2010
"""

import numpy as N
import scipy as S
import scipy.stats as stats
from PIL import Image
import pyx

printextra = 0

class PlotVariable(object):
    """
    A variable that we might want to plot
    """
    drawzero = 0
    title = 'Variable name, units'
    shorttitle = 'Variable'
    n = None
    min = None
    max = None

    badvalue = 0.0

    def __init__(self, data, drawzero=0):
	self.drawzero = drawzero
	self.data = data
	self.mask = self.data != self.badvalue

    def settitle(self, title, shorttitle=None):
	self.title = title
	if shorttitle: self.shorttitle = shorttitle

    def setminmaxn(self, min=None, max=None, n=None, logarithmic=0):
	"n is number of bins to use for the images"
	if not n is None: self.n = n
	self.min = min
	self.max = max
	self.logarithmic = logarithmic
	if self.logarithmic:
	    self.min = N.log10(self.min)
	    self.max = N.log10(self.max)
	    # this will fail if setminmaxn has already been called
	    self.data = N.log10(self.data)
	self.setmask((self.data > self.min) & (self.data <= self.max))

    def setmask(self, mask):
	"Adds mask to the data mask"
	self.mask = self.mask & mask

# end class PlotVariable

class Graph(pyx.graph.graphxy):

    epsilon = 1.e-4		# misses of axis extremes
    labelsonx = True
    labelsony = True
    figwidth = 10
    figheight = 10

    def __init__(self, xvar, yvar, weights=[None, None, None], gamma=1.0, statslevel=1):
	"""
	Create an individual pyx graph of one variable vs another

	Arguments:
		   xvar - x variable (instance of PlotVariable)
		   yvar - y variable (instance of PlotVariable)

	Derived from graphxy, so can be inserted in a canvas
	"""

	# Mask out any outlying points
	# This makes all the difference!
	xvar.data = xvar.data[xvar.mask & yvar.mask]
	yvar.data = yvar.data[xvar.mask & yvar.mask]
        # We now have three separate weights for the RGB channels
        self.weights = []
        for w in weights:
            if w is None:
                self.weights.append(None)
            else:
                self.weights.append(w[xvar.mask & yvar.mask])

	if statslevel > 0:
	    for var in (xvar, yvar):
		if printextra: print "Range of %s is [%0.4f, %0.4f]" % (var.shorttitle,
									var.data.min(),
									var.data.max())
		print "Mean +/- sigma of %s is %0.4f +/- %0.4f" % (var.shorttitle,
								   var.data.mean(),
								   var.data.std())
		if printextra: print "Median of %s is %0.4f" % (var.shorttitle,
								N.median(var.data))

	    print "Correlation coefficient of %s with %s is %0.4f" % (xvar.shorttitle,
								      yvar.shorttitle,
								      N.corrcoef(xvar.data, yvar.data)[0,1])
	    if printextra: 
                a, b = S.polyfit(xvar.data, yvar.data, 1)
                print "Linear fit is %s = %0.4f %s + %0.4f" % (yvar.shorttitle, a,
                                                               xvar.shorttitle, b)

	    print "Mean +/- sigma of (%s - %s) is %0.4f +/- %0.4f" % (yvar.shorttitle,
								      xvar.shorttitle,
								      (yvar.data-xvar.data).mean(),
								      (yvar.data-xvar.data).std())

	    if printextra: 
                quart25 = stats.scoreatpercentile(yvar.data - xvar.data, 25)
                quart50 = stats.scoreatpercentile(yvar.data - xvar.data, 50)
                quart75 = stats.scoreatpercentile(yvar.data - xvar.data, 75)
                print "Median of (%s - %s) is %0.4f or %0.4f" % (xvar.shorttitle,
                                                                 yvar.shorttitle,
                                                                 N.median(yvar.data - xvar.data),
                                                                 quart50)
                print "1st & 3rd quartiles of (%s - %s) are %0.4f and %0.4f (diff  %0.4f)" % (
                    xvar.shorttitle, yvar.shorttitle,
                    quart25, quart75, quart75 - quart25)

	    print

	# Make a grid containing histogram of how many pixels there are
	# with each combination of variables
        # Now one for each RGB channel 24 May 2010
        channels = []
        for w in self.weights:
            xygrid, xedges, yedges = N.histogram2d(
                xvar.data.flatten(), yvar.data.flatten(),
                bins=[xvar.n, yvar.n], range=[[xvar.min, xvar.max], [yvar.min, yvar.max]],
                normed=True, weights=w)
            # Turn it into an image, using PIL
            im = Image.new(mode='L', size=(xvar.n, yvar.n))
            # Make positive image
            xyscale = xygrid.flatten()/xygrid.max()
            if gamma != 1.0: xyscale = xyscale**(1./gamma)
            im.putdata(xyscale, scale=255.0)
            # note that transpose does not work in-place!
            channels.append(im.transpose(Image.ROTATE_90))

        ## Now combine the 3 RGB channels
        rgbim = Image.merge('RGB', channels)


	# Now make a graph to return
	if self.labelsonx:
	    xpainter = pyx.graph.axis.painter.regular()
	else:
	    xpainter = pyx.graph.axis.painter.linked()
	if self.labelsony:
	    ypainter = pyx.graph.axis.painter.regular()
	else:
	    ypainter = pyx.graph.axis.painter.linked()

	if xvar.logarithmic:
	    xaxis = pyx.graph.axis.logarithmic(title=xvar.title, painter=xpainter,
					       min=10**xvar.min+self.epsilon,
					       max=10**xvar.max-self.epsilon)
	else:
	    xaxis = pyx.graph.axis.linear(title=xvar.title, painter=xpainter,
					  min=xvar.min+self.epsilon, max=xvar.max-self.epsilon)
	if yvar.logarithmic:
	    yaxis = pyx.graph.axis.logarithmic(title=yvar.title, painter=ypainter,
					       min=10**yvar.min+self.epsilon,
					       max=10**yvar.max-self.epsilon)
	else:
	    yaxis = pyx.graph.axis.linear(title=yvar.title, painter=ypainter,
					  min=yvar.min+self.epsilon, max=yvar.max-self.epsilon)

	pyx.graph.graphxy.__init__(self, width=self.figwidth, height=self.figheight,
				   x=xaxis, y=yaxis)
	self.insert(pyx.bitmap.bitmap(0, 0, rgbim, width=self.figwidth, height=self.figheight))
	# add the zero-velocity line where needed
	self.dolayout()
	zerolinestyle = [pyx.color.transparency(0.5),
			 pyx.style.linestyle.solid,
			 pyx.style.linewidth.Thin]
	if xvar.drawzero: self.stroke(self.xgridpath(0), zerolinestyle)
	if yvar.drawzero: self.stroke(self.ygridpath(0), zerolinestyle)

# End class Graph
