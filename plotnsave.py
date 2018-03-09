
# coding: utf-8

import matplotlib as matlib
matlib.use('Agg')                    # This allows the standalone application to plot and save figs.
import matplotlib.pyplot as mplot
import numpy as np

from datetime import datetime

def plotnsave(y, plotargs = '', filename = '', dirname = '', savetext = True, savepic = False):
    """ Takes x and y as mandatory arguments and args and filename as optional. Plots y as a function of x. Passes plotargs to matlibplot.pyplot.plot. Saves as a .png file if filename is specified.
    """
    if dirname:
        filename = dirname + '/' + filename
    print(str(datetime.now())+': Beginning plot.')
    if plotargs:
        mplot.plot(y, plotargs)
    else:
        mplot.plot(y)
    if filename:
        if savepic:
            mplot.savefig(filename + '.png')
        if savetext:
            np.savetxt(filename + '.csv', y, delimiter=',')
        print(str(datetime.now())+ ': Plot saved.')
    else:
        print(str(datetime.now())+': Plot complete.')
