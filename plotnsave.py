
# coding: utf-8

import matplotlib as matlib
matlib.use('Agg')                    # This allows the standalone application to plot and save figs.
import matplotlib.pyplot as mplot
import numpy as np

from datetime import datetime

def plotnsave(x, y, plotargs = '', filename = '', savetext = True, savepic = False):
    """ Takes x and y as mandatory arguments and args and filename as optional. Plots y as a function of x. Passes plotargs to matlibplot.pyplot.plot. Saves as a .png file if filename is specified.
    """
    print(str(datetime.now())+': Beginning plot.')
    if plotargs:
        mplot.plot(x, y, plotargs)
    else:
        mplot.plot(x, y)
    if filename:
        if savepic:
            mplot.savefig(filename + '.png')
        if savetext:
            np.savetxt(filename + '_y.csv', y, delimiter=',')
            np.savetxt(filename + '_x.csv', x, delimiter=',')
        print(str(datetime.now())+ ': Plot saved.')
    else:
        print(str(datetime.now())+': Plot complete.')
