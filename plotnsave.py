
# coding: utf-8

import matplotlib as matlib
matlib.use('Agg')                    # This allows the standalone application to plot and save figs. Don't move.
import matplotlib.pyplot as mplot
import numpy as np
import csv
from plog import plog

from datetime import datetime

def plotnsave(y, plotargs = '', filename = '', dirname = '', savetext = True, savepic = False, inputfile = '', axis = 0, ylab = '', xlab = ''):
    """ Takes x and y as mandatory arguments and args and filename as optional. Plots y as a function of x. Passes plotargs to matlibplot.pyplot.plot. Saves as a .png file if filename is specified.
    """
    #print(y)
    if dirname:
        filename = dirname + '/' + filename
    print(str(datetime.now())+': Beginning plot.')
    if plotargs:
        mplot.plot(y, plotargs)
    else:
        if len(y) < 10:
            for i in range(len(y)):
                mplot.plot(y[i])
        else:
            mplot.plot(y)
    if axis:
        if len(axis) == 4:
            mplot.axis(axis)
        else:
            plog('Axis has the wrong dimensions (len(axis) = ' + str(len(axis)) + '), ignoring.')
    if xlab:
        mplot.xlabel(xlab)
    if ylab:
        mplot.ylabel(ylab)
            
    if filename:
        if savepic:
            plog('Saving ' + filename + '.png.')
            mplot.savefig(filename + '.png')
            plog('Save successful.')
            mplot.clf()
        if savetext:
            plog('Saving ' + filename + '.csv.')
            np.savetxt(filename + '.csv', y, delimiter=',')
            plog('Save successful.')
        print(str(datetime.now())+ ': Plot saved.')
    else:
        print(str(datetime.now())+': Plot complete.')

def collectdata():
    filename = input('Enter csv file to load (sans extension).\n')
    try:
        with open(filename + '.csv') as f:
            reader = csv.reader(f)
            lister = list(reader)
            row_count = sum(1 for row in lister)
            plog('File has {0} rows. Which row to plot?'.format(row_count))
            row = input('')
            #print(lister)
            try:
                row = int(row)
                plotnsave(np.loadtxt(lister[row], delimiter = ','), filename = filename, savepic = True, savetext = False)
            except ValueError:
                plog('Row should be an integer. Terminating.')
                return 0
    except IOError:
        plog('Couldn\'t find or open file.')
        
if __name__ == '__main__':
    collectdata()
