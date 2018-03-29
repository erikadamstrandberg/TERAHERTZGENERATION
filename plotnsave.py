
# coding: utf-8

import matplotlib as matlib
matlib.use('Agg')                    # This allows the standalone application to plot and save figs. Don't move.
import matplotlib.pyplot as mplot
import numpy as np
import csv
from plog import plog

from datetime import datetime

def plotnsave(x, y, plotargs = '', filename = '', dirname = '', savetext = True, savepic = False, inputfile = '', axis = np.array([]), ylab =  np.array([]), xlab =  np.array([]), xticks =  np.array([]), yticks =  np.array([]), cols = False):
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
            if x.any():
                mplot.plot(x, y)
            else:
                mplot.plot(y)
    if axis.any():
        if len(axis) == 4:
            mplot.axis(axis)
        else:
            plog('Axis has the wrong dimensions (len(axis) = ' + str(len(axis)) + '), ignoring.')
    if xlab.any():
        mplot.xlabel(xlab)
    if ylab.any():
        mplot.ylabel(ylab)
    if xticks.any():
        mplot.xticks(xticks)
    if yticks.any():
        mplot.yticks(yticks)
    mplot.yscale('log')
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
        data = np.genfromtxt(filename + '.csv')
        transpose = input('Transpose input? y/n')
        rowcol = 'row'
        if transpose.lower() == 'y':
            rowcol = 'col'
            data = np.transpose(data)
        if len(np.shape(data)) != 1:
            plog('File has {0} {1}. Which row to plot?'.format(np.shape(data)[0], rowcol))
            row = input('')
        try:
            row = int(row)
            plotnsave(data[row], filename = filename, savepic = True, savetext = False)
        except ValueError:
            plog('Row should be an integer. Terminating.')
            return 0
    except IOError:
        plog('Couldn\'t find or open file.')
        
if __name__ == '__main__':
    collectdata()
