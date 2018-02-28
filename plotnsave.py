# coding: utf-8

from matlibplot.plot import plot as mplot

def plotnsave(x, y, plotargs, filename):
    """ Takes x and y as mandatory arguments and args and filename as optional. Plots y as a function of x. Passes args to matlibplot.plot.plot. Saves as a .png file if filename is specified.
    """
    print(str(datetime.now())+': Beginning plot.')
    if args:
        mplot.plot(x, y, plotargs)
    else:
        mplot.plot(x, y)
    if filename:
        mplot.savefig(filename)
        np.savetxt(filename + '_y.csv', y, delimiter=',')
        np.savetext(filename + '_x.csv', x, delimiter=',')
        print(str(datetime.now())+ ': Plot saved.')
    else:
        print(str(datetime.now())+': Plot complete.')
