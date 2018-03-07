from sim1d import runsim, plog
from sys import argv
from ast import literal_eval   # Unsure if this does anything, try commenting out.
import numpy as np

def main():
    ''' Run as sweep.py dt [param] [start] [stop] [-dx or -steps] [step size or steps]
    '''
    optionalargs = [
        'dt', 'dz', 'time', 'size',
        'pulsestart', 'pulselength', 'nu', 'wavelength',
        'ramplength', 'rampdamp', 'plasmastopp', 'plottime',
        'fname']
    param = {}
    
    '''argv contains the arguments entered after python3 runsimulation.py.
    The first will always be "runsimulation.py" itself, so this one is uninteresting
    for our purposes.'''
    args = argv[1:]

    sweep = {}

    if args[0] in optionalargs: 
        sweep['param'] = args[0] # Doesn't need to check if it's a string since it'll imply the correct type by checking if the argument is in optionalargs.
    else:
        plog('Parameter invalid. Valid parameters are')
        print(string(optionalargs))
        return 0                 # Exits the program without proceeding.
    # Parameter is confirmed valid

    try:
        args[1] = float(args[1]) # Make sure the start is a number.
    except ValueError:
        plog('Start invalid.')
        return 0
    sweep['start'] = args[1]
    # Start is confirmed valid
    
    try:
        args[2] = float(args[2])
    except ValueError:
        plog('End invalid.')
        return 0
    sweep['stop'] = args[2]
    # Stop is confirmed valid

    if stop < start:
        plog('Stop is greater than start. Please reverse order.')
        return 0
    
    if args[3] == '-dx':
        try:
            sweep['stepsize'] = float(args[4])
        except ValueError:
            plog('Step size is not a float.')
            return 0
    elif args[3] == '-steps':
        try:
            sweep['steps'] = int(args[4])
        except ValueError:
            plog('Steps is not an integer')
            return 0
    else:
        plog('Specify either -dx or -i followed by stepsize or amount of steps.')
        return 0
    # Steps or step size is confirmed valid. At this point, no more checks need to be done.

    # Generates the vector with all the parameter values to be tested..
    msg = 'Sweeping parameter ' + sweep['param'] + " from " + str(sweep['start']) + " to " + str(sweep['stop'])

    if 'stepsize' in sweep:
        msg +=  " with a step size of " + str(sweep['stepsize']) + "."
        sweepparam = np.arange(sweep['start'], sweep['stop'], sweep['stepsize'])
        plog(msg)
        
    if 'steps' in sweep:
        sweepparam = np.linspace(sweep['start'], sweep['stop'], num = sweep['steps'])
        msg += " in " + str(sweep['steps']) + "steps." 
        plog(msg)

    if len(sweepparam) == 0:
        print('Somehow the program failed right before looping. This is the parameters that failed.')
        print(param)
        
    for i in range(len(sweepparam)):
        arguments = {str(sweep['param']): sweepparam[i]} 
        #runsim(**arguments)
        print(arguments)
    
if __name__ == '__main__':
    main()