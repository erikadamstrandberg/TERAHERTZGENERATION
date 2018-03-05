from sim1d import runsim
from sys import argv
from ast import literal_eval

def main():
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
    etrigger = False
    while args:
        if args[0][0] == '-':                      # Checks the first letter of the first word. If it's a dash, we might have a parameter.
            if args[0][1:] in optionalargs:        # Checks if the parameter is valid.
                param[args[0][1:]] = int(args[1])
            else:
                print(args[0][1:] + ' is not a valid argument, ignoring. ')
                etrigger = True
        args = args[1:]                            # Move to the next parameter.
    if etrigger:
        print('Valid parameters are:\n')
        print(string(optionargs))
    runsim(**param)

if __name__ == '__main__':
    main()
