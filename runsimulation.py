from sim1d import runsim
from sys import argv

def main():
    args = {}
    if argv[1:]:
        args = argv[1:]
        print(args)
    
    runsim(dt = float(1), dz = float(1))

if __name__ == '__main__':
    main()
