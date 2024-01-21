# scsc.py - command line interface


import os
import sys
from sys import stdout, stderr

from .app import APP, VERSION
from .pileup import pileup_main
from .simulation import simu_main


def usage(fp = stderr):
    s =  "\n" 
    s += "Version: %s\n" % (VERSION, )
    s += "Usage: %s [commands|options]\n" % (APP, )  
    s += "\n" 
    s += "Commands:\n"
    #s += "  main              Run all steps.\n"
    s += "  plp               Pileup allele-specific UMIs.\n"
    s += "  simu              CNV Simulation.\n"
    s += "\n" 
    s += "Options:\n"
    s += "  -h, --help        Print this message and exit.\n"
    s += "  -V, --version     Print version and exit.\n"
    s += "\n"

    fp.write(s)


def main():
    func = "main"

    if len(sys.argv) <= 1:
        usage(stderr)
        sys.exit(1)

    cmd = sys.argv[1]
    #if cmd == "main": usage(stderr); sys.exit(1)
    if cmd == "plp": pileup_main(sys.argv)
    elif cmd == "simu": simu_main(sys.argv)
    elif cmd in ("-V", "--version"): stderr.write("%s\n" % VERSION); sys.exit(1)
    elif op in ("-h", "--help"): usage(); sys.exit(1)
    else:
        stderr.write("[E::%s] invalid command: '%s'.\n" % (func, cmd))
        return(-1) 


if __name__ == "__main__":
    main()

