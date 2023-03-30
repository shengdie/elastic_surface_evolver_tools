#!/usr/bin/env python3

import sys
import argparse
from mesh import Ring

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--innerradius', '-ir', dest='in_r', type=float, help='inner radius')
    parser.add_argument('--outerradius', '-or', dest='out_r', type=float, help='outer radius')
    parser.add_argument('--in_free', '-if', dest='in_free', action='store_true', help='inner boundary is free')
    parser.add_argument('--gbend', '-gb', dest='gbend', action='store_true', help='add guassian energy')
    parser.add_argument('--egbend', '-egb', dest='egbend', action='store_true', help='add guassian energy except boundary')
    parser.add_argument('--twopoint', '-tp', dest='twopoint', action='store_true', help='inner boundary is free')
    parser.add_argument('--edgesize', '-es', dest='edgesize', type=float, help='Edge size')
    parser.add_argument('--dynamicsize', '-ds', dest='dynamicsize', type=str, default=None,
                        help='dynamicsize, ie, lambda p: abs(ic.dist(p))/5 + [edgesize]')
    parser.add_argument('--kargs', '-kw', dest='kargs', type=str, default='{}',
                        help='kargs for generate')                    
    parser.add_argument('filename', nargs='?', metavar='filename', default=None, help='"File_Name_Without_ext"')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    pargs = parser.parse_args()
    
    ring = Ring(pargs.in_r, pargs.out_r, pargs.edgesize, oname=pargs.filename, dynamicsize=pargs.dynamicsize, 
                in_free=pargs.in_free, twopoint=pargs.twopoint, gbend=pargs.gbend, egbend=pargs.egbend, **eval(pargs.kargs))
    # load_cmd = False
    # if pargs.in_free:
    ring.write_input(load_cmd=pargs.in_free, twopoint=pargs.twopoint)
