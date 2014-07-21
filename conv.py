#! /usr/bin/python

import re
import argparse
from utils import *

__author__="ph4r05"
__date__ ="$Jul 3, 2014 3:58:14 PM$"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Keccak1600 GF(2) conversion script.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-m','--multiply',  help='Perform multiplication in the equation. 0=no 1=yes', required=False, default=0, type=int)
    parser.add_argument('-f','--format',    help='Format of the output. 0=no change, 1=coordinates, 2=variable indexes', required=False, default=0, type=int)
    parser.add_argument('-v','--verbose',   help='Writes output to the standard output', required=False, default=0, type=int)
    parser.add_argument('file')
    args = parser.parse_args()
    
    print " [-] Processing file: %s" % (args.file)
    with open(args.file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = str(line).strip()
            if not line: 
                print line
                continue
            if line[0:2]=="//": 
                print line
                continue
            
            # Strip for left hand side and right hand side
            (lhs, rhs) = [str(x).strip() for x in line.split('=', 2)]
            
            # RHS follows the pattern: XX + (YY) * (ZZ)
            match = re.match(r'^(.+)\((.+)\)\*\((.+)\)\s*(.*)$', rhs)
            if not match: continue
            xx = match.group(1)
            yy = match.group(2)
            zz = match.group(3)
            cs = match.group(4)
            
            xxl=addition2termList(xx)
            yyl=addition2termList(yy)
            zzl=addition2termList(zz)
            cs =addition2termList(cs) 
            totalList = []
            
            if args.multiply>0:
                mList = multiplyLists(yyl, zzl)
                totalList = gf2ize(cs + xxl + mList)
            
            if args.format==3:
                lst = sorted(totalList, key=cmp_to_key(bitnameCmp))
                orderArr = [[],[],[]]
                for e in lst:
                    curOrder = getTermOrder(e)
                    curVar = bitName2variableStr(e, "", 4)
                    orderArr[curOrder].append(curVar)
                
                quadsA = sorted(yyl, key=cmp_to_key(bitnameCmp))
                quadsA = [bitName2variableStr(x, "", 4) for x in quadsA if x!= '1']
                
                quadsB = sorted(zzl, key=cmp_to_key(bitnameCmp))
                quadsB = [bitName2variableStr(x, "", 4) for x in quadsB if x!= '1']
                
                res  = bitName2variableStr(lhs, "", 4)+"="
                res +=     (",".join(quadsA))
                res += ";"+(",".join(quadsB))+";"
                res += ";".join([",".join(x) for x in orderArr])
                
                print res
                continue
            
            if args.multiply>0:
                print listDump([lhs], args.format, True, True) + " = " + listDump(totalList, args.format, True, False)
                continue
            
            print listDump([lhs], args.format, True, True) + " = " + listDump(cs + xxl, args.format, True, False) \
                + " + (" + listDump(yyl, args.format, True, False) + ")*(" + listDump(zzl, args.format, True, False) + ")"
        pass
    pass
pass


