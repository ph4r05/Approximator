#! /usr/bin/python

import re
import argparse

__author__="ph4r05"
__date__ ="$Jul 3, 2014 3:58:14 PM$"

tranTableRound={'A':'0', 'B':'1'}
tranTableY={'b':0, 'g':1, 'k':2, 'm':3, 's':4}
tranTableX={'a':0, 'e':1, 'i':2, 'o':3, 'u':4}

def cmp_to_key(mycmp):
    'Convert a cmp= function into a key= function'
    class K(object):
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0
        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K

def bitName2bitCoords(bitname):
    '''Translates bitname to bit coordinates, Aba13 -> A0013'''
    bitname=str(bitname).strip()
    if bitname=="1": return "1"
    return bitname[0] + str(tranTableY[bitname[1]]) + str(tranTableX[bitname[2]]) + bitname[3:]

def bitName2variable(bitname):
    '''Translates bitname to variable index.'''
    bitname=str(bitname).strip()
    if bitname=="1": return -1
    x = tranTableX[bitname[2]]
    y = tranTableY[bitname[1]]
    z = int(bitname[3:])
    return (64*(5*y+x)+z)

def addition2termList(input):
    '''Converts addition of the variables to the list representation'''
    return [str(x).strip() for x in input.split('+') if (str(x).strip())]

def bitnameCmp(x, y):
    '''Lexicographic comparison of the terms'''
    if x==y:   return 0
    if x=="1": return -1
    if y=="1": return 1
    
    quadX = "*" in x
    quadY = "*" in y
    if quadX and not quadY: return 1
    if not quadX and quadY: return -1
    if quadX and quadY:
        qx = x.split('*', 2)
        qy = y.split('*', 2)
        cmpFst = bitnameCmp(qx[0], qy[0])
        if cmpFst!=0: 
            return cmpFst
        else: 
            return bitnameCmp(qx[1], qy[1])
    
    xx = bitName2variable(x)    
    yy = bitName2variable(y)
    return xx-yy

def gf2ize(lst):
    '''Removes even number of elements from the addition list - GF2 addition'''
    res = set()
    for elem in lst:
        if elem in res:
            res.remove(elem)
        else:
            res.add(elem)
    return list(res)
    
def multiplyLists(a,b):
    '''Multiplies lists of the variables in GF(2)'''
    linList = []
    res = set()
    for ea in a:
        for eb in b:
            tmpRes=""
            if ea=="1":
                linList.append(eb)
                continue
            elif eb=="1":
                linList.append(ea)
                continue
            elif bitnameCmp(ea, eb) <= 0:
                tmpRes=ea+"*"+eb
            else:
                tmpRes=eb+"*"+ea
            
            # In GF(2), add only if not in the set. 
            # If present, remove from set.
            if tmpRes in res:
                res.remove(tmpRes)
            else:
                res.add(tmpRes)
                
    resList = list(res)
    resList.sort()
    linList = gf2ize(linList)
    linList.sort()
    return linList + resList

def listDump(lst, fmt, doSort=True):
    '''Dumps variable list in given format as string'''
    if doSort:
        lst = sorted(lst, key=cmp_to_key(bitnameCmp))
    res=""
    if fmt==0:
        res = " + ".join([x for x in lst])
    elif fmt==1:
        res = " + ".join([bitName2bitCoords(x) for x in lst])
    elif fmt==2:
        bitIdx = [bitName2variable(x) for x in lst]
        res = " + ".join("x_"+str(x) if x>=0 else "1" for x in bitIdx)
    return res

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
            
            if args.multiply>0:
                mList = multiplyLists(yyl, zzl)
                totalList = gf2ize(cs + xxl + mList)
                
                print listDump([lhs], args.format) + " = " + listDump(totalList, args.format)
                continue
            
            
            print listDump([lhs], args.format) + " = " + listDump(cs + xxl, args.format) \
                + " + (" + listDump(yyl, args.format) + ")*(" + listDump(zzl, args.format) + ")"
        pass
    pass
pass


