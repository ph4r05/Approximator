#! /usr/bin/python

import sys
import argparse
import itertools
from utils import *

__author__="ph4r05"
__date__ ="$Jul 3, 2014 3:58:14 PM$"

def pruneHighTerms(tlist):
    '''Remove high terms not reachable by classical attack, 256 and higher'''
    (l1,l2) = tlist
    l1 = [x for x in l1 if x<256]
    l2 = [x for x in l2 if x<256]
    return [l1, l2]

def pruneKeyTerms(tlist):
    '''Remove terms with key variables'''
    (l1,l2) = tlist
    l1 = [x for x in l1 if x<128]
    l2 = [x for x in l2 if x<128]
    return [l1, l2]

def prunePlainTerms(tlist):
    '''Remove terms with plaintext variables'''
    (l1,l2) = tlist
    l1 = [x for x in l1 if (x>=128 and x<256)]
    l2 = [x for x in l2 if (x>=128 and x<256)]
    return [l1, l2]

fidDepMap = []   # [term] -> set of functions that have this term in quadratic term
fidQuadMap = []
fidTermMap = [[],[],[],[],[]]  # quadratic term map

def expandTermList(fidx, level):
    '''Expands function to the list of the usefull terms for cube attack'''
    prod = itertools.product(fidQuadMap[fidx][0], fidQuadMap[fidx][1])
    tmpSet = set([])
    for term in prod:
        tA = fidTermMap[level-1][term[0]]
        tB = fidTermMap[level-1][term[1]]
        #print "ta=",tA,"; tb=",tB
        prod2 = itertools.product(tA, tB)
        for tmpTerm in prod2:
            tmpTup = () + tmpTerm[0] + tmpTerm[1]
            tmpTerm = list(set(tmpTup))
            tmpTerm.sort()
            newTerm = tuple(tmpTerm)
            
            # Duplicity reduction.
            #if len(tmpTup) != len(newTerm):
            #    continue
                
            # Allow maximally one key variable.
            #keyVars=len(filter(lambda y: y >= 128 and y < 256, newTerm))
            #if keyVars > 1:
            #    continue
            
            tmpSet.add(newTerm)
    fidTermMap[level][fidx] += list(tmpSet)
    fidTermMap[level][fidx].sort()
    #print "x"
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Keccak1600 GF(2) bottom-up script.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-m','--multiply',  help='Perform multiplication in the equation. 0=no 1=yes', required=False, default=0, type=int)
    parser.add_argument('-f','--format',    help='Format of the output. 0=no change, 1=coordinates, 2=variable indexes', required=False, default=0, type=int)
    parser.add_argument('-v','--verbose',   help='Writes output to the standard output', required=False, default=0, type=int)
    parser.add_argument('file')
    args = parser.parse_args()
    
    #print " [-] Processing file: %s" % (args.file)
    for i in range(0, 1600):
        fidDepMap.append(set([]))
        fidQuadMap.append([])
        fidTermMap[0].append([])
        fidTermMap[1].append([])
        fidTermMap[2].append([])
    
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
            
            # Strip for left hand side and right hand side.
            (lhs, rhs) = [str(x).strip() for x in line.split('=', 2)]
            (q1,q2,const,linear,quads) = rhs.split(";")
            q1 = [int(str(x).strip()) for x in q1.split(",")]
            q2 = [int(str(x).strip()) for x in q2.split(",")]
            qTerms = set(q1+q2)
            
            # Build dependencies of the functions to the next round.
            fid = int(lhs)
            for curQuad in (qTerms):
                fidDepMap[curQuad].add(fid)
                
            # Build quad map list
            curQuadTermList = [q1, q2]
            fidQuadMap[fid] = curQuadTermList
            
            q1l = len(filter(lambda y: y >= 128 and y < 256, q1))
            q2l = len(filter(lambda y: y >= 128 and y < 256, q2))
            
            fidTermMap[0][fid] = []
            pruneList = pruneHighTerms(fidQuadMap[fid])
            
            # Init all term posibilities expandable from this function
            prod = itertools.product(pruneList[0], pruneList[1])
            for term in prod:
                fidTermMap[0][fid].append(term)
        pass
    pass

    # Preprocessing - list gen
    listDepMap = []
    setQuadMap = []
    for i in range(0,1600):
        lpi = list(fidDepMap[i])
        lpi.sort()
        listDepMap.append(lpi)
        
        # Set-ize quad lists
        setQuadMap.append([set(fidQuadMap[i][0]), set(fidQuadMap[i][1])])
    
    # Find such pi and pj such that the quadratic term is not able to construct
    # for any functions.
    padBytes=[]#range(128,256)#[]#[256, 1023]
    
    piRels = [set([]) for x in range(0,1600)]
    for pi in (range(0,256)+padBytes):        
        for pj in (range(0,256)+padBytes):
            
            # Find set of the functions such that pi and pj are both in the
            # function contained in some quadratic terms.
            fIntersect = fidDepMap[pi] & fidDepMap[pj]
            
            # If the set is not null, verify if it is really possible to construct
            # quadratic term pi*pj in that function.
            if len(fIntersect)>0:
                usablePj=True
                for fidx in fIntersect:
                    if not usablePj: break
                    quads = setQuadMap[fidx]
                    
                    # Test if product pi*pj is not possible.
                    if (pi in quads[0] and pj in quads[1]) or (pi in quads[1] and pj in quads[0]):
                        usablePj=False
                        break
                if not usablePj:
                    continue
                    
            piRels[pi].add(pj)
        pass
    pass
    
    def dotDraw(rng):
        '''Draws DOT diagram for the relations'''
        graph = "graph finite_state_machine {\n"
        #graph += "    layout=\"circo\";\n"
        for i in rng:
            color = '00ff005f' if i < 128 else '0000ff5f'
            graph += '    node [shape=point, fillcolor="#'+color+'", label=""]N_'+str(i)+";\n" 
        
        for i in rng:
            for c in rng:
                if c < i: continue
                if c in piRels[i]: continue
                
                graph += 'N_' + str(i) + ' -- N_'+ str(c) + ";\n" #+'[penwidth="2.0", arrowsize="2.5"]'
            
        graph += "label=\"inRelation\"\n"
        graph += "fontsize=32;}\n"
        return graph
        pass
    
    
    def recFind(level, interList, curList):
        #Fuck yeah, end of the recursion
        if len(interList)==0 or level >= 8+len(padBytes):
            curListL = list(curList)
            curListL.sort()
            iListL = list(interList)
            iListL.sort()
            compl = list(set(range(0,256)) - set(iListL))
            compl.sort()
            
            # Convert to decimal number representation;
            dec=0
            for i in curListL:
                if i >= 128 or (i in padBytes): continue
                dec += 2**i
                
            keys=len(filter(lambda y: y>=128 and y not in curList, compl))
            if keys<=60:
                print "!!!level:", level, "; curList:", curListL , "; dec:", dec, "; interList: ", iListL, "; complement: ", compl, " keys left:", keys
            
            #for i in filter(lambda y: y<128 and y not in curList, compl):
            #    print "  ",i,": ", filter(lambda y: y in curListL, piRels[i])
            #    pass
            return
        
        # Choose pi
        reducedCurList = filter(lambda y: y<128, curList)
        mm = max(reducedCurList) if len(reducedCurList)>0 else 0
        for i in interList:
            # pi != pj
            if i in curList or i >= 128 or i < mm: continue
            
            niList = ((piRels[i] & interList)-curList)-set([i])
            #print "level:", level, "i:", i, "; curList:", curList, "; niList: ", niList
            
            recFind(level+1, niList, curList | set([i]))
        pass
    pass

    print dotDraw(range(0,255))
    sys.exit(3)
    
    # Find pairwise disjoint.
    if len(padBytes)>0:
        padIntersection = set(range(0,1600))
        for i in padBytes:
            padIntersection = padIntersection & piRels[i]
            
        if len(padIntersection)==0:
            print "No such solution exists for padding as defined.", padBytes
            sys.exit(2)
            
        curStart = 0
        recFind(len(padBytes), padIntersection, set(padBytes))
    else:
        curStart = 0
        recFind(1, piRels[curStart], set([curStart]))        
            
    sys.exit(2)
    
    
    
    
    
    
    
    
    piList = []
    curfidx = 0
    curFList = fidQuadMap[curfidx][0] + fidQuadMap[curfidx][1]
    for pi in range(0,128):
        piList.append([])
        lpi = listDepMap[pi]
        
        for pj in range(pi+1,256):
            lpj = listDepMap[pj]
            
            # If there is an intersection bewteen functions containing
            # pi and pj in quadratic terms, they cannot be used since
            # pi and pj can be multiplied in the lowest level of the 
            # multiplication tree
            fIntersect = fidDepMap[pi] & fidDepMap[pj]
            if len(fIntersect)!=0:
                usablePj=True
                for fidx in fIntersect:
                    if not usablePj: break
                    quads = setQuadMap[fidx]
                    
                    # Test if product pi*pj is not possible.
                    if (pi in quads[0] and pj in quads[1]) or (pi in quads[1] and pj in quads[0]):
                        usablePj=False
                        break
                if not usablePj:
                    continue
            
            # For each function
            fidxCtr=0
            for fidx in curFList:
                quads = setQuadMap[fidx]
                totalQuads = len(quads[0]) * len(quads[1])
                
                ffail=False # if there is intersection -> ffail=True
                for a in lpi:
                    if ffail: break
                    for b in lpj:
                        if (a in quads[0] and b in quads[1]) or (a in quads[1] and b in quads[0]):
                            ffail=True
                            break
                        
                if not ffail:
                    #print "#", fidxCtr, " No intersection here, pi=",pi, "; pj=",pj,"; fidx=", fidx
                    fidxCtr+=1
            
            if (fidxCtr==len(curFList)):
                print "pi=",pi, "; pj=",pj,"; #nema prienik=", fidxCtr, "!!!"
                piList[pi].append(pj)
    
    for (i,pi) in enumerate(piList):
        print "pi=",i,"; list=", pi
    
    
    #l38 = list(fidDepMap[38])
    #l38.sort()
    #l39 = list(fidDepMap[39])
    #l39.sort()

    #print "DepMap 38: ", l38
    #print "DepMap 39: ", l39
    #print "Priamo 38: ", fidQuadMap[38]
    #print "Priamo 58: ", fidQuadMap[58]
    #print "Priamo 123: ", fidQuadMap[123]
    #print ""
    #print "Priamo 39: ", fidQuadMap[39]
    
    #for i in l38:
    #    print "cur f=",i,fidQuadMap[i]
    #    if i in fidQuadMap[38][0] or i in fidQuadMap[38][1]:
    #        print " HA!, i=", i
    #print "\n39"
    #for i in l39:
    #    print "cur f=",i,fidQuadMap[i]
    #    if i in fidQuadMap[39][0] or i in fidQuadMap[39][1]:
    #        print " HA!, i=", i
    
    #for val in range(0, 1600):
    #    valCn=0
    #    for x in fidDepMap[val]:
    #        #print "x=",x
    #        (q1, q2) = fidQuadMap[x]
    #        q1l = len(filter(lambda y: y >= 128 and y < 256, q1))
    #        q2l = len(filter(lambda y: y >= 128 and y < 256, q2))
    #        if (val in q1 and q2l==0) or (val in q2 and q1l==0):
    #            valCn+=1
    #            #print "Mam to vole!", q1l, " ", q2l
    #            #print fid, " q1=", q1, " q2=", q2
    #    if valCn==len(fidDepMap[val]):
    #        print "WIN!, val=",val 
    #    else:
    #        print "valCn=", valCn, "; len=", len(fidDepMap[val])
                
    sys.exit(0)
    
    #pruneList = pruneHighTerms(fidQuadMap[0])
    #print "f00=",pruneList
    #print "f20=",pruneHighTerms(fidQuadMap[20])
    
    # Expand terms up to the second round.
    # Start at f0.
    for fidx in range(0,1600):
        expandTermList(fidx, 1)
        print fidx
    
    for k in range(0,1600):
        for x in fidTermMap[1][k]:
            if 0 in x and 49 in x:
                print "k=",k," zasa: ", x
    
    
     
    #print "\n".join(fidTermMap[1][0])
    
    #for fidx in range(0,1600):
    #    expandTermList(fidx, 2)
    #    print fidx
    #fidTermMap[2][0]
    
    
    olist = {0:[],1:[]}
    for (clid, clist) in enumerate(fidQuadMap[0]):
        # descent at lower level, level zero, thus prunning can be done.
        for (cmid, cmono) in enumerate(clist):
            cmonoPruneTerms = pruneHighTerms(fidQuadMap[cmono])
            olist[clid].append(cmonoPruneTerms)
    #print olist
    #print fidDepMap
    #print fidQuadMap
    
    
    
pass


