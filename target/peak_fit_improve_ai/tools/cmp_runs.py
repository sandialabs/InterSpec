#!/usr/bin/env python3
"""Presence on a FIXED truth denominator: matched counts per class, plus the peaks that flipped."""
import sys, csv, collections, os
def load(d):
    out=collections.defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(d,'per_peak.tsv')),delimiter='\t'):
        if r['set']!='truth': continue
        out[r['id']][round(float(r['energy']),1)]=r
    return out
def cls(z):
    z=float(z); return 'strong' if z>=8 else ('moderate' if z>=3 else 'weak')
A,B=load(sys.argv[1]),load(sys.argv[2])
tab=collections.defaultdict(lambda:[0,0,0]); flips=[]
for pid in sorted(set(A)&set(B)):
    for e in set(A[pid])&set(B[pid]):
        ra,rb=A[pid][e],B[pid][e]
        if float(ra['energy'])<20: continue
        if ra['verdict'].startswith('dontcare') or rb['verdict'].startswith('dontcare'): continue
        c=cls(ra['z_det']); fa=ra['verdict'].startswith('matched'); fb=rb['verdict'].startswith('matched')
        t=tab[c]; t[0]+=1; t[1]+=fa; t[2]+=fb
        if fa!=fb: flips.append((pid,float(ra['energy']),float(ra['z_det']),c,ra['verdict'],rb['verdict']))
print('%-9s %6s %16s %16s'%('class','n','found A','found B'))
for c in ('strong','moderate','weak'):
    n,fa,fb=tab[c]
    if n: print('%-9s %6d %8d (%5.1f%%) %8d (%5.1f%%)'%(c,n,fa,100.*fa/n,fb,100.*fb/n))
lost=[f for f in flips if f[4].startswith('matched')]; gained=[f for f in flips if not f[4].startswith('matched')]
print('\ngained in B: %d   lost in B: %d'%(len(gained),len(lost)))
for title,rows in (('LOST in B',lost),('GAINED in B',gained)):
    rows=sorted(rows,key=lambda r:-r[2])[:30]
    if rows:
        print('\n%s:'%title)
        for pid,e,z,c,va,vb in rows: print('  %-22s %8.1f keV z=%7.1f %-9s %s -> %s'%(pid,e,z,c,va,vb))
