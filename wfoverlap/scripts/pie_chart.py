#!/usr/bin/env python2

#******************************************
#
#    SHARC Program Suite
#
#    Copyright (c) 2018 University of Vienna
#
#    This file is part of SHARC.
#
#    SHARC is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SHARC is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    inside the SHARC manual.  If not, see <http://www.gnu.org/licenses/>.
#
#******************************************

#!/usr/bin/env python2
"""
Plot pie charts.
"""

import numpy
try:
    import matplotlib
    matplotlib.use('Agg')
    import pylab
except:
    print "pylab/matplotlib not installed - plotting not possible"
    raise

class pie_chart:
    def __init__(self):
        self.Smat = None
        self.labels = None
#        self.transpose = False
        self.Astates = None
        self.Bstates = None
        self.oformat = 'png'
        self.fsize=25
        self.colors = None
        self.ststr = 'Overlap matrix <P'
    
    def setup(self):
        pylab.ax = pylab.axes([0.1, 0.1, 0.8, 0.8])
        matplotlib.rc('font', size=self.fsize)
        
    def read(self, fname):
        ifile = open(fname, 'r')
        
        while(True):
            line = ifile.next()
            if self.ststr in line: break
            
        line = ifile.next()
        Slist = []
        while(True):
            try:
                words = ifile.next().split()
            except StopIteration:
                print "End of file reached"
                break
            if len(words) == 0: break
            
            Slist.append(words[2:])
            
        self.Smat = numpy.array(Slist, float)
        
        ifile.close()
        
        if self.Astates == None:
            self.Astates = range(len(self.Smat))
        else:
            self.Astates = [st-1 for st in self.Astates]
        if self.Bstates == None:
            self.Bstates = range(len(self.Smat[0]))
        else:
            self.Bstates = [st-1 for st in self.Bstates]
    
    def pie(self):
        self.pie_core(self.Smat, self.Astates, self.Bstates, 'PsiA')
        self.pie_core(self.Smat.transpose(), self.Bstates, self.Astates, 'PsiB')
        
    def pie_core(self, S, refstates, pstates, state_pre=''):
        for i in refstates:
            pylab.figure(figsize=(6,6))
            
            vals = [S[i, j]**2. for j in pstates]

            #explode = len(self.pstates) * [0.05]
            #explode[self.refstates.index(i)] = 0.0
            
            pylab.pie(vals, labels=self.labels, autopct=self.mk_autpct, explode=None, colors=self.colors)
            #pylab.pie(vals, labels=self.labels, autopct='%1.1f%%')
            
            #pylab.legend()
            #pylab.title(r'$x^x$')
            pylab.savefig("S_%s_%i.%s"%(state_pre,i+1,self.oformat))
    
    def mk_autpct(self, val, minprt=10., fmt='%1.0f%%'):
        if val >= minprt:
            return fmt%val
        else:
            return ''
            
if __name__ == '__main__':
    import sys
  
    pie = pie_chart()
    
    fname = sys.argv[1]
    
    args = sys.argv[2:]
    while(len(args)>0):
        arg = args.pop(0)
        if arg == '-l' or arg == '--labels':
            #pie.labels = [r'$%s$'%lab for lab in args.pop(0).split()] # does not work like this
            pie.labels = args.pop(0).split()
        elif arg == '-A' or arg == '--Astates':
            pie.Astates = [int(word) for word in args.pop(0).split()]
        elif arg == '-B' or arg == '--Bstates':
            pie.Bstates = [int(word) for word in args.pop(0).split()]
        elif arg == '-f' or arg == '--fsize':
            pie.fsize = int(args.pop(0))
        elif arg == '-o' or arg == '--oformat':
            pie.oformat = args.pop(0)
        elif arg == '-c' or arg == '--colors':
            pie.colors = args.pop(0).split()
        elif arg == '-r' or arg == '--renormalized':
            pie.ststr = 'Renormalized overlap matrix'
        elif arg == '-l' or arg == '--lowdin':
            pie.ststr = 'Orthonormalized overlap matrix'
        else:
            print "Unsupported option: '%s'"%arg
            exit(1)
    
    pie.setup()
    pie.read(fname)
    pie.pie()