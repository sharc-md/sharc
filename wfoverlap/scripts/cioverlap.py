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
Short python script to test the main Fortran executable.
If the bra and ket MOs are the same, use --same_mos for a shortcut to the overlap computation.
"""
import numpy

print "cioverlap.py <dets_a> <dets_b> [S_mo]"

class cioverlap:
    def __init__(self):
        self.nmo_a = -1
        self.nmo_b = -1
        self.smo = None

    def read_smo(self, smo_file):
        print "Reading MO-overlap matrix ..."
        sf = open(smo_file, 'r')
        line = sf.next()
        words = line.split()
        
        self.nmo_a = int(words[0])
        self.nmo_b = int(words[1])
        self.smo = numpy.zeros([self.nmo_a, self.nmo_b])
        
        for imo in xrange(self.nmo_a):
            line = sf.next()
            vals = [float(word) for word in line.split()]
            for jmo in xrange(self.nmo_b):
                self.smo[imo, jmo] = vals[jmo]        
        sf.close()
        
    def overlap(self, dets_a, dets_b):
        ovl = numpy.zeros([dets_a.nstate, dets_b.nstate])
        
        nel = dets_a.nel
        assert(nel == dets_b.nel)
        
        dets_a.set_orbinds()
        dets_b.set_orbinds()
        
        for det_a in dets_a.detlist:
            orb_ind_a = det_a[1]
            for det_b in dets_b.detlist:
                orb_ind_b = det_b[1]
                det_mat = numpy.zeros([nel, nel])
                for i_a in xrange(nel):
                    oi_a = orb_ind_a[i_a]
                    for i_b in xrange(nel):
                        oi_b = orb_ind_b[i_b]
                        if oi_a * oi_b < 0:
                            det_mat[i_a, i_b] = 0
                        else:
                            det_mat[i_a, i_b] = self.smo[abs(oi_a)-1, abs(oi_b)-1]
                
                det = numpy.linalg.det(det_mat)            
                for astate in xrange(dets_a.nstate):
                    for bstate in xrange(dets_b.nstate):
                        ovl[astate, bstate] += det_a[astate + 2] * det_b[bstate + 2] * det
                        
        self.prt_ovl(ovl, dets_a, dets_b)
        
    def prt_ovl(self, ovl, dets_a, dets_b):
        print "Overlap matrix <PsiA_i|PsiB_j>"
        self.prt_ovl_core(ovl)
        
        # Renormalize
        inorm_a = 1. / numpy.sqrt(dets_a.sqnorm)
        inorm_b = 1. / numpy.sqrt(dets_b.sqnorm)
        ovl_ren = ovl
        for i in xrange(len(ovl)):
            ovl_ren[i] *= inorm_a[i] * inorm_b
        print " Renormalized overlap matrix <PsiA_i|PsiB_j>"
        self.prt_ovl_core(ovl)
        
    def prt_ovl_core(self, ovl):
        print 11 * " ",
        for i in xrange(len(ovl[0])):
            print "|PsiB%3i>    "%(i+1),
        print
        
        for i,ovl_a in enumerate(ovl):
            print "<PsiA%3i|"%(i+1),
            for ovl_ab in ovl_a:
                print "% .10f"%ovl_ab,
            print
        print

class cioverlap_same(cioverlap):
    def read_smo(self, smo_file):
        print "Skipping MO-overlap input"
    
    def overlap(self, dets_a, dets_b):
        dets_a.sort()
        dets_b.sort()
        
        aind = 0
        bind = 0
        ovl = numpy.zeros([dets_a.nstate, dets_b.nstate])
        
        while(aind<dets_a.ndet and bind<dets_b.ndet):
            det_a = dets_a.detlist[aind]
            det_b = dets_b.detlist[bind]
            astring = det_a[0]
            bstring = det_b[0]            
            
            if astring == bstring:
                for astate in xrange(dets_a.nstate):
                    for bstate in xrange(dets_b.nstate):
                        ovl[astate, bstate] += det_a[astate + 2] * det_b[bstate + 2]
                aind += 1
                bind += 1
            elif astring < bstring:
                aind += 1
            elif bstring < astring:
                bind += 1
                
        self.prt_ovl(ovl, dets_a, dets_b)
        
class dets:
    def __init__(self):
        self.nstate  = -1
        self.norb = -1
        self.ndet    = -1
        self.nel   = -1
        self.detlist = None
        self.sqnorm = None
        
    def read_detfile(self, fname):
        detf = open(fname, 'r')
        line = detf.next()
        words = line.split()
        
        self.nstate = int(words[0])
        self.norb = int(words[1])
        self.ndet   = int(words[2])
        self.detlist= []
        self.sqnorm = numpy.zeros(self.nstate)
        
        for idet in xrange(self.ndet):
            line = detf.next()
            words = line.split()
            # Read: [determinant string, , coeff1, coeff2, ...]
            vals = [float(word) for word in words[1:]]
            self.detlist.append([words[0], []] + vals)
            valarr = numpy.array(vals)
            self.sqnorm += valarr * valarr
            
        detf.close()
        
        self.nel = len(self.orb_inds(self.detlist[0][0]))
        print "nel:", self.nel
        print "sqnorm:", self.sqnorm
        assert(self.nel == len(self.orb_inds(self.detlist[-1][0])))
        
    def orb_inds(self, detstr):
        """
        Extract the orbital indices.
        Positive: alpha
        Negative: beta
        """
        ret_list = []
        
        for i, char in enumerate(detstr):
            if char   == 'd':
                ret_list += [i+1, -i-1]
            elif char == 'a':
                ret_list += [i+1]
            elif char == 'b':
                ret_list += [-i-1]
                
        return ret_list
    
    def sort(self):
        self.detlist.sort()
    
    def set_orbinds(self):
        for idet in xrange(self.ndet):
            self.detlist[idet][1] = self.orb_inds(self.detlist[idet][0])

if __name__ == '__main__':
    import sys, time
    (tc, tt) = (time.clock(), time.time())
    fdets_a = sys.argv[1]
    fdets_b = sys.argv[2]
    
    if len(sys.argv) >= 4:
        S_mo = sys.argv[3]
        cio = cioverlap()    
    else:
        S_mo = None
        cio = cioverlap_same()

    cio.read_smo(S_mo)
    
    da = dets()
    da.read_detfile(fdets_a)
    
    db = dets()
    db.read_detfile(fdets_b)
    print "CPU time: % .1f s, walltime: %.1f s (reading files)"%(time.clock() - tc, time.time() - tt)
    
    cio.overlap(da, db)
    print "CPU time: % .1f s, walltime: %.1f s"%(time.clock() - tc, time.time() - tt)
