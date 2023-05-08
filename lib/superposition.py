"""
version 1.0
author: Felix Plasser
description: Superposition of two sets of weighted coordinates using a quaternion fit, C.F.F. Karney, J Mol Graphics and Modelling 25 (2007) p.595.
"""

import sys
try:
    import numpy
except:
    print('numpy not installed')
    sys.exit()

class superposition:
    """
    Class for superimposing two sets of weighted coordinates.
    """    
    def superimpose(self, ref_points, mv_points, weights):
        """
        The superposition is carried out.
        <ref_points> and <mv_points> are numpy 3xN matrices. <Weights> is a numpy N vector that contains the statistical weighing (typically the atomic masses.
        """
        self.ref_points = ref_points
        self.mv_points = mv_points
        self.weights = weights
        self.W = numpy.add.reduce(weights)
        
        self.ref_av = self.average(self.ref_points) # weighted average
        self.mv_av = self.average(self.mv_points)
        
        ref_c = self.ref_points - self.ref_av # centered coordinates
        mv_c = self.mv_points - self.mv_av
        
        Ak_T_list = [] # list with all the Ak        
        for i in range(len(ref_c)):
            Ak = self.A(ref_c[i]+mv_c[i],ref_c[i]-mv_c[i])
            Ak_T_list += [numpy.dot(Ak.transpose(),Ak)] # list of Ak's multiplied with transposed matrix
        Ak_T_array = numpy.array(Ak_T_list, float)
        
        #print 'Ak_T_array'
        #print Ak_T_array
        
        B = self.average(Ak_T_array)
        #print 'B', B
        
        evals, temp_evecs = numpy.linalg.eig(B)
        evecs = temp_evecs.transpose() # eigenvectors are initially written into columns
        
        #print evals, evecs
    
        i = numpy.argmin(evals)
        self.msd = evals[i]
        self.rot_quat = evecs[i] / numpy.sqrt(numpy.dot(evecs[i],evecs[i])) # vector is typically already normalised
        
    def average(self, array):
        """
        Returns the weighted average of <array>.
        """ 
        #print array
        #print array.shape[1:]
        res_array = numpy.zeros(array.shape[1:], float) # the resulting array has the shape of the inner part of <array>
        
        for i,sub_array in enumerate(array):
            res_array += sub_array * self.weights[i]
        
        return 1/self.W*res_array
    
    def A(self, a, b):
        """
        Create the matrix a as shown in Karney (2007).
        <a> and <b> are 3x1 numpy arrays.
        """
        A_list =  [[    0,-b[0],-b[1],-b[2]], \
                   [ b[0],    0,-a[2], a[1]], \
                   [ b[1], a[2],    0,-a[0]], \
                   [ b[2],-a[1], a[0],    0]]
        
        return numpy.array(A_list)
    
    # setting of vales
    def set_rotation(self, theta, vec):
        vec = vec / numpy.sqrt(numpy.dot(vec,vec))
        self.rot_quat = numpy.array([numpy.cos(theta/2), numpy.sin(theta/2)*vec[0],numpy.sin(theta/2)*vec[1],numpy.sin(theta/2)*vec[2]])
    
    # returning of values
    
    def ret_rmsd(self):
        """
        Returns the root mean squared deviation, i.e. the weighted norm of the difference vector. The result is in Angstrom.
        """
        if self.msd < 0.0:
            rmsd = 0.0
        else:
            rmsd = numpy.sqrt(self.msd)
        
        return rmsd
        
    def ret_ref_av(self):
        """
        Returns the average (weighted center) of the reference structure.
        """
        return self.ref_av
    
    def ret_mv_av(self):
        """
        Returns the average (weighted center) of the moving structure.
        """
        return self.mv_av
    
    def ret_rotation_matrix(self):
        """
        Returns the rotation matrix.
        """
        q = self.rot_quat
        
        Rq_list = [[  1-2*q[2]**2-2*q[3]**2, 2*q[1]*q[2]-2*q[0]*q[3], 2*q[1]*q[3]+2*q[0]*q[2]], \
                   [2*q[2]*q[1]+2*q[0]*q[3],   1-2*q[3]**2-2*q[1]**2, 2*q[2]*q[3]-2*q[0]*q[1]], \
                   [2*q[3]*q[1]-2*q[0]*q[2], 2*q[3]*q[2]+2*q[0]*q[1],   1-2*q[1]**2-2*q[2]**2]]
        
        return numpy.array(Rq_list)
    
    def ret_rotation_angle(self):
        """
        Returns the value of the rotation angle.
        """
        return numpy.arccos(self.rot_quat[0])*2
    
    def ret_rotation_axis(self):
        """
        Returns the rotation axis.
        """
        return self.rot_quat[1:4] / numpy.sin(self.ret_rotation_angle() / 2)