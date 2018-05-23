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
author: Felix Plasser, Andrew Atkins
version: 2.0
descr: Script for doing essential dynamics analysis.
"""

# runs on hawk4,5,6,11,14

import os, sys, datetime, readline, re, shutil
sys.path.insert(0, os.environ['SHARC']+'/../lib')
try:
    import numpy
except ImportError:
    print 'numpy package not installed'
    sys.exit()
try:
    import file_handler, vib_molden, traj_manip, struc_linalg
except ImportError:
    print 'file_handler, vib_molden, traj_manip or struc_linalg not found. They should be part of this package. Check the installation and if $SHARC/../lib is part of the PYTHONPATH environment variable.'
    sys.exit()

version='2.0'
versiondate=datetime.date(2018,2,1)

# ======================================================================= #
def centerstring(string,n,pad=' '):
  l=len(string)
  if l>=n:
    return string
  else:
    return  pad*((n-l+1)/2)+string+pad*((n-l)/2)

# ======================================================================= #
def displaywelcome():
  print 'Script for performing essential dynamics analysis started ...\n'
  string='\n'
  string+='  '+'='*80+'\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Essential dynamics analysis for SHARC dynamics',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Author: Felix Plasser, Andrew Atkins',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Version:'+version,80)+'||\n'
  string+='||'+centerstring(versiondate.strftime("%d.%m.%y"),80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  string+='''
This script reads output.xyz files and calculates the essential dynamics 
(i.e., Shows you which are the most important motions).
  '''
  print string

# ======================================================================= #
def readfile(filename):
  try:
    f=open(filename)
    out=f.readlines()
    f.close()
  except IOError:
    print 'File %s does not exist!' % (filename)
    sys.exit(12)
  return out

# ======================================================================= #
def writefile(filename,content):
  # content can be either a string or a list of strings
  try:
    f=open(filename,'w')
    if isinstance(content,list):
      for line in content:
        f.write(line)
    elif isinstance(content,str):
      f.write(content)
    else:
      print 'Content %s cannot be written to file!' % (content)
    f.close()
  except IOError:
    print 'Could not write to file %s!' % (filename)
    sys.exit(13)

# ======================================================================= #
def open_keystrokes():
  global KEYSTROKES
  KEYSTROKES=open('KEYSTROKES.tmp','w')

# ======================================================================= #
def close_keystrokes():
  KEYSTROKES.close()
  shutil.move('KEYSTROKES.tmp','KEYSTROKES.trajana_essdyn')

# ======================================================================= #
def question(question,typefunc,default=None,autocomplete=True,ranges=False):
  if typefunc==int or typefunc==float:
    if not default==None and not isinstance(default,list):
      print 'Default to int or float question must be list!'
      quit(1)
  if typefunc==str and autocomplete:
    readline.set_completer_delims(' \t\n;')
    readline.parse_and_bind("tab: complete")    # activate autocomplete
  else:
    readline.parse_and_bind("tab: ")            # deactivate autocomplete

  while True:
    s=question
    if default!=None:
      if typefunc==bool or typefunc==str:
        s+= ' [%s]' % (str(default))
      elif typefunc==int or typefunc==float:
        s+= ' ['
        for i in default:
          s+=str(i)+' '
        s=s[:-1]+']'
    if typefunc==str and autocomplete:
      s+=' (autocomplete enabled)'
    if typefunc==int and ranges:
      s+=' (range comprehension enabled)'
    s+=' '

    line=raw_input(s)
    line=re.sub('#.*$','',line).strip()
    if not typefunc==str:
      line=line.lower()

    if line=='' or line=='\n':
      if default!=None:
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return default
      else:
        continue

    if typefunc==bool:
      posresponse=['y','yes','true', 't', 'ja',  'si','yea','yeah','aye','sure','definitely']
      negresponse=['n','no', 'false', 'f', 'nein', 'nope']
      if line in posresponse:
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return True
      elif line in negresponse:
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return False
      else:
        print 'I didn''t understand you.'
        continue

    if typefunc==str:
      KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
      return line

    if typefunc==float:
      # float will be returned as a list
      f=line.split()
      try:
        for i in range(len(f)):
          f[i]=typefunc(f[i])
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return f
      except ValueError:
        print 'Please enter floats!'
        continue

    if typefunc==int:
      # int will be returned as a list
      f=line.split()
      out=[]
      try:
        for i in f:
          if ranges and '~' in i:
            q=i.split('~')
            for j in range(int(q[0]),int(q[1])+1):
              out.append(j)
          else:
            out.append(int(i))
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return out
      except ValueError:
        if ranges:
          print 'Please enter integers or ranges of integers (e.g. "-3~-1  2  5~7")!'
        else:
          print 'Please enter integers!'
        continue

### input
# ref_struc_file # file with a reference structure
# ref_struc_type # type of that file
# first_traj = 1
# last_traj = 50
# num_steps = 601 # maximum number of timesteps
# ana_ints = [[0,201],[0,601]] # time intervals to be analysed
# descr '' # descr added to the filenames
###

 ##  read in variable assignments from ess_dyn.inp, this is always done when the program is called
 # first default definitions
#dt = .5
# mass_wt_pw = 0.# mass weighting does not make sense as the variance does not depend on the mass
    # the variance, proportional to the amplitude of a normal mode, does not depend on mass with a given energy
#descr = ''
 
 # read from file
#execfile('ess_dyn.inp')


# ======================================================================= #
# ======================================================================= #
# ======================================================================= #
def get_general():

  INFOS={}

  print centerstring('Paths to trajectories',60,'-')
  print '\nPlease enter the paths to all directories containing the "TRAJ_0XXXX" directories.\nE.g. Sing_2/ and Sing_3/. \nPlease enter one path at a time, and type "end" to finish the list.'
  count=0
  paths=[]
  while True:
    path=question('Path: ',str,'end')
    if path=='end':
      if len(paths)==0:
        print 'No path yet!'
        continue
      print ''
      break
    path=os.path.expanduser(os.path.expandvars(path))
    if not os.path.isdir(path):
      print 'Does not exist or is not a directory: %s' % (path)
      continue
    if path in paths:
      print 'Already included.'
      continue
    ls=os.listdir(path)
    print ls
    for i in ls:
      if 'TRAJ' in i:
        count+=1
    print 'Found %i subdirectories in total.\n' % count
    paths.append(path)
  INFOS['paths']=paths
  print 'Total number of subdirectories: %i\n' % (count)

  # try to obtain the maximum number of time steps in the trajectories
  maxlen=0
  dt=0.0
  forbidden=['crashed','running','dead','dont_analyze']
  for idir in INFOS['paths']:
    ls=os.listdir(idir)
    for itraj in ls:
      if not 'TRAJ_' in itraj:
        continue
      path=idir+'/'+itraj
      pathfile=path+'/output.lis'
      if not os.path.isfile(pathfile):
        continue
      lstraj=os.listdir(path)
      valid=True
      for i in lstraj:
        if i.lower() in forbidden:
          valid=False
          break
      if not valid:
        continue
      f=readfile(pathfile)
      for line in f:
        if '#' in line:
          continue
        s=line.split()
        step=int(s[0])
        if step>maxlen:
          maxlen=step
        if dt==0.:
          dt=float(s[1])


  print centerstring('Path to reference structure',60,'-')
  print '\nPlease enter the path to the equilibrium structure of your system (in the same atomic order as that given in the dynamics output)'
  print ''
  refpath=question('Path: ',str,'ref.xyz')
  refpath=os.path.expanduser(os.path.expandvars(refpath))
  print ''
  reftype=question('Please give the type of coordinate file',str,os.path.splitext(refpath)[1][1:])
  INFOS['refstruc']=refpath
  INFOS['reftype']=reftype
  print ''



  massweightques=question('Do you wish to use mass weighted coordinates?',bool,True)
  INFOS['massweight']=massweightques
  print ''


  print centerstring('Number of total steps in your trajectories',60,'-')
  #print '\n total simulation time *2 and +1 if timestep is 0.5 fs'
  print ''
  while True:
    numsteps=question('Number of time steps: ',int,[maxlen+1])[0]
    if numsteps <=0:
      print 'Number of steps must be positive!'
      continue
    break
  INFOS['numsteps']=numsteps
  print ''

  print centerstring('The time step of your calculation',60,'-')
  print ''
  while True:
    timestep=question('Length of time step: ',float,[dt])[0]
    if timestep <=0.:
      print 'Time step must be positive!'
      continue
    break
  INFOS['timestep']=timestep
  print ''

  intervallist=[]
  intervallist=[]
  print centerstring('Time steps to be analysed',60,'-')
  print '\nPlease enter the time step intervals for which the statistical analysis should be carried out. '
  print ''
  while True:
    interval=question('Time step interval: ',int,[0,maxlen])
    #endtime=question('End time of interval: ',int,[maxlen])[0]
    print ''
    #interval=[starttime,endtime]
    intervallist.append(interval)
    moreinterval=question('Do you want to add another time interval for analysis?',bool,False)
    if not moreinterval:
      break
  INFOS['interval']=intervallist
  print ''

  print centerstring('Results directory',60,'-')
  print 'Please give the name of the subdirectory to be used for the results (use to save similar analysis in separate subdirectories).'
  destin=question('Name for subdirectory?',str,'essdyn')
  INFOS['descr']=destin
  print ''
  
  return INFOS


def ess_dyn(INFOS):
    """
    Essential dynamics analysis analysis. Typically this procedure is carried out.
    """
    print 'Preparing essential dynamics analysis ...'
    
    num_steps=INFOS['numsteps']
    descr=INFOS['descr']
    ref_struc_file=INFOS['refstruc']
    ref_struc_type=INFOS['reftype']
    ana_ints=INFOS['interval']
    dt=INFOS['timestep']
    mawe=INFOS['massweight']

    try:
        os.makedirs('ESS_DYN/'+descr+'/total_cov')
    except:
        print 'Output directory could not be created. It either already exists or you have no writing access.'
        
    try:
        os.makedirs('ESS_DYN/'+descr+'/cross_av')
    except:
        print 'Output directory could not be created. It either already exists or you have no writing access.'
    
    ref_struc = struc_linalg.structure('ref_struc') # define the structure that all the time step structures are superimposed onto
    ref_struc.read_file(ref_struc_file, ref_struc_type)
    num_at = ref_struc.ret_num_at()
    mol_calc = struc_linalg.mol_calc(def_file_path=ref_struc_file, file_type=ref_struc_type)
    
    # used for computing the covariance for each pair of coordinates over all trajectories and timesteps
    num_points = numpy.zeros(len(ana_ints)) # number of all timesteps in all the trajectories for each interval analysed
    X_sum_array = numpy.zeros([len(ana_ints), num_at*3], float) # sum for every coordinate in the specified time interval
    XY_sum_array = numpy.zeros([len(ana_ints), num_at*3, num_at*3], float) # a number for every time interval analysed and pair of coordinates

    # used for computing time resolved mean and variance
    cross_num_array = numpy.zeros(num_steps)
    cross_sum_array = numpy.zeros([num_steps,num_at*3], float) # a number for every time step and coordinate; sum, has to be divided by num_array

    forbidden=['crashed','running','dead','dont_analyze']
    width=30
    files=[]
    ntraj=0
    print 'Checking the directories...'
    for idir in INFOS['paths']:
      ls=os.listdir(idir)
      for itraj in ls:
        if not 'TRAJ_' in itraj:
          continue
        path=idir+'/'+itraj
        s=path+' '*(width-len(path))
        pathfile=path+'/output.xyz'
        if not os.path.isfile(pathfile):
          s+='%s NOT FOUND' % (pathfile)
          print s
          continue
        lstraj=os.listdir(path)
        valid=True 
        for i in lstraj:
          if i.lower() in forbidden:
            s+='DETECTED FILE %s' % (i.lower())
            print s
            valid=False
            break
        if not valid:
          continue
        s+='OK'
        print s
        ntraj+=1
        files.append(pathfile)
    print 'Number of trajectories: %i' % (ntraj)
    if ntraj==0:
      print 'No valid trajectories found, exiting...'
      sys.exit(0)

#    Numtraj=(last_traj+1)-first_traj 
    for i in xrange(ntraj):
         
#        ls=os.listdir(os.getcwd())
#        numfile=len(ls)
#        k=i+first_traj
#        string=str(k).rjust(5, '0')
#        trajfolder=None
#        j=0
#        for j in range(numfile):
#            trajfolder=re.search(string,str(ls[j]))
#            if trajfolder!=None:
#               break
#        trajcheck=None
#        if trajfolder !=None:
#           trajcheck=re.search('TRAJ',str(ls[j]))
#        if trajcheck !=None:
           print 'Reading trajectory ' + str(files[i]) + ' ...'
           folder_name = str(files[i])[:-10]
           trajectory = traj_manip.trajectory(folder_name, ref_struc, dt=dt)
    
           coor_matrix = trajectory.ret_coor_matrix()
           # addition for total variance
           for i,interv in enumerate(ana_ints):
               part_mat = coor_matrix[interv[0]:interv[1]]
               num_points[i] += part_mat.shape[0] # not st - en if the matrix does not go until en
               X_sum_array[i] += numpy.add.reduce(part_mat) # add the values in one column vector
               XY_sum_array[i] += numpy.dot(part_mat.transpose(), part_mat)
    
           # addition for time resolved results
           for nr, tstep in enumerate(coor_matrix):
               try:
                   cross_num_array[nr] += 1
                   cross_sum_array[nr] += tstep
               except:
                   print 'num_steps has to be at least as large as the maximum number of time steps in any trajectory!'
                   sys.exit()
                
            
    print 'Processing data ...'
    for ind,num in enumerate(cross_num_array):
        if num == 0:   # if num_steps was set larger than needed
            cross_num_array = cross_num_array[0:ind]
            cross_sum_array = cross_sum_array[0:ind]
            num_steps = ind # num_steps has to be passed as an argument. so it can be changed here.
            break
            
    cross_mean_array = numpy.zeros([num_steps,num_at*3], float) # cross_sum_array / num_array        
    #print len(cross_num_array),len(cross_mean_array)
    #mass_mat = mol_calc.ret_mass_matrix(power=mass_wt_pw/2.)
            
    # total covariance
    for i,interv in enumerate(ana_ints):
        cov_mat_i = numpy.zeros([num_at*3, num_at*3], float)
        exp_X_i = X_sum_array[i] / num_points[i] # expected values of x and x*y
        exp_XY_i = XY_sum_array[i] / num_points[i]
        
        av_struc = mol_calc.make_structure(exp_X_i)
        
        for ii in xrange(3*num_at):
            for iii in xrange(3*num_at):
                cov_mat_i[ii,iii] = exp_XY_i[ii,iii] - exp_X_i[ii]*exp_X_i[iii]
                cov_mat_i[iii,ii] = exp_XY_i[ii,iii] - exp_X_i[ii]*exp_X_i[iii]
        
        if mawe == True:
           mass_mat = mol_calc.ret_mass_matrix(power=0.5)
           cov_mat_mawe=numpy.dot(mass_mat, numpy.dot(cov_mat_i, mass_mat))        
           cov_mat_i=cov_mat_mawe
        #cov_mat_i = numpy.dot(mass_mat, cov_mat_i) # mass weighting
        #cov_mat_i = numpy.dot(cov_mat_i, mass_mat)
        
        cov_eigvals, t_eigvects = numpy.linalg.eigh(cov_mat_i)
        if mawe== True:
           mass_mat_inv=numpy.linalg.inv(mass_mat)
#mass_mat**-1
           unity = numpy.dot(mass_mat,mass_mat_inv)
           t_eigvects_nonmawe = numpy.dot(mass_mat_inv,t_eigvects)
           t_eigvects=t_eigvects_nonmawe
         
        cov_eigvects = t_eigvects.transpose()
        #print cov_eigvals[0]
        #print cov_eigvects[0]
        vib_molden.make_molden_file(struc=av_struc, freqs=cov_eigvals, vibs=cov_eigvects, out_file='ESS_DYN/'+descr+'/total_cov/'+      str(interv[0])+'-'+str(interv[1])+'.molden')
    
    # covariance of time averaged structures
    for i in xrange(num_steps):
        cross_mean_array[i] = cross_sum_array[i] / cross_num_array[i]
        
    for i,interv in enumerate(ana_ints):
        st,en = interv
        cov_mat_i = numpy.zeros([num_at*3, num_at*3], float)
        exp_X_i = numpy.add.reduce(cross_mean_array[st:en]) / (en-st) # expected values of x and x*y
        exp_XY_i = numpy.dot(cross_mean_array[st:en].transpose(), cross_mean_array[st:en]) / (en-st)
        
        av_struc = mol_calc.make_structure(exp_X_i)
        
        for ii in xrange(3*num_at):
            for iii in xrange(3*num_at): 
                cov_mat_i[ii,iii] = exp_XY_i[ii,iii] - exp_X_i[ii]*exp_X_i[iii]
                cov_mat_i[iii,ii] = exp_XY_i[ii,iii] - exp_X_i[ii]*exp_X_i[iii]

        if mawe == True:
           mass_mat = mol_calc.ret_mass_matrix(power=0.5)
           cov_mat_mawe=numpy.dot(mass_mat, numpy.dot(cov_mat_i, mass_mat))
           cov_mat_i=cov_mat_mawe                
#         cov_mat_i = numpy.dot(mass_mat, cov_mat_i) # mass weighting
#         cov_mat_i = numpy.dot(cov_mat_i, mass_mat)
        
        cov_eigvals, t_eigvects = numpy.linalg.eigh(cov_mat_i)
        if mawe==True:
           mass_mat_inv=numpy.linalg.inv(mass_mat)
           unity = numpy.dot(mass_mat,mass_mat_inv)
           t_eigvects_nonmawe = numpy.dot(mass_mat_inv,t_eigvects)
           t_eigvects=t_eigvects_nonmawe
        cov_eigvects = t_eigvects.transpose()
         
        vib_molden.make_molden_file(struc=av_struc, freqs=cov_eigvals, vibs=cov_eigvects, out_file='ESS_DYN/'+descr+'/cross_av/'+str(interv[0])+'-'+str(interv[1])+'.molden')
    

def main():
    '''Main routine'''

    displaywelcome()
    open_keystrokes()

    INFOS=get_general()
    ess_dyn(INFOS)
        
    close_keystrokes()

if __name__=='__main__':
  try:
    main()
  except KeyboardInterrupt:
    print '\nCtrl+C makes me a sad SHARC ;-(\n'
    quit(0)

#    ess_dyn(num_steps=num_steps)
