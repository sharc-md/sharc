"""
version 1.0
author: Felix Plasser
description: Procedures for plotting using the Matplotlib/pylab interface.
"""

import os, sys
import numpy
import pylab
import file_handler

def mean_std_from_files(mean_file, out_dir, xlist, std_file=None, col_list=None):
    """
    Makes a number of mean/std plots from a file containing means and a file containing standard deviations and saves them to directory out_dir.
    Optionally <std_file> can be left out and just the means are plotted.
    If <col_list> is given, a joined plot with the modes contained in <col_list> is given.
    """
    figwidth = 2*8.5/2.54

    try:
        os.makedirs(out_dir)
    except OSError:
        pass

    if not std_file == None:
        with_std = True
    else:
        with_std = False

    mean_lines = open(mean_file, 'r').readlines()
    if with_std: std_lines = open(std_file, 'r').readlines()

    header = [file_handler.line_to_words(line[:-1]) for line in mean_lines[0:3]]

    #print mean_lines[1:2]

    mean_list = []
    if with_std: std_list = []
    for line in mean_lines[3:]:
        mean_list += [[eval(word) for word in file_handler.line_to_words(line[:-1])]]
        
    if with_std:
        for line in std_lines[3:]:
            std_list += [[eval(word) for word in file_handler.line_to_words(line[:-1])]]

    mean_array = numpy.array(mean_list, float).transpose()
    if with_std: std_array = numpy.array(std_list, float).transpose()
    
    #print mean_array
    #print std_array

    if col_list == None:
        for i in xrange(numpy.size(mean_array, 0)):
            pylab.figure(figsize=(figwidth, .75*figwidth))
            pylab.plot(xlist[:len(mean_array[i].tolist())], mean_array[i].tolist()[:len(xlist)], color='k', linewidth=1)
            if with_std: plot_std(mean_array[i].tolist()[:len(xlist)], std_array[i].tolist()[:len(xlist)], xlist[:len(mean_array[i].tolist())])
            
            #pylab.title('20ag'+' - '+str(header[1][i])+'/cm - '+str(header[2][i]+'fs'))
            pylab.title(str(header[0][i])+' - '+str(header[1][i])+'/cm - '+str(header[2][i]+'fs'))

            pylab.xlabel("time (fs)")
            
            pylab.savefig(out_dir + '/mode_' + str(i+1).rjust(3,'0') + '.png', orientation='landscape')
    else:
        linestyles = ['-']#,'--',':','-.']
        colors = ['k', 'b', 'r', 'g']
        
        pylab.figure(figsize=(figwidth, .75*figwidth)) # - temp
        for run_ind, col_ind in enumerate(col_list):
            i = col_ind - 1
            pylab.plot(xlist[:len(mean_array[i].tolist())], mean_array[i].tolist()[:len(xlist)], color=colors[run_ind%len(colors)], linewidth=1, label=str(col_ind), linestyle=linestyles[run_ind%len(linestyles)])
            if with_std: plot_std(mean_array[i].tolist()[:len(xlist)], std_array[i].tolist()[:len(xlist)], xlist[:len(mean_array[i].tolist())])
    
            pylab.legend() # - temp
            pylab.xlabel("time (fs)")
            pylab.ylabel("Displacement (A)")
            
        out_file = 'modes_' + str(col_list[0])
        for col_ind in col_list[1:]:
            out_file += '_' + str(col_ind)
        #out_file += '_new' # + temp
        out_file += '.png' 
            
        pylab.savefig(os.path.join(out_dir, out_file), orientation='landscape')
        #print os.path.join(out_dir, out_file)
        #pylab.show()
            
def plot_std(mean_list, std_list, xlist):
    """
    Plot the standard deviation.
    """
    # mean_list and std_list can be lists or numpy.arrays
    rev_xlist = [x for x in xlist] # reverse list for plotting std
    rev_xlist.reverse()

    #print xlist
    #print rev_xlist
    
    plot_list = [] # list with the values to be plotted
    temp_list = []
    #print mean_list
    for nr,mean in enumerate(mean_list):
        plot_list += [mean + std_list[nr]]
    for nr,mean in enumerate(mean_list):
        #print mean, std_list[nr]
        temp_list += [mean - std_list[nr]]

    temp_list.reverse()
    plot_list += temp_list

    
    pylab.fill(xlist + rev_xlist, plot_list, alpha=.3, facecolor='grey')

def bars_from_file(in_file, out_dir):
    """
    Create a bar graph from data in a file.
    """
    figwidth = 2*8.5/2.54

    try:
        os.makedirs(out_dir)
    except OSError:
        pass

    bar_lines = open(in_file, 'r').readlines()
    
    header = [file_handler.line_to_words(line[:-1]) for line in bar_lines[0:3]]
    
    bar_list = []
    bar_titles = []
    for line in bar_lines[3:]:
        t_list = file_handler.line_to_words(line[:-1])
        bar_list += [[eval(word) for word in t_list[1:]]]
        bar_titles += [t_list[0]]

    #print mean_lines[1:2]

    for i,bars in enumerate(bar_list):
        pylab.figure(figsize=(2*figwidth, .75*figwidth))
        xlist = [eval(num) for num in header[0][1:]]
        pylab.bar(xlist, bars)
        
        pylab.title(bar_titles[i])
        
        pylab.savefig(out_dir + '/' + bar_titles[i] + '.png', orientation='landscape')

def mean_std_from_file(file, out_dir, col_list=None, num_steps=None):
    """
    Makes a number of mean/std plots from one file containing means standard deviations and saves them to directory out_dir.
    mean_value.3 can be plotted with this routine.
    If <col_list> is given, a joined plot with the modes contained in <col_list> is given.
    <num_steps> gives the maximum number of plotted steps.
    """
    figwidth = 2*8.5/2.54

    try:
        os.makedirs(out_dir)
    except OSError:
        pass

    lines = open(file, 'r').readlines()
    
    read_list = []
    for line in lines[:num_steps]:
        words = file_handler.line_to_words(line)
        if not words[1] == '0':
            read_list += [[eval(word) for word in words]]
        else:
            print 'No trajectories at timestep ' + str(words[0])
        
    plot_array = numpy.array(read_list, float).transpose()
    
    if col_list == None:
        for i in xrange(numpy.size(plot_array, 0)/2-1):
            pylab.figure(figsize=(figwidth, .75*figwidth))
            pylab.plot(plot_array[0].tolist(), plot_array[2+2*i].tolist(), color=[0.,0.,0.], linewidth=1)
            plot_std(plot_array[2+2*i].tolist(), plot_array[3+2*i].tolist(), plot_array[0].tolist())
    
            pylab.title(str(i+1))
            pylab.xlabel("time (fs)")
            
            pylab.savefig(out_dir + '/int_coor_' + str(i+1).rjust(3,'0') + '.png', orientation='landscape')
    else:
        linestyles = ['-','--',':','-.']
        
        pylab.figure(figsize=(figwidth, .75*figwidth))
        for run_ind, col_ind in enumerate(col_list):
            i = col_ind - 1
            pylab.plot(plot_array[0].tolist(), plot_array[2+2*i].tolist(), color=[0.,0.,0.], linewidth=1, label=str(col_ind), linestyle=linestyles[run_ind%4])
            plot_std(plot_array[2+2*i].tolist(), plot_array[3+2*i].tolist(), plot_array[0].tolist())
            
            pylab.legend()
            pylab.xlabel("time (fs)")
            
        out_file = 'int_coors_' + str(col_list[0])
        for col_ind in col_list[1:]:
            out_file += '_' + str(col_ind)
        
            
        pylab.savefig(os.path.join(out_dir, out_file), orientation='landscape')        
        
def nma_positions(xind, yind, i_tst_list, symbol='ko', file_name='nma.txt'):
    """
    Plot time steps specified in nested list <i_t_list> on the surface of normal coordinates with indices <xind> and <yind>.
    """
    
    # not efficiently written but time is not critical
    
    for i,tstep in i_tst_list:
        lines = open('TRAJ' + str(i) + '/RESULTS/' + file_name).readlines()
        vec = [eval(word) for word in file_handler.line_to_words(lines[tstep+3])]
        pylab.plot([vec[xind + 1]], [vec[yind + 1]], symbol)
       
if __name__ == '__main__':
    nma_positions(10, 14, [[4, 200], [4, 201], [4, 502]], file_name='no_abs_nma.txt')
    #pylab.show()
    
if False:
    mean_std_from_files(mean_file='NMA/mean_against_time.txt',std_file='NMA/std_against_time.txt',out_dir='NMA/time_plots',xlist=[.5 * i for i in xrange(601)])
