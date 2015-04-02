#!/usr/bin/env python2
#
# Read and format the output of a timers.dat file
#
# TODO: make it capable of processing timers from MPI runs
__author__ = 'marchdf'
def usage(): 
    print '\nUsage: {0:s} [options ...]\nOptions:\n  -d, --dir\tdirectory to post-process (default to cwd)\n  -h, --help\tshow this message and exit\n'.format(sys.argv[0])

import sys, getopt
from numpy import *
import matplotlib.axis as axis
from pylab import*
import os

# Parse command line arguments using getopt (argparse is better but require python 1.7 at least)
# from: http://www.cyberciti.biz/faq/python-command-line-arguments-argv-example/
fdir = './' # post-process current dir by default
try:
    myopts, args = getopt.getopt(sys.argv[1:],"hd:",["help","dir="])
except getopt.GetoptError as e:
    print (str(e))
    usage()
    sys.exit(2)
 
for o, a in myopts:
    if o in ('-h', '--help'):
        usage()
        sys.exit()
    elif o == '-d':
        fdir=a+'/'
print 'Parsing directory', fdir

rc('text', usetex=True)
rc('font', family='serif', serif='Times')

#cmap =['#F15A60','#7AC36A','#5A9BD4','#FAA75B','#9E67AB','#CE7058','#D77FB4','#737373']
cmap =['#EE2E2F','#008C48','#185AA9','#F47D23','#662C91','#A21D21','#B43894','#010202']
markertype = ['s','d','o','p','h']
linetype = ['-','--','-.',':']
dashseq = [[1,0],[10,5],[10, 4, 3, 4],[3, 3],[10, 4, 3, 4, 3, 4],[3, 3],[3, 3]];

def get_timer(fname,name):
    # Returns the time and count of a timer called name in file fname
    with open(fdir+fname,"r") as f:
        for line in f: 
            if name in line:
                info = next(f).split()
                return float(info[0]), int(info[1])
    return 0,0

def print_children(fname,parent_time,children,other_name):
    # Format the output for the subroutines in a given list (children)
    # wrt a parent time. Returns the sorted times, percentages, and
    # the name of the children associated to those times.

    # Get the children timers
    times = []
    for k,child in enumerate(children):
        t,_ = get_timer(fname,child);
        times.append(t)

    # Append the leftover times + a name in children
    times.append(parent_time-sum(times))
    children.append(other_name)

    # Sort in descending order
    idx   = argsort(times)[::-1]
    times = sort(times)[::-1]
    times_percent = times/parent_time*100

    # Output
    for k in range(0,len(children)):
        print '\t{0:6.2f}% time (or {1:6.2f}s) spent in {2:s}'.format(times_percent[k],times[k],children[idx[k]])

    return times,times_percent,array(children)[idx]


# Get the timer files 
tfiles = [fname for fname in os.listdir(fdir) if fname.startswith("timers.dat")]

# Loop over them all
for cnt, fname in enumerate(tfiles):
    
    print '\n================================================================================'
    print 'Parsing the timers for processor {0:d}'.format(cnt)

    # Main program
    main,_ = get_timer(fname,'main')
    print 'Main program ({0:.2f}s):'.format(main) 
    print_children(fname,main,['RK'],'other (setup,...)')

    # RK portion
    rk,_   = get_timer(fname,'RK')
    rk_children = ['output','solver','DG','new_dt','limiting','communication','advect_particles']
    print 'RK subroutines ({0:.2f}s):'.format(rk) 
    print_children(fname,rk,rk_children,'other (setup,copies,axpys,...)')
    
    # DG portion
    dg,_   = get_timer(fname,'DG')
    dg_children = ['mapToFace','rflctiveBoundary','collocationU','collocationdU','collocationUF','evaluate_sf','evaluate_q','redistribute_sf','gemm_s','gemm_f','redistribute_q','gemm_q','mapToElement','addSFQ']
    print 'DG subroutines ({0:.2f}s):'.format(dg) 
    print_children(fname,dg,dg_children,'other')
    
    # Limiting portion
    limiting,_   = get_timer(fname,'limiting')
    limiting_children = ['sensors']
    print 'limiting subroutines ({0:.2f}s):'.format(limiting) 
    print_children(fname,limiting,limiting_children,'actual limiting')

    # Communication
    communication,_   = get_timer(fname,'communication')
    communication_children = ['communication','comm_packaging','comm_unpackaging','comm_memcpy','comm_sendrecv']
    print 'communication subroutines ({0:.2f}s):'.format(communication) 
    print_children(fname,communication,communication_children,'other (barriers, waits,...)')

    # Output
    output,_   = get_timer(fname,'output')
    output_children = ['format_output','write_output','format_sensor','write_sensor','write_particles']
    print 'output subroutines ({0:.2f}s):'.format(output) 
    print_children(fname,output,output_children,'other')


# # Plot
# figure(0);
# p = plot((t-offset[cnt])*tConvert,A*gfConvert,color=cmap[cnt],lw=3,ls='-')
# p[0].set_dashes(dashseq[cnt])

# # Experiment
# figure(0);
# plot(t,A,markertype[0],markeredgecolor=cmap[cnt],markerfacecolor=cmap[cnt],markersize=5);

# # Legends
# figure(0);
# text(10.1, 4.5,r"$\frac{h}{\lambda}=2$", ha='left', va='center',fontsize=22,fontweight='bold', color=cmap[0])
# text(10.1, 5.85,r"$\frac{h}{\lambda}=4$", ha='left', va='center',fontsize=22,fontweight='bold', color=cmap[1])
# text(10.1, 6.6,r"$\frac{h}{\lambda}=6$", ha='left', va='center',fontsize=22,fontweight='bold', color=cmap[2])

# figure(0);
# xlabel(r"$\mathbf{t~[ms]}$",fontsize=22,fontweight='bold')
# ylabel(r"$\mathbf{a(t)-a_0~[cm]}$",fontsize=22,fontweight='bold')
# setp(gca().get_xmajorticklabels(),fontsize=18,fontweight='bold');
# setp(gca().get_ymajorticklabels(),fontsize=18,fontweight='bold');
# #setp(gca(),xlim=[-0.1,11],ylim=[-0.5,6.8]) 
# setp(gca(),xlim=[-0,6.5],ylim=[0,2.5]) 
# #gca().set_xticks(arange(0,11,2))
# #gca().set_yticks(arange(0,7,2))
# #gca().spines['right'].set_color('none')
# #gca().spines['top'].set_color('none')
# #gca().xaxis.set_ticks_position('bottom')
# #gca().yaxis.set_ticks_position('left')
# savefig('jacobsgf.pdf',format='pdf') 
# savefig('jacobsgf.png',format='png')
# savefig('jacobsgf.eps',format='eps')
# show()


