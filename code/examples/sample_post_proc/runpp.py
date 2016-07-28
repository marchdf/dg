#!/usr/bin/env python
#
# Run a bunch of gmsh and python post-processors on a set of problems
# 
__author__ = 'marchdf'
def usage():
    print '\nUsage: {0:s} [options ...]\nOptions:\n -d, --development\tsubmit to development queue\n -h, --help\tshow this message and exit\n'.format(sys.argv[0])

#================================================================================
#
# Imports
#
#================================================================================
import sys,getopt,os,time,shutil
from numpy import *
from subprocess import call
import subprocess as sp
import time
import aux_functions as auxf

#================================================================================
#
# Functions
#
#================================================================================

#================================================================================
#
# Parse arguments
#
#================================================================================
development = False
try:
    myopts, args = getopt.getopt(sys.argv[1:],"hd",["help","development"])
except getopt.GetoptError as e:
    print (str(e))
    usage()
    sys.exit(2)

for o, a in myopts:
    if o in ('-h', '--help'):
        usage()
        sys.exit()
    elif o in ('-d', '--development'):
        development = True

#================================================================================
#
# Directories to post-process
#
#================================================================================

# ppdirs = ['lores/simblst_1_0.1_1_0.03_0.33333333','lores/simblst_1_0.1_1.2_0.03_0.33333333','lores/simblst_1_0.1_1.5_0.03_0.33333333','lores/simblst_1_0.1_2_0.03_0.33333333','lores/simblst_1_0.1_2.5_0.03_0.33333333','lores/simblst_1_0.1_3_0.03_0.33333333','lores/simblst_1_0.1_1_0.03_3','lores/simblst_1_0.1_1.2_0.03_3','lores/simblst_1_0.1_1.5_0.03_3','lores/simblst_1_0.1_2_0.03_3','lores/simblst_1_0.1_2.5_0.03_3','lores/simblst_1_0.1_3_0.03_3']
# NUM_PROCS = '2'

#ppdirs = ['simblst_1_0.1_1_0.03_0.33333333','simblst_1_0.1_1.2_0.03_0.33333333','simblst_1_0.1_1.5_0.03_0.33333333','simblst_1_0.1_1_0.03_0.33333333','simblst_1_0.1_2.5_0.03_0.33333333','simblst_1_0.1_3_0.03_0.33333333','simblst_1_0.1_1_0.03_3','simblst_1_0.1_1.2_0.03_3','simblst_1_0.1_1.5_0.03_3','simblst_1_0.1_1_0.03_3','simblst_1_0.1_2.5_0.03_3','simblst_1_0.1_3_0.03_3']
#ppdirs = ['simblst_2_0.1_1_0.03_0.33333333','simblst_2_0.1_1.2_0.03_0.33333333','simblst_2_0.1_1.5_0.03_0.33333333','simblst_2_0.1_2_0.03_0.33333333','simblst_2_0.1_2.5_0.03_0.33333333','simblst_2_0.1_3_0.03_0.33333333','simblst_2_0.1_1_0.03_3','simblst_2_0.1_1.2_0.03_3','simblst_2_0.1_1.5_0.03_3','simblst_2_0.1_2_0.03_3','simblst_2_0.1_2.5_0.03_3','simblst_2_0.1_3_0.03_3']

#ppdirs = ['simblst_1_0.05_1_0.03_0.33333333','simblst_1_0.05_1_0.03_3','simblst_1_0.05_1.2_0.03_0.33333333','simblst_1_0.05_1.2_0.03_3','simblst_1_0.05_1.5_0.03_0.33333333','simblst_1_0.05_1.5_0.03_3','simblst_1_0.05_2_0.03_0.33333333','simblst_1_0.05_2_0.03_3','simblst_1_0.05_2.5_0.03_0.33333333','simblst_1_0.05_2.5_0.03_3','simblst_1_0.05_3_0.03_0.33333333','simblst_1_0.05_3_0.03_3','simblst_1_0.1_1_0.03_0.33333333','simblst_1_0.1_1_0.03_3','simblst_1_0.1_1.2_0.03_0.33333333','simblst_1_0.1_1.2_0.03_3','simblst_1_0.1_1.5_0.03_0.33333333','simblst_1_0.1_1.5_0.03_3','simblst_1_0.1_2_0.03_0.33333333','simblst_1_0.1_2_0.03_3','simblst_1_0.1_2.5_0.03_0.33333333','simblst_1_0.1_2.5_0.03_3','simblst_1_0.1_3_0.03_0.33333333','simblst_1_0.1_3_0.03_3','simblst_1_0.3_1_0.03_0.33333333','simblst_1_0.3_1_0.03_3','simblst_1_0.3_1.2_0.03_0.33333333','simblst_1_0.3_1.2_0.03_3','simblst_1_0.3_1.5_0.03_0.33333333','simblst_1_0.3_1.5_0.03_3','simblst_1_0.3_2_0.03_0.33333333','simblst_1_0.3_2_0.03_3','simblst_1_0.3_2.5_0.03_0.33333333','simblst_1_0.3_2.5_0.03_3','simblst_1_0.3_3_0.03_0.33333333','simblst_1_0.3_3_0.03_3','simblst_1_0.5_1_0.03_0.33333333','simblst_1_0.5_1_0.03_3','simblst_1_0.5_1.2_0.03_0.33333333','simblst_1_0.5_1.2_0.03_3','simblst_1_0.5_1.5_0.03_0.33333333','simblst_1_0.5_1.5_0.03_3','simblst_1_0.5_2_0.03_0.33333333','simblst_1_0.5_2_0.03_3','simblst_1_0.5_2.5_0.03_0.33333333','simblst_1_0.5_2.5_0.03_3','simblst_1_0.5_3_0.03_0.33333333','simblst_1_0.5_3_0.03_3']
#ppdirs = ['simblst_2_0.05_1_0.03_0.33333333','simblst_2_0.05_1_0.03_3','simblst_2_0.05_1.2_0.03_0.33333333','simblst_2_0.05_1.2_0.03_3','simblst_2_0.05_1.5_0.03_0.33333333','simblst_2_0.05_1.5_0.03_3','simblst_2_0.05_2_0.03_0.33333333','simblst_2_0.05_2_0.03_3','simblst_2_0.05_2.5_0.03_0.33333333','simblst_2_0.05_2.5_0.03_3','simblst_2_0.05_3_0.03_0.33333333','simblst_2_0.05_3_0.03_3','simblst_2_0.1_1_0.03_0.33333333','simblst_2_0.1_1_0.03_3','simblst_2_0.1_1.2_0.03_0.33333333','simblst_2_0.1_1.2_0.03_3','simblst_2_0.1_1.5_0.03_0.33333333','simblst_2_0.1_1.5_0.03_3','simblst_2_0.1_2_0.03_0.33333333','simblst_2_0.1_2_0.03_3','simblst_2_0.1_2.5_0.03_0.33333333','simblst_2_0.1_2.5_0.03_3','simblst_2_0.1_3_0.03_0.33333333','simblst_2_0.1_3_0.03_3','simblst_2_0.3_1_0.03_0.33333333','simblst_2_0.3_1_0.03_3','simblst_2_0.3_1.2_0.03_0.33333333','simblst_2_0.3_1.2_0.03_3','simblst_2_0.3_1.5_0.03_0.33333333','simblst_2_0.3_1.5_0.03_3','simblst_2_0.3_2_0.03_0.33333333','simblst_2_0.3_2_0.03_3','simblst_2_0.3_2.5_0.03_0.33333333','simblst_2_0.3_2.5_0.03_3','simblst_2_0.3_3_0.03_0.33333333','simblst_2_0.3_3_0.03_3','simblst_2_0.5_1_0.03_0.33333333','simblst_2_0.5_1_0.03_3','simblst_2_0.5_1.2_0.03_0.33333333','simblst_2_0.5_1.2_0.03_3','simblst_2_0.5_1.5_0.03_0.33333333','simblst_2_0.5_1.5_0.03_3','simblst_2_0.5_2_0.03_0.33333333','simblst_2_0.5_2_0.03_3','simblst_2_0.5_2.5_0.03_0.33333333','simblst_2_0.5_2.5_0.03_3','simblst_2_0.5_3_0.03_0.33333333','simblst_2_0.5_3_0.03_3']
#ppdirs = ['simblst_3_0.05_1_0.03_0.33333333','simblst_3_0.05_1_0.03_3','simblst_3_0.05_1.2_0.03_0.33333333','simblst_3_0.05_1.2_0.03_3','simblst_3_0.05_1.5_0.03_0.33333333','simblst_3_0.05_1.5_0.03_3','simblst_3_0.05_2_0.03_0.33333333','simblst_3_0.05_2_0.03_3','simblst_3_0.05_2.5_0.03_0.33333333','simblst_3_0.05_2.5_0.03_3','simblst_3_0.05_3_0.03_0.33333333','simblst_3_0.05_3_0.03_3','simblst_3_0.1_1_0.03_0.33333333','simblst_3_0.1_1_0.03_3','simblst_3_0.1_1.2_0.03_0.33333333','simblst_3_0.1_1.2_0.03_3','simblst_3_0.1_1.5_0.03_0.33333333','simblst_3_0.1_1.5_0.03_3','simblst_3_0.1_2_0.03_0.33333333','simblst_3_0.1_2_0.03_3','simblst_3_0.1_2.5_0.03_0.33333333','simblst_3_0.1_2.5_0.03_3','simblst_3_0.1_3_0.03_0.33333333','simblst_3_0.1_3_0.03_3','simblst_3_0.3_1_0.03_0.33333333','simblst_3_0.3_1_0.03_3','simblst_3_0.3_1.2_0.03_0.33333333','simblst_3_0.3_1.2_0.03_3','simblst_3_0.3_1.5_0.03_0.33333333','simblst_3_0.3_1.5_0.03_3','simblst_3_0.3_2_0.03_0.33333333','simblst_3_0.3_2_0.03_3','simblst_3_0.3_2.5_0.03_0.33333333','simblst_3_0.3_2.5_0.03_3','simblst_3_0.3_3_0.03_0.33333333','simblst_3_0.3_3_0.03_3','simblst_3_0.5_1_0.03_0.33333333','simblst_3_0.5_1_0.03_3','simblst_3_0.5_1.2_0.03_0.33333333','simblst_3_0.5_1.2_0.03_3','simblst_3_0.5_1.5_0.03_0.33333333','simblst_3_0.5_1.5_0.03_3','simblst_3_0.5_2_0.03_0.33333333','simblst_3_0.5_2_0.03_3','simblst_3_0.5_2.5_0.03_0.33333333','simblst_3_0.5_2.5_0.03_3','simblst_3_0.5_3_0.03_0.33333333','simblst_3_0.5_3_0.03_3']
# ppdirs = ['simblst_3_0.1_2_0.03_0.33333333']
# NUM_PROCS = '4'

#ppdirs = ['simblst_3_0.1_1_0.03_0.33333333','simblst_3_0.1_1.2_0.03_0.33333333','simblst_3_0.1_1.5_0.03_0.33333333','simblst_3_0.1_2_0.03_0.33333333','simblst_3_0.1_2.5_0.03_0.33333333','simblst_3_0.1_3_0.03_0.33333333','simblst_3_0.1_1_0.03_3','simblst_3_0.1_1.2_0.03_3','simblst_3_0.1_1.5_0.03_3','simblst_3_0.1_2_0.03_3','simblst_3_0.1_2.5_0.03_3','simblst_3_0.1_3_0.03_3']
#ppdirs = ['simblst_3_0.1_2_0.03_0.33333333','simblst_3_0.1_2_0.03_3']
#NUM_PROCS = '4'

# ppdirs = ['hires/simblst_1_0.1_1_0.03_0.33333333','hires/simblst_1_0.1_1.2_0.03_0.33333333','hires/simblst_1_0.1_1.5_0.03_0.33333333','hires/simblst_1_0.1_2_0.03_0.33333333','hires/simblst_1_0.1_2.5_0.03_0.33333333','hires/simblst_1_0.1_3_0.03_0.33333333','hires/simblst_1_0.1_1_0.03_3','hires/simblst_1_0.1_1.2_0.03_3','hires/simblst_1_0.1_1.5_0.03_3','hires/simblst_1_0.1_2_0.03_3','hires/simblst_1_0.1_2.5_0.03_3','hires/simblst_1_0.1_3_0.03_3']
# NUM_PROCS = '8'

# ppdirs = ['hhres/simblst_1_0.1_1_0.03_0.33333333','hhres/simblst_1_0.1_1.2_0.03_0.33333333','hhres/simblst_1_0.1_1.5_0.03_0.33333333','hhres/simblst_1_0.1_2_0.03_0.33333333','hhres/simblst_1_0.1_2.5_0.03_0.33333333','hhres/simblst_1_0.1_3_0.03_0.33333333','hhres/simblst_1_0.1_1_0.03_3','hhres/simblst_1_0.1_1.2_0.03_3','hhres/simblst_1_0.1_1.5_0.03_3','hhres/simblst_1_0.1_2_0.03_3','hhres/simblst_1_0.1_2.5_0.03_3','hhres/simblst_1_0.1_3_0.03_3']
# NUM_PROCS = '16'


# For the convergence runs
#ppdirs = ['simblst_convergence/'+s for s in ['simblst_100_1', 'simblst_100_2', 'simblst_25_1', 'simblst_25_2', 'simblst_50_1', 'simblst_50_2'] ] # 'simblst_200_1', 'simblst_200_2', 
#NUM_PROCS = '4'
# ppdirs = ['simblst_convergence/'+s for s in ['simblst_50_0'] ]
# NUM_PROCS = '4'
# ppdirs = ['simblst_convergence/'+s for s in ['simblst_25_2','simblst_50_2','simblst_100_2','simblst_25_1','simblst_50_1','simblst_100_1'] ]
# NUM_PROCS = '4'
# NUM_T='100'

# ppdirs = ['simblst_convergence/'+s for s in ['simblst_200_1']]
# NUM_PROCS = '8'
# NUM_T='100'

# ppdirs = ['simblst_convergence/'+s for s in ['simblst_200_2']]
# NUM_PROCS = '16'
# NUM_T='100'

# ppdirs = ['simblst_convergence/'+s for s in ['simblst_25_2','simblst_50_2','simblst_100_2','simblst_200_2','simblst_25_1','simblst_50_1','simblst_100_1','simblst_200_1'] ]

# ppdirs = ['rminstb_convergence/'+s for s in ['rminstb_16_1', 'rminstb_16_2', 'rminstb_32_1', 'rminstb_32_2', 'rminstb_64_1', 'rminstb_64_2', 'rminstb_128_1', 'rminstb_128_2', 'rminstb_256_1', 'rminstb_256_2'] ]
ppdirs=['rminstb_convergence/rminstb_128_2']
NUM_PROCS = '1'
NUM_T='200'


# A0 = 0.03
# configurations = [0.33333333]
# lengths   = [0.5,1,2,3,4];
# strengths = [0.5,0.3,0.1,0.05];
# machs     = [1.2];
# dirnames = []
# for config in configurations:
#     for length in lengths:
#         for strength in strengths: 
#             for mach in machs:  
#                 dirnames.extend(['simblst_'+str(length)+'_'+str(strength)+'_'+str(mach)+'_'+str(A0)+'_'+str(config)])
# ppdirs = dirnames
# NUM_PROCS = '8'
# NUM_T='100'


# directories containing the postprocessing files
ppfdir = 'simblst/'
gmsh_ppfiles = []
python_ppfiles = []
gmsh_ppfiles = ['get_interface.geo','circ_enst_full.geo','interface_edgemid.geo','interface_minmax.geo']
#gmsh_ppfiles = ['ddt_circ_full.geo']
#python_ppfiles = ['circ_dGdt_full.py']

# gmsh_ppfiles = ['get_interface.geo','circ_enst_half.geo','interface_edgemid.geo','interface_minmax.geo','ddt_circ_half.geo']
# python_ppfiles = ['circ_dGdt_half.py']

#================================================================================
#
# Default directories
#
#================================================================================
BASEDIR=os.getcwd()
DATADIR=os.environ['SCRATCH']+'/'

#================================================================================
#
# Define the resources to use
#
#================================================================================
RUNTIME="24:00:00"
NP = 1
QUEUE="normal"
if development:
    QUEUE="development"
RESOURCES='-n '+str(NP)

#================================================================================
#
# Commands that we will need
#
#================================================================================
GMSH = '/home1/02366/marchdf/gmsh-2.11.0-Linux/bin/gmsh - '

# counters
num2submit = len(ppdirs)
cnt = 1
print 'Total number of post-processing jobs to be submitted = {0:.0f}'.format(num2submit)

#================================================================================
#
# Loop on all the directories to postprocess
#
#================================================================================
for k,ppdir in enumerate(ppdirs):

    print '\nSubmitting post-processing job [{0:.0f}/{1:.0f}]'.format(cnt,num2submit); cnt += 1;

    JOBNAME = 'PP_'+ppdir.replace('/','_')
    SUBMIT=[JOBNAME,RUNTIME,QUEUE,RESOURCES]

    # Go to the work directory
    WORKDIR = DATADIR+ppdir
    os.chdir(WORKDIR)
    print 'Entering', WORKDIR

    # initialize
    CMD = ''

    # loop on gmsh post-processing files
    for gmsh_ppfile in gmsh_ppfiles:
        
        # create command to run
        CMD += GMSH+gmsh_ppfile+';\n'

        # copy the post-proc file into our directory
        try:
            shutil.copyfile(BASEDIR+'/'+ppfdir+gmsh_ppfile, gmsh_ppfile)
        except IOError, e:
            print "Unable to copy gmsh pp file. %s" % e
            sys.exit()

        # replace the number of processors that were used for the run
        auxf.replace_in_file(gmsh_ppfile,{'NUM_PROCS':NUM_PROCS,'NUM_T':NUM_T})

    # loop on python post-processing files
    for python_ppfile in python_ppfiles:
        
        # create command to run
        CMD += 'python '+python_ppfile+';\n'

        # copy the post-proc file into our directory
        try:
            shutil.copyfile(BASEDIR+'/'+ppfdir+python_ppfile, python_ppfile)
        except IOError, e:
            print "Unable to copy python pp file. %s" % e
            sys.exit()

    # write the submit file in our directory
    submitname = auxf.write_submit(WORKDIR,SUBMIT,CMD);

    # submit the run to the queue
    call('sbatch '+submitname, shell=True)

    # go back to base directory
    os.chdir(BASEDIR)


    
