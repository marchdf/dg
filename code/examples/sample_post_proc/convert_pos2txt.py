#!/usr/bin/env python
#
# Convert some data to the gmsh txt format
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
order = '2'
A0 = 0.03
configurations = [0.33333333]
lengths   = [0.5,1,2,3,4];
strengths = [0.5,0.3,0.1,0.05];
machs     = [1.2];
dirnames = []
for config in configurations:
    for length in lengths:
        for strength in strengths: 
            for mach in machs:  
                dirnames.extend(['simblst_'+str(length)+'_'+str(strength)+'_'+str(mach)+'_'+str(A0)+'_'+str(config)])

p2tdirs = dirnames
NUM_PROCS = '8'
NUM_T='100'

fields = ['rho','p','ux','uy']

# Directory containing sample pos2txt.geo
p2tfdir = 'simblst/'
p2tfile = 'pos2txt.geo'

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
num2submit = len(p2tdirs)
cnt = 1
print 'Total number of pos2txt jobs to be submitted = {0:.0f}'.format(num2submit)

#================================================================================
#
# Loop on directories and do the conversions
#
#================================================================================
for k,p2tdir in enumerate(p2tdirs):

    print '\nSubmitting pos2txt job [{0:.0f}/{1:.0f}]'.format(cnt,num2submit); cnt += 1;

    JOBNAME = 'P2T_'+p2tdir.replace('/','_')
    SUBMIT=[JOBNAME,RUNTIME,QUEUE,RESOURCES]

    # Go to the work directory
    WORKDIR = DATADIR+p2tdir
    os.chdir(WORKDIR)
    print 'Entering', WORKDIR

    # initialize
    CMD = ''

    for field in fields:

        # Field specific file that will do the conversion
        p2tfile_field = p2tfile+'_'+field

        # create command to run
        CMD += GMSH+p2tfile_field+';\n'        

        # copy the post-proc file into our directory
        try:
            shutil.copyfile(BASEDIR+'/'+p2tfdir+p2tfile, p2tfile_field)
        except IOError, e:
            print "Unable to copy p2t file. %s" % e
            sys.exit()

        # replace the number of processors that were used for the run and
        # the field that we want to convert
        auxf.replace_in_file(p2tfile_field,{'NUM_PROCS':NUM_PROCS,'FIELD':field,'ORDER':order,'NUM_T':NUM_T})

    # write the submit file in our directory
    submitname = auxf.write_submit(WORKDIR,SUBMIT,CMD);
    
    # submit the run to the queue
    call('sbatch '+submitname, shell=True)

    # go back to base directory
    os.chdir(BASEDIR)
