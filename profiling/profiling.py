#!/usr/bin/env python
#
#
#
__author__ = 'marchdf'
def usage():
    print '\nUsage: {0:s} [options ...]\nOptions:\n -f, --force\tforce recompile\n -h, --help\tshow this message and exit\n'.format(sys.argv[0]) 

#================================================================================
#
# Imports
#
#================================================================================
import sys,getopt,os,time,shutil
from numpy import *
from subprocess import call

#================================================================================
#
# Parse arguments
#
#================================================================================
recompile = False
try:
    myopts, args = getopt.getopt(sys.argv[1:],"hf",["help","force"])
except getopt.GetoptError as e:
    print (str(e))
    usage()
    sys.exit(2)

for o, a in myopts:
    if o in ('-h', '--help'):
        usage()
        sys.exit()
    elif o in ('-f', '--force'):
        recompile = True

#================================================================================
#
# Define some functions
#
#================================================================================
def write_deck(WORKDIR,preprocessors,defs):
    deckname = WORKDIR+'/deck.inp'
    with open(deckname, "w") as f:
        f.write('#Make options\n')
        f.write('LIEU='+preprocessors[0]+'\n');
        f.write('ARCH='+preprocessors[1]+'\n');
        f.write('PROF='+preprocessors[2]+'\n');
        f.write('PARA='+preprocessors[3]+'\n');
        f.write('PREC='+preprocessors[4]+'\n');
        f.write('DIMS='+preprocessors[5]+'\n');
        f.write('PROB='+defs[0]+'\n');
        f.write('NFLD='+defs[1]+'\n');
        f.write('GAMM='+defs[2]+'\n');
        f.write('FLUX='+defs[3]+'\n');
        f.write('MASS='+defs[4]+'\n');
        f.write('#Code options\n')
        f.write('#time integration\nrk4\n')
        f.write('#output time step size\n'+defs[8]+'\n')
        f.write('#final time\n'+defs[9]+'\n')
        f.write('#Courant-Friedrichs-Lewy condition\n'+defs[10]+'\n')
        f.write('#order\n'+defs[5]+'\n')
        f.write('#mesh file\nmesh.msh\n'+defs[11]+'\n')
        f.write('#limiter\n'+defs[7]+'\n')
        f.write('#initial condition\n'+defs[6]+'\n')
    return deckname


#================================================================================
def add_lagrange_deck(deckname,lagrange_particles):
    with open(deckname, "a") as f:
        f.write('#lagrange particles\n'+lagrange_particles+'\n')

#================================================================================
def add_sensor_deck(deckname,sensor_thresholds):
    with open(deckname, "a") as f:
        f.write('#sensor thresholds\n'+sensor_thresholds+'\n')
        
#================================================================================
def get_make_command(deckname):
    makecmd='make '
    with open(deckname) as f: 
        for line in f:
            if '#Code options' in line:
                return makecmd
            if '#Make options' not in line:
                makecmd = makecmd + line.rstrip('\n') + ' '

#================================================================================
def test_recompile(CODEDIR,recompile,makecmd):

    if recompile:
        return True
    else: 
        with open(CODEDIR+'last_makecmd') as f:
            last_makecmd = f.readline()
        if makecmd == last_makecmd:
            return False
        else:
            return True
        
#================================================================================
def compile_exec(recompile,CODEDIR,WORKDIR,makecmd):
    if recompile:
        os.chdir(CODEDIR)
        return_code = call('make clean', shell=True)
        return_code = call(makecmd, shell=True)
        if return_code != 0: 
            print '\nFailure at compilation\n'
            shutil.rmtree(WORKDIR, ignore_errors=True)
            sys.exit()
    print 'Compiled successfully!'
    with open(CODEDIR+'last_makecmd', "w") as f:
        f.write(makecmd)
    return True

#================================================================================
def write_submit(WORKDIR,SUBMIT,CMD):
    submitname = WORKDIR+'/submit.batch'
    with open(submitname, "w") as f:
        f.write('#!/bin/bash\n\n')
        f.write('# Taken from https://www.tacc.utexas.edu/user-services/user-guides/stampede-user-guidef.write#running\n\n')
        f.write('#SBATCH -J '+SUBMIT[0]+'         # job name\n')
        f.write('#SBATCH -o '+SUBMIT[0]+'.o%j     # output and error file name (%j expands to jobID)\n')
        f.write('#SBATCH -p '+SUBMIT[2]+'\n')
        f.write('#SBATCH '+SUBMIT[3]+'\n')
        f.write('#SBATCH -A TG-CTS130005\n')
        f.write('#SBATCH -t '+SUBMIT[1]+'       # run time\n')
        f.write('#SBATCH --mail-user=marchdf@umich.edu\n')
        f.write('#SBATCH --mail-type=begin  # email me when the job starts\n')
        f.write('#SBATCH --mail-type=end    # email me when the job finishes\n\n')
        f.write('module load cuda\n\n')
        f.write('#  Put your job commands after this line\n')
        f.write(CMD+'\n')
    return submitname

#================================================================================
#
# PARAMETERS to play with 
#
#================================================================================

# resolution
resolution = 128;

# method order (number from 1 to 5)
order      = 1;

# number of processors to use
NP         = 1;

# type of problem to run, possible values: 1 or 2
# if using problem 2, you should probaly have limiting != none
problem_type = 1

# type of limiting to use, possible values: none, hrl, m2l, hri, m2i
limiting   = 'none'

# architecture to use, possible values: USE_CPU, USE_GPU
ARCH       = 'USE_CPU'  


#================================================================================
# 
# Defaults preprocessors
#
#================================================================================
LIEU='HOME'
PROF='OFF'
PARA='NO_MPI'
PREC='USE_DOUBLE'
DIMS='TWOD'
preprocessors=[LIEU,ARCH,PROF,PARA,PREC,DIMS]

#================================================================================
#
# Problem definitions
#
#================================================================================
MESHTYPE='qua'
MESHFILE='square_'+MESHTYPE+'_'+str(resolution)+'_'+str(order)+'_np'+str(NP)+'.msh';
print MESHFILE
JOBNAME='square'
RUNTIME="03:00:00"
if ARCH == 'USE_GPU':
    QUEUE="gpu"
else: 
    QUEUE="normal"

#================================================================================
#
# Directory setup
#
#================================================================================
BASEDIR=os.environ['HOME']+'/dg/'
CODEDIR=BASEDIR+'code/'
BINDIR =BASEDIR+'code/bin/'
MESHDIR=BASEDIR+'profiling/mesh/'
DATADIR=BASEDIR+'profiling/'


#================================================================================
#
# Define the resources to use
#
#================================================================================
RESOURCES='-n '+str(NP)
if ARCH == 'USE_GPU':
    RESOURCES += ' -N '+str(NP)
SUBMIT=[JOBNAME,RUNTIME,QUEUE,RESOURCES]

#================================================================================
#
# Commands that we will need
#
#================================================================================
DG='./dgexec -d deck.inp';
if PARA == 'USE_MPI':
    DG='mpirun '+DG

#================================================================================
#
# Setup and run the code
#
#================================================================================

# problem definition
if problem_type == 1:
    defs=['MULTIFLUID','5','GAMNCON','ROE','0',str(order),'sinegam',limiting,'5','5','0.5',MESHTYPE]
elif problem_type == 2:
    defs=['MULTIFLUID','5','GAMNCON','ROE','0',str(order),'sodcirc',limiting,'0.2','0.2','0.5',MESHTYPE]

# Create directory
WORKDIR=DATADIR+defs[6]+'_'+str(resolution)+'_'+str(order)+'_'+str(NP)+'_'+limiting+'_'+ARCH
print 'Creating directory',WORKDIR;
shutil.rmtree(WORKDIR, ignore_errors=True)
os.makedirs(WORKDIR)

# Write the deck
deckname = write_deck(WORKDIR,preprocessors,defs);

# read the deck to get the make command
makecmd = get_make_command(deckname)

# should I recompile the code
recompile = test_recompile(CODEDIR,recompile,makecmd)

# compile to appropriate executable
success = compile_exec(recompile,CODEDIR,WORKDIR,makecmd)

# Go to the work directory
os.chdir(WORKDIR)
print 'Entering', WORKDIR

# write the submit file in our directory
submitname = write_submit(WORKDIR,SUBMIT,DG);

# copy the executable into our directory
try:
    shutil.copy2(BINDIR+'dgexec', 'dgexec')
except IOError, e:
    print "Unable to copy executable file. %s" % e
    sys.exit()

# copy the mesh into our directory
try:
    shutil.copyfile(MESHDIR+MESHFILE, 'mesh.msh')
except IOError, e:
    print "Unable to copy mesh file. %s" % e
    sys.exit()

# submit the run to the queue
call(DG, shell=True)

# go back to base directory
os.chdir(BASEDIR)
