#!/usr/bin/env python2
#
__author__ = 'marchdf'


def usage():
    print '\nUsage: {0:s} [options ...]\nOptions:\n -f, --force\tforce recompile\n -h, --help\tshow this message and exit\n'.format(sys.argv[0])


#=========================================================================
#
# Imports
#
#=========================================================================
import sys
import getopt
import os
import time
import shutil
from numpy import *
from subprocess import call
import subprocess as sp
import time

#=========================================================================
#
# Functions
#
#=========================================================================


def write_deck(WORKDIR, preprocessors, defs):
    deckname = WORKDIR + '/deck.inp'
    with open(deckname, "w") as f:
        f.write('#Make options\n')
        f.write('LIEU=' + preprocessors[0] + '\n')
        f.write('ARCH=' + preprocessors[1] + '\n')
        f.write('PROF=' + preprocessors[2] + '\n')
        f.write('PARA=' + preprocessors[3] + '\n')
        f.write('PREC=' + preprocessors[4] + '\n')
        f.write('DIMS=' + preprocessors[5] + '\n')
	f.write('ERRR=' + preprocessors[6] + '\n') #PEJ 07/05/2017
	f.write('VISC=' + preprocessors[7] + '\n') #PEJ 08/17/2017
	f.write('FSOL=' + preprocessors[8] + '\n') #PEJ 08/22/2017
	f.write('CONS=' + preprocessors[9] + '\n') #PEJ 10/01/2017
	f.write('TGVS=' + preprocessors[10] + '\n') #PEJ 10/01/2017
        f.write('PROB=' + defs[0] + '\n')
        f.write('NFLD=' + defs[1] + '\n')
        f.write('GAMM=' + defs[2] + '\n')
        f.write('FLUX=' + defs[3] + '\n')
        f.write('MASS=' + defs[4] + '\n')
        f.write('#Code options\n')
        f.write('#time integration\nrk4\n')
        f.write('#output time step size\n' + defs[8] + '\n')
        f.write('#final time\n' + defs[9] + '\n')
        f.write('#Courant-Friedrichs-Lewy condition\n' + defs[10] + '\n')
        f.write('#Von Neumann condition\n' + defs[11] + '\n') #PEJ 06/01/2017
	f.write('#sensor thresholds\n' + defs[12] + '\n') #PEJ 10/11/2017
        f.write('#order\n' + defs[5] + '\n')
        f.write('#mesh file\nmesh.msh\nhex\n') #need to specify element type here
        f.write('#limiter\n' + defs[7] + '\n')
        f.write('#initial condition\n' + defs[6] + '\n')
    return deckname
        

#=========================================================================


def add_lagrange_deck(deckname, lagrange_particles):
    with open(deckname, "a") as f:
        f.write('#lagrange particles\n' + lagrange_particles + '\n')

#=========================================================================


def add_sensor_deck(deckname, sensor_thresholds):
    with open(deckname, "a") as f:
        f.write('#sensor thresholds\n' + sensor_thresholds + '\n')

#=========================================================================


def get_make_command(deckname):
    makecmd = 'make '
    with open(deckname) as f:
        for line in f:
            if '#Code options' in line:
                return makecmd
            if '#Make options' not in line:
                makecmd = makecmd + line.rstrip('\n') + ' '


#=========================================================================


def test_recompile(CODEDIR, recompile, makecmd):

    if recompile:
        return True
    else:
        try:
            with open(CODEDIR + 'last_makecmd') as f:
                last_makecmd = f.readline()
        except:  # if the last make command file does not exist, recompile
            return True
        if makecmd == last_makecmd:
            return False
        else:
            return True


#=========================================================================


def compile_exec(recompile, CODEDIR, WORKDIR, makecmd):
    if recompile:
        os.chdir(CODEDIR)
        return_code = call('make clean', shell=True)
        return_code = call(makecmd, shell=True)
        if return_code != 0:
            print '\nFailure at compilation\n'
            shutil.rmtree(WORKDIR, ignore_errors=True)
            sys.exit()
    print 'Compiled successfully!'
    with open(CODEDIR + 'last_makecmd', "w") as f:
        f.write(makecmd)
    return True

#=========================================================================


def write_submit(WORKDIR, SUBMIT, CMD):
    submitname = WORKDIR + '/jobscript.lsf'
    with open(submitname, "w") as f:
        f.write('#!/bin/bash\n')
        f.write('#\n')
        f.write('#BSUB -a openmpi\n')
        f.write('#BSUB -c 35000000\n')
        f.write('#BSUB -n 72\n') #maybe should be number of processors. Not sure yet
        f.write('#BSUB -R "span[ptile=18]"\n')
#        f.write('#BSUB -R rusage[mem=48000]\n')
 #       f.write('#BSUB -M 750\n')
        f.write('#BSUB -J Job_p2hex_64x85\n')
        f.write('#BSUB -e errors.p2hex_64x85\n')
        f.write('#BSUB -o output.p2hex_64x85\n')
        f.write('#BSUB -q gpu_p100\n')
        f.write('#BSUB -B\n')
        f.write('#BSUB -N\n')
        f.write('#BSUB -u philzconfluxrelay@gmail.com\n')
        f.write('#\n')
        f.write(CMD + '\n')
    return submitname

#=========================================================================


def execute_dg_code(DG, fname=None):
    # Run the DG code on the local machine

    print 'Running simulation'
    if fname:
        f = open(fname, "w")
        proc = sp.Popen(DG, shell=True, stdout=f,
                        stderr=sp.PIPE)  # send output to file
    else:
        # send output to screen
        proc = sp.Popen(DG, shell=True, stderr=sp.PIPE)

    # Wait for it to end and get the output
    out, err = proc.communicate()
    errcode = proc.returncode

    if errcode == 0:
        print 'Success running the executable'
    else:
        print 'Failed to run the executable'
        print err

#=========================================================================
#
# Default arguments
#
#=========================================================================
recompile = False

#=========================================================================
#
# Parse arguments
#
#=========================================================================
try:
    myopts, args = getopt.getopt(sys.argv[1:], "hf", ["help", "force"])
except getopt.GetoptError as e:
    print(str(e))
    usage()
    sys.exit(2)

for o, a in myopts:
    if o in ('-h', '--help'):
        usage()
        sys.exit()
    elif o in ('-f', '--force'):
        recompile = True


#=========================================================================
#
# Defaults preprocessors
#
#=========================================================================
LIEU = 'CONFLUX'
ARCH = 'USE_CPU'
PROF = 'OFF'
PARA = 'USE_MPI'
PREC = 'USE_DOUBLE'
DIMS = 'THREED'
ERRR = 'NOERROR'
VISC = 'CONSTANTVIS'
FSOL = 'ICBN'
CONS = 'NULL'
TGVS = 'TGVSTATS'
preprocessors = [LIEU, ARCH, PROF, PARA, PREC, DIMS, ERRR, VISC, FSOL, CONS, TGVS]

#=========================================================================
#
# Problem definitions (most of these are just for the cluster)
#
#=========================================================================
RUNTIME = "20:00:00"
NP = '64'
JOBNAME = 'example'
MESHPFX = 'sample_mesh'
MESHFILE = MESHPFX + '_1.msh'
if ARCH == 'USE_GPU':
    QUEUE = "gpu"
else:
    QUEUE = "normal"

#=========================================================================
#
# Directories
#
#=========================================================================
BASEDIR = os.getcwd()
BASEDIR = BASEDIR.split("/code/")[0] + '/'
CODEDIR = BASEDIR + 'code/'
MESHDIR = os.getcwd() + '/mesh/'
DATADIR = os.getcwd() + '/'

#=========================================================================
#
# Define the resources to use (also just for the cluster, mostly)
#
#=========================================================================
RESOURCES = '-n ' + str(NP)
if ARCH == 'USE_GPU':
    RESOURCES += ' -N ' + str(NP)
SUBMIT = [JOBNAME, RUNTIME, QUEUE, RESOURCES]

#=========================================================================
#
# Commands that we will need
#
#=========================================================================
DG = './dgexec -d deck.inp'
if PARA == 'USE_MPI':
#    DG = 'mpirun ' + DG
#    DG = 'mpirun -np 3 ' + DG	    	
    DG = 'mpirun -np ' + NP + ' ' + DG

# Define the problem
PDG="2"
#Syntax:[[0]physics, [1]number of fields, [2]eos, [3]Riemann, [4]mass fraction, [5]p, [6]IC, [7]limiter, [8]output increment, [9]t_final, [10]CFL, [11]VNN, [12]sensor threshholds]
# Define the problem

defs = ['SINGLEFLUID', '5', 'GAMNCON', 'SLAU', '0', PDG,
        'tgrvrtx', 'null', '0.5', '10.5', '0.3','1.4','0.1 0.1 0.1']
#running to one extra output increment to make sure I at least make it to t=6.25

# Create directory
WORKDIR = DATADIR + defs[6].replace(" ", "_")
print 'Creating directory', WORKDIR
shutil.rmtree(WORKDIR, ignore_errors=True)
os.makedirs(WORKDIR)

# Write the deck
deckname = write_deck(WORKDIR, preprocessors, defs)

# read the deck to get the make command
makecmd = get_make_command(deckname)

# should I recompile the code
recompile = test_recompile(CODEDIR, recompile, makecmd)

# compile to appropriate executable
success = compile_exec(recompile, CODEDIR, WORKDIR, makecmd)

# Go to the work directory
os.chdir(WORKDIR)
print 'Entering', WORKDIR

# write the submit file in our directory
submitname = write_submit(WORKDIR, SUBMIT, DG)

# copy the executable into our directory
try:
    shutil.copy2(CODEDIR + 'bin/dgexec', 'dgexec')
except IOError, e:
    print "Unable to copy executable file. %s" % e
    sys.exit()

# copy the mesh into our directory
try:
    shutil.copyfile(MESHDIR + MESHFILE, 'mesh.msh')
except IOError, e:
    print "Unable to copy mesh file. %s" % e
    sys.exit()

# Run the code
#execute_dg_code(DG) //Marc's original run command
#os.chdir(WORKDIR)
print "Now run the jobscript with {bsub <jobscript.lsf} in the working directory"

# go back to base directory
os.chdir(BASEDIR)