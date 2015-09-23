#!/usr/bin/env python
#
__author__ = 'marchdf, modified by: Brandon Patterson'
def usage():
    print '\nUsage: {0:s} [options ...]\nOptions:\n -f, --force\tforce recompile\n -h, --help\tshow this message and exit\n'.format(sys.argv[0]) 


#
# Imports
#
import sys,getopt,os,time,shutil
from numpy import *
from subprocess import call
import subprocess as sp
import time

#
# Functions
#
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
        f.write('#mesh file\nmesh.msh\nqua\n')
        f.write('#limiter\n'+defs[7]+'\n')
        f.write('#initial condition\n'+defs[6]+'\n')
    return deckname

def add_lagrange_deck(deckname,lagrange_particles):
    with open(deckname, "a") as f:
        f.write('#lagrange particles\n'+lagrange_particles+'\n')

def add_sensor_deck(deckname,sensor_thresholds):
    with open(deckname, "a") as f:
        f.write('#sensor thresholds\n'+sensor_thresholds+'\n')

def get_make_command(deckname):
    makecmd='make '
    with open(deckname) as f: 
        for line in f:
            if '#Code options' in line:
                return makecmd
            if '#Make options' not in line:
                makecmd = makecmd + line.rstrip('\n') + ' '

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
        f.write('#SBATCH --mail-user=awesome@umich.edu\n')
        f.write('#SBATCH --mail-type=begin  # email me when the job starts\n')
        f.write('#SBATCH --mail-type=end    # email me when the job finishes\n\n')
        f.write('module load cuda\n\n')
        f.write('#  Put your job commands after this line\n')
        f.write(CMD+'\n')
    return submitname


# def execute_dg_code(DG,fname=None):
#     # Run the DG code

#     print 'Running simulation'
#     if fname: 
#         f = open(fname, "w")
#         proc = sp.Popen(DG, shell=True, stdout=f,stderr=sp.PIPE) # send output to file
#     else:
#         proc = sp.Popen(DG, shell=True, stderr=sp.PIPE) # send output to screen

#     # Wait for it to end and get the output
#     out,err = proc.communicate()
#     errcode = proc.returncode

#     if errcode == 0:
#         print 'Success running the executable'
#     else:
#         print 'Failed to run the executable'
#         print err


#
# Default arguments
#
recompile = False

#
# Parse arguments
#
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
        


#
# Length, strength and mach numbers
#
pa_ = [1e7]

# 
# Defaults preprocessors
#
LIEU='TACC'
ARCH='USE_GPU'
PROF='OFF'
PARA='NO_MPI'
PREC='USE_DOUBLE'
DIMS='TWOD'
preprocessors=[LIEU,ARCH,PROF,PARA,PREC,DIMS]

#
#Problem definitions
#
RUNTIME="00:20:00"
NP=1;
RES='25'
JOBNAME='rmawave_'+RES
MESHPFX='rmawave_'+RES+'_';
MESHFILE=MESHPFX+'1.msh'
if ARCH == 'USE_GPU':
    QUEUE="gpu"
else: 
    QUEUE="normal"

#
# Directories
#
BASEDIR='/home1/03766/tg830659/dg/'
CODEDIR=BASEDIR+'code/'
BINDIR =BASEDIR+'code/bin/'
MESHDIR=BASEDIR+'mesh/'
DATADIR=os.environ['SCRATCH']+'/'

#
# Define the resources to use
#
RESOURCES='-n '+str(NP)
if ARCH == 'USE_GPU':
    RESOURCES += ' -N '+str(NP)
SUBMIT=[JOBNAME,RUNTIME,QUEUE,RESOURCES]

#
# Commands that we will need
#
DG='./dgexec -d deck.inp';
if PARA == 'USE_MPI':
    DG='ibrun '+DG

# counters
num2submit = len(pa_)
cnt = 1
print 'Total number of simulations to be submitted = {0:.0f}'.format(num2submit)
    
for pa in pa_:
                
    print '\nSubmitting simulation [{0:.0f}/{1:.0f}]'.format(cnt,num2submit); cnt += 1;
                
    # Define the problem
    # defs=['STIFFENED','6=number_of_unknowns','GAMNCON'=nonconservative_gamma_model,'RUS'=Rusinov_flux_solver,'0'=number_of_mass_fracs,'1'=p(dg_order),'rmawave '+str(pa)+' 0.03 5.0 0.0'=(IC_ string_with_name_and_inputs),'m2l'=marcs_limiter,'0.1'=output_timestep,'4'=final_time,'0.5'=cfl_number]
    defs=['STIFFENED','7','GAMNCON','ROE','1','1','rmawave '+str(pa)+' 0.03 1.0 0.0','m2l','0.1','50','0.5']

    # Create directory
    WORKDIR=DATADIR+defs[6].replace(" ","_")+'_'+RES
    print 'Creating directory',WORKDIR;
    shutil.rmtree(WORKDIR, ignore_errors=True)
    os.makedirs(WORKDIR)
                
    # Write the deck
    deckname = write_deck(WORKDIR,preprocessors,defs);

    # Add sensors for runs with discontinuities
    #sensor_thresholds = '0.01 0.001 0.01'
    #add_sensor_deck(deckname,sensor_thresholds);

    # Add lagrange particles
    lagrange_particles=' '+str(int(RES)+1)+' '
    for ii in range(0, (int(RES)+1)):
        xx = 0.03*(sin(float(ii)/float(RES)*2*pi-pi/2)+1)
        yy = float(ii)/float(RES)
        part0 = str(xx)
        part1 = str(yy)    
        lagrange_particles = lagrange_particles+part0+' '+part1+' '
    print lagrange_particles

    add_lagrange_deck(deckname,lagrange_particles);

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

    # Run the code
    call('sbatch submit.batch', shell=True)

    # go back to base directory
    os.chdir(BASEDIR)
