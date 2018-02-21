#!/usr/bin/env python
#
# Run all the tests in this directory.
#
__author__ = 'marchdf'


def usage():
    print('\nUsage: {0:s} [options ...]\nOptions:\n -f, --force\tforce recompile\n -h, --help\tshow this message and exit\n'.format(sys.argv[0]))


#=========================================================================
#
# Imports
#
#=========================================================================
# Ignore deprecation warnings
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import sys
import getopt
import os
import time
import shutil
from numpy import *
import subprocess as sp

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))
from scripts.mytermcolor import colored

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
# Parse command line arguments using getopt (argparse is better but require python 1.7 at least)
# from:
# http://www.cyberciti.biz/faq/python-command-line-arguments-argv-example/
try:
    myopts, args = getopt.getopt(sys.argv[1:], "hf", ["help", "force"])
except getopt.GetoptError as e:
    print((str(e)))
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
# Some defaults variables
#
#=========================================================================
cmap_med = ['#F15A60', '#7AC36A', '#5A9BD4', '#FAA75B',
            '#9E67AB', '#CE7058', '#D77FB4', '#737373']
cmap = ['#EE2E2F', '#008C48', '#185AA9', '#F47D23',
        '#662C91', '#A21D21', '#B43894', '#010202']
successcolor = cmap[1]
errcolor = cmap[0]

#=========================================================================
#
# Define some functions
#
#=========================================================================


def get_make_command(deckname):
    makecmd = 'make '
    with open(deckname) as f:
        for line in f:
            if '#Code options' in line:
                return makecmd
            if '#Make options' not in line:
                makecmd = makecmd + line.rstrip('\n') + ' '


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


def compile_exec(recompile, CODEDIR, CWD, makecmd):
    if recompile:
        os.chdir(CODEDIR)
        return_code = sp.call('make clean', shell=True)
        return_code = sp.call(makecmd, shell=True)
        os.chdir(CWD)
        if return_code != 0:
            print('\nFailure at compilation\n')
            sys.exit()
    with open(CODEDIR + 'last_makecmd', "w") as f:
        f.write(makecmd)
    return True


def execute_dg_code(DG):
    # Run the DG code but do so only if there are less than X active
    # simulations running. R

    # Launch the process. Pipe the error (but not the output since we
    # want it to go directly to a file)
    proc = sp.Popen(DG + ' > out', shell=True, stderr=sp.PIPE)

    # Wait for it to end and get the output
    err = proc.communicate()
    errcode = proc.returncode

    if errcode == 0:
        print(colored('TEST SUCCESS\n', successcolor))
        return 0
    else:
        print(colored('TEST FAIL', errcolor))
        print(colored('     with error code:' + str(errcode), errcolor))
        print(colored('     with error:' + err[1], errcolor))
        return 1


def get_immediate_subdirectories(a_dir):
    # From
    # http://stackoverflow.com/questions/800197/get-all-of-the-immediate-subdirectories-in-python
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]


def clean_data(dir, patterns):
    # remove all files in directory ending with any pattern in patterns
    for f in os.listdir(dir):
        if f.endswith(tuple(patterns)):
            os.remove(os.path.join(dir, f))


#=========================================================================
#
# Basic information/setup
#
#=========================================================================
CODEDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) + '/'
BINDIR = CODEDIR + 'bin/'
BINNAME = 'dgexec'
TESTDIR = CODEDIR + 'tests/'
subdirs = get_immediate_subdirectories(TESTDIR)
subdirs.sort()


# submit counters
num2test = len(subdirs)
print('Total number of simulations to be submitted = {0:.0f}'.format(num2test))
fails = 0

#=========================================================================
#
# Loop on the tests directories
#
#=========================================================================
for cnt, subdir in enumerate(subdirs):

    this_test_dir = TESTDIR + subdir
    print('=' * 80, '\n')
    print(' Running test [{0:.0f}/{1:.0f}]'.format(cnt +
                                                   1, num2test), 'in directory', this_test_dir, '\n')
    print('=' * 80)

    # go to the specific test directory
    os.chdir(this_test_dir)

    # clean the directory of old data files
    clean_data(this_test_dir, ['.pos', '.dat', 'out'])

    # get the name of the deck
    deckname = [deckname for deckname in os.listdir(
        this_test_dir) if deckname.endswith('.inp')]
    deckname = deckname[0]

    # read the deck to get the make command
    print('\t- Reading the make command')
    makecmd = get_make_command(deckname)

    # Executable commands
    DG = './' + BINNAME + ' -d ' + deckname

    # ammend the exec command if using MPI
    if 'USE_MPI' in makecmd:
        DG = 'mpirun ' + DG

    # should I recompile the code
    print('\t- Test if recompile is necessary')
    recompile = test_recompile(CODEDIR, recompile, makecmd)

    # compile to appropriate executable
    print('\t- Start compile')
    success = compile_exec(recompile, CODEDIR, os.getcwd(), makecmd)
    if success:
        print('\t- Compiled successfully!')

    # copy the executable into our directory
    try:
        shutil.copy2(BINDIR + BINNAME, BINNAME)
    except IOError as e:
        print("Unable to copy executable file. %s" % e)
        sys.exit()

    # Run the code
    print('\t- Running test')
    fails += execute_dg_code(DG)

    # go back to main test directory
    os.chdir(TESTDIR)


print(colored('Number of succesful tests = {0:2.0f}'.format(
    num2test - fails), successcolor))
print(colored('Number of failed tests = {0:2.0f}'.format(fails), errcolor))
