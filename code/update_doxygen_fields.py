#!/usr/bin/env python
#
# Loop through all the necessary code files and update the copyright
# and license doxygen information
#
__author__    = 'marchdf'
__copyright__ = "Copyright (C) 2012-2015, Regents of the University of Michigan"
__license__   = "This project is released under the GNU Public License. See LICENSE."
__email__     = "marchdf@umich.edu"
def usage():
    print '\nUsage: {0:s} [options ...]\nOptions:\n  -h, --help\tshow this message and exit\n'.format(sys.argv[0])

#================================================================================
#
# Imports
#
#================================================================================
# Ignore deprecation warnings
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import sys, getopt, os, re


#================================================================================
#
# Parse arguments
#
#================================================================================
try:
    myopts, args = getopt.getopt(sys.argv[1:],"h",["help"])
except getopt.GetoptError as e:
    print (str(e))
    usage()
    sys.exit(2)

for o, a in myopts:
    if o in ('-h', '--help'):
        usage()
        sys.exit()

#================================================================================
#
# Some defaults variables
#
#================================================================================
file_suffixes = ['.h','.cc','.cu','.mk']
ignored_directories = ['doc','matlab_functions','dependencies','quadratures','objects']
new_copyright_line = 'copyright Copyright (C) 2012-2015, Regents of the University of Michigan\n'
new_license_line = 'license This project is released under the GNU Public License. See LICENSE.\n'

#================================================================================
#
# Define some functions
#
#================================================================================


#================================================================================
#
# Basic information/setup
#
#================================================================================
cwd = os.getcwd()


#================================================================================
#
# Get all the files which end in one of the suffixes
#
#================================================================================
pp_files = []

# traverse the directory
for root, dirs, files in os.walk(cwd):

    # Make sure to ignore some directories
    if not any( fdir in root for fdir in ignored_directories):

        # loop on all the files
        for file in files:

            # loop on all the suffixes
            for suffix in file_suffixes:

                # check to see if the file ends with the suffix
                if file.endswith(suffix):
                    pp_files.append(os.path.join(root, file))


#================================================================================
#
# For each of these files, update the copyright and license line
#
#================================================================================
for pp_file in pp_files:
    # read all the lines in the file
    with open(pp_file, "r") as sources:
        lines = sources.readlines()

    # write all the lines back to the file but update the copyright line
    with open(pp_file, "w") as sources:
        for line in lines:
                
            # if this is the copyright line, replace it with the new one
            if '\\copyright' in line and 'Gmsh' not in line:
                prefix = re.split('copyright',line)[0]
                line = prefix + new_copyright_line


            # if this is the license line, replace it with the new one
            if '\\license' in line:
                prefix = re.split('license',line)[0]
                line = prefix + new_license_line

            # write the line to the file
            sources.write(line)


