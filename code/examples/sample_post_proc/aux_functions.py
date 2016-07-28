#================================================================================
def write_submit(WORKDIR,SUBMIT,CMD):
    submitname = WORKDIR+'/pp_submit.batch'
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
        f.write('#  Put your job commands after this line\n')
        f.write(CMD+'\n')
    return submitname

#================================================================================
def replace_in_file(fname,replacements):
    # Replace using a dictionary of replacements (={old_string:new_string} in a file
    lines = []
    with open(fname) as fin:
        for line in fin:
            for old_string, new_string in replacements.iteritems():
                line = line.replace(old_string, new_string)
            lines.append(line)
    with open(fname, 'w') as fout:
        for line in lines:
            fout.write(line)
