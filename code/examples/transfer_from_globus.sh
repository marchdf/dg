#!/bin/bash
#
# Transfer data.tar.gz from stampede to a local computer using Globus
#

dirs=('bblwedg_10_2' 'bblwedg_10_3')
fnames=('data.tar.gz')

for ((i=0; i<${#dirs[@]}; i++)); do #C style loop with the array length

    # make source and target directory names
    sname='xsede#stampede'
    tname='u_j3rlh7guyyi6lfghlpe4fmrlgm#066d2bae-d4dd-11e5-975d-22000b9da45e'

    # Source directory
    sdir='/~/scratch/'${dirs[i]}

    # Target directory
    tdir='/mnt/datadrive/marchdf/hpcdata/bblwedg/20160711/'${dirs[i]}

    # output names and make the necessary target directory
    echo 'Transfering from' $sname 'to' $tname
    echo '    source folder is' $sdir 'and target folder is' $tdir
    mkdir -p $tdir

    # loop on the files to be transfered
    for ((k=0; k<${#fnames[@]}; k++)); do #C style loop with the array length
	echo '        transfer file' ${fnames[k]}
	ssh marchdf@cli.globusonline.org transfer -- $sname$sdir$'/'${fnames[k]} $tname$tdir$'/'${fnames[k]} 
    done
done
