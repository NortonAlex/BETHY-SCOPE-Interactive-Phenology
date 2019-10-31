#! /usr/bin/python

import sys
import csv
import os

## Read me:
# This script produces the control files for block parallel runs. Ensure the output 
# directories specified in these control files actually exist, otherwise errors occur.
# Example run from the command line:
#
# python write_block_control_files.py 20 control_script_default_hires.nml control_script_block
#
## Command line arguments:
## first argument = nsplits (i.e. nblocks) 
## second argument = file name of default control file for which to duplicate
## third argument = prefix of block control file names

## Change to working directory where input files are located
wd = '/home/563/ajn563/short/fluoro2/control_blockruns/'
os.chdir(wd)

# SET A COUPLE OF THINGS
# Command line arguments 

if len(sys.argv) > 1:
   for arg in sys.argv:
      print 'Command line argument:', arg
else:
   print 'Require 3 inputs as command line argument... exiting'
   sys.exit()

# SET SOME VARIABLES FROM COMMAND LINE ARGUMENTS

nblocks = int(sys.argv[1])
control_default = str(sys.argv[2])    # name of control file to duplicate 
blockfilename = str(sys.argv[3])   # start of control file name

# IMPORT DEFAULT CONTROL FILE

results = []
with open('%s' % control_default) as inputfile:
    for row in csv.reader(inputfile):
        results.append(row)


# WRITE OUT BLOCK CONTROL FILES

filelist = []

for i in range(nblocks):
    filelist.append('%s_%02d.nml' % (blockfilename,i+1))

for i in range(nblocks):
    ifile = '%s' % filelist[i]
    print 'Index', i, '  Writing File:',ifile
    ifile = open('%s' % ifile,'w')
    iline = 0
    for item in results:
        print 'Line:',iline+1,'  item:',item
        if ((iline == 0) | (item[0][0] == "/")):
            # We skip any lines without the "=" char in them (i.e. 1st and last lines)
            ifile.write("%s\n" % item[0])
            iline += 1
            continue
        ind = item[0].index("=")
        if item[0][:ind].strip() == "outdir":
            # Outdir
            outdir = " outdir       = './output/job%02d/'" % (i+1) 
            ifile.write("%s\n" % outdir)
        elif item[0][:ind].strip() == "nblocks":
            # number of blocks to parallelise over (i.e. number of servers)
            string = " nblocks      = %02d" % nblocks
            ifile.write("%s\n" % string)
        elif item[0][:ind].strip() == "iblock":
            # which block is this
            string = " iblock       = %d" % (i+1)
            ifile.write("%s\n" % string)
        else:
            ifile.write("%s\n" % item[0])
        iline += 1
    ifile.close()



