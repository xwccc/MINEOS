# MINEOS
Fortran code for normal modes calculation and mode summation

This is the non-parallel version of MINEOS. The original MINEOS program, mineos.f, is the 
product of many people's efforts, including Freeman Gilbert, John Woodhouse and Guy Masters. 

Although 98% similar, there are perhaps as many versions of mineos.f as there are users. 
This version was originally obtained in 1991 by Jeroen Tromp at Princeton from Jeffery Park 
of Yale. Parts of it have been modified later by Li Zhao and Hsin-Ying Yang. 

The wrapper mineos_drv.f was originally written by Dimitri Komatitsch, and modified 
later by Li Zhao and Hsin-Ying Yang. 

A few programs were written by Li Zhao to facilitate checking for missing and problematic 
modes, and for the purpose of calculating synthetic seismograms by normal-mode summation.
They include:

checkmode.f
mineos_merge.f
mineos_strip.f
modesum7_nonmpi.f

Other utilities:

1. read_mineos.f can be used to read the .fun files from mineos.f for plotting.
2. eigenfun.gmt can be used to plot the eigenfunctions read by read_mineos.f

To compile: 

All programs can be compiled by Intel compiler with the command

ifort -132 -O4 -save -zero -o name_of_executable name_of_program.f

The include file mineos.inc is needed to compile mineos.f.

Procedures for computing normal modes and synthetics: 

1. Prepare the model file (e.g. ak135.card).
2. Edit run.mineos.drv to set running parameters for mineos (l and f limits, etc.)
3. Run run.mineos.drv to compute the modes and check for missing and problematic ones. 
4. For each file in CHECKMODE_OUTPUT_FILES and CHECKMODE_FREQERR_FILES: 
   4.1: Run mineos manually to find missing modes or improve the accuracy. 
        For each run, the corresponding input file in MINEOS_INPUT_FILES must be edited.
        NOTE: make sure to use different output file names for .fre and .fun in 
        MINEOS_INPUT_FILES so that results from 4.1 do not over-write those from Step 3.
   4.2: Run mineos_merge to merge the multiple sets of .fre and .fun from Step 3 and 
        Step 4.1 into a single set of .fre and .fun files.
5. Run mineos_strip for toroidal, spheroidal and radial modes separately to generate the 
   mode table files used for calculating synthetics by mode summation. 
   For each type of mode, an input file is needed (e.g. mineos_strip_tor.inp).
6. Run modesum7_nonmpi to compute synthetics. 
   This run can calculate multiple synthetics. Each synthetic needs an input file 
   (e.g. 20070228_231315_BOSA.inp). All the input files are listed in input_file_list. 
   modesum7_nonmpi will loop over all input files in input_file_list. 
