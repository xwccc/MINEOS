from glob import glob
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('Agg')
# read fre files
fmax = 50
lmax = 545
filelist = glob('/home/xuwc/global_seismology/MINEOS/tmp/new_sph_00*.fre')
#filelist = glob('/home/2021gs/global_seismology_2021/test_xwc/tmp/premANIC_sph_00*.fre')
#filelist = glob('/home/2021gs/global_seismology_2021/pengjx/Do_not_modify_anything_in_this_folder/MINEOS/tm1/premANIC.card_tor_00*.fre')
filelist.sort()

nlf = np.array([])

for file in filelist:
    with open(file,'r') as f:
        lines = f.readlines()
    n = np.array([int(line.split()[0]) for line in lines])
    l = np.array([int(line.split()[2]) for line in lines])
    f = np.array([float(line.split()[4]) for line in lines])
    nlf_tmp = np.vstack((n,l,f))
    if len(nlf) == 0:
        nlf = nlf_tmp
    else:
        nlf = np.hstack((nlf,nlf_tmp))
  
# sort by n

nlf = nlf[:,np.argsort(nlf[0,:])]

# plot dispersion curve

nmax = int(nlf[0].max())
plt.figure(figsize=[16,9])
for n in range(nmax):
    # sort by l
    lf = nlf[1:,nlf[0]==n]
    lf = lf[:,np.argsort(lf[0,:])]
    plt.plot(lf[0],lf[1],'ko-',linewidth=0.4,markersize=0.6)
plt.xlabel('l')
plt.ylabel('f (mHz)')
plt.title('Toridal Mode Dispersion Curve')
plt.ylim([0,fmax])
plt.xlim([0,lmax+10])
plt.savefig('Disp_curve.png')
plt.savefig('Disp_curve.pdf')
plt.close()
