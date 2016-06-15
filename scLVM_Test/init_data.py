import sys
import scipy as SP
import pylab as PL
from matplotlib import cm
import h5py
#make sure your paths point to limix and scLVM directories
limix_path = '/Users/simonsteiger/anaconda/pkgs/limix-0.8.0-py27_0/lib/python2.7/site-packages/limix'
sclvm_path = '~/RNA_Seq/scLVM_Test'
sys.path.append(limix_path)
sys.path.append(sclvm_path)
#import scLVM
sys.path.append('./../scLVM')
from scLVM import scLVM
