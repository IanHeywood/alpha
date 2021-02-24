#!/usr/bin/env python
# ian.heywood@physics.ox.ac.uk


import glob
import os
import sys


imsize = 6000
prefix = sys.argv[1]


models = sorted(glob.glob(prefix+'*00*-model.fits'))
residuals = sorted(glob.glob(prefix+'*00*-residual.fits'))


for model in models:
	os.system('fitstool.py -z '+str(imsize)+' '+model)

for residual in residuals:
	os.system('fitstool.py -z '+str(imsize)+' '+residual)