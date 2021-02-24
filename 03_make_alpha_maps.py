import os
import sys
import numpy
#import pylab
import glob
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from scipy.optimize import curve_fit


def fitfunc(x, A, B):
    return A*x + B


def getfreqs(fitslist):
	freqs = []
	for fitsfile in fitslist:
		hdu = fits.open(fitsfile)
		hdr = hdu[0].header
		freq = hdr['CRVAL3']
		freqs.append(freq)
	freqs = numpy.array(freqs)
	return freqs


def getImage(fitsfile):
        input_hdu = fits.open(fitsfile)[0]
        if len(input_hdu.data.shape) == 2:
                image = numpy.array(input_hdu.data[:,:])
        elif len(input_hdu.data.shape) == 3:
                image = numpy.array(input_hdu.data[0,:,:])
        else:
                image = numpy.array(input_hdu.data[0,0,:,:])
        return image


def flushFits(newimage,fitsfile):
	# Write numpy array newimage to fitsfile
	# Dimensions must match (obv)
	f = fits.open(fitsfile,mode='update')
	input_hdu = f[0]
	if len(input_hdu.data.shape) == 2:
		input_hdu.data[:,:] = newimage
	elif len(input_hdu.data.shape) == 3:
		input_hdu.data[0,:,:] = newimage
	else:
		input_hdu.data[0,0,:,:] = newimage
	f.flush()


def makecube(fitslist):
	temp = []
	for fitsfile in fitslist:
		img = getImage(fitsfile)
		temp.append(img)
	cube = numpy.dstack(temp)
	return cube


prefix = sys.argv[1]
threshold = 1.0e-3

fitslist = sorted(glob.glob(prefix+'*00*pbcor.fits'))
freqlist = sorted(glob.glob(prefix+'*00*pbcor.fits'))[:-1]

alphafits = prefix+'_alpha_'+str(threshold).replace('.','p')+'.fits'
alphaerrorfits = prefix+'_alphaerror_'+str(threshold).replace('.','p')+'.fits'
f0fits = prefix+'_f0_'+str(threshold).replace('.','p')+'.fits'
nfluxfits = prefix+'_nflux_'+str(threshold).replace('.','p')+'.fits'
overwrite = False

os.system('cp '+fitslist[0]+' '+alphafits)
os.system('cp '+fitslist[0]+' '+alphaerrorfits)
os.system('cp '+fitslist[0]+' '+f0fits)
os.system('cp '+fitslist[0]+' '+nfluxfits)

freqs = getfreqs(freqlist)
logf = numpy.log10(freqs)

print(freqs)
print(logf)

cube = makecube(fitslist)
mask = cube > threshold
cube[~mask] = numpy.nan

flattened_mask = numpy.mean(mask,axis=2)

skymask = flattened_mask > 0.5

idx = numpy.column_stack(numpy.where(skymask==True))

alphaimage = numpy.empty((cube.shape[0],cube.shape[1]))
alphaerrorimage = numpy.empty((cube.shape[0],cube.shape[1]))
f0image = numpy.empty((cube.shape[0],cube.shape[1]))
nfluximage = numpy.empty((cube.shape[0],cube.shape[1]))

alphaimage[:] = numpy.nan
alphaerrorimage[:] = numpy.nan
f0image[:] = numpy.nan
nfluximage[:] = numpy.nan

nspec = len(idx)
tenpc = int(nspec/10.0)
count = 0
pcount = 0

for i,j in idx:
	if count == 0:
		print(str(pcount)+'%...')
	spec = numpy.array(cube[i,j,:])
	if numpy.isnan(spec.sum()):
		specmask = ~numpy.isnan(spec)
	else:
		specmask = [True]*len(freqs)
	num_points = len(spec[specmask])
	f0 = numpy.mean(freqs[specmask])
	logspec = numpy.log10(spec[specmask])
	popt,pcov = curve_fit(fitfunc,logf[specmask],logspec)
	A = popt[0]
        alpha_err = numpy.sqrt(numpy.diag(pcov))[0]
	alphaimage[i,j] = A
	alphaerrorimage[i,j] = alpha_err
	f0image[i,j] = f0
	nfluximage[i,j] = num_points
	count+=1
	if count == tenpc:
		count = 0
		pcount += 10


flushFits(alphaimage,alphafits)
flushFits(alphaerrorimage,alphaerrorfits)
flushFits(f0image,f0fits)
flushFits(nfluximage,nfluxfits)
