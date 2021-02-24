#!/usr/bin/env python
# ian.heywood@physics.ox.ac.uk


import numpy
import os
import glob
import pylab
import sys
import scipy.signal
from shutil import copyfile
from astropy.convolution import Gaussian2DKernel
from astropy.io import fits
from astropy.wcs import WCS
from optparse import OptionParser
from scipy import ndimage,stats



def get_image(fitsfile):
    input_hdu = fits.open(fitsfile)[0]
    if len(input_hdu.data.shape) == 2:
            image = numpy.array(input_hdu.data[:,:])
    elif len(input_hdu.data.shape) == 3:
            image = numpy.array(input_hdu.data[0,:,:])
    elif len(input_hdu.data.shape) == 4:
            image = numpy.array(input_hdu.data[0,0,:,:])
    else:
            image = numpy.array(input_hdu.data[0,0,0,:,:])
    return image


def deg2rad(xx):
    return numpy.pi*xx/180.0


def get_beam(fitsfile):
    input_hdu = fits.open(fitsfile)[0]
    hdr = input_hdu.header
    bmaj = hdr.get('BMAJ')
    bmin = hdr.get('BMIN')
    bpa = hdr.get('BPA')
    pixscale = hdr.get('CDELT2')
    return bmaj,bmin,bpa,pixscale


def flush_fits(newimage,fitsfile):
    f = fits.open(fitsfile,mode='update')
    input_hdu = f[0]
    if len(input_hdu.data.shape) == 2:
            input_hdu.data[:,:] = newimage
    elif len(input_hdu.data.shape) == 3:
            input_hdu.data[0,:,:] = newimage
    elif len(input_hdu.data.shape) == 4:
            input_hdu.data[0,0,:,:] = newimage
    else:
            input_hdu.data[0,0,0,:,:] = newimage
    f.flush()


def bb(xx):
    return str(round((xx*3600.0),3))


def get_target_beam(residuals):
    majors = []
    minors = []
    pas = []
    for residual in residuals:
        bmaj,bmin,bpa,pixscale = get_beam(residual)
        majors.append(bmaj)
        minors.append(bmin)
        pas.append(bpa)
    target_bmaj = numpy.max(majors)
    target_bmin = numpy.max(minors)
    target_bpa = numpy.median(numpy.array(pas))
    print('       | Target beam is: '+bb(target_bmaj)+','+bb(target_bmin)+','+str(round(target_bpa,3)))
    return target_bmaj,target_bmin,target_bpa


def drop_deg_axes(fitsfile):
    input_hdu = fits.open(fitsfile,mode='update')
    data = input_hdu[0].data[0,0,:,:]
    input_hdu[0].data = data
    input_hdu.flush()
    hdr = input_hdu[0].header


def fix_beam_header(fitsfile,bmaj,bmin,bpa):
    hdu = fits.open(fitsfile,mode='update')
    hdr = hdu[0].header
    hdr.set('BMAJ',bmaj)
    hdr.set('BMIN',bmin)
    hdr.set('BPA',bpa)
    hdu.flush()  


def pbcor(pbdir,infits,threshold=0.3):
    idx = infits.index('-00')
    chan = infits[idx+1:idx+4]
    pb_fits = glob.glob(pbdir.rstrip('/')+'/*'+chan+'*fits')[0]
    pb_image = get_image(pb_fits)
    mask = pb_image < threshold
    pb_image[mask] = numpy.nan
    in_image = get_image(infits)
    pbcor_image = in_image / pb_image
    return pbcor_image


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

prefix = sys.argv[1]

print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
print('       | '+prefix)

kernel_size = 51
pbdir = '/users/ianh/pbimages'

residuals = sorted(glob.glob(prefix+'*00*-residual.fits'))
print('       | Found '+str(len(residuals))+' residual sub-bands')
target_bmaj,target_bmin,target_bpa = get_target_beam(residuals)

for residual_fits in residuals:
    print('       | '+residual_fits)
    # Get restoring beam
    bmaj,bmin,bpa,pixscale = get_beam(residual_fits)

    # Setup FITS files for beams and kernels
    template_fits = residual_fits.replace('.fits','_template.fits')
    target_beam_fits = template_fits.replace('.fits','_target_beam.fits')
    restoring_beam_fits = template_fits.replace('.fits','_restoring_beam.fits')
    kernel_fits = template_fits.replace('.fits','_kernel.fits')

    restored_fits = residual_fits.replace('residual','image-conv')
    pbcor_fits = restored_fits.replace('.fits','_pbcor.fits')

    print(' ----> | Creating kernel template image')
    os.system('fitstool.py -f -z '+str(kernel_size)+' -o '+template_fits+' '+residuals[0])
    drop_deg_axes(template_fits)

    # Target beam image
    print (' ----> | Creating target beam image')
    copyfile(template_fits,target_beam_fits)
    target_xstd = target_bmin/(2.3548*pixscale)
    target_ystd = target_bmaj/(2.3548*pixscale)
    target_theta = deg2rad(target_bpa)
    target_gaussian = Gaussian2DKernel(x_stddev=target_xstd,y_stddev=target_ystd,theta=target_theta,x_size=kernel_size,y_size=kernel_size,mode='center')
    target_image = target_gaussian.array
    target_image = target_image / numpy.max(target_image)
    flush_fits(target_image,target_beam_fits)

    # Restoring beam image
    print(' ----> | Creating restoring beam image')
    copyfile(template_fits,restoring_beam_fits)
    xstd = bmin/(2.3548*pixscale)
    ystd = bmaj/(2.3548*pixscale)
    theta = deg2rad(bpa)
    restoring_gaussian = Gaussian2DKernel(x_stddev=xstd,y_stddev=ystd,theta=theta,x_size=kernel_size,y_size=kernel_size,mode='center')
    restoring_image = restoring_gaussian.array
    restoring_image = restoring_image / numpy.max(restoring_image)
    flush_fits(restoring_image,restoring_beam_fits)

    # Run PyPHER to get homogenising kernel
    os.system('pypher '+restoring_beam_fits+' '+target_beam_fits+' '+kernel_fits)

    print(' <---- | Reading kernel image')
    kernel_image = get_image(kernel_fits)
    print(' <---- | Reading residual image')
    residual_image = get_image(residual_fits)
    print('  (X)  | Convolving residual with homogenising kernel')
    residual_conv_image = scipy.signal.fftconvolve(residual_image, kernel_image, mode='same')

    print(' <---- | Reading model image')
    model_image = get_image(residual_fits.replace('residual','model'))
    print('  (X)  | Convolving model with target beam')
    model_conv_image = scipy.signal.fftconvolve(model_image, target_image, mode='same')

    print('   +   | Restoring model to residual')
    restored_image = residual_conv_image + model_conv_image
    print(' ----> | Creating restored image')
    copyfile(residual_fits,restored_fits)
    flush_fits(restored_image,restored_fits)
    fix_beam_header(restored_fits,target_bmaj,target_bmin,target_bpa)

    print('   /   | Applying primary beam correction')
    pbcor_image = pbcor(pbdir,restored_image)
    print(' ----> | Creating primary beam corrected image')
    copyfile(restored_fits,pbcor_fits)
    flush_fits(pbcor_image,pbcor_fits)

    print('       | Done')
    print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
