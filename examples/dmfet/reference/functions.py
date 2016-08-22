#!/usr/bin/python2.7
from math import *
from numpy import reshape
from numpy import fft
from numpy import sum
from numpy import meshgrid
from scipy import dot
from scipy import array
from scipy import average
from scipy import optimize
from scipy import random

lambd = 0.0e-10

def laplacian3d(field, nx, ny, nz):
    kx = fft.fftfreq(nx,1.0)*2*pi
    ky = fft.fftfreq(ny,1.0)*2*pi
    kz = fft.fftfreq(nz,1.0)*2*pi
    KX, KY, KZ = meshgrid(kx, ky, kz, indexing='ij')
    return fft.ifft(-KX**2*fft.fft(field, axis = 0), axis = 0).real + \
        fft.ifft(-KY**2*fft.fft(field, axis = 1), axis = 1).real + \
        fft.ifft(-KZ**2*fft.fft(field, axis = 2), axis = 2).real

def penalty(extpot, nx, ny, nz):
    field = reshape(array(extpot),(nx,ny,nz),order='F')
    laplacian = laplacian3d(field,nx,ny,nz)
    pen = sum(field * laplacian) * lambd # \int V*\nambla V
    pen_grad = reshape(laplacian,(nx*ny*nz), order='F') * lambd*2.0 # gradient of penalty
    return pen, pen_grad
