#George Emanuel
#emanuega0@gmail.com
#Copyright Presidents and Fellows of Harvard College, 2017.


import numpy as np
import scipy as sp
import scipy.ndimage
from scipy import fftpack

from core.analysis import analysis

def cart2polar(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return r, theta

def polar2cart(r, theta):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

def index_coords(data, origin=None):
    ny, nx = data.shape[:2]
    if origin is None:
        origin_x, origin_y = nx // 2, ny // 2
    else:
        origin_x, origin_y = origin
    x, y = np.meshgrid(np.arange(nx), np.arange(ny))
    x -= origin_x
    y -= origin_y
    return x, y

def reproject_image_into_polar(data, origin=None):
    ny, nx = data.shape[:2]
    if origin is None:
        origin = (nx//2, ny//2)

    x, y = index_coords(data, origin=origin)
    r, theta = cart2polar(x, y)

    r_i = np.linspace(r.min(), r.max(), nx)
    theta_i = np.linspace(theta.min(), theta.max(), ny)
    theta_grid, r_grid = np.meshgrid(theta_i, r_i)

    xi, yi = polar2cart(r_grid, theta_grid)
    xi += origin[0] 
    yi += origin[1]
    xi, yi = xi.flatten(), yi.flatten()
    coords = np.vstack((xi, yi)) 

    zi = sp.ndimage.map_coordinates(data, coords, order=1)
    mapped = zi.reshape((nx,ny))

    return mapped, r_i, theta_i

def psd(data):
    return np.abs(fftpack.fftshift(fftpack.fft2(data)))**2

def radial_psd(data):
    return np.sum(reproject_image_into_polar(psd(data))[0], axis=1)

def deviation_point(data):
    logPSD = np.log(radial_psd(data))
    baseline = np.mean(logPSD[250:400])

    index = 250
    deviationFound = False
    deviationCount = 0
    while not deviationFound and index > 0:
        index -= 1
        if logPSD[index] > baseline:
            deviationCount += 1
        else:
            deviationCount = 0
        if deviationCount >= 10:
            deviationFound = True

    return index 



def is_focused(data):
    return deviation_point(data) > 150


    '''
    radialPower = radial_psd(data[50:-50, 50:-50]) \
            if all(x >= 100 for x in data.shape) \
            else radial_psd(data)
    return sum(radialPower[25:65]) < sum(radialPower[70:140])
    '''

