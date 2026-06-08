from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from radio_beam import Beam
import astropy.units as u
from scipy.interpolate import RegularGridInterpolator as rgi

class Plot:
    """Class for plotting polarization plots with RA as x-axis and DEC as y-axis"""

    def __init__(self, ra, dec):
        self.ra = ra
        self.dec = dec

    def set_beam(self, bmaj, bmin, bpa):
        self.bmaj = bmaj
        self.bmin = bmin
        self.bpa = bpa

    def plot (self, value, value_name, unit, beam=False,contour=None,ax=None, color='Viridis'):
        if ax is None:
            ax = plt.gca()
        if contour is not None:
            im = ax.contourf(self.ra, self.dec, value, cmap=color, levels=50)
        else:
            #using convention of east to left, north to up
            extent = [
                self.ra.max(), self.ra.min(),
                self.dec.min(), self.dec.max()
            ]
            im = ax.imshow( value,cmap=color, extent= extent)

        ax.set_xlabel(r'$\Delta$RA (arcsec)')
        ax.set_ylabel(r'$\Delta$Dec (arcsec)')
        ax.set_aspect('equal')

        if beam == True:


        cbar = plt.colorbar(im)
        cbar.set_label(f'{value_name} [{unit}]')

        #change later - but base limit and axis switch
        plt.xlim(3, -3)
        plt.ylim(-3, 3)
        plt.show()

    def plot_beam(self,ax):
        my_beam = Beam(self.bmaj * u.arcsec, self.bmin * u.arcsec, self.bpa * u.degree)
        y_cen_pix, x_cen_pix = -2, 2
        pixscale = u.arcsec
        ellipse_artist = my_beam.ellipse_to_plot(x_cen_pix, y_cen_pix, pixscale, facecolor='black')
        ax.add_artist(ellipse_artist)