import numpy as np
import matplotlib.pyplot as plt



def image_plot_pf(pf, extent):
    ax = plt.gca()
    im = ax.imshow( pf, cmap='viridis', extent = extent)
    ax.set_xlabel(r'$\Delta$RA (arcsec)')
    ax.set_ylabel(r'$\Delta$Dec (arcsec)')
    ax.set_aspect('equal')

    cbar = plt.colorbar(im)
    cbar.set_label('Polarization Fraction [%]')

    plt.xlim(3, -3) #manually plots x in other direction
    plt.ylim(-3, 3)
    plt.show()

def cont_plot_pf(pf, extent, ra, dec):
    ax = plt.gca()
    im = ax.contourf(ra, dec ,pf, cmap='viridis', extent = extent, levels = 50 ) #RA plotted on x-axis
    ax.set_xlabel(r'$\Delta$RA (arcsec)')
    ax.set_ylabel(r'$\Delta$Dec (arcsec)')
    ax.set_aspect('equal')

    cbar = plt.colorbar(im)
    cbar.set_label('Polarization Fraction [%]')

    plt.xlim(3, -3) #manually plots x in other direction
    plt.ylim(-3, 3)
    plt.show()

def im_plot_p(P,extent, noise_I, I):
    ax = plt.gca()
    im = ax.imshow( P,cmap='viridis', extent = extent)
    ax.set_xlabel(r'$\Delta$RA (arcsec)')
    ax.set_ylabel(r'$\Delta$Dec (arcsec)')
    ax.set_aspect('equal')

    contours = ax.contour(I, colors='white', linewidths=.5, extent=extent,      #plots I contour levels
                          levels=[3 * noise_I, 10 * noise_I, 25 * noise_I, 50 * noise_I, 100 * noise_I, 200 * noise_I,
                                  325 * noise_I, 500 * noise_I, 1000 * noise_I])
    ax.clabel(contours, inline=1, fontsize=10)

    cbar = plt.colorbar(im)
    cbar.set_label('Linear Polarization Intensity [Jy/beam]')

    plt.xlim(3,-3)
    plt.ylim(-3,3)
    plt.show()

def cont_plot_p(P, I, noise_I, extent, ra, dec):
    ax = plt.gca()
    im = ax.contourf(ra, dec, P, cmap='viridis', extent = extent, levels=50)
    ax.set_xlabel(r'$\Delta$RA (arcsec)')
    ax.set_ylabel(r'$\Delta$Dec (arcsec)')
    ax.set_aspect('equal')

    contours = ax.contour(ra, dec, I, colors='white', linewidths=.5, extent=extent,
                          levels=[3 * noise_I, 10 * noise_I, 25 * noise_I, 50 * noise_I, 100 * noise_I, 200 * noise_I,
                                  325 * noise_I, 500 * noise_I, 1000 * noise_I])
    ax.clabel(contours, inline=1, fontsize=10)

    cbar = plt.colorbar(im)
    cbar.set_label('P [Jy/beam]')

    plt.xlim(3,-3)
    plt.ylim(-3, 3)
    plt.show()