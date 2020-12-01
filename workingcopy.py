import numpy as np
from physunits import *

class GaussianBeam():

    def __init__(self, wavelength, min_waist, medium_index=1.0, **kwargs):
        self.wavelength = wavelength
        self.min_waist = min_waist
        self.q = 1.0j * np.pi * self.min_waist ** 2 / (self.wavelength * medium_index)

    def _abcd_transform(self, abcd_matrix=np.array([[1., 0.], [0., 1.]])):
        A, B, C, D = abcd_matrix.flatten()
        # print(self.q, self.min_waist)
        self.q = (self.q * A + B) / (self.q * C + D)
        self.min_waist = np.sqrt(self.wavelength * self.q.imag / np.pi)
        # print(self.q, self.min_waist)

    def render(self, z_grid):
        self.waist = self.min_waist * np.sqrt(1 + z_grid**2 / self.q.imag**2)

    def translate(self, distance, n_ref=1.):
        self._abcd_transform(np.array([[1., distance / n_ref], [0., 1.]]))

    def refraction(self, index_medium_1, index_medium_2):
        """ Flat interface refraction abcd matrix.
        
        Arguments:
        index_medium_1 -- index of refraction of medium 1
        index_medium_2 -- index of refraction of medium 2
        """
        self._abcd_transform(np.array([[1., 0.], [0., index_medium_1/index_medium_2]]))

    def thin_lens(self, focal_length):
        """ Thin lens abcd matrix.
        
        Arguments:
        f -- focal length, f > 0 positive lens
        """
        self._abcd_transform(np.array([[1., 0.], [-1/focal_length, 1.]]))


    def profile(self, axis, z_grid, color='salmon', alpha=0.1):
        self.render(z_grid)
        axis.plot(z_grid/self.q.imag, self.waist/self.min_waist, c=color)
        axis.plot(z_grid/self.q.imag, -self.waist/self.min_waist, c=color)
        axis.fill_between(
            z_grid/self.q.imag,
            self.waist/self.min_waist,
            -self.waist/self.min_waist, 
            color=color,
            alpha=alpha,
            )


if __name__ in "__main__":
    import matplotlib.pyplot as plt

    plt.figure(1)
    ax = plt.subplot(111)

    beam = GaussianBeam(wavelength=1064.1*nm, min_waist=200*um,)
    beam.translate(-100*mm)
    z0 = np.linspace(-200*mm, 200*mm, 64)
    beam.profile(ax, z0, color='goldenrod')

    # z1 = np.linspace(200*mm, 400*mm, 64)
    # beam.translate(z0[-1])
    # beam.thin_lens(focal_length=100*mm)
    # beam.profile(ax, z1, color='orange')
    
    plt.show()