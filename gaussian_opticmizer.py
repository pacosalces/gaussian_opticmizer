#! /usr/bin/python
#
# __author__ = [pacosalces,]

import numpy as np
import matplotlib.pyplot as plt

# Update matplotlib style
plt.rcParams.update(
    {
        "text.usetex": True,
        "scatter.edgecolors": "k",
        "image.cmap": "RdBu_r",
        "image.origin": "lower",
        "font.family": "serif",
        "font.size": 14,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "savefig.dpi": 100,
        "pdf.compression": 9,
        "figure.figsize": (16, 4),
    }
)

from physunits import *
from catalog import loader


class GaussianBeam:
    def __init__(self, wavelength, w0, z0):
        self.wavelength = wavelength
        self.w0 = w0
        self.z0 = z0

        # Rayleigh range
        self.zR = np.pi * self.w0 ** 2 / self.wavelength

        # q-parameter at min waist
        self.q0 = self.z0 + 1j * self.zR

    def q(self, z):
        return self.q0 - z

    def waist(self, z):
        """Evaluate waist over z"""
        return self.w0 * np.sqrt(
            1 + (self.q(z).real ** 2 / self.q(z).imag ** 2)
        )


class OpticalSystem:
    def __init__(self, z, seed):
        # Initially empty over full span
        self.segments = [z]
        self.beams = [seed]

    def segmentate(self, z0, abcd_matrix):
        """Each time a new element is added, break segments and update
        beams"""
        A, B, C, D = abcd_matrix.flatten()
        beam_counter = -1
        new_beams, new_segments = [], []
        for segment in list(self.segments):
            beam_counter += 1
            if z0 > segment[0] and z0 <= segment[-1]:
                input_beam = self.beams[beam_counter]
                new_segments.append(segment[segment <= z0])
                new_segments.append(segment[segment > z0])
                q_z0 = self.beams[beam_counter].q(z0)
                q_new = (q_z0 * A + B) / (q_z0 * C + D)
                output_beam = GaussianBeam(
                    wavelength=input_beam.wavelength,
                    w0=np.sqrt(input_beam.wavelength * q_new.imag / np.pi),
                    z0=q_new.real + z0,
                )
                new_beams.append(input_beam)
                new_beams.append(output_beam)
            else:
                new_segments.append(segment)
                new_beams.append(self.beams[beam_counter])

        self.beams = new_beams
        self.segments = new_segments

    def thin_lens(self, z0, f):
        """Thin lens abcd matrix.

        Arguments:
        z0 -- lens position
        f -- lens focal length (f > 0 positive convention)
        """
        self.segmentate(z0, np.array([[1.0, 0.0], [1 / f, 1.0]]))

    def normal_refraction(self, z0, n1, n2):
        """Flat surface refraction abcd matrix at normal incidence.

        Arguments:
        z0 -- surface position
        n1 -- index of refraction of input
        n2 -- index of refraction of output
        """
        self.curved_refraction(z0, np.inf, n1, n2)

    def curved_refraction(self, z0, R, n1, n2):
        """Spherical surface with R radius of curvature separating
        two media with refractive indices n1, n2

        Args:
            z0 (float): surface position along optical axis
            R (float): radius of curvature (np.inf for flat surface)
            n1 (float): index of refraction of input medium
            n2 (float): index of refraction of output medium
        """
        if R == np.inf:
            C = 0.0
        elif R == 0.0:
            raise ValueError("Unphysical curvature!")
        else:
            C = (n1 - n2) / (-R * n2)
        self.segmentate(z0, np.array([1.0, 0.0], [C, n1 / n2]))

    def thick_lens(self, z0, R1, R2, n1, n2, d, **other):
        """Thick lens abcd matrix

        Args:
            z0 (float): lens center along optical axis
            R1 (float): input radius of curvature (np.inf for flat surface)
            R2 (float): output radius of curvature (np.inf for flat surface)
            n1 (float): index of refraction of input medium
            n2 (float): index of refraction of output medium
            d (float): lens thickness along optical axis
        """
        if R1 == np.inf:
            C1 = 0.0
        elif R1 == 0.0:
            raise ValueError("Unphysical curvature!")
        else:
            C1 = (n1 - n2) / (-R1 * n2)

        if R2 == np.inf:
            C2 = 0.0
        elif R2 == 0.0:
            raise ValueError("Unphysical curvature!")
        else:
            C2 = (n2 - n1) / (-R2 * n1)

        # Construct full transfer matrix
        Mi = np.array([[1.0, 0.0], [C1, n1 / n2]])
        Mt = np.array([[1.0, d], [0.0, 1.0]])
        Mo = np.array([[1.0, 0.0], [C2, n2 / n1]])

        self.segmentate(z0, Mo @ Mt @ Mi)

    def catalog_lens(self, z0, part_no):
        """Find a lens in catalog files

        Args:
            z0 (float): lens center along optical axis
            part_no (str): Name of file (must exist as yaml file)
        """
        catalog_params = loader(part_no)
        self.thick_lens(z0, **catalog_params)

    def draw(self, **plot_kwargs):
        """Draw all beam profiles in their segments"""
        plt.figure()
        ax = plt.subplot(111)
        for j, segment in enumerate(self.segments):
            beam = self.beams[j]
            min_waist = np.sqrt(beam.wavelength * beam.q0.imag / np.pi)
            if "color" not in plot_kwargs:
                color = f"C{j}"
            else:
                color = plot_kwargs["color"]
            wj = beam.waist(segment)
            plt.plot(segment / mm, wj / um, c=color)
            plt.plot(segment / mm, -wj / um, c=color)
            if (
                beam.q0.real > self.segments[0][0]
                and beam.q0.real < self.segments[-1][-1]
            ):
                plt.text(
                    s=fR"{min_waist / um:.2f} um",
                    x=beam.q0.real / mm,
                    y=1.6 * min_waist / um,
                    c=color,
                )
                plt.scatter(
                    beam.q0.real / mm,
                    min_waist / um,
                    color=color,
                    marker=7,
                    lw=3,
                )
                plt.scatter(
                    beam.q0.real / mm,
                    -min_waist / um,
                    color=color,
                    marker=6,
                    linewidth=3,
                )
            plt.fill_between(
                segment / mm,
                wj / um,
                -wj / um,
                **plot_kwargs,
            )

        plt.xlim((self.segments[0][0] / mm, self.segments[-1][-1] / mm))
        plt.ylabel(fR"Waist ($\mu \rm m$)")
        plt.xlabel(fR"$z \,(\rm mm$)")
        plt.tight_layout()
        plt.show()


class Opticsmizer:
    def __init__(self, optical_system):
        self.optical_system = optical_system

    def minimize_waist(
        self,
    ):
        pass


if __name__ == "__main__":
    # First create an optical system and seed with an input gaussian beam
    z_span = np.linspace(-200 * mm, 800 * mm, 2 ** 11)
    seed = GaussianBeam(wavelength=1064.1 * nm, z0=-40 * cm, w0=500 * um)
    telescope = OpticalSystem(z=z_span, seed=seed)

    # Add a couple of thin lenses
    telescope.catalog_lens(z0=0 * mm, part_no="./catalog/LA4380")
    telescope.thin_lens(z0=207 * mm, f=100 * mm)
    telescope.catalog_lens(z0=260 * mm, part_no="./catalog/LBF254050C")
    telescope.thin_lens(z0=467 * mm, f=157 * mm)
    telescope.draw(alpha=0.1)
    w_input = telescope.beams[0].waist(z_span[0])
    w_output = telescope.beams[-1].waist(z_span[-1])
    print(f"The transverse magnification is {w_output/w_input:.2f}")
