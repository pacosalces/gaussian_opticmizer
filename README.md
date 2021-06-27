# Gaussian Opticmizer

A(nother) simple, one-dimensional gaussian beam optimizer using ray-transfer matrix methods and q-beam. The goal is to eventually populate a "catalog" from the publicly available commercial lens parameters from mainstream manufacturers. 

## Installation
Clone the repository into a local file directory. Will hopefully package it in the near future... but no promises.

## Usage 
If you want to quickly estimate what a paraxial optical system will do to a monochromatic gaussian beam, the code in this repository should be helpful.

As a quick example consider the following script fragments to simulate a magnifying telescope.

```python
from physunits.length import *
from gaussian_opticmizer import GaussianBeam, OpticalSystem
```
We then create an optical system spanning some distance "z" by instantiating the optica system class. For this we also use a "seed" gaussian beam whose wavelength, minimum waist, and center (minimum waist position) are specified.
```python
z_span = np.linspace(-200 * mm, 800 * mm, 2 ** 11)
seed = GaussianBeam(wavelength=1064.1 * nm, z0=-40 * cm, w0=500 * um)
telescope = OpticalSystem(z=z_span, seed=seed)
```
Add a couple of different lenses at locations "z0". Every time the optical system gets a new element, the seed beam is updated. The elements can be inserted in an arbitrary order.
```python
telescope.catalog_lens(z0=0 * mm, part_no="./catalog/LA4380")
telescope.thin_lens(z0=207 * mm, f=100 * mm)
telescope.catalog_lens(z0=260 * mm, part_no="./catalog/LBF254050C")
telescope.thin_lens(z0=467 * mm, f=157 * mm)
```
Here we may use the built-in "draw" method to visualize the system.
```python
telescope.draw(alpha=0.1)
```

![alt text](/examples/example_1.png?raw=true)

Finally, the "telescope" object contains some useful attributes, so we can access information at any time; for example we can estimate the magnification of the telescope by issuing:
```python
w_input = telescope.beams[0].waist(z_span[0])
w_output = telescope.beams[-1].waist(z_span[-1])
print(f"The transverse magnification is {w_output/w_input:.2f}")

>> The transverse magnification is 3.15

```
