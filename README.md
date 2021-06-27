# Gaussian Opticmizer

A(nother) simple, one-dimensional gaussian beam optimizer using ray-transfer matrix methods and q-beam. The goal is to eventually populate a "catalog" from the publicly available commercial lens parameters from mainstream manufacturers. 

## Installation
Clone the repository into a local file directory. Will hopefully package it in the near future... but no promises.

## Usage 
If you want to quickly estimate what a paraxial optical system will do to a monochromatic gaussian beam, this should be helpful. As a quick example consider the following script:

```python
>> import numpy as np
>> from physunits.length import nm, um, mm, cm
>> from gaussian_opticmizer import GaussianBeam, OpticalSystem

# Create an optical system instance over some distance "z" and 
# a seed gaussian beam that serves as the input beam 
>> z_span = np.linspace(-200 * mm, 800 * mm, 2 ** 11)
>> seed = GaussianBeam(wavelength=1064.1 * nm, z0=-40 * cm, w0=500 * um)
>> telescope = OpticalSystem(z=z_span, seed=seed)

# Add a couple of lenses
>> telescope.catalog_lens(z0=0 * mm, part_no="./catalog/LA4380")
>> telescope.thin_lens(z0=207 * mm, f=100 * mm)
>> telescope.catalog_lens(z0=260 * mm, part_no="./catalog/LBF254050C")
>> telescope.thin_lens(z0=467 * mm, f=157 * mm)

# All beams are automatically updated, so we can just print information,
# for example below we can estimate the magnification of the telescope
>> w_input = telescope.beams[0].waist(z_span[0])
>> w_output = telescope.beams[-1].waist(z_span[-1])
>> print(f"The transverse magnification is {w_output/w_input:.2f}")

The transverse magnification is 3.15

# We can also use the built-in "draw" method to visualize the system
>> telescope.draw(alpha=0.1)
```
![alt text](/examples/example_1.png?raw=true)