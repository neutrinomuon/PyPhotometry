### PyPhotometry
####  Fully based on a Fortran legacy package to easily compute the photometric fluxes and magnitudes in different systems
email: [antineutrinomuon@gmail.com](mailto:antineutrinomuon@gmail.com), [jean@astro.up.pt](mailto:jean@astro.up.pt)

github repository: <a href="https://github.com/neutrinomuon/PyPhotometry">PyPhotometry</a>

last stable version: 0.0.7

© Copyright ®

J.G. - Jean Gomes @ 2023

<hr>

<img src="https://skillicons.dev/icons?i=python,fortran,c,numpy&theme=light"><br>
<img src="https://img.shields.io/pypi/pyversions/PyPhotometry"><img src="https://anaconda.org/neutrinomuon/PyPhotometry/badges/license.svg">

<hr>

<div align="center">
<img src='https://github.com/neutrinomuon/PyPhotometry/raw/main/tutorials/PyPhotometry.png' width='85%'>
</div>

<hr>

#### Requirements

The following packages are required to run this project:

- astropy>=5.0.4
- matplotlib>=3.7.1
- setuptools>=61.2.0
- SQLAlchemy>=1.4.32

and pyphot may be used for comparison and tests:

- pyphot>=1.4.4

You can install all the required packages by running the following command:

<pre>
pip install -r requirements.txt
</pre>

Additionally, you may optionally install pyphot for testing or comparison
purposes. 

Please note that pyphot is not a mandatory requirement for running this
project and is only recommended if you intend to test or compare with it.

#### <b>RESUME</b>

<img src='https://github.com/neutrinomuon/PyPhotometry/raw/main/figures/PyPhotometry.png' width=120>

<strong>PyPhotometry</strong> is a Python package that builds upon a
collection of Fortran 2003+ routines originally developed between 2003 and
2004. These routines are the foundation of the package and can be traced back
to that time period. The licensing details for the Fortran routines can be
found in the LICENSE.txt file included with the package.

The main purpose of PyPhotometry is to enable the computation of photometric
fluxes and magnitudes in various photometric systems. It offers support for
multiple magnitude systems, including the VEGA standard, the VEGA system
proposed by Bohlin and Gilland in 2004, the AB system, the TG standard system
(Thuan & Gunn), the WFPC2 system, the FOCA system at 2000, and also provides
an option without any calibration.

It's important to note that PyPhotometry does not include the Pyphot package
developed by M. Fouesneau, but it can be used for comparison purposes.

However, it is not mandatory to install Pyphot in order to use PyPhotometry. The PyPhotometry package comes with its own set of accompanying routines that provide the necessary functionality.

Original Fortran 2003+ routines date back to 2003-2004. Read the <a
href='https://github.com/neutrinomuon/PyPhotometry/blob/main/LICENSE.txt'>LICENSE.txt</a>
file.

PyPhotometry is a Python package based on a Fortran legacy package that allows
you to compute photometric fluxes and magnitudes in various photometric
systems. The package provides different magnitude systems, such as VEGA
standard, VEGA proposed by Bohlin and Gilland 2004, AB system, TG standard
system (Thuan & Gunn), WFPC2 system, FOCA at 2000, and without any
calibration.

Pyphot from M. Fouesneau is *NOT* part of the distribution, but used as a
comparison: <a
href='https://mfouesneau.github.io/pyphot/index.html#package-main-content'>https://mfouesneau.github.io/pyphot/index.html#package-main-content</a>. If
you want to install for comparison then:

<pre>
pip install pyphot
</pre>

However, it is not necessary for the usage of this package. This package is
meant for a comparison, but PyPhotometry legacy routines are more
general. Accompanying there are several other routines.

<hr>

#### <b>Brief Tutorial</b>

A brief tutorial can be found at <a
href='https://github.com/neutrinomuon/PyPhotometry/blob/main/tutorials/PyPhotometry%20-%20Example%201.ipynb'>PyPhotometry
Example1.ipynb</a>

<hr>

#### <b>INSTALLATION</b>

You can easily install <a
href=https://pypi.org/project/PyPhotometry/>PyPhotometry</a> by using pip -
<a href='https://pypi.org/'>PyPI - The Python Package Index</a>:

<pre>
pip install PyPhotometry
</pre>

<br>or by using a generated conda repository <a
href='https://anaconda.org/neutrinomuon/PyPhotometry'>https://anaconda.org/neutrinomuon/PyPhotometry</a>:

<img src="https://anaconda.org/neutrinomuon/PyPhotometry/badges/version.svg"><img src="https://anaconda.org/neutrinomuon/PyPhotometry/badges/latest_release_date.svg"><img src="https://anaconda.org/neutrinomuon/PyPhotometry/badges/platforms.svg">

<pre>
conda install -c neutrinomuon pyphotometry
</pre>

<br>OBS.: Linux, OS-X and Windows pre-compilations available in conda.

You can also clone the repository and install by yourself in your machine:

<pre>
git clone https://github.com/neutrinomuon/PyPhotometry
python setup.py install
</pre>

<hr>

#### <b>METHOD & REFERENCES</b>

##### Magnitude Systems
The following magnitude systems are supported by PyPhotometry:

- VEGA standard: Based on the Bessel (2005), Cousins & Jones (1976), and Kitchin (2003) references.
- VEGA proposed by Bohlin and Gilland 2004;
- AB standard system: Based on Oke (1974) reference;
- TG standard system (Thuan & Gunn): Based on Oke & Gunn (1983), Schild (1984), Schneider et al. (1983), Thuan & Gunn (1976), and Wade et al. (1979) references;
- WFPC2 system: Based on the Stone (1996) reference;
- FOCA at 2000 system;
- Without any calibration.

##### Calibration Stars

PyPhotometry provides calibration stars used in the magnitude systems:

- VEGA spectrum: Intrinsic Flux - [erg/s/cm2/A].
- SUN spectrum: Intrinsic Flux - [erg/s/A].
- F subdwarf: Used to calibrate the Thuan & Gunn system.

For more details on the usage and options, please refer to the PyPhotometry
GitHub repository.

<hr>

#### <b>STRUCTURE</b>

The main structure of the directories and files are:

<pre>
PyPhotometry
├── Literature
│   └── Bohlin and Gilland (2004) - Absolute Flux Distribution of the SDSS Standard BD +17_4708.pdf
├── README.md
├── __pycache__
│   └── Filters.cpython-39.pyc
├── showdown.min.js
├── scripts
│   └── update_readme.py
├── index.html
├── LICENSE.txt
├── data
│   ├── IRAS.12mu.txt
│   ├── Herschel_SPIRE.PSW_ext.txt
│   ├── GalexNUV.txt
│   ├── SDSSg.txt
│   ├── 2MASSJ.txt
│   ├── WISE4.txt
│   ├── 2MASSKs.txt
│   ├── Herschel_SPIRE.PLW_ext.txt
│   ├── GalexFUV.txt
│   ├── filters.db
│   ├── Herschel_SPIRE.PMW_ext.txt
│   ├── IRAS.60mu.txt
│   ├── SDSSr.txt
│   ├── calibration_stars
│   │   ├── sun_reference_stis_001.fits
│   │   ├── Filters_ReadMe.txt
│   │   ├── kp00_6000.ascii
│   │   ├── Sun_LR.dat
│   │   ├── bd17d4708_stisnic_001.fits
│   │   ├── BD+17o4708.dat
│   │   ├── BD+17d4708.dat
│   │   ├── VegaLR_OLD.dat
│   │   ├── VegaLR.dat
│   │   └── Sun.dat
│   ├── 2MASSH.txt
│   ├── IRAS.100mu.txt
│   ├── WISE3.txt
│   ├── ListFilters.txt
│   ├── SDSSi.txt
│   ├── SDSSz.txt
│   ├── IRAS.25mu.txt
│   ├── SDSSu.txt
│   ├── WISE2.txt
│   └── WISE1.txt
├── setup.py
├── tutorials
│   └── PyPhotometry.png
├── src
│   ├── python
│   │   ├── __pycache__
│   │   ├── __init__.py
│   │   └── Photometry.py
│   └── fortran
│       ├── EvalFilters.compile
│       ├── IntegralALL.f90
│       ├── PropFilters.compile
│       ├── DataTypes.f90
│       ├── PropFilters.f90
│       ├── EvalFilters.cpython-39-x86_64-linux-gnu.so
│       ├── ReadFilters.f90
│       ├── GaussLegendreQuadrature.f90
│       ├── PropFilters.cpython-39-x86_64-linux-gnu.so
│       ├── makefile
│       ├── LINinterpol.f90
│       └── EvalFilters.f90
├── version.txt
├── PyPhotometry.egg-info
│   ├── PKG-INFO
│   ├── dependency_links.txt
│   ├── SOURCES.txt
│   ├── top_level.txt
│   └── requires.txt
└── build
    ├── lib.linux-x86_64-3.9
    │   └── PyPhotometry
    ├── src.linux-x86_64-3.9
    │   ├── PyPhotometry
    │   ├── build
    │   └── numpy
    └── temp.linux-x86_64-3.9
        ├── ccompiler_opt_cache_ext.py
        ├── src
        ├── PyPhotometry
        ├── .libs
        └── build

23 directories, 63 files
</pre>

<br>PyPhotometry.py is a python wrapper to the library in fortran called
PyPhotometry.flib. The fortran directory can be compiled separately for
each individual subroutine.

<hr>

#### ISSUES AND CONTRIBUTIONS

If you encounter any issues with this project, please feel free to submit an
issue on the GitHub repository. We appreciate your feedback and are committed
to improving the quality of our codebase.

If you'd like to contribute to this project, we welcome pull requests from the
community. Before submitting a pull request, please make sure to fork the
repository and create a new branch for your changes. Once your changes are
complete, submit a pull request and we'll review your code as soon as
possible.

For any questions or concerns about contributing, please contact the project
maintainer at antineutrinomuon@gmail.com. Thank you for your interest in
contributing to our project!

<hr>

#### <b>LICENSE</b>

This software is provided "AS IS" (see DISCLAIMER below). Permission to use,
for non-commercial purposes is granted. Permission to modify for personal or
internal use is granted, provided this copyright and disclaimer are included
in ALL copies of the software. All other rights are reserved. In particular,
redistribution of the code is not allowed without explicit permission by the
author.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
