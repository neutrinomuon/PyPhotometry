### PyPhotometry
####  Fully based on a Fortran legacy package to easily compute the photometric fluxes and magnitudes in different systems
email: [antineutrinomuon@gmail.com](mailto:antineutrinomuon@gmail.com), [jean@astro.up.pt](mailto:jean@astro.up.pt)

github repository: <a href="https://github.com/neutrinomuon/PyPhotometry">PyPhotometry</a>

last stable version: 0.0.9

© Copyright ®

J.G. - Jean Gomes
<!-- https://zenodo.org/badge/doi/10.5281/zenodo.10527286.svg -->
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.10527286.svg)](https://zenodo.org/badge/doi/10.5281/zenodo.10527286.svg)

<hr>

<img src="https://skillicons.dev/icons?i=python,fortran,c,numpy&theme=light"><br>
<img src="https://img.shields.io/pypi/pyversions/PyPhotometry"><img src="https://anaconda.org/neutrinomuon/PyPhotometry/badges/license.svg">

<hr>

<div align="center">
<img src='https://github.com/neutrinomuon/PyPhotometry/raw/main/figures/PyPhotometry.png' width=85%>
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

<img src='https://github.com/neutrinomuon/PyPhotometry/raw/main/figures/PyPhotometryIcon.png' width=120>

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

Now, the package PyPhotometry is in agreement with PEP 8 guidelines:

<br> --------------------------------------------------------------------
<br> Your code has been rated at 10.00/10 (previous run: 10.00/10, +0.00)

<hr>

#### <b>Brief Tutorial</b>

<div align="center">
<img src='https://github.com/neutrinomuon/PyPhotometry/raw/main/figures/FigurePyPhotometry.png' width='85%'>
</div>

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
- F subdwarf: Used to calibrate the Thuan & Gunn system. Only used for backward compatibility.

For more details on the usage and options, please refer to the PyPhotometry
GitHub repository.

<hr>

#### <b>STRUCTURE</b>
<pre>
#################################################
workspace
├── setup.py
├── data
│   ├── PLANCK_HFI.217GHz.txt
│   ├── PLANCK_LFI.030GHz.txt
│   ├── Spitzer_MIPS.160mu.txt
│   ├── SDSSu.txt
│   ├── SDSSr.txt
│   ├── Herschel_Pacs.green.txt
│   ├── Herschel_SPIRE.PMW.txt
│   ├── WISE1.txt
│   ├── ListFilters.txt
│   ├── Herschel_SPIRE.PLW.txt
│   ├── PLANCK_HFI.545GHz.txt
│   ├── Spitzer_MIPS.24mu.txt
│   ├── Herschel_SPIRE.PSW.txt
│   ├── Spitzer_IRAC.I4.txt
│   ├── PLANCK_LFI.070GHz.txt
│   ├── Spitzer_IRAC.I3.txt
│   ├── templates
│   │   ├── bc2003_hr_m62_chab_ssp_190.spec
│   │   ├── bc2003_hr_m62_salp_ssp_001.spec
│   │   ├── bc2003_hr_m62_chab_ssp_160.spec
│   │   └── bc2003_hr_m62_chab_ssp_001.spec
│   ├── SDSSz.txt
│   ├── IRAS.12mu.txt
│   ├── PLANCK_HFI.100GHz.txt
│   ├── Herschel_Pacs.red.txt
│   ├── PLANCK_LFI.044GHz.txt
│   ├── Herschel_SPIRE.PLW_ext.txt
│   ├── IRAS.100mu.txt
│   ├── Spitzer_MIPS.70mu.txt
│   ├── GalexFUV.txt
│   ├── WISE2.txt
│   ├── SDSSi.txt
│   ├── SDSSg.txt
│   ├── IRAS.25mu.txt
│   ├── Herschel_SPIRE.PSW_ext.txt
│   ├── Herschel_SPIRE.PMW_ext.txt
│   ├── PLANCK_HFI.353GHz.txt
│   ├── 2MASSKs.txt
│   ├── WISE4.txt
│   ├── Herschel_Pacs.blue.txt
│   ├── PLANCK_HFI.143GHz.txt
│   ├── calibration_stars
│   │   ├── bd17d4708_stisnic_001.fits
│   │   ├── Sun.dat
│   │   ├── Sun_LR.dat
│   │   ├── BD+17d4708.dat
│   │   ├── BD+17o4708.dat
│   │   ├── Vega.dat
│   │   ├── Filters_ReadMe.txt
│   │   ├── sun_reference_stis_001.fits
│   │   ├── VegaLR.dat
│   │   ├── kp00_6000.ascii
│   │   └── VegaLR_OLD.dat
│   ├── PLANCK_HFI.857GHz.txt
│   ├── README.md
│   ├── 2MASSH.txt
│   ├── IRAS.60mu.txt
│   ├── Spitzer_IRAC.I2.txt
│   ├── 2MASSJ.txt
│   ├── Spitzer_IRAC.I1.txt
│   ├── WISE3.txt
│   └── GalexNUV.txt
├── LICENSE.txt
├── Notes.txt
├── README_setup.txt
├── scripts
│   └── update_readme.py
├── __pycache__
│   └── Filters.cpython-39.pyc
├── tutorials
│   ├── PyPhotometry - Example 1.ipynb
│   ├── .ipynb_checkpoints
│   │   └── PyPhotometry - Example 1-checkpoint.ipynb
│   ├── .jupyter_ystore.db
│   └── PyPhotometry.png
├── .git
│   ├── HEAD
│   ├── objects
│   │   ├── pack
│   │   │   ├── pack-45351592350c6110dc6d5e83234328cfd06b823f.idx
│   │   │   ├── pack-45351592350c6110dc6d5e83234328cfd06b823f.pack
│   │   │   └── pack-45351592350c6110dc6d5e83234328cfd06b823f.rev
│   │   └── info
│   ├── config
│   ├── FETCH_HEAD
│   ├── info
│   │   └── exclude
│   ├── hooks
│   │   ├── push-to-checkout.sample
│   │   ├── prepare-commit-msg.sample
│   │   ├── pre-rebase.sample
│   │   ├── fsmonitor-watchman.sample
│   │   ├── post-update.sample
│   │   ├── pre-push.sample
│   │   ├── pre-commit.sample
│   │   ├── applypatch-msg.sample
│   │   ├── pre-merge-commit.sample
│   │   ├── pre-receive.sample
│   │   ├── pre-applypatch.sample
│   │   ├── commit-msg.sample
│   │   ├── sendemail-validate.sample
│   │   └── update.sample
│   ├── shallow
│   ├── refs
│   │   ├── heads
│   │   │   └── main
│   │   ├── remotes
│   │   │   └── origin
│   │   │       └── main
│   │   └── tags
│   ├── description
│   ├── branches
│   ├── logs
│   │   ├── HEAD
│   │   └── refs
│   │       ├── heads
│   │       │   └── main
│   │       └── remotes
│   │           └── origin
│   │               └── main
│   └── index
├── Literature
│   ├── Bohlin, Gordon, Tremblay (2014) - Techniques and Review of Absolute Flux Calibration from the Ultraviolet to the Mid-Infrared.pdf
│   └── Bohlin and Gilland (2004) - Absolute Flux Distribution of the SDSS Standard BD +17_4708.pdf
├── PyPhotometry
│   ├── win-arm64
│   │   ├── pyphotometry-0.0.6-py39_0.tar.bz2
│   │   └── pyphotometry-0.0.5-py39_0.tar.bz2
│   ├── osx-64
│   │   ├── pyphotometry-0.0.6-py39_0.tar.bz2
│   │   └── pyphotometry-0.0.5-py39_0.tar.bz2
│   ├── README.txt
│   ├── linux-s390x
│   │   ├── pyphotometry-0.0.6-py39_0.tar.bz2
│   │   └── pyphotometry-0.0.5-py39_0.tar.bz2
│   ├── meta.yaml
│   ├── linux-armv7l
│   │   ├── pyphotometry-0.0.6-py39_0.tar.bz2
│   │   └── pyphotometry-0.0.5-py39_0.tar.bz2
│   ├── linux-armv6l
│   │   ├── pyphotometry-0.0.6-py39_0.tar.bz2
│   │   └── pyphotometry-0.0.5-py39_0.tar.bz2
│   ├── osx-arm64
│   │   ├── pyphotometry-0.0.6-py39_0.tar.bz2
│   │   └── pyphotometry-0.0.5-py39_0.tar.bz2
│   ├── linux-64
│   │   ├── pyphotometry-0.0.6-py39_0.tar.bz2
│   │   └── pyphotometry-0.0.5-py39_0.tar.bz2
│   ├── linux-ppc64
│   │   ├── pyphotometry-0.0.6-py39_0.tar.bz2
│   │   └── pyphotometry-0.0.5-py39_0.tar.bz2
│   ├── linux-aarch64
│   │   ├── pyphotometry-0.0.6-py39_0.tar.bz2
│   │   └── pyphotometry-0.0.5-py39_0.tar.bz2
│   ├── win-32
│   │   ├── pyphotometry-0.0.6-py39_0.tar.bz2
│   │   └── pyphotometry-0.0.5-py39_0.tar.bz2
│   ├── linux-32
│   │   ├── pyphotometry-0.0.6-py39_0.tar.bz2
│   │   └── pyphotometry-0.0.5-py39_0.tar.bz2
│   ├── win-64
│   │   ├── pyphotometry-0.0.6-py39_0.tar.bz2
│   │   └── pyphotometry-0.0.5-py39_0.tar.bz2
│   └── linux-ppc64le
│       ├── pyphotometry-0.0.6-py39_0.tar.bz2
│       └── pyphotometry-0.0.5-py39_0.tar.bz2
├── showdown.min.js
├── version.txt
├── requirements.txt
├── index.html
├── figures
│   ├── PyPhotometryIcon.png
│   ├── logs.jpg
│   ├── FigurePyPhotometry.png
│   └── PyPhotometry.png
├── README.md
├── src
│   ├── fortran
│   │   ├── IntegralALL.f90
│   │   ├── PropFilters.f90
│   │   ├── moddatatype.mod
│   │   ├── makefile
│   │   ├── ReadFilters.f90
│   │   ├── GaussLegendreQuadrature.f90
│   │   ├── PropFilters.compile
│   │   ├── EvalFilters.f90
│   │   ├── EvalFilters.compile
│   │   ├── LINinterpol.f90
│   │   ├── DataTypes.f90
│   │   └── eval.exe
│   └── python
│       ├── __pycache__
│       │   └── PyPhotometry.cpython-39.pyc
│       ├── __init__.py
│       └── photometry.py
└── .github
    └── workflows
        └── update_readme.yml

128 directories, 533 files
#################################################
Generated with tree_colored @ 2023 - © Jean Gomes
#################################################
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


#### <b>LOGS</b>

<table>

<tr><td align='center'><img width=50 src="https://github.com/neutrinomuon/PyPhotometry/blob/main/figures/logs.jpg?raw=true"></td><td>LOGS</td></tr>

<tr><td>05/10/2023</td><td>J.G. changed some parts of the code to be compatible with Herschel filters. Also, added log dates explicitly in README file.</td></tr>

<tr><td>07/12/2023</td><td>J.G. changed some parts of the code to extend the calibration stars. Also, many changes in the interpolation scheme. Verification of the SQL tables. Also, added log dates explicitly in README file.</td></tr>

</table>

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
