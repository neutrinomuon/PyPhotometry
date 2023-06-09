### PyPhotometry
####  Fully based on a Fortran legacy package to easily compute the photometric fluxes and magnitudes in different systems
email: [antineutrinomuon@gmail.com](mailto:antineutrinomuon@gmail.com), [jean@astro.up.pt](mailto:jean@astro.up.pt)

github repository: <a href="https://github.com/neutrinomuon/PyPhotometry">PyPhotometry</a>

last stable version: 0.0.12

© Copyright ®

J.G. - Jean Gomes @ 2023

<hr>

<img src="https://skillicons.dev/icons?i=python,fortran,c,numpy&theme=light"><br>
<img src="https://img.shields.io/pypi/pyversions/PyPhotometry"><img src="https://anaconda.org/neutrinomuon/PyPhotometry/badges/license.svg">

<hr>

<div align="center">
<img src='https://raw.githubusercontent.com/neutrinomuon/PyPhotometry/main/tutorials/PyPhotometry.png' width='60%'>
</div>

<hr>

#### <b>RESUME</b>

Original Fortran 2003+ routines date back to 2003-2004. Read the <a
href='https://github.com/neutrinomuon/PyPhotometry/blob/main/LICENSE.txt'>LICENSE.txt</a>
file.

PyPhotometry is a Python package based on a Fortran legacy package that allows
you to compute photometric fluxes and magnitudes in various photometric
systems. The package provides different magnitude systems, such as VEGA
standard, VEGA proposed by Bohlin and Gilland 2004, M_AB standard system, M_TG
standard system (Thuan & Gunn), WFPC2 system, FOCA at 2000, and without any
calibration.



Pyphot from M. FouesneauA is *NOT* part of the distribution, but used as a
comparison: <a
href='https://mfouesneau.github.io/pyphot/index.html#package-main-content'>https://mfouesneau.github.io/pyphot/index.html#package-main-content</a>. If
you want to install for comparison then:

<pre>
pip install pyphot
</pre>

However, it is not necessary for the usage of this package. Accompanying there
are several routines.

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
conda install -c neutrinomuon PyPhotometry
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
- VEGA proposed by Bohlin and Gilland 2004.
- M_AB standard system: Based on Oke (1974) reference.
- M_TG standard system (Thuan & Gunn): Based on Oke & Gunn (1983), Schild (1984), Schneider et al. (1983), Thuan & Gunn (1976), and Wade et al. (1979) references.
- WFPC2 system: Based on the Stone (1996) reference.
- FOCA at 2000 system.
- Without any calibration.

##### Calibration Stars

PyPhotometry provides calibration stars used in the magnitude systems:

- VEGA spectrum: Intrinsic Flux - erg/s/cm2/A.
- SUN spectrum: Intrinsic Flux - erg/s/A.
- F subdwarf: Used to calibrate the Thuan & Gunn system.

For more details on the usage and options, please refer to the PyPhotometry
GitHub repository.

<hr>

#### <b>STRUCTURE</b>

The main structure of the directories and files are:

<pre>
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
