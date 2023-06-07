### PyPhotometry
####  Fully based on a Fortran legacy package to easily compute the photometric fluxes and magnitudes in different systems
Obs.: A Fortran legacy Interpolation routines also furnished <br>
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
file. When analyzing astronomical spectra, astronomers often bin the data to
increase the signal-to-noise ratio and reduce the effects of noise in the
data. Binning refers to the process of averaging the intensity of adjacent
spectral channels, or pixels, to produce a new, coarser set of data.

In the process of binning, it is important to ensure that the principle of
flux density conservation is maintained. This means that the total energy
emitted by the object, and hence its flux density, must remain constant after
binning.

To conserve flux density, the intensity of each binned pixel should be scaled
by the number of pixels it represents. For example, if two adjacent pixels are
binned together, the intensity of the resulting bin should be the sum of the
intensities of the two original pixels, divided by two. This ensures that the
total energy in the bin is conserved, and that the flux density of the object
remains the same.

It's worth noting that binning can introduce errors in the spectral data,
especially if the signal-to-noise ratio is low or if the binning is too
coarse. In general, astronomers choose a binning size that balances the need
for a high signal-to-noise ratio with the desire to maintain the spectral
resolution and avoid introducing significant errors in the data.

In summary, the principle of flux density conservation is important to
consider when binning astronomical spectra, and astronomers need to scale the
intensity of each binned pixel to ensure that the total energy emitted by the
object is conserved. SpectRes from A. Carnall is *NOT* part of the distribution,
but used as a comparison: <a
href='https://github.com/ACCarnall/SpectRes'>https://github.com/ACCarnall/SpectRes</a>. If you want to install for comparison then:

<pre>
pip install spectres
</pre>

However, it is not necessary for the usage of this package. Accompanying there are several routines for interpolations.

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

Here, the method used is with the Cumulative function to produce a new
flux-conserved, some options can be chosen for the interpolation:

<table border="1">

<tr><td width='20%'>Integer Number</td><td width='33%'>Option: Interpolation Schemes</td><td>Brief Description</td></tr>

<tr><td>0)</td><td>AkimaSpline</td><td>Akima Spline interpolation. The Akima spline is a C1 differentiable function (that is, has a continuous first
derivative) but, in general, will have a discontinuous second derivative at the knot points.<p>

<ul>
<il>Akima, H. (1970). A New Method of Interpolation and Smooth Curve Fitting Based on Local Procedures. <i>Journal of the ACM</i>, 17(4),
589-602. Association for Computing Machinery. DOI: 10.1145/321607.321609. Link: <a href="https://dl.acm.org/doi/10.1145/321607.321609">https://dl.acm.org/doi/10.1145/321607.321609</a></il><p>

<il>De Boor, C. (1978). A Practical Guide to Splines. Springer-Verlag. ISBN: 978-0-387-95366-3. Link: <a href="https://www.springer.com/gp/book/9780387953663">https://www.springer.com/gp/book/9780387953663</a></il><p>

<il>Forsythe, G.E. (1979) <i>Computer Methods for Mathematical Computations</i>. Prentice-Hall, Inc. DOI: 10.1002/zamm.19790590235. Link: <a href="https://onlinelibrary.wiley.com/doi/10.1002/zamm.19790590235">https://onlinelibrary.wiley.com/doi/10.1002/zamm.19790590235</a></il><p>

<il>Bartels, R.H., Beatty, J.C., and Barsky, B.A. (1987). An Introduction to Splines for Use in Computer Graphics and Geometric Modeling. Morgan Kaufmann
Publishers. Link: <a href="https://www.osti.gov/biblio/5545263">https://www.osti.gov/biblio/5545263</a></il><p>

<il>Press, W. H., Teukolsky, S. A., Vetterling, W. T., &amp; Flannery, B. P. (2007). Numerical Recipes: The Art of Scientific Computing (3rd
ed.). Cambridge University Press. ISBN: 978-0521880688. Link: <a href="http://numerical.recipes/">http://numerical.recipes/</a></il><p>

<il>Cheney, W. and Kincaid, D. (2008). Numerical Mathematics and Computing (6th ed.). Brooks/Cole. ISBN: 978-0495114758. Link: <a href="https://www.amazon.com/Numerical-Mathematics-Computing-Ward-Cheney/dp/0495114758">https://www.amazon.com/Numerical-Mathematics-Computing-Ward-Cheney/dp/0495114758</a></il><p>

<il>Burden, R. L., &amp; Faires, J. D. (2010). Numerical Analysis (9th ed.). Brooks/Cole. ISBN: 978-0538733519. Link: <a href="https://www.amazon.com/Numerical-Analysis-Richard-L-Burden/dp/0538733519">https://www.amazon.com/Numerical-Analysis-Richard-L-Burden/dp/0538733519</a></il><p>

</ul>

</td></tr>

<tr><td>1)</td><td>Interpolado</td><td>Based on a linear interpolation within
a table of pair values.<p>

<ul>

<il>Atkinson, K. E. (1991). <i>An introduction to numerical analysis</i> (2nd ed.). John Wiley & Sons. ISBN: 978-0-471-62489-9. Link: <a href="https://www.wiley.com/en-us/An+Introduction+to+Numerical+Analysis%2C+2nd+Edition-p-9780471624899">https://www.wiley.com/en-us/An+Introduction+to+Numerical+Analysis%2C+2nd+Edition-p-9780471624899</a></il><p>

<il>Press, W. H., Teukolsky, S. A., Vetterling, W. T., &amp; Flannery, B. P. (2007). <em>Numerical Recipes: The Art of Scientific Computing</em> (3rd
ed.). Cambridge University Press. ISBN: 978-0521880688. Link: <a href="http://numerical.recipes/">http://numerical.recipes/</a></il><p>

<il>Chapra, S. C., & Canale, R. P. (2010). <i>Numerical methods for engineers</i> (6th ed.). McGraw-Hill. ISBN: 978-0073401065. Link: <a href="https://www.amazon.com/Numerical-Methods-Engineers-Steven-Chapra/dp/0073401064">https://www.amazon.com/Numerical-Methods-Engineers-Steven-Chapra/dp/0073401064</a></il><p>

<il>Burden, R. L., Faires, J. D. & A.M. Burden (2015). <i>Numerical analysis (10th ed.)</i> (9th ed.). Cengage Learning. ISBN: 978-1305253667. Link: <a href="https://www.amazon.com/Numerical-Analysis-Richard-L-Burden/dp/1305253663">https://www.amazon.com/Numerical-Analysis-Richard-L-Burden/dp/1305253663</a></il><p>

</ul>

</td></tr>

<tr><td>2)</td><td>LINdexerpol</td><td>Based on a linear interpolation within
a table of pair values using indexing. The references are the same as in
1).</td></tr>

<tr><td>3)</td><td>LINinterpol</td><td>Based on a linear interpolation within
a table of pair values. The references are the same as in 1).</td></tr>

<tr><td>4)</td><td>SPLINE1DArr</td><td>This is a Fortran 2003 subroutine
called SPLINE1DArr that takes an array of values x to interpolate from the
arrays t and y. It has ten input arguments, six output arguments, and two
optional arguments. The interpolation is linear.<p>

<ul>
<il>Forsythe, G.E. (1977) <em>Computer Methods For Mathematical Computations</em>. Ed. Prentice-Hall, Inc. <a href="https://onlinelibrary.wiley.com/doi/abs/10.1002/zamm.19790590235">https://onlinelibrary.wiley.com/doi/abs/10.1002/zamm.19790590235</a></il><p>

<il>De Boor, C. (1978). <em>A Practical Guide to Splines</em>. Springer-Verlag. ISBN: 978-0-387-95366-3. Link: <a href="https://www.springer.com/gp/book/9780387953663">https://www.springer.com/gp/book/9780387953663</a></il><p>

<il>Bartels, R.H., Beatty, J.C., and Barsky, B.A. (1998). <em>An Introductionto Splines for Use in Computer Graphics and Geometric Modeling</em>. Morgan
Kaufmann Publishers. ISBN: 978-1558604000. Link: <a href="https://www.amazon.com/Introduction-Computer-Graphics-Geometric-Modeling/dp/1558604006">https://www.amazon.com/Introduction-Computer-Graphics-Geometric-Modeling/dp/1558604006</a></il><p>

<il>Press, W. H., Teukolsky, S. A., Vetterling, W. T., &amp; Flannery, B. P. (2007). <em>Numerical Recipes: The Art of Scientific Computing</em> (3rd
ed.). Cambridge University Press. ISBN: 978-0521880688. Link: <a href="http://numerical.recipes/">http://numerical.recipes/</a></il><p>

<il>Cheney, W. and Kincaid, D. (2007). <em>Numerical Mathematics and Computing</em> (6th ed.). Brooks/Cole. ISBN: 978-0495114758. Link: <a href="https://www.amazon.com/Numerical-Mathematics-Computing-Ward-Cheney/dp/0495114758">https://www.amazon.com/Numerical-Mathematics-Computing-Ward-Cheney/dp/0495114758</a></il><p>

</td></tr>

<tr><td>5)</td><td>SPLINE3DFor</td><td>This function evaluates the cubic
spline interpolation. The same references as in 4).</td></tr>

</table>

<hr>

#### <b>STRUCTURE</b>

The main structure of the directories and files are:

<pre>
PyPhotometry
├── dist
│   ├── PyPhotometry-0.0.4.tar.gz
│   ├── PyPhotometry-0.0.7.tar.gz
│   ├── PyPhotometry-0.0.3.tar.gz
│   ├── PyPhotometry-0.0.6.tar.gz
│   ├── PyPhotometry-0.0.8.tar.gz
│   ├── PyPhotometry-0.0.9.tar.gz
│   ├── PyPhotometry-0.0.1.tar.gz
│   ├── PyPhotometry-0.0.2.tar.gz
│   ├── PyPhotometry-0.0.10.tar.gz
│   ├── PyPhotometry-0.0.11.tar.gz
│   └── PyPhotometry-0.0.5.tar.gz
├── README.md
├── .#README.md
├── showdown.min.js
├── scripts
│   └── update_readme.py
├── index.html
├── LICENSE.txt
├── setup.py
├── PyPhotometry
│   ├── calculate_sha256.py
│   ├── meta.yaml
│   └── README.txt
├── tutorials
│   ├── PyPhotometry.png
│   ├── .ipynb_checkpoints
│   │   └── PyPhotometry Example-checkpoint.ipynb
│   └── PyPhotometry Example.ipynb
├── PyPhotometry.egg-info
│   ├── PKG-INFO
│   ├── dependency_links.txt
│   ├── SOURCES.txt
│   ├── top_level.txt
│   └── requires.txt
├── src
│   ├── python
│   │   ├── PyFluxConSpec.py
│   │   ├── califa_cmap_alternative.py
│   │   ├── PyLinear__int.py
│   │   ├── PyLINinterpol.py
│   │   ├── fluxconserve.py
│   │   ├── __init__.py
│   │   ├── PySPLINE1DArr.py
│   │   ├── PySPLINE3DFor.py
│   │   ├── fluxconserving.png
│   │   ├── PyAkimaSpline.py
│   │   ├── PyInterpolado.py
│   │   └── PyLINdexerpol.py
│   └── fortran
│       ├── SPLINE3DFor.cpython-39-darwin.so
│       ├── LINdexerpol.cpython-311-darwin.so
│       ├── SPLINE3DFor.cpython-310-darwin.so
│       ├── LINdexerpol.f90
│       ├── LINinterpol.cpython-39-darwin.so
│       ├── SPLINE3DFor.cpython-39-x86_64-linux-gnu.so
│       ├── SPLINE1DFlt.cpython-39-x86_64-linux-gnu.so
│       ├── SPLINE3DFor.compile
│       ├── FluxConSpec.compile
│       ├── FluxConSpec.cpython-39-darwin.so
│       ├── SPLINE1DArr.compile
│       ├── SPLINE1DFlt.cpython-39-darwin.so
│       ├── AkimaSpline.cpython-310-darwin.so
│       ├── Interpolado.cpython-39-x86_64-linux-gnu.so
│       ├── Interpolado.cpython-38-x86_64-linux-gnu.so
│       ├── SPLINE1DFlt.f90
│       ├── FluxConSpec.cpython-39-x86_64-linux-gnu.so
│       ├── SPLINE1DFlt.cpython-310-darwin.so
│       ├── FluxConSpec.cpython-311-darwin.so
│       ├── AkimaSpline.f90
│       ├── SPLINE1DArr.f90
│       ├── SPLINE1DFlt.cpython-311-darwin.so
│       ├── FluxConSpec.cpython-38-x86_64-linux-gnu.so
│       ├── AkimaSpline.cpython-39-darwin.so
│       ├── SPLINE1DFlt.cpython-38-x86_64-linux-gnu.so
│       ├── SPLINE3DFor.cpython-38-x86_64-linux-gnu.so
│       ├── AkimaSpline.compile
│       ├── DataTypes.f90
│       ├── SPLINE1DArr.cpython-39-darwin.so
│       ├── LINdexerpol.cpython-310-darwin.so
│       ├── LINdexerpol.compile
│       ├── SPLINE3DFor.f90
│       ├── SPLINE1DFlt.compile
│       ├── Interpolado.f90
│       ├── LINinterpol.cpython-311-darwin.so
│       ├── SPLINE1DArr.cpython-310-darwin.so
│       ├── FluxConSpec.cpython-310-x86_64-linux-gnu.so
│       ├── LINinterpol.cpython-310-darwin.so
│       ├── SPLINE1DFlt.cpython-310-x86_64-linux-gnu.so
│       ├── LINinterpol.compile
│       ├── SPLINE1DArr.cpython-311-darwin.so
│       ├── SPLINE3DFor.cpython-311-darwin.so
│       ├── LINinterpol.compile~
│       ├── FluxConSpec.cpython-310-darwin.so
│       ├── FluxConSpec.f90
│       ├── Interpolado.cpython-310-x86_64-linux-gnu.so
│       ├── LINdexerpol.cpython-310-x86_64-linux-gnu.so
│       ├── README.txt
│       ├── Interpolado.compile
│       ├── SPLINE3DFor.cpython-310-x86_64-linux-gnu.so
│       ├── LINdexerpol.cpython-39-darwin.so
│       ├── LINinterpol.f90
│       ├── LINdexerpol.cpython-39-x86_64-linux-gnu.so
│       ├── LINdexerpol.cpython-38-x86_64-linux-gnu.so
│       └── AkimaSpline.cpython-311-darwin.so
├── version.txt
└── build
    ├── lib.linux-x86_64-3.9
    │   └── PyPhotometry
    ├── lib.macosx-11.1-arm64-cpython-39
    │   └── PyPhotometry
    ├── temp.macosx-11.1-arm64-cpython-39
    │   ├── __pycache__
    │   ├── PyPhotometry
    │   ├── ccompiler_opt_cache_ext.py
    │   ├── src
    │   └── build
    ├── src.linux-x86_64-3.9
    │   ├── PyPhotometry
    │   ├── build
    │   └── numpy
    ├── temp.linux-x86_64-3.9
    │   ├── __pycache__
    │   ├── PyPhotometry
    │   ├── ccompiler_opt_cache_ext.py
    │   ├── src
    │   └── build
    └── src.macosx-11.1-arm64-3.9
        ├── PyPhotometry
        ├── build
        └── numpy

32 directories, 99 files
</pre>

<br>PyPhotometry.py is a python wrapper to the library in fortran called
PyPhotometry.flib. The fortran directory can be compiled separately for
each individual subroutine.

<hr>

#### ISSUES AND CONTRIBUTIONS

If you encounter any issues with this project, please feel free to submit an issue on the GitHub repository. We appreciate your feedback and are committed to improving the quality of our codebase.

If you'd like to contribute to this project, we welcome pull requests from the community. Before submitting a pull request, please make sure to fork the repository and create a new branch for your changes. Once your changes are complete, submit a pull request and we'll review your code as soon as possible.

For any questions or concerns about contributing, please contact the project maintainer at antineutrinomuon@gmail.com. Thank you for your interest in contributing to our project!

<hr>

#### <b>LICENSE</b>

This software is provided "AS IS" (see DISCLAIMER below). Permission to
use, for non-commercial purposes is granted. Permission to modify for personal
or internal use is granted, provided this copyright and disclaimer are
included in ALL copies of the software. All other rights are reserved. In
particular, redistribution of the code is not allowed without explicit
permission by the author.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
