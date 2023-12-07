'''Setup toinstall PyPhotometry'''
# --------------------------------------------------------------------
# Your code has been rated at 10.00/10 (previous run: 10.00/10, +0.00)

# Import libraries used in the script
import os

from numpy.distutils.core import Extension
from numpy.distutils.core import setup

with open("README.md","r",encoding="utf-8") as fh:
    long_description = fh.read()

with open("version.txt","r",encoding="utf-8") as vh:
    version_description = vh.read()

# Function to get all files in a directory and its subdirectories
def get_files(directory):
    '''Find all files within a directory to add'''
    files = []
    for root, _, filenames in os.walk(directory):
        for filename in filenames:
            files.append(os.path.join(root, filename))
    return files

# Define the directories you want to include
DATA_DIRECTORY = 'data'
CALIBRATION_DIRECTORY = 'data/calibration_stars'
TEMPLATES_DIRECTORY = 'data/templates'

# Get all files within the data directory and its subdirectories
data_files = get_files(DATA_DIRECTORY)

# Get all files within the calibration directory
calibration_files = get_files(CALIBRATION_DIRECTORY)

# Get all files within the templates directory
template_files = get_files(TEMPLATES_DIRECTORY)

ext1 = Extension(  name='pyphotometry.flib',
                   sources=['src/fortran/DataTypes.f90',
                            'src/fortran/LINinterpol.f90',
                            'src/fortran/GaussLegendreQuadrature.f90',
                            'src/fortran/IntegralALL.f90',
                            'src/fortran/PropFilters.f90',
                            'src/fortran/EvalFilters.f90'],
                 )

setup( name='pyphotometry',
       version=version_description,
       ext_modules=[ ext1 ],
       extra_compile_args=['-O3'],
       description=('pyphotometry is a Python package based on a '
                    'Fortran legacy package that allows you to '
                    'compute photometric fluxes and magnitudes '
                    'in various photometric systems.'),
       long_description=long_description,      # Long description read from the the readme file
       long_description_content_type="text/markdown",
       author='Jean Gomes',
       author_email='antineutrinomuon@gmail.com',
       maintainer='Jean Gomes',
       maintainer_email='antineutrinomuon@gmail.com',
       keywords='photometry,stars,galaxies,magnitude,systems',
       url='https://github.com/neutrinomuon/PyPhotometry',
       docs_url='https://github.com/neutrinomuon/PyPhotometry',
       download_url='https://github.com/neutrinomuon/PyPhotometry',
       install_requires=[ 'numpy','matplotlib','sqlalchemy','astropy' ],
       requires_python='>=3.9',
       classifiers=[
           "Programming Language :: Python :: 3",
           "Programming Language :: Fortran",
           "Operating System :: OS Independent",
                   ],
       package_dir={"pyphotometry": "src/python"},
       packages=['pyphotometry'],
       data_files=[
           ('pyphotometry/data', data_files),
           ('pyphotometry/data/calibration_stars', calibration_files),
           ('pyphotometry/data/templates', template_files)
           ],
       include_package_data=True,
      )
