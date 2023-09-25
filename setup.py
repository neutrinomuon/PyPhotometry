from numpy.distutils.core import Extension
from numpy.distutils.core import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("version.txt", "r") as vh:
    version_description = vh.read()

# Function to get all files in a directory and its subdirectories
def get_files(directory):
    import os
    
    files = []
    for root, _, filenames in os.walk(directory):
        for filename in filenames:
            files.append(os.path.join(root, filename))
    return files

# Define the directories you want to include
data_directory = 'data'
calibration_directory = 'data/calibration_stars'
templates_directory = 'data/templates'

# Get all files within the data directory and its subdirectories
data_files = get_files(data_directory)

# Get all files within the calibration directory
calibration_files = get_files(calibration_directory)

# Get all files within the templates directory
template_files = get_files(templates_directory)

ext1 = Extension(  name='PyPhotometry.flib',
                   sources=['src/fortran/DataTypes.f90',
                            'src/fortran/LINinterpol.f90',
                            'src/fortran/GaussLegendreQuadrature.f90',
                            'src/fortran/IntegralALL.f90',
                            'src/fortran/PropFilters.f90',
                            'src/fortran/EvalFilters.f90'],
                 )
    
setup( name='PyPhotometry',
       version=version_description,
       ext_modules=[ ext1 ],
       extra_compile_args=['-O3'],
       description='PyPhotometry is a Python package based on a Fortran legacy package that allows you to compute photometric fluxes and magnitudes in various photometric systems.',
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
       package_dir={"PyPhotometry": "src/python"},
       packages=['PyPhotometry'],
       data_files=[('PyPhotometry/data', data_files), ('PyPhotometry/data/calibration_stars', calibration_files), ('PyPhotometry/data/templates', template_files)],
       include_package_data=True,
       #data_files=[ ('', ['version.txt','LICENSE.txt']),
       #             ('data', data_files),
       #             ('data/calibration', calibration_files) ],
       #package_data={'PyPhotometry': [
       #                               'data/*',
       #                               'data/calibration/*',
       #                              ],
       #             },
       #package_data={'PyPhotometry': [                                                                                                                                                                                               
       #                               'data/*',                                                                                                                                                                                      
       #                               'data/calibration/*',                                                                                                                                                                          
       #                              ],                                                                                                                                               
       #             },                          
       #include_package_data=True,
       #data_files=[('', ['version.txt', 'LICENSE.txt'])]
      )
    
