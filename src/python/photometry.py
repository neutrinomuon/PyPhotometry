#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Revised interface on Sat Jan 30 12:07:21 2021
Revised interface on Sex  1 Dez 2023 19:43:22 WET

@author: Jean Gomes

RESUME :  Filters

Version: v0.0.9

PYTHON : Python compatibility using f2py revised. Better usage  with numpy.  

Written: Jean Michel Gomes © Copyright
Created on Sat Jan 30 12:07:21 2021
"""

# Import necessary libraries
import os
import sys
import time
# sys.path.insert(1, 'Filters/')

# import sqlalchemy for the database
import sqlalchemy
from sqlalchemy import create_engine, exc, Column, String,  Float, Integer, PickleType
from sqlalchemy.pool import NullPool

# Deprecation of call, so changed - 25/09/2023
#from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import declarative_base

from sqlalchemy.orm import class_mapper, sessionmaker
from sqlalchemy import and_, or_
from sqlalchemy import select
from sqlalchemy.ext.automap import automap_base
from sqlalchemy import MetaData, Table

# import numpy
import numpy as np
import pkg_resources # Need to verify this package in the near future

# import graphical interface
import matplotlib as mpl
# from pylab import *
from pylab import subplot2grid, FormatStrFormatter
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter

# import astropy libraries
import astropy.io.fits as fits
import astropy.constants as const
import astropy.units as units

# import fortran legacy routines in python
# import PropFilters as prop
# import EvalFilters as evalf
from pyphotometry.flib import propfilters as prop
from pyphotometry.flib import evalfilters as evalf

# import pyphot to check the results
import pyphot
from pyphot import Sun, Vega, unit

global LSun
LSun = 3.839e33

# ! *** Open ListFilters ******************************************************
# !    arq_fil1 = file_dir(1:ilastnum)//'ListFilters.txt'
# !    open  (unit=21,file=arq_fil1,status='old',ERR=23)
# !    read  (21,*,ERR=23) Nfilters
# !
# !    if ( IsShowOn == 1_IB ) then
# !        write (*,'(4x,a)')   '[PropFilters]'
# !        write (W1aux,'(i15)') Nfilters
# !        write (*,'(4x,a,a)') '... Nfilters: ',trim(adjustl(W1aux))
# !    end if
# !
# !    do i_filter=1,Nfilters
# !        read  (21,*) axfilter
# !        name_fil(i_filter) = axfilter
# !
# !        arq_fil1 = file_dir(1:ilastnum)//axfilter
# !        open  (unit=22,file=arq_fil1,status='old',ERR=23)
# !        do j_filter=1,15
# !            read  (22,*,END=20,ERR=23)
# !        end do
# !        do j_filter=1,NVecsize
# !           read  (22,*,END=20,ERR=23) T_lambda(j_filter,i_filter),          &
# !                                      T_fluxes(j_filter,i_filter)
# !        end do
# !20      Numb_lbd(i_filter) = j_filter-1
# !        close (22)
# ! *** Open ListFilters ******************************************************

# sqlalchemy scheme ***************************************************************

# General Database Class
class create_database_filters():

    global Base
    Base = declarative_base()

    def __init__( self,database_path=None,database_filename='filters.db',verbose=False ):
        self.store_filters = 0
        self.readfilters = 0
        self.read_calibration_stars = 0

        # database_path by default is define as None, if so, look for installation of PyPhotometry
        # Specify the package name
        if database_path == None:
            package_name = 'pyphotometry'

            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)

            # Get the base directory path of the package
            database_path = package_dist.location \
                          + '/' \
                          + package_name \
                          + '/data/'
            if verbose:
                print("... Package directory:", database_path)

        self.database_path = database_path
        self.database_filename = database_filename
        self.database_filters = os.path.join(database_path, database_filename)

        # By using pkg_resources.resource_filename(__name__, 'filters.db'),
        # the code retrieves the absolute file path of 'filters.db' within
        # the package or module. This path is then assigned to
        # self.database_filters for further use, such as
        # creating an SQLAlchemy engine with the SQLite database located
        # at that path. However, it does not work if where you're running
        # it does not contain the filters.db file. So, a better approach
        # is from above.
        # self.database_filters = pkg_resources.resource_filename(__name__, 'filters.db')

        self.Engine = create_engine(
            'sqlite:///' + self.database_filters, 
            poolclass=NullPool,
            echo=False,
            )

        self.Session = sessionmaker(bind=self.Engine)
        self.session = self.Session()

    # Create the table information object class
    class _filters_class( Base ):
        """Information of Columns to create the filters SQL table"""

        __tablename__ = 'filters'

        filterid = Column(Integer, primary_key=True)
        name_filter = Column(String, primary_key=True)
        detector = Column(String)
        units = Column(String)

        wavelength = Column(PickleType)
        transmission = Column(PickleType)

        N_lambda = Column(Integer)
        t_l_area = Column(Float)
        t_n_area = Column(Float)

        lamb_eff = Column(Float)
        widtheff = Column(Float)
        vegaflux = Column(Float)

        magabsys = Column(Float)
        magtgsys = Column(Float)

    class _spectra_class( Base ):
        """Information of calibration stars to create the SQL table"""

        __tablename__ = 'spectra'

        specid = Column(Integer, primary_key=True)
        name_spec = Column(String, primary_key=True)
        flux_units = Column(String)
        lbd_units = Column(String)

        wavelength = Column(PickleType)
        fluxes = Column(PickleType)

        N_lambda = Column(Integer)

    def read_spectra(
            self,
            path_Vega=None,
            path1Vega=None,
            path_Sun=None,
            path1Sun=None,
            path2Sun=None,
            path_BD=None,
            verbose=False
            ):

        # Read Vega *******************************************************************
        if path_Vega == None:
            package_name = 'pyphotometry'

            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)

            # Get the base directory path of the package
            path_Vega = package_dist.location \
                      + '/' \
                      + package_name \
                      + '/data/calibration_stars/VegaLR.dat'
            if verbose:
                print("... path_Vega directory:", path_Vega)
        # Read Vega *******************************************************************

        # Read Vega ALTERNATIVE *******************************************************
        if path_Vega == None:
            package_name = 'pyphotometry'

            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)

            # Get the base directory path of the package
            path1Vega = package_dist.location \
                      + '/' \
                      + package_name + \
                      '/data/calibration_stars/Vega.dat'
            if verbose:
                print("... path1Vega directory:", path1Vega)
        # Read Vega ALTERNATIVE *******************************************************

        # Read Sun ********************************************************************
        if path_Sun == None:
            package_name = 'pyphotometry'

            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)

            # Get the base directory path of the package
            path_Sun = package_dist.location \
                     + '/' \
                     + package_name \
                     + '/data/calibration_stars/Sun_LR.dat'
            if verbose:
                print("... path_Sun directory:", path_Sun)
        # Read Sun ********************************************************************

        # Read Sun_1 ******************************************************************
        if path1Sun == None:
            package_name = 'pyphotometry'

            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)

            # Get the base directory path of the package
            path1Sun = package_dist.location \
                     + '/' \
                     + package_name \
                     + '/data/calibration_stars/Sun.dat'
            if verbose:
                print("... path1Sun directory:", path1Sun)
        # Read Sun_1 ******************************************************************

        # Read Sun_2 ******************************************************************
        if path2Sun == None:
            package_name = 'pyphotometry'

            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)

            # Get the base directory path of the package
            path2Sun = package_dist.location \
                     + '/' \
                     + package_name \
                     + '/data/calibration_stars/sun_reference_stis_001.fits'
            if verbose:
                print("... path2Sun directory:", path2Sun)
        # Read Sun_2 ******************************************************************

        # Read BD *********************************************************************
        if path_BD == None:
            package_name = 'pyphotometry'

            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)

            # Get the base directory path of the package
            path_BD = package_dist.location \
                    + '/' \
                    + package_name \
                    + '/data/calibration_stars/BD+17d4708.dat'
            if verbose:
                print("... path_BD directory:", path2Sun)
        # Read BD *********************************************************************

        o = Filters()
        o.ReadCalibrationStars(
            path_Vega=path_Vega,
            path1Vega=path1Vega,
            path_Sun=path_Sun,
            path1Sun=path1Sun,
            path2Sun=path2Sun,
            path_BD=path_BD
            )

        self.lambVega = o.lambVega
        self.fluxVega = o.fluxVega
        self.lambvega = self.lambVega
        self.fluxvega = self.fluxVega

        self.lambVega_1 = o.lambVega_1
        self.fluxVega_1 =o.fluxVega_1
        self.lambvega_1 = self.lambVega_1
        self.fluxvega_1 = self.fluxVega_1

        self.lamb_Sun = o.lamb_Sun
        self.flux_Sun = o.flux_Sun
        self.lamb_sun = self.lamb_Sun
        self.flux_sun  =self.flux_Sun

        self.lamb1Sun = o.lamb1Sun
        self.flux1Sun = o.flux1Sun
        self.lamb1sun = self.lamb1Sun
        self.flux1sun  =self.flux1Sun
        #self.lamb2sun = self.lamb2Sun
        #self.flux2sun  =self.flux2Sun

        self.lamb_FBD = o.lamb_FBD
        self.flux_FBD = o.flux_FBD
        self.lamb_fbd  = self.lamb_FBD
        self.flux_fbd = self.flux_FBD

        self.read_calibration_stars = 1

        if verbose:
            print("[read_spectra]")
            print("... Read spectra of calibration stars")
            print("[read_spectra]")

        return

    def read_filters( self,path_data=None,N_lambda=5000,verbose=False ):

        if verbose:
            print("[read_filters]")
            print("... Reading filters")

        # database_path by default is define as None, if so, look for installation of PyPhotometry
        # Specify the package name
        if path_data == None:
            package_name = 'pyphotometry'

            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)

            # Get the base directory path of the package
            path_data = package_dist.location \
                      + '/' \
                      + package_name \
                      + '/data/'
            if verbose:
                print("... path_data directory:", path_data)

        # Verify if database contains models for Photometry
        try:
            read = self.session.query(self._filters_class).filter(
                np.size(self._filters_class.name_filter) > 0 
                )
        except:
            read = 'ERROR' #np.zeros()

        # print('###########################################')
        # print("Filter Name")
        # for name_filter, wavelength in self.session.query(self._filters_class.name_filter,self._filters_class.wavelength):
        #     print(" {0:<20s} {1:<10d} ".format(name_filter,wavelength.size))

        # print('###########################################')
        
        # How to search
        #query = select( self._filters_class.name_filter ).filter( self._filters_class.name_filter == 'WISE1' ).all()
        #print( query )

        #stmt = self.session.query( self._filters_class.detector ).all()
        #print(stmt)
        #Base_table = automap_base()
        #Base_table.prepare( self.Engine, reflect=True )

        #Table = Base_table.classes.filters

        #stmt = select('*').select_from(self._filters_class)
        #result = self.session.execute(stmt).fetchall()

        try:
            if verbose:
                print("... Verified: {0:} entries in filters.db".format(np.count_nonzero(read.all())))
        except:
            if verbose:
                print("... Failed")

        if verbose:
            print("... Start reading filters database")

        try:
            if np.count_nonzero(read.all()) > 0:
                if verbose:
                    print("... Filters database already exists")

                self.store_filters = 1

                meta = MetaData()

                filter_table = Table(
                    'filters', meta, 
                    Column('filterid', Integer, primary_key=True), 
                    Column('name_filter', String, primary_key=True),
                    Column('detector',String),
                    Column('units',String),
                    Column('N_lambda',Integer),
                    Column('wavelength',PickleType),
                    Column('transmission',PickleType),
                    Column('lamb_eff',Float),
                    Column('widtheff',Float),
                    Column('vegaflux',Float),
                    Column('magabsys',Float),
                    Column('magtgsys',Float),
                    Column('t_l_area',Float),
                    Column('t_n_area',Float),
                    )

                conn = self.Engine.connect() 
                tabl = filter_table.select()
                result = conn.execute(tabl)

                filters_in_database = []
                for row in result:
                    filters_in_database.append(row[1])

                #filters_in_database = nparray( filters_in_database, dtype=object )
                N_filters_in_database =  np.size(filters_in_database)

                #print(filters_in_database)

                # Add filter if is not in database                
                path = path_data #'../../data/'
                # Just an ascii file
                arq_fil1 = 'ListFilters.txt'

                # Store name of filters
                o_list = open(path + arq_fil1,'r')
                orlist = o_list.readlines()
                name_list = []
                N_filters = int( orlist[0].split()[0] )
                #print(N_filters)

                index_filters_in_database = np.zeros([N_filters],dtype=int)
                for i in enumerate(orlist[1:N_filters+1]):
                    filter_string = i[1].split('.txt')[0]
                    name_list.append( filter_string )
                    for j in filters_in_database:
                        if j == filter_string:
                            #print(j,filter_string)
                            index_filters_in_database[i[0]] += 1

                name_list = np.array( name_list, dtype=object )
                index_filters_in_database = np.array( index_filters_in_database, dtype=int )

                #print( index_filters_in_database )
                #print( name_list[ index_filters_in_database == 0 ] )

                o = Filters()
                path = path_data
                arq_fil1 = 'ListFilters.txt'

                for j in name_list[ index_filters_in_database == 0 ]:
                    #print(j)
                    string_ = '%{ }%'.format(j)
                    # print(string_)
                    tabl = filter_table.select().where( filter_table.c.name_filter.ilike(string_) )
                    result = conn.execute(tabl)
                    # print( string_ )
                    # print( result )

                    count = 0
                    for row in result:
                        # print (row[1])
                        count += 1

                    # print( count )

                    if count <= 0:
                        #ReadONEFilter( self,path,arq_fil1,N_lambda=5000 )
                        path = path_data #'../../data/'
                        arq_fil1 = j + '.txt'
                        #print(path + arq_fil1)
                        o.ReadONEFilter( path=path, arq_fil1=arq_fil1 )

                        if o.read_calibration_stars == 1 and self.read_calibration_stars == 0:
                            self.lambVega = o.lambVega
                            self.fluxVega = o.fluxVega 
                            self.lambvega = o.lambVega
                            self.fluxvega = o.fluxVega

                            self.lamb_Sun = o.lamb_Sun
                            self.flux_Sun = o.flux_Sun 
                            self.lamb_sun = o.lamb_Sun
                            self.flux_sun = o.flux_Sun
                        
                            self.lamb1Sun = o.lamb1Sun
                            self.flux1Sun = o.flux1Sun
                            self.lamb1sun = o.lamb1Sun
                            self.flux1sun = o.flux1Sun

                            self.lamb_FBD = o.lamb_FBD
                            self.flux_FBD = o.flux_FBD
                            self.lamb_fbd = o.lamb_FBD
                            self.flux_fbd = o.flux_FBD
                            
                            self.read_calibration_stars = 1

                        d = o.onefilter[0]
                        filter_object = self._filters_class( name_filter=d[0], filterid=N_filters_in_database, detector=d[1], units=d[4]['units'], wavelength=np.array(d[2], dtype=float), transmission=d[3], N_lambda=int(d[4]['N_lambda']), t_l_area=d[4]['t_l_area'], t_n_area=d[4]['t_n_area'], vegaflux=d[4]['standard'], lamb_eff=d[4]['lamb_eff'], widtheff=d[4]['widtheff'], magabsys=d[4]['magabsys'], magtgsys=d[4]['magtgsys'] )

                        #print(filter_object)
                        self.session.add( filter_object )

                        try:
                            self.session.commit()
                            if verbose:
                                print("... The filter {0:} was included in the database filters.db.".format(d[0]))
                        except:
                            self.session.rollback()

                self.store_filters = 1                                

                # Read database
                s = filter_table.select()
                result = conn.execute(s)

                rows = self.session.query( filter_table ).count()
                self.N_filters = rows / 2
                self.Nfilters = self.N_filters

                filterid = np.zeros( [self.N_filters],dtype=int )
                t_l_area = np.zeros( [self.N_filters],dtype=float )
                t_n_area = np.zeros( [self.N_filters],dtype=float )
                magabsys = np.zeros( [self.N_filters],dtype=float )
                magtgsys = np.zeros( [self.N_filters],dtype=float )
                standard = np.zeros( [self.N_filters],dtype=float )
                numb_lbd = np.zeros( [self.N_filters],dtype=int )
                lamb_eff = np.zeros( [self.N_filters],dtype=float )
                widtheff = np.zeros( [self.N_filters],dtype=float )
                f__units = np.zeros( [self.N_filters],dtype=object )
                detector = np.zeros( [self.N_filters],dtype=object )
                name_filter = np.zeros( [self.N_filters],dtype=object )
                
                wavelength = np.zeros( [N_lambda, self.N_filters],dtype=float )
                transmission = np.zeros( [N_lambda, self.N_filters],dtype=float )

                self.filters = {}

                # Column('filterid', Integer, primary_key = True), 
                # Column('name_filter', String, primary_key=True),
                # Column('detector',String),
                # Column('units',String),
                # Column('N_lambda',Integer),
                # Column('wavelength',PickleType),
                # Column('transmission',PickleType),
                # Column('lamb_eff',Float),
                # Column('widtheff',Float),
                # Column('standard',Float),
                # Column('magabsys',Float),
                # Column('magtgsys',Float),
                
                for row in result:
                    #print(row) #.append()
                    filterid[ row[0] ] = row[0]
                    name_filter[ row[0] ] = row[1]
                    detector[ row[0] ] = row[2]
                    f__units[ row[0] ] = row[3]
                    numb_lbd[ row[0] ] = row[4]
                    
                    #print( row[5] )
                    
                    wavelength[ 0:np.size(row[5]), row[0] ] = row[5]
                    transmission[ 0:np.size(row[6]), row[0] ] = row[6]

                    lamb_eff[ row[0] ] = row[7]
                    widtheff[ row[0] ] = row[8]
                    standard[ row[0] ] = row[9]
                    magabsys[ row[0] ] = row[10]
                    magtgsys[ row[0] ] = row[11]

                    t_l_area[ row[0] ] = row[12]
                    t_n_area[ row[0] ] = row[13]

                    i = row[0]
                    Ntlambda = numb_lbd[ i ]
                    d = {}
                    d['N_lambda'] = Ntlambda
                    d['t_l_area'] = t_l_area[ i ]
                    d['t_n_area'] = t_n_area[ i ]
                    d['magabsys'] = magabsys[ i ]
                    d['magtgsys'] = magtgsys[ i ]
                    d['standard'] = standard[ i ]
                    d['lamb_eff'] = lamb_eff[ i ]
                    d['widtheff'] = widtheff[ i ]
                    d['units'] = f__units[ i ]

                    self.filters[ i ] = [ str(name_filter[i]),detector[i],wavelength[0:Ntlambda,i], transmission[0:Ntlambda,i], d ]
                    string_ =  str(name_filter[i]) + '.txt'
                    self.filters[ string_ ] = self.filters[ i ]
                    self.filters[ name_filter[i] ] = self.filters[ i ]

                self.session.close()
                
                self.T_lambda = wavelength
                self.T_fluxes = transmission
                self.Numblbd = numb_lbd
                self.numb_lbd = self.Numblbd
                self.name_fil = name_filter
                self.detector = detector
                self.t_l_area = t_l_area
                self.magabsys = magabsys
                self.magtgsys = magtgsys
                self.standard = standard
                self.lamb_eff = lamb_eff
                self.widtheff = widtheff

                self.readfilters = 1
                return
        
            else:
                raise ValueError("Force exception")
            
        except:
            if verbose:
                print("... Need to construct filters database")
        
            metadata = Base.metadata.create_all(self.Engine)
        
            o = Filters()
            #print(o)
            path = path_data #'../../data/'
            arq_fil1 = 'ListFilters.txt'
            IsKeepOn = o.ReadFilters( path,arq_fil1,verbose=verbose )

            if o.read_calibration_stars == 1 and self.read_calibration_stars == 0:
                self.lambVega = o.lambVega
                self.fluxVega = o.fluxVega 
                self.lambvega = o.lambVega
                self.fluxvega = o.fluxVega

                self.lamb_Sun = o.lamb_Sun
                self.flux_Sun = o.flux_Sun 
                self.lamb_sun = o.lamb_Sun
                self.flux_sun = o.flux_Sun
            
                self.lamb1Sun = o.lamb1Sun
                self.flux1Sun = o.flux1Sun
                self.lamb1sun = o.lamb1Sun
                self.flux1sun = o.flux1Sun

                self.lamb_FBD = o.lamb_FBD
                self.flux_FBD = o.flux_FBD
                self.lamb_fbd = o.lamb_FBD
                self.flux_fbd = o.flux_FBD

                self.read_calibration_stars = 1

            if IsKeepOn != 1:
                print("... Problem running filters database")
                return
        
            # 4 different entry for each filter
            N_entries = int( len(o.filters) / 4 )
            # print( len(o.filters) )
            if verbose:
                print("... N_entries: {0:}".format(N_entries))

            self.N_filters = N_entries
            self.Nfilters = self.N_filters
            
            # Add to SQL table
            # Create object
            self.T_lambda = np.zeros( [N_lambda,N_entries], dtype=float ) 
            self.T_fluxes = np.zeros( [N_lambda,N_entries], dtype=float )
            
            self.t_l_area =  np.zeros( [N_entries], dtype=float )
            self.magabsys = np.zeros( [N_entries], dtype=float )
            self.magtgsys =  np.zeros( [N_entries], dtype=float )
            self.standard =  np.zeros( [N_entries], dtype=float )

            self.Numblbd = np.zeros( [N_entries], dtype=int )
            self.name_fil = np.zeros( [N_entries], dtype=object )
            self.detector = np.zeros( [N_entries], dtype=object )

            self.lamb_eff = np.zeros( [N_entries], dtype=float )
            self.widtheff = np.zeros( [N_entries], dtype=float )

            for i in range(N_entries):
                #print(o.filters[i])
                d =  o.filters[i]
                # print( d )
                #print( d[0],d[1] )
                # filter_object = self._filters_class( name_filter=d[0], filterid=int(i), detector=d[1], units=d[4]['units'], wavelength=np.array(d[2], dtype=float), transmission=d[3], N_lambda=int(d[4]['N_lambda']), t_l_area=d[4]['t_l_area'], t_n_area=d[4]['t_n_area'], vegaflux=d[4]['standard'], lamb_eff=d[4]['lamb_eff'], widtheff=d[4]['widtheff'], magabsys=d[4]['magabsys'], magtgsys=d[4]['magtgsys'] )
               
                # Changed structure of flter dictionairy to include name of filter and not only name of file of the filter
                filter_object = self._filters_class( name_filter=d[0], filterid=int(i), detector=d[2], units=d[5]['units'], wavelength=np.array(d[3], dtype=float), transmission=d[4], N_lambda=int(d[5]['N_lambda']), t_l_area=d[5]['t_l_area'], t_n_area=d[5]['t_n_area'], vegaflux=d[5]['standard'], lamb_eff=d[5]['lamb_eff'], widtheff=d[5]['widtheff'], magabsys=d[5]['magabsys'], magtgsys=d[5]['magtgsys'] )

                # DEBUG
                # print(np.array(d[3], dtype=float)[0],np.array(d[3], dtype=float)[-1])
                # print(d[0],d[5]['lamb_eff'])   
                # print("lamb_eff ==> ",d[5]['lamb_eff'])

                #print(filter_object.name_filter)
                self.session.add( filter_object )
                
                # wave = np.array(d[2], dtype=float)
                # fluxes = d[3]
                wave = np.array(d[3], dtype=float)
                fluxes = d[4]
                j = wave.size
                self.T_lambda[ 0:j,i ] = wave[ 0:j ] 
                self.T_fluxes[ 0:j,i ] = fluxes[ 0:j ] 
                                
                # self.Numblbd [ i ] = int(d[4]['N_lambda'])
                # self.name_fil[ i ] = d[0]
                # self.detector[ i ] = d[1]
                # self.t_l_area[ i ]  = d[4]['t_l_area']
                # self.magabsys[ i ] = d[4]['magabsys']
                # self.magtgsys[ i ] = d[4]['magtgsys']
                # self.standard[ i ] = d[4]['standard']
                # self.lamb_eff[ i ] = d[4]['lamb_eff']
                # self.widtheff[ i ] =d[4]['widtheff']
                
                self.Numblbd [ i ] = int(d[5]['N_lambda'])
                self.name_fil[ i ] = d[0]
                self.detector[ i ] = d[2]
                self.t_l_area[ i ] = d[5]['t_l_area']
                self.magabsys[ i ] = d[5]['magabsys']
                self.magtgsys[ i ] = d[5]['magtgsys']
                self.standard[ i ] = d[5]['standard']
                self.lamb_eff[ i ] = d[5]['lamb_eff']
                self.widtheff[ i ] = d[5]['widtheff']
            
                # print('AQUI',self.name_fil[ i ],i,N_entries)    
            
            self.numb_lbd = self.Numblbd

            # Debug
            # print( 'PASSOU',N_entries )            
            #print(self.session.new)
            try:
                self.session.commit()
                if verbose:
                    print("... The filters were included in the database filters.db.")
                self.store_filters = 1
            
            except:
                self.session.rollback()
                if verbose:
                    print("... The filters are already included in the database filters.db.")

            self.filters = o.filters
            self.session.close()

            self.readfilters = 1

            return
    
    def FindFilters( self,word_search=None,verbose=False ):
    
        if self.readfilters == 0:
            print("... You need to first read the filters")
            return

        # print(self.filters.keys())

        # Filters in the database
        if word_search == None:
            #N_filters = self.N_filters
            d = list( self.filters.keys() )
            # print(d)
            if verbose:
                print("... List all filters - N_filters: {}".format(self.N_filters))
            
                for i in enumerate(d[::4]):
                    print( "    {0:>04d} {1:}".format(i[0],self.filters[ i[1] ][0]) )
            i = d[::4]
            r = d[3::4]
        else:               
            r,i = self.search( self.filters,word_search )
            #print( r )

        self.ind_names = r
        self.ind_filters = i

        # print(self.name_fil[i])
        return r,i

    def search( self,values,searchFor ):

        isfilterinlist = 0
        keys = list( values.keys() )
        # Debug
        # print(keys)

        if type( searchFor ) == tuple:
            searchFor = list( searchFor )
        # Otherwise transform in list if not list
        if type( searchFor ) != list:
            searchFor = list( [searchFor] )

        searchFor = list(searchFor)
        listsearched = []
        indsearched = []
        for j in searchFor:
            substring = str(j).lower()
            #print(substring)
            
            # l = 0
            for k in enumerate(keys[3::4]):
                fullstring = str(k[1]).lower()

                # search block 1
                keys_block_1 = keys[1::4][k[0]]
                fullstring_block_1 = str(keys_block_1).lower()
                # search block 2
                keys_block_2 = keys[2::4][k[0]]
                fullstring_block_2 = str(keys_block_2).lower()

                if substring in fullstring:
                    # print( fullstring,k )
                    # print( self.filters([ k[1] ]) )
                    # print(k)
                    # indsearched.append( k[0] )
                    indsearched.append( k[0] )
                    listsearched.append( k[1] )
                    # l += 1
                    isfilterinlist = 1    
                elif substring in fullstring_block_1:
                    indsearched.append( k[0] )
                    listsearched.append( k[1] )
                    isfilterinlist = 1
                elif substring in fullstring_block_2:
                    indsearched.append( k[0] )
                    listsearched.append( k[1] )
                    isfilterinlist = 1

        if isfilterinlist:
            indsearched = np.array( indsearched,dtype=int )
            # listsearched = listsearched
            return listsearched,indsearched
        else:
            return None, None
                    
    def plotfilter( self,lbd_norm=5000.,ind_filters=None,verbose=False ):

        if self.readfilters == 0:
            if verbose:
                print("[plotfilter]")
                print("... You need to first read the filters")
            return

        if ind_filters is not None:
            select_filters = np.array( ind_filters )
        else:
            select_filters = np.arange(0,self.N_filters,1) 

        o_lambdavega = self.lambvega
        o_fluxesvega = self.fluxvega
        o_lambda_sun = self.lamb_sun
        o_fluxes_sun = self.flux_sun
        o_lambda_fbd = self.lamb_fbd
        o_fluxes_fbd = self.flux_fbd

        N_filters = select_filters.size
        t_lambda = self.T_lambda[:,select_filters]
        t_fluxes = self.T_fluxes[:,select_filters]
        N_lambda = self.T_fluxes[:,0].size

        # Reshape
        t_lambda = t_lambda.reshape(N_lambda,N_filters)
        t_fluxes = t_fluxes.reshape(N_lambda,N_filters)
        # print(t_lambda.shape)

        # t_l_area = self.t_l_area[select_filters].reshape(N_filters)
        # magabsys = self.magabsys[select_filters].reshape(N_filters)
        # magtgsys = self.magtgsys[select_filters].reshape(N_filters)
        # standard = self.standard[select_filters].reshape(N_filters)
        lamb_eff = self.lamb_eff[select_filters].reshape(N_filters)
        numb_lbd = self.numb_lbd[select_filters].reshape(N_filters)
        detector = np.array(self.detector[select_filters]).reshape(N_filters)
        name_fil = np.array(self.name_fil[select_filters]).reshape(N_filters)

        # print(name_fil)

        norm_sun = np.interp( lbd_norm, self.lamb_sun, self.flux_sun )
        # norm1sun = np.interp( lbd_norm, self.lamb1sun, self.flux1sun )
        #norm2sun = np.interp( lbd_norm, self.lamb2sun, self.flux2sun )
        normvega = np.interp( lbd_norm, self.lambvega,self.fluxvega )
        norm_fbd = np.interp( lbd_norm, self.lamb_fbd,self.flux_fbd )

        # print(normvega)

        # int_type = 2
        # aux_flx = o_fluxes_fbd / norm_fbd
        # print(np.min(aux_flx[aux_flx>0.0]))
        # print(np.min(self.flux_sun[self.flux_sun>0.0]))
        # print(np.min(self.fluxvega[self.fluxvega>0.0]))

        # mag_spec,photfluxvega_,mcalibra,iskeepon = evalf( o_lambdavega,o_fluxesvega/normvega,t_lambda,t_fluxes,t_l_area,magabsys,magtgsys,standard,numb_lbd,int_type=int_type,verbosity=0 )
        # mag_spec,photflux_sun_,mcalibra,iskeepon = evalf( o_lambda_sun,o_fluxes_sun/norm_sun,t_lambda,t_fluxes,t_l_area,magabsys,magtgsys,standard,numb_lbd,int_type=int_type,verbosity=0 )
        # mag_spec,photflux_fbd_,mcalibra,iskeepon = evalf( o_lambda_fbd,o_fluxes_fbd/norm_fbd,t_lambda,t_fluxes,t_l_area,magabsys,magtgsys,standard,numb_lbd,int_type=int_type,verbosity=0 )

        photfluxvega = np.zeros( [N_filters], dtype=float )
        photflux_sun = np.zeros( [N_filters], dtype=float )
        photflux_fbd = np.zeros( [N_filters], dtype=float )
        for i in enumerate(range(N_filters)):
            # Debug
            # print(i[0],self.N_filters)
            # print(i,numb_lbd[i[0]],i[0])
            # print(t_lambda[1:numb_lbd[i[0]],i[0]])
            dtlambda = np.min( t_lambda[1:numb_lbd[i[0]],i[0]] - t_lambda[0:numb_lbd[i[0]]-1,i[0]] )
            
            flag = (o_lambda_fbd >= t_lambda[0,i[0]]-100.0) & (o_lambda_fbd <= t_lambda[numb_lbd[i[0]]-1,i[0]]+100.0)

            # Ensure at least ~20 points are selected
            j = 1
            while np.count_nonzero(flag) < 20 and j < 20:
                # Expand the range
                j += 1
                flag = (o_lambda_fbd >= t_lambda[0, i[0]] - j * 100.0) & (o_lambda_fbd <= t_lambda[numb_lbd[i[0]] - 1, i[0]] + j * 100.0)

            nolambda = o_lambda_fbd[flag].size
            dolambda = np.min( o_lambda_fbd[flag][1:nolambda] - o_lambda_fbd[flag][0:nolambda-1] )
        
            if dtlambda > dolambda:
                dtlambda = dolambda / 2.
            
            trlambda = np.arange( t_lambda[0, i[0]]-2.*dtlambda,t_lambda[numb_lbd[i[0]]-1, i[0]]+2.*dtlambda,dtlambda )
            trfluxes = np.interp( trlambda,t_lambda[0:numb_lbd[i[0]],i[0]],t_fluxes[0:numb_lbd[i[0]],i[0]] )
            trfluxes = np.where( trfluxes < 0.0, 0.0 ,trfluxes )
            
            orlambdavega = trlambda
            orlambda_sun = trlambda
            orlambda_fbd = trlambda
            orfluxesvega = np.interp( orlambdavega,o_lambdavega,o_fluxesvega )
            orfluxes_sun = np.interp( orlambda_sun,o_lambda_sun,o_fluxes_sun )
            orfluxes_fbd = np.interp( orlambda_fbd,o_lambda_fbd,o_fluxes_fbd )
            # ntlambda = trlambda.size
        
            # mag_spec,photflux,mcalibra,iskeepon = evalf( orlambda_fbd,orfluxes_fbd/norm_fbd,trlambda,trfluxes,t_l_area[i[0]],magabsys[i[0]],magtgsys[i[0]],standard[i[0]],ntlambda,int_type=int_type,verbosity=0 )
            photfluxvega[i[0]] = np.trapz(trfluxes*orfluxesvega) / np.trapz(trfluxes) / normvega
            photflux_sun[i[0]] = np.trapz(trfluxes*orfluxes_sun) / np.trapz(trfluxes) / norm_sun
            photflux_fbd[i[0]] = np.trapz(trfluxes*orfluxes_fbd) / np.trapz(trfluxes) / norm_fbd
            # photflux_fbd[i[0]] = photflux

        # print(photflux_vega,normvega)
        # print(photflux_sun,norm_sun)
        # print(photflux_fbd,norm_fbd)

        # Plot start here
        fig = plt.figure( figsize=(12,10),dpi=120,facecolor='w',edgecolor='w' )
        plt.subplots_adjust(bottom=.02, left=.06, right=.95, top=.98, wspace=0.0, hspace=0.0) 
        
        def define_size_cm(desired_size_cm):
            dpi = fig.get_dpi()
            if desired_size_cm < 0.0:
                desired_size_cm = 0.5
            # Convert the desired size to inches
            desired_size_inches = desired_size_cm / 2.54
            # Calculate the corresponding area in points^2
            desired_area = (desired_size_inches * dpi) ** 2
            return desired_area
        
        # fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
        def colorFader(c1,c2,mix=0):
            c1=np.array(mpl.colors.to_rgb(c1))
            c2=np.array(mpl.colors.to_rgb(c2))
            return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

        def colorFaderThree(c1, c2, c3, mix1=0, mix2=0, half=None):
            c1 = np.array(mpl.colors.to_rgb(c1))
            c2 = np.array(mpl.colors.to_rgb(c2))
            c3 = np.array(mpl.colors.to_rgb(c3))

            if half:                
                # Adjust the mixing parameters dynamically to emphasize green in the middle
                #mix1 = half / 2 if half < 0.5 else 0.5  # Increase mix1 for the first half
                #mix2 = 0.5 if half < 0.5 else (half - 0.5) / 0.5  # Increase mix2 for the second half
                mix = half
                if half < 0.5:                
                    final_color = (1-mix)*c1 + mix*c2
                else:
                    final_color = (1-mix)*c2 + mix*c3
            else:
                # Linear interpolation between c1 and c2
                color_intermediate = (1 - mix1) * c1 + mix1 * c2

                # Linear interpolation between c2 and c3
                final_color = (1 - mix2) * color_intermediate + mix2 * c3
    
            return mpl.colors.to_hex(final_color)

        c1='#1f77b4' # blue
        c2='green' # green
        c3='red' # red
        n=N_filters

        #ax = plt.subplot(111)

        # Top plot ###########################################################
        ax1_top = subplot2grid( (20,20), (0,0), colspan=20, rowspan=10 )
                                        
        # Sets the position and size of the panel for Plot #01
        # ax1_top.axis('on')
        # ax1_top.axes.get_xaxis().set_visible(False)
        # ax1_top.set_xscale('log')

        # ax1_top.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        # ax1_top.set_xticks( np.geomspace(1000, 42000 ,20).round() )
        # ax1_top.axis.set_minor_formatter(NullFormatter())

        # Scatter plots for photometry
        desired_size_cm = 0.3
        size_pt = define_size_cm(desired_size_cm)
        
        sorted_index = np.argsort(np.argsort(lamb_eff))
        colors = [] #np.array[ (N_filters), dtype=object ]
        for i_ind in range(N_filters):
            # Color fader with two or three colors
            # colors.append( colorFader(c1,c3,sorted_index[i]/n) )
            colors.append( colorFaderThree(c1,c2,c3,half=sorted_index[i_ind]/n) )

        # Debug color
        # print(colors)
        colors = np.array( colors, dtype=object )

        # Décalage
        decalage = [0.0,2.5,4.5]
        zorder = 0
        ax1_top.plot( np.log10(self.lambvega[ self.fluxvega>0.0 ]),np.log10(self.fluxvega[self.fluxvega>0.0]/normvega)-decalage[0],linewidth=5,color='darkgreen',label='Vega',zorder=zorder)
        # ax1_top.scatter( np.log10(self.lambvega[ self.fluxvega>0.0 ]),np.log10(self.fluxvega[self.fluxvega>0.0]/normvega),s=10,color='darkgreen',label='Vega')
        zorder += 1
        ax1_top.plot( np.log10(self.lamb_fbd[ self.flux_fbd>0.0 ]),np.log10(self.flux_fbd[self.flux_fbd>0.0]/norm_fbd)-decalage[1],linewidth=5,color='brown',label='BD+17d4708',zorder=zorder)
        # ax1_top.scatter( np.log10(self.lamb_fbd[ self.flux_fbd>0.0 ]),np.log10(self.flux_fbd[self.flux_fbd>0.0]/norm_fbd),color='brown',label='BD+17d4708')
        zorder += 1
        ax1_top.plot( np.log10(self.lamb_sun[ self.flux_sun>0.0 ]),np.log10(self.flux_sun[self.flux_sun>0.0]/norm_sun)-decalage[2],linewidth=5,color='darkorange',label='Sun',zorder=zorder)
        # ax1_top.scatter( np.log10(self.lamb_sun[ self.flux_sun>0.0 ]),np.log10(self.flux_sun[self.flux_sun>0.0]/norm_sun+1e-30),color='darkorange',label='Sun')

        # Photometric points - Photometric fluxes
        zorder += 1
        ax1_top.scatter( np.log10(lamb_eff[photfluxvega>0.0]),np.log10(photfluxvega[photfluxvega>0.0]-decalage[0]),s=size_pt*2,marker='s',color='grey',edgecolors='black',alpha=0.5,zorder=zorder )
        ax1_top.scatter( np.log10(lamb_eff[photfluxvega>0.0]),np.log10(photfluxvega[photfluxvega>0.0]-decalage[0]),s=size_pt,color=colors[ photfluxvega>0.0 ],alpha=1.0,zorder=zorder )
        zorder += 1
        ax1_top.scatter( np.log10(lamb_eff[photflux_fbd>0.0]),np.log10(photflux_fbd[photflux_fbd>0.0])-decalage[1],s=size_pt*2,marker='s',color='grey',edgecolors='black',alpha=0.5,zorder=zorder )
        ax1_top.scatter( np.log10(lamb_eff[photflux_fbd>0.0]),np.log10(photflux_fbd[photflux_fbd>0.0])-decalage[1],s=size_pt,color=colors,alpha=1.0,zorder=zorder )
        zorder += 1
        ax1_top.scatter( np.log10(lamb_eff[photflux_sun>0.0]),np.log10(photflux_sun[photflux_sun>0.0])-decalage[2],s=size_pt*2,marker='s',color='grey',edgecolors='black',alpha=0.5,zorder=zorder )
        ax1_top.scatter( np.log10(lamb_eff[photflux_sun>0.0]),np.log10(photflux_sun[photflux_sun>0.0])-decalage[2],s=size_pt,color=colors,alpha=1.0,zorder=zorder )

        # Limits for showing plot
        positive_values = t_lambda[np.where((t_lambda > 0) & (t_fluxes > 0))]

        min_xvec = np.round( np.log10( np.min(positive_values) ), 2 )
        max_xvec = np.round( np.log10( np.max(positive_values) ), 2 ) 
        
        j = 0
        delta_xvec = 0.0
        while delta_xvec <= 0.0 and j < 20:
            j += 1
            delta_xvec = np.round((max_xvec - min_xvec) / 12, j)
            # print(delta_xvec)

        if delta_xvec <= 0.0:
            delta_xvec = 1.0

        # min_xvec = np.round( np.log10( np.min(positive_values) ), 2 )
        # max_xvec = np.round( np.log10( np.max(positive_values) ), 2 ) 

        min_xvec = min_xvec-delta_xvec
        max_xvec = max_xvec+delta_xvec

        # Debug for checking min, max and delta xvec
        # print(min_xvec)
        # print(max_xvec)
        # print(delta_xvec)

        x_vec = np.arange( min_xvec,max_xvec,delta_xvec )
        ax1_top.xaxis.set_minor_formatter(FormatStrFormatter("%.2f"))

        ax1_top.set_xlim(x_vec[0],x_vec[-1])

        flag     = (np.log10(self.lambvega) >= x_vec[0]) & (np.log10(self.lambvega) <= x_vec[-1]) & (self.fluxvega > 0.0)
        flagvega = flag
        flag_sun = (np.log10(self.lamb_sun) >= x_vec[0]) & (np.log10(self.lamb_sun) <= x_vec[-1]) & (self.flux_sun > 0.0)
        flag_fbd = (np.log10(self.lamb_fbd) >= x_vec[0]) & (np.log10(self.lamb_fbd) <= x_vec[-1]) & (self.flux_fbd > 0.0)

        yvec = np.concatenate( [np.log10(self.fluxvega[flagvega]/normvega),
                                np.log10(self.flux_sun[flag_sun]/norm_sun),
                                np.log10(self.flux_fbd[flag_fbd]/norm_fbd)] )

        ymin = np.min( yvec ) - np.max(decalage)
        ymax = np.max( yvec ) + 1.0

        # Debug for checking min y
        # print(ymin,x_vec[0],x_vec[-1],np.log10(self.lambvega[0]),np.log10(self.lambvega[-1]))
        ax1_top.set_ylim(ymin,ymax)

        # Debug for number of major tickmarks
        # Ny = 8
        # Nx = 12
        # ax1_top.xaxis.set_major_locator(plt.MaxNLocator(Nx))
        # ax1_top.xaxis.set_minor_locator(plt.MaxNLocator(Nx*2))
        # ax1_top.yaxis.set_major_locator(plt.MaxNLocator(Ny))
        # ax1_top.yaxis.set_minor_locator(plt.MaxNLocator(Ny*2))

        for axis in [ax1_top.xaxis,ax1_top.yaxis]:
            axis.set_major_formatter(ScalarFormatter())
            axis.set_minor_formatter(NullFormatter())

        ax1_top.set_xticks(x_vec)
        #ax1_top.set_yticks(np.arange(0,1.1,0.1)) 

        ax1_top.legend(loc=1, prop={'size': 16}, title='Calibration Stars', title_fontsize='18')

        # Debug - not necessary for upper plot
        # ax1_top.set_xlabel("log Wavelength [$\AA$]", fontsize=16)
        ax1_top.set_ylabel(f"log F$_\lambda$ [Normalized at {lbd_norm:<3.0f} $\AA$]", fontsize=18)
        # Top plot ###########################################################

        # Bottom plot ########################################################
        ax1_bot = subplot2grid( (20,20), (11,0), colspan=20, rowspan=10 )                                           
        # Sets the position and size of the panel for Plot #01
        # ax1_top.axis('on')
        # ax1_top.axes.get_xaxis().set_visible(False)
        
        ax1_bot.set_xlim(x_vec[0],x_vec[-1])
        ax1_top.xaxis.set_minor_formatter(FormatStrFormatter("%.2f"))

        # Debug for major tickmarks
        # Ny = 8
        # Nx = 12
        # ax1_top.xaxis.set_major_locator(plt.MaxNLocator(Nx))
        # ax1_top.xaxis.set_minor_locator(plt.MaxNLocator(Nx*2))
        # ax1_top.yaxis.set_major_locator(plt.MaxNLocator(Ny))
        # ax1_top.yaxis.set_minor_locator(plt.MaxNLocator(Ny*2))

        for axis in [ax1_bot.xaxis,ax1_bot.yaxis]:
            axis.set_major_formatter(ScalarFormatter())
            axis.set_minor_formatter(NullFormatter())

        ax1_bot.set_xticks(x_vec)

        # Transmission curve goes from 0 to 1
        ax1_bot.set_yticks(np.arange(0,1.1,0.1))

        order_index = np.argsort(lamb_eff)

        # Debug - Checking the order of filters
        # sorted_index = np.argsort( order_index )
        # print(order_index)
        # print(self.name_fil[order_index])

        for i_ind in enumerate(t_lambda[0,0:N_filters]):
            j_ind = order_index[i_ind[0]]
            k_ind = order_index[j_ind]
            N_lambda = numb_lbd[j_ind]

            # Debug - Used in ColorFaderThree
            # Calculate mixing parameters based on the index and total number of filters
            # mix1 = i_ind[0] / n
            # mix2 = (i_ind[0] + 1) / (1.3*n)  # Assuming you want a smooth transition to the next color

            mix = i_ind[0] / n
            zorder = 100 + k_ind

            if detector[j_ind] == 'energy':
                ax1_bot.plot(np.log10(t_lambda[0:N_lambda,j_ind]),t_fluxes[0:N_lambda,j_ind]*t_lambda[0:N_lambda,j_ind],label=str(name_fil[j_ind]).split('.txt')[0],color=colorFaderThree(c1,c2,c3,half=mix),zorder=zorder)
                ax1_bot.fill_between(np.log10(t_lambda[0:N_lambda,j_ind]),t_fluxes[0:N_lambda,j_ind]*t_lambda[0:N_lambda,j_ind],color=colorFaderThree(c1,c2,c3,half=mix),alpha=0.5)
                
            else:
                ax1_bot.plot(np.log10(t_lambda[0:N_lambda,j_ind]),t_fluxes[0:N_lambda,j_ind],label=str(name_fil[j_ind]).split('.txt')[0],color=colorFaderThree(c1,c2,c3,half=mix),zorder=zorder)
                ax1_bot.fill_between(np.log10(t_lambda[0:N_lambda,j_ind]),t_fluxes[0:N_lambda,j_ind],color=colorFaderThree(c1,c2,c3,half=mix),alpha=0.5)

            #print(N_lambda,str(self.name_fil[i[0]]).split('.txt')[0])
            #print(self.T_fluxes[0:N_lambda,i[0]])

        # It is possible to handle the legends
        # handles, labels = ax1_bot.get_legend_handles_labels()
        # handles = np.array(handles)[sorted_index]
        # labels = np.array(labels)[sorted_index]
        # sorted_handles = [handles[i] for i in sorted_index]
        # sorted_labels = [labels[i] for i in sorted_index]

        ax1_bot.legend(loc=1, prop={'size': 10.5}, bbox_to_anchor=(1.30, 2.23), title='Filters', title_fontsize='16')
        ax1_bot.set_xlabel("log Wavelength [$\AA$]", fontsize=18)
        ax1_bot.set_ylabel("Transmission", fontsize=18)
        # Bottom plot ########################################################

        plt.show()

        if verbose:
            print("[plotfilter]")

        return

    def evaluate_photometry( self,wave,flux,ind_filters=None,plot=False  ):        

        if self.readfilters == 0:
            if verbose:
                print("[evaluate_photometry]")
                print("... You need to first read the filters")
            return

        if ind_filters is not None:
            select_filters = np.array( ind_filters )
        else:
            select_filters = np.arange(0,self.N_filters,1) 

        o_lambda = wave
        o_fluxes = flux

        #print(o_lambda.shape,o_fluxes.shape)
        N_filters = select_filters.size
        N_lambda = self.T_lambda[:,0].size
        t_lambda = self.T_lambda[:,select_filters].reshape(N_lambda,N_filters)
        t_fluxes = self.T_fluxes[:,select_filters].reshape(N_lambda,N_filters)
        
        # t_l_area = self.t_l_area[select_filters].reshape(N_filters)
        # magabsys = self.magabsys[select_filters].reshape(N_filters)
        # magtgsys = self.magtgsys[select_filters].reshape(N_filters)
        # standard = self.standard[select_filters].reshape(N_filters)
        numb_lbd = self.numb_lbd[select_filters].reshape(N_filters)

        # int_type = 2
        # mag_spec,photflux,mcalibra,iskeepon = evalfilters(o_lambda,o_fluxes,t_lambda,t_fluxes,t_l_area,magabsys,magtgsys,standard,numb_lbd,[int_type,verbosity])
        # mag_spec,photflux,mcalibra,iskeepon = evalf( o_lambda,o_fluxes,t_lambda,t_fluxes,t_l_area,magabsys,magtgsys,standard,numb_lbd,int_type=int_type,verbosity=0 )    
        # photflux_ = np.copy(photflux)
        # print("photflux all spectrum:", photflux)
        # print(numb_lbd)
        # print(standard.shape)

        # self.mag_spec = np.zeros([N_filters,7], dtype=float)
        self.photflux = np.zeros([N_filters], dtype=float)
        # self.lambmean = np.zeros([N_filters], dtype=float)
        photflux = self.photflux
        
        for i in enumerate(range(N_filters)):
            # Debug
            # for j in range(numb_lbd[i[0]]):
                # if t_lambda[j,i[0]] <= 0.0:
                    # print(i[0],'PROBLEM')
            flag = (o_lambda >= t_lambda[0,i[0]]-100.0) & (o_lambda <= t_lambda[numb_lbd[i[0]]-1,i[0]]+100.0)

            # Ensure at least ~20 points are selected
            j = 1
            while np.count_nonzero(flag) < 20 and j < 20:
                # Expand the range
                j += 1
                flag = (o_lambda >= t_lambda[0, i[0]] - j * 100.0) & (o_lambda <= t_lambda[numb_lbd[i[0]] - 1, i[0]] + j * 100.0)

            dtlambda = np.min(t_lambda[1:numb_lbd[i[0]]-2,i[0]] - t_lambda[0:numb_lbd[i[0]]-3,i[0]])
            nolambda = o_lambda[flag].size
            # if nolambda >= 5:
            dolambda = np.min(o_lambda[flag][1:nolambda-2] - o_lambda[flag][0:nolambda-3])
            # else:
                # dolambda = 1.0
                
            if dolambda < dtlambda:
                dtlambda = dolambda
                
            # if dtlambda <= dolambda:
            # dtlambda = dtlambda / 2.
            orlambda = np.arange(o_lambda[flag][0]-dtlambda,o_lambda[flag][-1]+dtlambda,dtlambda)
            orfluxes = np.interp(orlambda,o_lambda[flag],o_fluxes[flag])

            # print(t_lambda[0:numb_lbd[i[0]],i[0]].shape)
            # print(t_fluxes[0:numb_lbd[i[0]],i[0]].shape)

            # trlambda = orlambda
            trfluxes = np.interp(orlambda,t_lambda[0:numb_lbd[i[0]],i[0]],t_fluxes[0:numb_lbd[i[0]],i[0]])
            trfluxes = np.where( trfluxes < 0.0, 0.0, trfluxes )
            
            # self.lambmean[i[0]] = np.trapz(trfluxes*trlambda) / np.trapz(trfluxes)
            # else:
                # dolambda = dolambda / 2.

                # orlambda = np.arange(o_lambda[flag][0]-dolambda,o_lambda[flag][-1]+dolambda,dolambda)
                # orfluxes = np.interp(orlambda,o_lambda[flag],o_fluxes[flag])

                # trlambda = orlambda
                # trfluxes = np.interp(orlambda,t_lambda[0:numb_lbd[i[0]],i[0]],t_fluxes[0:numb_lbd[i[0]],i[0]])
                # trfluxes = np.where( trfluxes < 0.0, 0.0, trfluxes )

            photflux[i[0]] = np.trapz(trfluxes*orfluxes) / np.trapz(trfluxes)
            # nrlambda = orlambda.size

            # Already done when reading filter
            # if self.detector[ i[0] ] == 'energy':
                # plt.plot(trlambda,trfluxes*trlambda)
                # trfluxes = trfluxes / trlambda

            # print(t_lambda[0,i[0]]-100.0,t_lambda[numb_lbd[i[0]]-1,i[0]]+100.0)
            # print()
            # print("flag ======> ",flag,np.count_nonzero(flag[ flag==True ]))
            # count_true = np.sum(flag)
            # print(self.name_fil[i[0]],count_true)
            # print(t_lambda[0,i[0]],t_lambda[numb_lbd[i[0]]-1,i[0]])
            # print(dtlambda,dolambda)
            # if self.name_fil[i[0]] == '2MASS_H':
            #     # print(orlambda)
            #     # print(orfluxes)
            #     # print(trfluxes)
            #     print(t_l_area[i[0]],magabsys[i[0]],magtgsys[i[0]],standard[i[0]],numb_lbd[i[0]])
            #     print(nrlambda,trfluxes.size)
            # print(26000.0,40000.0,'WISE1.txt')
            # print(o_lambda[flag])
            # if np.count_nonzero(flag) > 5:
                # mag_spec,photflux,mcalibra,iskeepon = evalf( orlambda,orfluxes,t_lambda[0:numb_lbd[i[0]],i[0]],t_fluxes[0:numb_lbd[i[0]],[i[0]]],t_l_area[i[0]],magabsys[i[0]],magtgsys[i[0]],standard[i[0]],numb_lbd[i[0]],int_type=int_type,verbosity=0 )
                # mag_spec,photflux,mcalibra,iskeepon = evalf( orlambda,orfluxes,trlambda,trfluxes,t_l_area[i[0]],magabsys[i[0]],magtgsys[i[0]],standard[i[0]],nrlambda,int_type=int_type,verbosity=0 )
                # mag_spec = np.ones([1,7], dtype=float) * -999.0
                # photflux = np.ones(  [1], dtype=float) * -999.0
                # mcalibra = np.ones([1,7], dtype=float) * -999.0
                # iskeepon = 0
                # print("photon from interpolated: ",photflux[0])
                
                # print(np.trapz(trfluxes*orfluxes) / np.trapz(trfluxes))
                # photflux[0] = np.trapz(trfluxes*orfluxes) / np.trapz(trfluxes)
                # print(orfluxes)
                # print( np.trapz(trfluxes*orfluxes) )
            # else:
                # mag_spec = np.ones([1,7], dtype=float) * -999.0
                # photflux = np.ones(  [1], dtype=float) * -999.0
                # mcalibra = np.ones([1,7], dtype=float) * -999.0
                # iskeepon = 0

            # print(mag_spec.shape,photflux)

            # self.mag_spec[i[0],:] = mag_spec
            self.photflux[i[0]] = photflux[i[0]]

        # mag_spec,photflux,mcalibra,iskeepon = evalf( o_lambda,o_fluxes,t_lambda,t_fluxes,t_l_area,magabsys,magtgsys,standard,numb_lbd,int_type=int_type,lsun_con=LSun,verbosity=0 )

        #O_lambda,O_fluxes,NOlambda,T_lambda,T_fluxes,       &
         #               Ntlambda,Nfilters,MAG_spec,PhotFlux,Mcalibra,       &
          #              T_l_Area,magABsys,magTGsys,standard,Numb_lbd,       &
           #             IsKeepOn,Int_Type,LSun_con,verbosity
        #mag_spec,photflux,mcalibra,iskeepon = evalfilters(o_lambda,o_fluxes,t_lambda,t_fluxes,t_l_area,magabsys,magtgsys,standard,numb_lbd,[int_type,lsun_con,verbosity])

        #print(mag_spec)

        #print("")
        #print( iskeepon,mcalibra )
        #print( self.magtgsys )
        #print("")
        #print( self.magabsys )

        # self.mag_spec = mag_spec
        # self.photflux = photflux
        #self.iskeepon_photflux = iskeepon
        #self.lamb_eff =  self.lamb_eff

        #print(magabsys,mag_spec)
        # print(self.lamb_eff)

        if plot:
            # Adjust the figure size as needed
            fig, ax1 = plt.subplots(figsize=(6, 7.1))  

            lamb_eff = self.lamb_eff[select_filters].reshape(N_filters)
            detector = np.array(self.detector[select_filters]).reshape(N_filters)
            name_fil = np.array(self.name_fil[select_filters]).reshape(N_filters)
            
            zorder  = 0
            ind = np.argsort( lamb_eff )
            # ind = ind_filters[ind]

            #print( self.lamb_eff[ind].size )

            xlow = t_lambda[:,ind].flatten()
            xlow = xlow[ xlow>0.0 ].min() * 0.85
            xupp = t_lambda[:,ind].flatten().flatten()
            xupp = xupp[ xupp>0.0 ].max() * 1.15

            n = ind.size
            colors = pl.cm.jet(np.linspace(0,1,n))

            for i in enumerate(t_lambda[0,ind]):
                # print( ind[i[0]] , self.name_fil[ind[i[0]]], self.lamb_eff[ind[i[0]]], self.photflux[ind[i[0]]] )

                if detector[ ind[i[0]] ] == 'energy':
                    ax1.plot( t_lambda[0:numb_lbd[ind[i[0]]],ind[i[0]]],t_fluxes[0:numb_lbd[ind[i[0]]],ind[i[0]]]*t_lambda[0:numb_lbd[ind[i[0]]],ind[i[0]]],alpha=0.6,zorder=zorder,label=name_fil[ind[i[0]]],color=colors[i[0]] )
                else:
                    ax1.plot( t_lambda[0:numb_lbd[ind[i[0]]],ind[i[0]]],t_fluxes[0:numb_lbd[ind[i[0]]],ind[i[0]]],alpha=0.6,zorder=zorder,label=name_fil[ind[i[0]]],color=colors[i[0]] )

                zorder += 1
                if detector[ ind[i[0]] ] == 'energy':
                    ax1.fill_between(t_lambda[0:numb_lbd[ind[i[0]]],ind[i[0]]],t_fluxes[0:numb_lbd[ind[i[0]]],ind[i[0]]]*t_lambda[0:numb_lbd[ind[i[0]]],ind[i[0]]],color=colors[i[0]],alpha=0.4,zorder=zorder)
                else:
                    ax1.fill_between(t_lambda[0:numb_lbd[ind[i[0]]],ind[i[0]]],t_fluxes[0:numb_lbd[ind[i[0]]],ind[i[0]]],color=colors[i[0]],alpha=0.4,zorder=zorder)

                # plt.plot( self.T_lambda[0:self.Numblbd[ind[i[0]]],ind[i[0]]],self.T_fluxes[0:self.Numblbd[ind[i[0]]],ind[i[0]]]/max(self.T_fluxes[0:self.Numblbd[ind[i[0]]],ind[i[0]]])*max(o_fluxes[  (o_lambda >= xlow) & (o_lambda <= xupp) ]),alpha=0.6,zorder=zorder,label=self.name_fil[ind[i[0]]],color=colors[i[0]] )
                # zorder += 1
                # plt.fill_between(self.T_lambda[0:self.Numblbd[ind[i[0]]],ind[i[0]]],self.T_fluxes[0:self.Numblbd[ind[i[0]]],ind[i[0]]]/max(self.T_fluxes[0:self.Numblbd[ind[i[0]]],ind[i[0]]])*max(o_fluxes[  (o_lambda >= xlow) & (o_lambda <= xupp) ]),color=colors[i[0]],alpha=0.4,zorder=zorder)

            ax1.set_xscale('log')
            ax1.set_yscale('linear')
            ax1.set_ylabel('Transmission', color='blue')
            ax1.tick_params('y', colors='blue')
            ax1.set_ylim( 0.0,1.05 )
            ax1.set_xlabel(r"log Wavelength [$\AA$]")

            # plt.xlim(xlow,xupp)
            # plt.ylim( 0.0,1.2 )
            # plt.yscale('linear')

            # Create a second y-axis for the log-log curve
            ax2 = ax1.twinx()

            zorder += 1
            lbd =  lamb_eff[ind] # [ (self.lamb_eff >= xlow) & (self.lamb_eff <= xupp) ]
            jnd = np.argsort(lbd)

            lbd = lbd[jnd]
            pht = photflux[ind][jnd]

            #print(self.lamb_eff[ (self.lamb_eff >= xlow) & (self.lamb_eff <= xupp) ][ind]  )

            ax2.plot(lbd,pht,linewidth=2,color='black',zorder=zorder)
            # ax2.plot( np.log10(self.lamb_eff[ ind ][ jnd ]),np.log10(photflux[ ind ][ jnd ]),color='black',zorder=zorder)
            #print(self.lamb_eff[ind],self.name_fil[ind])

            zorder += 1
            # Debug - Vega spectrum
            # ax2.scatter(self.lambVega[ (self.lambVega >= xlow) & (self.lambVega <= xupp) ],self.fluxVega[  (self.lambVega >= xlow) & (self.lambVega <= xupp) ],color='magenta',zorder=zorder)

            ax2.plot(o_lambda[ (o_lambda >= xlow) & (o_lambda <= xupp) ],o_fluxes[  (o_lambda >= xlow) & (o_lambda <= xupp) ],linewidth=4,color='#1f77b4',zorder=zorder)
            # ax2.scatter(o_lambda[ (o_lambda >= xlow) & (o_lambda <= xupp) ],o_fluxes[  (o_lambda >= xlow) & (o_lambda <= xupp) ],linewidth=4,s=3,color='yellow',zorder=zorder)

            # ax2.scatter(o_lambda[ (o_lambda >= xlow) & (o_lambda <= xupp) ],o_fluxes[  (o_lambda >= xlow) & (o_lambda <= xupp) ],s=20,color='blue',zorder=zorder)

            zorder += 1
            
            dpi = fig.get_dpi()
            desired_size_cm = 0.6
            # Convert the desired size to inches
            desired_size_inches = desired_size_cm / 2.54
            # Calculate the corresponding area in points^2
            desired_area = (desired_size_inches * dpi) ** 2
            
            ax2.scatter(lamb_eff[ind],photflux[ind],color='black',s=desired_area,alpha=0.5,zorder=zorder)
            ax2.scatter(lamb_eff[ind],photflux[ind],color=colors,s=desired_area/4.,zorder=zorder)
            
            zorder += 1
            # ax2.scatter(self.lamb_eff[ind],photflux_[ind],color='black',s=20,zorder=zorder)

            ax2.set_xscale('log')
            ax2.set_yscale('log')

            # ax2.set_xlim(1e4,1e5)
            ax2.set_xlim(xlow,xupp)
            # ax2.set_ylim(1e-13,1.0e-3)

            legend = ax1.legend(loc=1, prop={'size': 6}, facecolor='white', fancybox=True, framealpha=0.4, edgecolor="black", title='Filters', bbox_to_anchor=(1.52, 1.02))
            legend.get_frame().set_alpha(None)
            legend.set_zorder( zorder + 102 ) 
            
            # ax2.set_xlabel("log Wavelength [$\AA$]")
            ax2.set_ylabel("F$_\lambda$ [L$_\odot \AA^{-1}$]" )
            
            plt.show()
            
    def evaluate_photometry_pyphot( self,l,f,plot=False  ):
        o_lambda = l #self.lamb_sun
        o_fluxes = f #self.flux_sun 
        
        # t_lambda = self.T_lambda
        # t_fluxes = self.T_fluxes
        
        # get the internal default library of passbands filters
        lib = pyphot.get_library()
        
        name_lib = [ 'GALEX_FUV',
                     'GALEX_NUV',
                     'SDSS_u',
                     'SDSS_g',
                     'SDSS_r',
                     'SDSS_i',
                     'SDSS_z',
                     '2MASS_H',
                     '2MASS_J',
                     '2MASS_Ks',
                     'WISE_RSR_W1',
                     'WISE_RSR_W2',
                     'WISE_RSR_W3',
                     'WISE_RSR_W4'     ]
        
        # Compute photometric flux with pyphot   
        n = len(name_lib) 
        lambd = np.zeros( [n], dtype=float )
        photflux = np.zeros( [n], dtype=float )
        mag_spec = np.zeros( [n], dtype=float )
        
        for i in range( n ):
            filters = lib[ name_lib[ i ] ]
            #print( filters.info() )
            
            # compute the integrated flux through the filter f
            # note that it work on many spectra at once
            val = filters.get_flux( l * unit['Angstrom'], f * unit['erg/s/cm**2/Angstrom'] , axis=1 )
            lbd = filters.leff.to("Angstrom")
            
            photflux[ i ] = val.value
            lambd[ i ] = lbd.value
            mag_spec[ i ] = -2.5 * np.log10( photflux[ i ] ) - filters.Vega_zero_mag
        
        self.mag_spec = mag_spec
        self.photflux = photflux
        self.lamb_eff = lambd
        #print(magabsys,mag_spec)
        
        if plot:
            xlow = self.T_lambda.flatten()
            xlow = xlow[ xlow>0.0 ].min() * 0.85
            xupp = self.T_lambda.flatten()
            xupp = xupp[ xupp>0.0 ].max() * 1.15
            
            zorder  = 0
            ind =  np.argsort(self.lamb_eff)
            n = ind.size
            colors  = pl.cm.jet(np.linspace(0,1,n))
            
            for i in enumerate(self.T_lambda[0,:]):
                plt.plot( self.T_lambda[0:self.Numblbd[ind[i[0]]],ind[i[0]]],self.T_fluxes[0:self.Numblbd[ind[i[0]]],ind[i[0]]]/max(self.T_fluxes[0:self.Numblbd[ind[i[0]]],ind[i[0]]])*max(o_fluxes[  (o_lambda >= xlow) & (o_lambda <= xupp) ]),alpha=0.6,zorder=zorder,label=self.name_fil[ind[i[0]]],color=colors[i[0]] )
                zorder += 1
                plt.fill_between(self.T_lambda[0:self.Numblbd[ind[i[0]]],ind[i[0]]],self.T_fluxes[0:self.Numblbd[ind[i[0]]],ind[i[0]]]/max(self.T_fluxes[0:self.Numblbd[ind[i[0]]],ind[i[0]]])*max(o_fluxes[  (o_lambda >= xlow) & (o_lambda <= xupp) ]),alpha=0.4,zorder=zorder)
                        
            zorder += 1
        
            lbd =  self.lamb_eff[ (self.lamb_eff >= xlow) & (self.lamb_eff <= xupp) ]
            ind = np.argsort(lbd)
        
            #print(self.lamb_eff[ (self.lamb_eff >= xlow) & (self.lamb_eff <= xupp) ][ind]  )
        
            plt.plot(self.lamb_eff[ (self.lamb_eff >= xlow) & (self.lamb_eff <= xupp) ][ind],photflux[ (self.lamb_eff >= xlow) & (self.lamb_eff <= xupp) ][ind],zorder=zorder)
            
            #print(self.lamb_eff,self.name_fil)
            
            zorder += 1
            plt.scatter(o_lambda[ (o_lambda >= xlow) & (o_lambda <= xupp) ],o_fluxes[  (o_lambda >= xlow) & (o_lambda <= xupp) ],color='blue',zorder=zorder)
            zorder += 1
            plt.scatter(self.lamb_eff[ (self.lamb_eff >= xlow) & (self.lamb_eff <= xupp) ] ,photflux[  (self.lamb_eff >= xlow) & (self.lamb_eff <= xupp) ],color='red',zorder=zorder)
              
            plt.xscale('log')
            plt.yscale('log')
    
            plt.xlim(xlow,xupp)
            # plt.xlim(1e4,1e5)

            plt.legend(loc=1, prop={'size': 7}, title='Filters')
            plt.xlabel("log Wavelength [$\AA$]")
            plt.ylabel("F$_\lambda$ [L$_\odot \AA^{-1}$]" )
            
            plt.show()

# path = '../../data/'
# arq_fil1 = 'ListFilters.txt'

# o.plotfilter()
# o.ReadFilters(path,arq_fil1)
# o.plotfilter()

# o.evaluate_photometry()

#print(o.T_lambda.shape)
#print(repr(o.name_fil))
        
# sqlalchemy scheme ***************************************************************

# Filters  **************************************************************************
class Filters():
    
    def __init__( self ):
        self.readfilters = 0
        self.filters = {}
        self.read_calibration_stars = 0

    def ReadFilters( self,path,arq_fil1,N_lambda=5000,verbose=False ):
        if verbose:
            print("[ReadFilters]")
            print( "... Reading Filters " )
            print( "... arq_fil1: {0:}".format(arq_fil1) ) #arq_fil1 = file_dir(1:ilastnum)//'ListFilters.txt'
        
        try:
            o_filter = open(path + arq_fil1,'r')
        except:
            print('... Filter List does not exist')
            o_filter.close()
            return -999
        
        r_filter = o_filter.readlines()
        Nfilters = int(r_filter[0].split()[0])
        self.Nfilters = Nfilters
        if verbose:
            print("... Nfilters: {0:d}".format(Nfilters))

        # Set filter arrays
        self.T_lambda = np.zeros([N_lambda,Nfilters], dtype=float)
        self.T_fluxes = np.zeros([N_lambda,Nfilters], dtype=float)

        #CH=30
        self.name_fil = np.zeros([Nfilters], dtype=object)
        #self.name_fil = np.empty((CH,Nfilters), dtype='c')
        self.Numblbd = np.zeros([Nfilters], dtype=int)

        self.detector_type = np.empty([Nfilters], dtype=object)
        self.units = np.empty([Nfilters], dtype=object)
        self.name_of_filter = np.empty([Nfilters], dtype=object)
        #print(self.name_fil.shape)

        # Test - Debug
        #print(r_filter)
        #print( len(r_filter) )
        
        # print("")
        for i in enumerate(r_filter[1:Nfilters+1]):
            i_split =str( i[1].split()[0] )
            #n = len(i_split)
            #self.name_fil[:,i[0]] = i_split + ' ' * ( CH - n  )
            self.name_fil[i[0]] = i_split
            # print("... name_fil: {0:004d} - {1:<10s}".format(i[0],i_split))
            
            o1filter = open(path + i_split,'r')
            r1filter = o1filter.readlines()
            
            # print(i_split)
            # print(self.name_fil[i[0]])
            # print(self.Numblbd[i[0]])
            for j in enumerate(r1filter):
                j_split = j[1].split()
                                
                if j_split[0][0] != '#':
                    k = self.Numblbd[i[0]]

                    self.T_lambda[k,i[0]] = j_split[0]
                    self.T_fluxes[k,i[0]] = j_split[1]
                    self.Numblbd[i[0]] += 1
                else:
                    #print(j_split)
                    if len(j_split) > 1:
                        det = j_split[1].split(':')[0].lower()

                        if det == 'detector':
                            if verbose:
                                print( "... Detector type: {0:}".format(j_split[2]) )
                            self.detector_type[i[0]] = j_split[2].lower()
                            
                        if det == 'units':
                            if verbose:
                                print( "... Units: {0:}".format(j_split[2]) )
                            self.units[i[0]] = j_split[2].lower()
                             
                        if det == 'name':
                            if verbose:
                                print( "... Name: {0:}".format( ''.join(j_split[2:])) )
                                print("... name_fil: {0:004d} - {1:<10s}".format(i[0],i_split))
                            self.name_of_filter[i[0]] = ''.join(j_split[2:])
                            
            o1filter.close()
            
            # Verify if repetitive lambdas or negative
            unique_lambdas, unique_indexes = np.unique(self.T_lambda[:,i[0]], return_index=True)
            
            # print(unique_lambdas)
            # print(unique_indexes)
            unique__fluxes = np.copy(self.T_fluxes[unique_indexes,i[0]])
            # print(unique__fluxes)
            
            # # Create a mask based on some conditions
            mask = (unique_lambdas > 0.0) & (unique__fluxes >= 0.0)

            # Apply the mask to get the desired unique elements and indexes
            filtered_unique_lambdas = unique_lambdas[mask]
            # filtered_unique_indexes = unique_indexes[mask]
            filtered_unique__fluxes = unique__fluxes[mask]

            self.Numblbd[i[0]] = np.size(filtered_unique_lambdas)

            self.T_lambda[:,i[0]] = 0.0
            self.T_fluxes[:,i[0]] = 0.0
            
            self.T_lambda[0:self.Numblbd[i[0]],i[0]] = filtered_unique_lambdas
            self.T_fluxes[0:self.Numblbd[i[0]],i[0]] = filtered_unique__fluxes

            #self.T_lambda[:,i[0]] = self.T_lambda[:,i[0]] * units.Angstrom

            # If detector type is energy then
            #a = np.trapz( ifT[ind] * _sflux, _slamb[ind], axis=axis )
            #b = np.trapz( ifT[ind], _slamb[ind])
            if self.detector_type[i[0]] == 'energy':
                # self.T_fluxes[0:self.Numblbd[i[0]],i[0]] /= self.T_lambda[0:self.Numblbd[i[0]],i[0]]
                self.T_fluxes[0:self.Numblbd[i[0]],i[0]] = np.where( self.T_lambda[0:self.Numblbd[i[0]],i[0]] > 0.0, self.T_fluxes[0:self.Numblbd[i[0]],i[0]]/self.T_lambda[0:self.Numblbd[i[0]],i[0]], 0.0 )

        o_filter.close()
        self.readfilters = 1
        
        #self.name_fil = np.array(self.name_fil, dtype='c')
        #print( self.name_fil )
        
        # Read Calibration Stars
        self.ReadCalibrationStars()
                
        # Compute properties of filter
        #file_dir = './'
        #int_type = 2
        #t_l_area,t_n_area,magabsys,magtgsys,standard,lamb_eff,widtheff,iskeepon = prop.propfilters( self.T_lambda,self.T_fluxes,self.Numblbd,file_dir,int_type=int_type,verbosity=0 )
        # t_l_area,t_n_area,magabsys,magtgsys,standard,lamb_eff,widtheff,iskeepon = prop.propfilters( self.T_lambda,self.T_fluxes,self.Numblbd,self.lambvega,self.fluxvega,self.lamb_sun,self.flux_sun,self.lamb_fbd,self.flux_fbd,verbosity=0 )

        self.numb_lbd = self.Numblbd        
        self.t_l_area = np.zeros([self.Numblbd.size], dtype=float)
        self.t_n_area = np.zeros([self.Numblbd.size], dtype=float)
        self.magabsys = np.zeros([self.Numblbd.size], dtype=float)
        self.magtgsys = np.zeros([self.Numblbd.size], dtype=float)
        self.standard = np.zeros([self.Numblbd.size], dtype=float)
        self.lamb_eff = np.zeros([self.Numblbd.size], dtype=float)
        self.widtheff = np.zeros([self.Numblbd.size], dtype=float)
        
        for i in enumerate(range(self.Nfilters)):
            # print(i)
            
            flag = (self.lambvega >= self.T_lambda[0,i[0]]) & (self.lambvega <= self.T_lambda[self.Numblbd[i[0]]-1,i[0]])
            
            dtlambda = np.min( self.T_lambda[1:self.Numblbd[i[0]],i[0]]-self.T_lambda[0:self.Numblbd[i[0]]-1,i[0]] )

            count_true = np.sum(flag)
            if count_true > 0:
                aux_lbd = np.copy( self.lambvega[flag] )
                # print(aux_lbd,self.T_lambda[0,i[0]],self.T_lambda[self.Numblbd[i[0]]-1,i[0]])
                d_lambda = np.min( aux_lbd[1:] - aux_lbd[0:-1] ) 
            else:
                d_lambda = dtlambda * 2.
                
            # print("dtlambda,d_lambda: ",dtlambda,d_lambda)
            
            if dtlambda <= d_lambda:
                aux_lbd = np.arange(self.T_lambda[0,i[0]]-2.*dtlambda,self.T_lambda[self.Numblbd[i[0]]-1,i[0]]+2.*dtlambda,dtlambda)
                aux_flx = np.interp(aux_lbd,self.lambvega,self.fluxvega)
            else:
                aux_lbd = np.arange(self.T_lambda[0,i[0]]-2.*d_lambda,self.T_lambda[self.Numblbd[i[0]]-1,i[0]]+2.*d_lambda,d_lambda)
                aux_flx = np.interp(aux_lbd,self.lambvega,self.fluxvega)
            
            aux_tra = np.interp(aux_lbd,self.T_lambda[0:self.Numblbd[i[0]],i[0]],self.T_fluxes[0:self.Numblbd[i[0]],i[0]])
            aux_tra = np.where( aux_tra < 0.0, 0.0, aux_tra )
            aux_flx = np.where( aux_flx < 0.0, 0.0, aux_flx )
            
            # Already done when reading filter
            # if self.detector_type[i[0]] == 'energy':
            #     aux_tra *= aux_lbd 
                
            # t_l_area,t_n_area,magabsys,magtgsys,standard,lamb_eff,widtheff,iskeepon = prop( self.T_lambda[0:self.Numblbd[i[0]],i[0]],self.T_fluxes[0:self.Numblbd[i[0]],i[0]],self.Numblbd[i[0]],self.lambvega,self.fluxvega,self.lamb_sun,self.flux_sun,self.lamb_fbd,self.flux_fbd,verbosity=0 )
            # t_l_area,t_n_area,magabsys,magtgsys,standard,lamb_eff,widtheff,iskeepon = prop( self.T_lambda[0:self.Numblbd[i[0]],i[0]],self.T_fluxes[0:self.Numblbd[i[0]],i[0]],self.Numblbd[i[0]],aux_lbd,aux_flx,self.lamb_sun,self.flux_sun,self.lamb_fbd,self.flux_fbd,verbosity=0 )
            num_tra = aux_tra.size
            t_l_area,t_n_area,magabsys,magtgsys,standard,lamb_eff,widtheff,iskeepon = prop( aux_lbd,aux_tra,num_tra,aux_lbd,aux_flx,self.lamb_sun,self.flux_sun,self.lamb_fbd,self.flux_fbd,verbosity=0 )

            # print(i,lamb_eff)
                        
            self.t_l_area[i[0]] = t_l_area[0]
            self.t_n_area[i[0]] = t_n_area[0]
            self.magabsys[i[0]] = magabsys[0]
            self.magtgsys[i[0]] = magtgsys[0]
            self.standard[i[0]] = standard[0]
            self.lamb_eff[i[0]] = lamb_eff[0]
            self.widtheff[i[0]] = widtheff[0]
 
        # t_l_area,t_n_area,magabsys,magtgsys,standard,lamb_eff,widtheff,iskeepon = prop( self.T_lambda,self.T_fluxes,self.Numblbd,self.lambvega,self.fluxvega,self.lamb_sun,self.flux_sun,self.lamb_fbd,self.flux_fbd,verbosity=0 )

        # print(iskeepon)
        # print( self.lambVega[0],self.lambVega[-1] )
        # print('lamb_eff ========>', lamb_eff)
        # print('widtheff ========>', widtheff)
        # print('standard ========>', standard)

        # self.t_l_area = t_l_area
        # self.t_n_area = t_n_area
        # self.magabsys = magabsys
        # self.magtgsys = magtgsys
        # self.standard = standard
        # self.numb_lbd = self.Numblbd
        # self.lamb_eff = lamb_eff
        # self.widtheff = widtheff

        for i in range(Nfilters):
            N_lambda = self.Numblbd[i]
            d = {}
            d['N_lambda'] = N_lambda
            d['t_l_area'] = self.t_l_area[i]
            d['t_n_area'] = self.t_n_area[i]
            d['magabsys'] = self.magabsys[i]
            d['magtgsys'] = self.magtgsys[i]
            d['standard'] = self.standard[i]
            d['lamb_eff'] = self.lamb_eff[i]
            d['widtheff'] = self.widtheff[i]
            d['units'] = self.units[i]
            
            # Added
            # self.name_of_filter[ i ]

            self.filters[ i ] = [ str(self.name_of_filter[i]),str(self.name_fil[i]).split('.txt')[0],self.detector_type[i],self.T_lambda[0:N_lambda,i], self.T_fluxes[0:N_lambda,i], d ]
            self.filters[ self.name_fil[i] ] = self.filters[ i ] #self.T_lambda[0:N_lambda,i], self.T_fluxes[0:N_lambda,i]
            self.filters[ str(self.name_fil[i]).split('.txt')[0] ] = self.filters[ i ]
            self.filters[ str(self.name_of_filter[i]) ] = self.filters[ i ]

            # Ok, potentially a problem if self.name_of_filter == self.name_fill.split('txt')[0]
            # fill_character = ' '  # Specify the character you want to use for filling

            # print( "{0:0>4} {1:<30} {2:<30} {3:<30}".format( i, self.name_fil[i].ljust(40, fill_character), str(self.name_of_filter[i]) , str(self.name_fil[i]).split('.txt')[0] ) )

        # print( "N_filters: ----- DEBUG: ",len(self.filters) )
        if verbose:
            print("[ReadFilters]")

        return 1
    
    def ReadONEFilter( self,path,arq_fil1,N_lambda=5000,verbose=False ):
        if verbose:
            print("[ReadONEFilter]")
            print( "... Reading One Filter " )
            print( "... arq_fil1: {0:}".format(arq_fil1) )
                
        Nfilters =1
        if verbose:
            print("... Nfilters: {0:d}".format(Nfilters))
                
        # Set filter arrays
        self.T1lambda = np.zeros([N_lambda,Nfilters], dtype=float)
        self.T1fluxes = np.zeros([N_lambda,Nfilters], dtype=float)

        self.name1fil = np.zeros([Nfilters], dtype=object)
        self.Numblbd1 = np.zeros([Nfilters], dtype=int)
        
        self.detector1type = np.empty([Nfilters], dtype=object)
        self.units1 = np.empty([Nfilters], dtype=object)
                
        print("")
        self.name1fil[0] = arq_fil1.split('.txt')[0]
        if verbose:
            print("... name1fil: {0:004d} - {1:<10s}".format(0,arq_fil1))
                    
        o1filter = open(path + arq_fil1,'r')
        r1filter = o1filter.readlines()
                    
        for j in enumerate(r1filter):
            j_split = j[1].split()
            
            if j_split[0][0] != '#':
                k = self.Numblbd1[0]

                self.T1lambda[k,0] = j_split[0]
                self.T1fluxes[k,0] = j_split[1]
                self.Numblbd1[0] += 1
            else:
                #print(j_split)
                if len(j_split) > 1:
                    det = j_split[1].split(':')[0].lower()

                    if det == 'detector':
                        if verbose:
                            print("... Detector type: {0:}".format(j_split[2]))
                        self.detector1type[0] = j_split[2].lower()
                            
                    if det == 'units':
                        if verbose:
                            print("... Units: {0:}".format(j_split[2]))
                        self.units1[0] = j_split[2].lower()

        o1filter.close()

        #self.T_lambda[:,i[0]] = self.T_lambda[:,i[0]] * units.Angstrom

        # If detector type is energy then
        #a = np.trapz( ifT[ind] * _sflux, _slamb[ind], axis=axis )
        #b = np.trapz( ifT[ind], _slamb[ind])
        if self.detector1type[0] == 'energy':
            self.T1fluxes[0:self.Numblbd1[0],0] /= self.T1lambda[0:self.Numblbd1[0],0]

        self.readonefilter = 1

        # Read Calibration Stars
        self.ReadCalibrationStars()

        # Compute properties of one filter
        #print( self.T1lambda,self.T1fluxes,self.Numblbd1 )
        # t_l_area,t_n_area,magabsys,magtgsys,standard,lamb_eff,widtheff,iskeepon = prop.propfilters( self.T1lambda,self.T1fluxes,self.Numblbd1,self.lambvega,self.fluxvega,self.lamb_sun,self.flux_sun,self.lamb_fbd,self.flux_fbd,verbosity=0 )
        t_l_area,t_n_area,magabsys,magtgsys,standard,lamb_eff,widtheff,iskeepon = prop( self.T1lambda,self.T1fluxes,self.Numblbd1,self.lambvega,self.fluxvega,self.lamb_sun,self.flux_sun,self.lamb_fbd,self.flux_fbd,verbosity=0 )

        # print("lamb_eff ==> ",lamb_eff)
    
        # self.t_l_area = t_l_area
        # self.t_n_area = t_n_area
        # self.magabsys = magabsys
        # self.magtgsys = magtgsys
        # self.standard = standard
        # self.numb_lbd = self.Numblbd
        # self.lamb_eff = lamb_eff

        self.onefilter = {}

        for i in range(Nfilters):
            N_lambda = self.Numblbd1[i]
            d = {}
            d['N_lambda'] = N_lambda
            d['t_l_area'] = t_l_area[ i ]
            d['t_n_area'] = t_n_area[ i ]
            d['magabsys'] = magabsys[ i ]
            d['magtgsys'] = magtgsys[ i ]
            d['standard'] = standard[ i ]
            d['lamb_eff'] = lamb_eff[ i ]
            d['widtheff'] = widtheff[ i ]
            d['units'] = self.units1[ i ]
            
            self.onefilter[ i ] = [ str(self.name1fil[i]).split('.txt')[0],self.detector1type[i],self.T1lambda[0:N_lambda,i], self.T1fluxes[0:N_lambda,i], d ]
            self.onefilter[ str(self.name1fil[i]) + '.txt' ] = self.onefilter[ i ]
            self.onefilter[ self.name1fil[i] ] = self.onefilter[ i ]

            #print( self.onefilter )
            return
            
    def ReadCalibrationStars( self,path_Vega=None,path1Vega=None,path_Sun=None,path1Sun=None,path2Sun=None,path_BD=None,verbose=False ):

        if self.read_calibration_stars == 1:
            return
        
# ! *** Read calibration stars ************************************************
# !     RESUME : VEGA spectrum.                                               !
# !              Intrinsic Flux: erg/s/cm2/A                                  !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     !arq_fil1 = file_dir(1:ilastnum)//'VegaLR.dat'
#     !open  (21,status='old',file=arq_fil1,ERR=22)
#     !read  (21,*,ERR=22) arq_lixo,NVegalbd
#     !
#     !allocate ( lambVega(NVegalbd) )
#     !allocate ( fluxVega(NVegalbd) )
#     !
#     !do i_lambda=1,NVegalbd
#     !    read  (21,*,ERR=22) lambVega(i_lambda),fluxVega(i_lambda)
#     !end do
#     !close (21)

# ! *** SUN spectrum **********************************************************
# !     Intrinsic flux: erg/s/A.                                              !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     !arq_fil1 = file_dir(1:ilastnum)//'Sun_LR.dat'
#     !open(21,status='old',file=arq_fil1,ERR=22)
#     !read(21,*,ERR=22) arq_lixo,NSun_lbd
#     !
#     !allocate ( lamb_Sun(NSun_lbd) )
#     !allocate ( flux_Sun(NSun_lbd) )
#     !
#     !do i_lambda=1,NSun_lbd
#     !    read(21,*,ERR=22) lamb_Sun(i_lambda),flux_Sun(i_lambda)
#     !end do
#     !close(21)

# ! *** Reading of the spectrum of BD+17o4708 *********************************
# !     RESUME : F subdwarf used to calibrate the Thuan & Gunn system.        !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     !arq_fil1 = file_dir(1:ilastnum)//'BD+17o4708.dat'
#     !arq_fil1 = file_dir(1:ilastnum)//'BD+17d4708.dat'
#     !open(21,status='old',file=arq_fil1,ERR=22)
#     !read(21,*,ERR=22) arq_lixo,NFBD_lbd
#     !
#     !allocate ( lamb_FBD(NFBD_lbd) )
#     !allocate ( flux_FBD(NFBD_lbd) )
#     !
#     !do i_lambda=1,NFBD_lbd
#     !    read(21,*,ERR=22) lamb_FBD(i_lambda),flux_FBD(i_lambda)
#     !end do
#     !close(21)
# ! *** Read calibration stars ************************************************
        
# https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec

# Read Vega *******************************************************************
        if path_Vega == None:
            package_name = 'pyphotometry'
        
            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)
        
            # Get the base directory path of the package
            path_Vega = package_dist.location \
                      + '/' \
                      + package_name \
                      + '/data/calibration_stars/VegaLR.dat'
            if verbose:
                print("... path_Vega directory:", path_Vega)
# Read Vega *******************************************************************

# Read Vega ALTERNATIVE *******************************************************
        if path1Vega == None:
            package_name = 'pyphotometry'
        
            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)
        
            # Get the base directory path of the package
            path1Vega = package_dist.location \
                      + '/' \
                      + package_name \
                      + '/data/calibration_stars/Vega.dat'
            if verbose:
                print("... path1Vega directory:", path1Vega)
# Read Vega ALTERNATIVE *******************************************************

# Read Sun ********************************************************************
        if path_Sun == None:
            package_name = 'pyphotometry'
        
            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)
        
            # Get the base directory path of the package
            path_Sun = package_dist.location + '/' + package_name + '/data/calibration_stars/Sun_LR.dat'
            if verbose:
                print("... path_Sun directory:", path_Sun)
# Read Sun ********************************************************************

# Read Sun_1 ******************************************************************
        if path1Sun == None:
            package_name = 'pyphotometry'
        
            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)
        
            # Get the base directory path of the package
            path1Sun = package_dist.location + '/' + package_name + '/data/calibration_stars/Sun.dat'
            if verbose:
                print("... path1Sun directory:", path1Sun)
# Read Sun_1 ******************************************************************

# Read Sun_2 ******************************************************************
        if path2Sun == None:
            package_name = 'pyphotometry'
        
            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)
        
            # Get the base directory path of the package
            path2Sun = package_dist.location \
                     + '/' \
                     + package_name \
                     + '/data/calibration_stars/sun_reference_stis_001.fits'
            if verbose:
                print("... path2Sun directory:", path2Sun)
# Read Sun_2 ******************************************************************

# Read BD *********************************************************************
        if path_BD == None:
            package_name = 'pyphotometry'
        
            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)
        
            # Get the base directory path of the package
            path_BD = package_dist.location \
                    + '/' \
                    + package_name \
                    + '/data/calibration_stars/BD+17d4708.dat'
            if verbose:
                print("... path_BD directory:", path2Sun)
# Read BD *********************************************************************

#  *** VEGA spectrum ***************************************************************
#         Intrinsic Flux: erg/s/cm2/A                                                                                                                  !
######################################################################
        # start_time = time.time()
        # o = open(path_Vega,'r')
        # r = o.readlines()
        # r_split  = r[0].split()[1]
        
        # self.NVegalbd = int(r_split)
        # if verbose:
        #     print("... NVegalbd: {0:}".format(self.NVegalbd))
        # self.lambVega = np.zeros([self.NVegalbd], dtype=float)
        # self.fluxVega = np.zeros([self.NVegalbd], dtype=float)
        
        # for i in enumerate(r[1:]):
        #     r_split = i[1].split()
        #     self.lambVega[i[0]] = r_split[0]
        #     self.fluxVega[i[0]] = r_split[1]
        # o.close()
        # end_time = time.time()
        # print(end_time-start_time)

        # start_time = time.time()
        data = np.loadtxt(path_Vega,skiprows=0)
        # print(data.shape)
        self.lambVega = np.copy(data[:,0])
        self.fluxVega = np.copy(data[:,1])
        self.NVegalbd = self.lambVega.size
        # end_time = time.time()
        # print(end_time-start_time)
#  *** VEGA spectrum ***************************************************************

#  *** VEGA spectrum - ALTERNATIVE *************************************************
#         Intrinsic Flux: erg/s/cm2/A                                                                                                                  !
#  *** VEGA spectrum - ALTERNATIVE *************************************************
        # start_time = time.time()
        # o = open(path1Vega,'r')
        # r = o.readlines()
        # r_split  = r[0].split()[1]
        
        # print(r_split)
        
        # self.NVegalbd_1 = int(r_split)
        # if verbose:
            # print("... NVegalbd_1: {0:}".format(self.NVegalbd_1))
        # self.lambVega_1 = np.zeros([self.NVegalbd_1], dtype=float)
        # self.fluxVega_1 = np.zeros([self.NVegalbd_1], dtype=float)
        
        # for i in enumerate(r[1:]):
            # r_split = i[1].split()
            # self.lambVega_1[i[0]] = r_split[0]
            # self.fluxVega_1[i[0]] = r_split[1]
        # o.close()
        # end_time = time.time()
        # print(end_time-start_time)

        # start_time = time.time()
        data = np.loadtxt(path1Vega,skiprows=0)

        self.lambVega_1 = np.copy(data[:,0])
        self.fluxVega_1 = np.copy(data[:,1])
        self.NVegalbd_1 = self.lambVega_1.size
        # end_time = time.time()
        # print(end_time-start_time)
#  *** VEGA spectrum ***************************************************************

#  *** SUN spectrum ****************************************************************
#         Intrinsic Flux: erg/s/cm2/A                                                                                                                  !
######################################################################     
        # o = open(path_Sun,'r')
        # r = o.readlines()
        # r_split  = r[0].split()[1]
        
        # self.NSun_lbd = int(r_split)
        # if verbose:
        #     print("... NSun_lbd: {0:}".format(self.NSun_lbd))
        # self.lamb_Sun = np.zeros([self.NSun_lbd], dtype=float)
        # self.flux_Sun = np.zeros([self.NSun_lbd], dtype=float)
        
        # for i in enumerate(r[1:]):
        #     r_split = i[1].split()
        #     self.lamb_Sun[i[0]] = r_split[0]
        #     self.flux_Sun[i[0]] = r_split[1]
        # o.close()
        
        # start_time = time.time()
        data = np.loadtxt(path1Sun,skiprows=0)

        self.lamb_Sun = np.copy(data[:,0])
        self.flux_Sun = np.copy(data[:,1])
        self.NSun_lbd = self.lamb_Sun.size
        # end_time = time.time()
        # print(end_time-start_time)

        # o = open(path1Sun,'r')
        # r = o.readlines()
        # r_split  = r[0].split()[1]

        # self.NSun1lbd = int(r_split)
        # if verbose:
        #     print("... NSun1lbd: {0:}".format(self.NSun1lbd))
        # self.lamb1Sun = np.zeros([self.NSun1lbd], dtype=float)
        # self.flux1Sun = np.zeros([self.NSun1lbd], dtype=float)

        # for i in enumerate(r[1:]):
        #     r_split = i[1].split()
        #     self.lamb1Sun[i[0]] = r_split[0]
        #     self.flux1Sun[i[0]] = r_split[1]
        # o.close()

        # start_time = time.time()
        data = np.loadtxt(path1Sun,skiprows=0)

        self.lamb1Sun = np.copy(data[:,0])
        self.flux1Sun = np.copy(data[:,1])
        self.NSun1lbd = self.lamb1Sun.size
        # end_time = time.time()
        # print(end_time-start_time)

        #fits_open = fits.open(path2Sun)
        #header = fits_open[0].data
        #data = fits_open[1].data

        #self.lamb2Sun = data['WAVELENGTH']
        #self.flux2Sun = data['FLUX']
        #self.NSun2lbd = self.lamb2Sun.size

        #fits_open.close()
#  *** SUN spectrum ****************************************************************

# *** Reading of the spectrum of BD+17d4708 ****************************************
#        RESUME : F subdwarf used to calibrate the Thuan & Gunn system.                                          !
# Coordinates : 22:11:31.37 +18:05:34.1  ±   0.001                                                               !
##################################################################### 
        # o = open(path_BD,'r')
        # r = o.readlines()
        # r_split  = r[0].split()[1]
        
        # self.NFBD_lbd = int(r_split)
        # if verbose:
        #     print("... NFBD_lbd: {0:}".format(self.NFBD_lbd))
        # self.lamb_FBD = np.zeros([self.NFBD_lbd], dtype=float)
        # self.flux_FBD = np.zeros([self.NFBD_lbd], dtype=float)
        
        # for i in enumerate(r[1:]):
        #     r_split = i[1].split()
        #     self.lamb_FBD[i[0]] = r_split[0]
        #     self.flux_FBD[i[0]] = r_split[1]
        # o.close()

        # start_time = time.time()
        data = np.loadtxt(path_BD,skiprows=0)

        self.lamb_FBD = np.copy(data[:,0])
        self.flux_FBD = np.copy(data[:,1])
        self.NFBD_lbd = self.lamb_FBD.size
        # end_time = time.time()
        # print(end_time-start_time)
# *** Reading of the spectrum of BD+17o4708 ******************************************

# *** Interpolate/Extrapolate boundaries *********************************************
        def interp_stars(lbd,flx):
            '''Log interpolate spectrum of calibrated stars upper part'''
            if np.log10(lbd[-1]) <= 4.5:
                delta_lbd = 0.1e5
            else:
                delta_lbd = 1e6
                
            lamb_extended = np.arange( lbd[-1]+delta_lbd/2.,50.1e7,delta_lbd )

            # Calculate the slope for extrapolation in log space
            slope = (np.median(np.log10(flx[-21:-1])) - np.median(np.log10(flx[-41:-21]))) / \
            (np.median(np.log10(lbd[-21:-1])) - np.median(np.log10(lbd[-41:-21])))

            flux_linear_extended_log = np.log10(flx[-1]) + slope * (np.log10(lamb_extended) - np.log10(lbd[-1]))
            flux_linear_extended = 10.**(flux_linear_extended_log)

            lamb_extended = np.concatenate( [lbd,lamb_extended] )
            flux_extended = np.concatenate( [flx,flux_linear_extended] )

            index_ = np.argsort(lamb_extended)
            lamb_extended = lamb_extended[index_]
            flux_extended = flux_extended[index_]

            flux_extended = np.where( np.isnan(flux_extended), 0.0, flux_extended )
            flux_extended = np.where( np.isinf(flux_extended), 0.0, flux_extended )

            lamb_extended = lamb_extended[ flux_extended > 0.0 ]
            flux_extended = flux_extended[ flux_extended > 0.0 ]

            return lamb_extended,flux_extended
        
        def lower_interp_stars(lbd,flx):
            '''Log interpolate spectrum of calibrated stars lower part'''
            delta_lbd = 1.0
            
            if lbd[0] > 1000.0:
                # print('EXTENDED')
                lamb_extended = np.arange( 1000.0,lbd[0],delta_lbd )

                # Calculate the slope for extrapolation in log space
                slope = (np.median(np.log10(flx[21:41])) - np.median(np.log10(flx[0:21]))) / \
                        (np.median(np.log10(lbd[21:41])) - np.median(np.log10(lbd[0:21])))

                flux_linear_extended_log = np.log10(flx[0]) + slope * (np.log10(lamb_extended) - np.log10(lbd[0]))
                flux_linear_extended = 10.**(flux_linear_extended_log)
                
                # flux_linear_extended = flux_linear_extended*0+flx[5]
                
                lamb_extended = np.concatenate( [lamb_extended,lbd] )
                flux_extended = np.concatenate( [flux_linear_extended,flx] )

                index_ = np.argsort(lamb_extended)
                lamb_extended = lamb_extended[index_]
                flux_extended = flux_extended[index_]

                flux_extended = np.where( np.isnan(flux_extended), 0.0, flux_extended )
                flux_extended = np.where( np.isinf(flux_extended), 0.0, flux_extended )

                lamb_extended = lamb_extended[ flux_extended > 0.0 ]
                flux_extended = flux_extended[ flux_extended > 0.0 ]
            else:
                lamb_extended = lbd
                flux_extended = flx
            return lamb_extended,flux_extended

        self.lambVega,self.fluxVega = lower_interp_stars( self.lambVega,self.fluxVega )
        self.lambVega,self.fluxVega = interp_stars( self.lambVega,self.fluxVega )
        self.lambvega = self.lambVega
        self.fluxvega = self.fluxVega
        
        self.lamb_Sun,self.flux_Sun = lower_interp_stars( self.lamb_Sun,self.flux_Sun )
        self.lamb_Sun,self.flux_Sun = interp_stars( self.lamb_Sun,self.flux_Sun )
        self.lamb_sun = self.lamb_Sun
        self.flux_sun = self.flux_Sun
        
        self.lamb1Sun,self.flux1Sun = lower_interp_stars( self.lamb1Sun,self.flux1Sun )
        self.lamb1Sun,self.flux1Sun = interp_stars( self.lamb1Sun,self.flux1Sun )
        self.lamb1sun = self.lamb1Sun
        self.flux1sun = self.flux1Sun

        # Debug
        #self.lamb2sun = self.lamb2Sun
        #self.flux2sun = self.flux2Sun

        self.lamb_FBD,self.flux_FBD = lower_interp_stars( self.lamb_FBD,self.flux_FBD )
        self.lamb_FBD,self.flux_FBD = interp_stars( self.lamb_FBD,self.flux_FBD )
        self.lamb_fbd = self.lamb_FBD
        self.flux_fbd = self.flux_FBD

        self.read_calibration_stars = 1

        return

    def evaluate_photometry( self ):
        o_lambda = self.lamb_sun
        o_fluxes = self.flux_sun 

        t_lambda = self.T_lambda
        t_fluxes = self.T_fluxes

        t_l_area = self.t_l_area
        magabsys = self.magabsys
        magtgsys = self.magtgsys
        standard = self.standard
        numb_lbd = self.numb_lbd
        int_type = 2
        # mag_spec,photflux,mcalibra,iskeepon = evalfilters(o_lambda,o_fluxes,t_lambda,t_fluxes,t_l_area,magabsys,magtgsys,standard,numb_lbd,[int_type,verbosity])
        # mag_spec,photflux,mcalibra,iskeepon = evalf.evalfilters( o_lambda,o_fluxes,t_lambda,t_fluxes,t_l_area,magabsys,magtgsys,standard,numb_lbd,int_type=int_type,verbosity=0 )

        # mag_spec,photflux,mcalibra,iskeepon = evalf( o_lambda,o_fluxes,t_lambda,t_fluxes,t_l_area,magabsys,magtgsys,standard,numb_lbd,int_type=int_type,verbosity=0 )

        self.mag_spec = np.zeros([numb_lbd.size], dtype=float)
        self.photflux = np.zeros([numb_lbd.size], dtype=float)

        for i in enumerate(range(self.Nfilters)):
            # Debug
            # for j in range(numb_lbd[i[0]]):
                # if t_lambda[j,i[0]] <= 0.0:
                    # print(i[0],'PROBLEM')
            flag = (o_lambda >= t_lambda[0,i[0]]) & (o_lambda >= t_lambda[numb_lbd[i[0]],i[0]])
            print("flag ======> ",flag)
            mag_spec,photflux,mcalibra,iskeepon = evalf( o_lambda,o_fluxes,t_lambda[0:numb_lbd[i[0]],i[0]],t_fluxes[0:numb_lbd[i[0]],[i[0]]],t_l_area[i[0]],magabsys[i[0]],magtgsys[i[0]],standard[i[0]],numb_lbd[i[0]],int_type=int_type,verbosity=0 )
            self.mag_spec[i[0]] = mag_spec[0]
            self.photflux[i[0]] = photflux[0]

        # self.mag_spec = mag_spec
        # self.photflux = photflux
        #print(magabsys,mag_spec)

        xlow = self.T_lambda.flatten()
        xlow = xlow[ xlow>0.0 ].min() * 0.85
        xupp = self.T_lambda.flatten()
        xupp = xupp[ xupp>0.0 ].max() * 1.15

        zorder  = 0
        for i in enumerate(self.T_lambda[0,:]):
            plt.plot(self.T_lambda[0:self.Numblbd[i[0]],i[0]],self.T_fluxes[0:self.Numblbd[i[0]],i[0]]/max(self.T_fluxes[0:self.Numblbd[i[0]],i[0]])*max(o_fluxes[  (o_lambda >= xlow) & (o_lambda <= xupp) ]),alpha=0.5,zorder=zorder)
            zorder += 1
            plt.fill_between(self.T_lambda[0:self.Numblbd[i[0]],i[0]],0,self.T_fluxes[0:self.Numblbd[i[0]],i[0]]/max(self.T_fluxes[0:self.Numblbd[i[0]],i[0]])*max(o_fluxes[  (o_lambda >= xlow) & (o_lambda <= xupp) ]),alpha=0.2,zorder=zorder)
                
        zorder += 1
        plt.plot(self.lamb_eff[ (self.lamb_eff >= xlow) & (self.lamb_eff <= xupp) ],photflux[ (self.lamb_eff >= xlow) & (self.lamb_eff <= xupp) ],zorder=zorder)
        zorder += 1
        plt.plot(o_lambda[ (o_lambda >= xlow) & (o_lambda <= xupp) ],o_fluxes[  (o_lambda >= xlow) & (o_lambda <= xupp) ],zorder=zorder)
        # plt.scatter(o_lambda[ (o_lambda >= xlow) & (o_lambda <= xupp) ],o_fluxes[  (o_lambda >= xlow) & (o_lambda <= xupp) ],color='black',s=30,zorder=zorder)
        zorder += 1
        plt.scatter(self.lamb_eff[ (self.lamb_eff >= xlow) & (self.lamb_eff <= xupp) ] ,photflux[  (self.lamb_eff >= xlow) & (self.lamb_eff <= xupp) ],color='red',zorder=zorder)
            
        plt.xscale('log')
        plt.yscale('log')

        plt.xlim(xlow,xupp)
        
        #plt.legend(loc=1, prop={'size': 8}, title='Filters')
        legend = plt.legend(loc=1, prop={'size': 6}, facecolor='grey', fancybox=True, framealpha=1, edgecolor="black", title='Filters')
        legend.get_frame().set_alpha(None)
        legend.set_zorder( zorder + 102 )
        
        plt.xlabel("log Wavelength")
        plt.ylabel("F$_\lambda$ [L$_\odot \AA^{-1}$]" )
        
        plt.show()
        
    def plotfilter( self ):

        if self.readfilters == 0:
            print("... You need to first read the filters")
            return

        o_lambda = self.lambVega
        o_fluxes = self.fluxVega
        t_lambda = self.T_lambda
        t_fluxes = self.T_fluxes
        t_l_area = self.t_l_area
        magabsys = self.magabsys
        magtgsys = self.magtgsys
        standard = self.standard
        lamb_eff = self.lamb_eff

        numb_lbd = self.numb_lbd

        int_type = 2 
        mag_spec,photflux,mcalibra,iskeepon = evalf( o_lambda,o_fluxes,t_lambda,t_fluxes,t_l_area,magabsys,magtgsys,standard,numb_lbd,int_type=int_type,verbosity=0 )

        Figure = plt.figure( figsize=(12,10),dpi=120,facecolor='w',edgecolor='w' )
        plt.subplots_adjust(bottom=.02, left=.06, right=.95, top=.98, wspace=0.0, hspace=0.0) 
    
        ax = plt.subplot(111)

        # Top plot ###########################################################
        ax1_top = subplot2grid( (20,20), (0,0), colspan=20, rowspan=10 )                                           
        # Sets the position and size of the panel for Plot #01
        #ax1_top.axis('on')
        #ax1_top.axes.get_xaxis().set_visible(False)

        #ax1_top.set_xscale('log')

        #ax1_top.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        #ax1_top.set_xticks( np.geomspace(1000, 42000 ,20).round() )
        #ax1_top.axis.set_minor_formatter(NullFormatter())
        lbd_norm = 5000.

        norm_sun = np.interp( lbd_norm, self.lamb_sun, self.flux_sun )
        norm1sun = np.interp( lbd_norm, self.lamb1sun, self.flux1sun )
        #norm2sun = np.interp( lbd_norm, self.lamb2sun, self.flux2sun )

        normvega = np.interp( lbd_norm, self.lambvega,self.fluxvega )
        norm_fbd = np.interp( lbd_norm, self.lamb_fbd,self.flux_fbd )

        ax1_top.plot( np.log10(self.lamb_sun[ self.flux_sun>0.0 ]), np.log10(self.flux_sun[self.flux_sun>0.0]/norm_sun+1e-30),label='Sun')
        ax1_top.plot( np.log10(self.lamb_sun[ self.flux1sun>0.0 ]), np.log10(self.flux1sun[self.flux_sun>0.0]/norm1sun+1e-30),label='Sun Check')
        #ax1_top.plot( np.log10(self.lamb2sun[ self.flux2sun>0.0 ]), np.log10(self.flux2sun[self.flux2sun>0.0]/norm2sun+1e-30),linestyle='--',label='Sun FITS')

        #ax1_top.plot( np.log10(self.lambvega),np.log10(self.fluxvega/normvega),label='Vega')
        #ax1_top.plot( np.log10(self.lamb_fbd),np.log10(self.flux_fbd/norm_fbd),label='BD+17d4708')

        ax1_top.scatter( lamb_eff,photflux,color='red' )

        x_vec = np.arange(2,6.7,0.2)
        #x_vec  =np.arange(1000,60000,3000) #10.0**x_vec
        #ax1_top.set_xticks( 10**x_vec ) #np.geomspace(1000, 20000 ,15).round() )
        ax1_top.xaxis.set_minor_formatter(FormatStrFormatter("%.2f"))

        ax1_top.set_xlim(x_vec[0],x_vec[-1])
        
        flag = (np.log10(self.lambvega) >= x_vec[0]) & (np.log10(self.lambvega) <= x_vec[-1])
        ymin = np.log10(self.fluxvega[ flag ])
        ymax = 1.1
        ax1_top.set_xlim(ymin,ymax)
        
        #Ny = 8
        #Nx = 12
        #ax1_top.xaxis.set_major_locator(plt.MaxNLocator(Nx))
        #ax1_top.xaxis.set_minor_locator(plt.MaxNLocator(Nx*2))
        #ax1_top.yaxis.set_major_locator(plt.MaxNLocator(Ny))
        #ax1_top.yaxis.set_minor_locator(plt.MaxNLocator(Ny*2))
        
        for axis in [ax1_top.xaxis,ax1_top.yaxis]:
            axis.set_major_formatter(ScalarFormatter())
            axis.set_minor_formatter(NullFormatter())

        ax1_top.set_xticks(x_vec)
        #ax1_top.set_yticks(np.arange(0,1.1,0.1)) 

        ax1_top.legend(loc=1, prop={'size': 8}, title='Calibration Stars')
        ax1_top.set_xlabel("log Wavelength [$\AA$]")
        ax1_top.set_ylabel(f"log F$_\lambda$ [Normalized at {lbd_norm:<3.0f} $\AA$]")
        # Top plot ###########################################################

        # Bottom plot ########################################################
        ax1_bot = subplot2grid( (20,20), (11,0), colspan=20, rowspan=10 )                                           
        # Sets the position and size of the panel for Plot #01
        #ax1_top.axis('on')
        #ax1_top.axes.get_xaxis().set_visible(False)
        
        x_vec = np.arange(2,6.7,0.2)
        ax1_bot.set_xlim(x_vec[0],x_vec[-1])

        #x_vec = np.arange(3,4.8,0.1)
        #x_vec  =np.arange(1000,60000,3000) #10.0**x_vec
        #ax1_top.set_xticks( 10**x_vec ) #np.geomspace(1000, 20000 ,15).round() )
        ax1_top.xaxis.set_minor_formatter(FormatStrFormatter("%.2f"))

        # Ny = 8
        # Nx = 12
        #ax1_top.xaxis.set_major_locator(plt.MaxNLocator(Nx))
        #ax1_top.xaxis.set_minor_locator(plt.MaxNLocator(Nx*2))
        #ax1_top.yaxis.set_major_locator(plt.MaxNLocator(Ny))
        #ax1_top.yaxis.set_minor_locator(plt.MaxNLocator(Ny*2))
        
        for axis in [ax1_bot.xaxis,ax1_bot.yaxis]:
            axis.set_major_formatter(ScalarFormatter())
            axis.set_minor_formatter(NullFormatter())

        ax1_bot.set_xticks(x_vec)
        ax1_bot.set_yticks(np.arange(0,1.1,0.1)) 

        for i in enumerate(self.T_lambda[0,0:self.Nfilters]):
            N_lambda = self.Numblbd[i[0]]
            ax1_bot.plot(np.log10(self.T_lambda[0:N_lambda,i[0]]),self.T_fluxes[0:N_lambda,i[0]],label=str(self.name_fil[i[0]]).split('.txt')[0])
                        
        ax1_bot.legend(loc=1, prop={'size': 8}, title='Filters')
        ax1_bot.set_xlabel("log Wavelength")
        ax1_bot.set_ylabel("Transmission")
        # Bottom plot ########################################################

        plt.show()

        return

# Main
if __name__ == '__main__':
    print("... Building filters database")
    print("... sqlalchemy version: {}".format(sqlalchemy.__version__))
    
    database_filters = create_database_filters()
    database_filters._filters_class()

    # read filters
    verbose = False
    database_filters.read_filters( verbose=verbose )
    # print( "... N_entries: ",len(database_filters.filters)/4 )

    # read spectra
    # database_filters.read_spectra( verbose=verbose ) 
    
    # search filters
    word_search = 'planck','wise','iras','herschel'
    fil,ind = database_filters.FindFilters(word_search)
    # print(fil,ind)
    
    # plot
    start0time = time.time()
    database_filters.plotfilter(ind_filters=ind)    
    end__0time = time.time()
    print(f"... Time to plot photometric fluxes of stars: {end__0time-start0time:<010.5} s")

    # Read a given spectrum test
    lbd = []
    flx = []
    
    # Read template test ******************************************************
    package_name = 'pyphotometry'
            
    # Get the distribution object for the package
    package_dist = pkg_resources.get_distribution(package_name)
            
    # Get the base directory path of the package
    path = package_dist.location \
         + '/' \
         + package_name \
         + '/data/templates/bc2003_hr_m62_chab_ssp_190.spec'
            
    print("... path directory:", path)
    # Read template test ******************************************************
    
    file = path
        
    o = open(file,'r')
    r = o.readlines()
    for i in r:
        i_split = i.split()
        if i_split[0][0] != '#':
            lbd.append(i_split[0])
            flx.append(i_split[1])
            
    lbd = np.array(lbd, dtype=float)
    flx = np.array(flx, dtype=float)

    lbd_ = lbd[ lbd <= 1e4 ]
    
    delta_l_log = 0.01
    
    log_l_aux = np.arange(4.0,6.2,delta_l_log, dtype=float)
    log_l_ = np.concatenate( [np.log10(lbd_),log_l_aux] )
    
    # interpolate in log l
    f_ = np.interp(log_l_,np.log10(lbd),np.log10(flx))
    l_ = 10**(log_l_)
    
    # log_l_ = np.arange(1.9,6.2,0.000001)
    # print(l_[0],l_[-1])
    # f_ = np.interp(log_l_,np.log(l),f)
    l = l_
    f = 10.**f_

    # print(l.shape,f.shape)
    o.close()
    
    # print(l[-1],f[-1])
    # lamb_extended = np.arange( 1.59e6,200.1e7,0.1e3 )
    log_lamb_extended = np.arange( 6.2+delta_l_log,9.3,delta_l_log )
    
    # # Calculate the slope for extrapolation in log space
    slope = (np.mean(np.log10(f[-2:])) - np.mean(np.log10(f[-4:-2]))) / \
            (np.mean(np.log10(l[-2:])) - np.mean(np.log10(l[-4:-2])))

    # print(slope)

    flux_linear_extended_log = np.log10(f[-1]) \
                             + slope \
                             * (log_lamb_extended - np.log10(l[-1]))
    flux_linear_extended = 10.**(flux_linear_extended_log)
    lamb_extended = np.concatenate( [l,10.**log_lamb_extended] )
    flux_extended = np.concatenate( [f,flux_linear_extended] )

    l = lamb_extended
    f = flux_extended

    # flag = (f > 0.0) & ~np.isnan(f) & ~np.isinf(f)
    # l = l[ flag ]
    # f = f[ flag ]

    # i= np.argsort(l)

    # l = l[i]
    # f = f[i]
    
    # log_l_ = np.log10( np.arange(90,200e7,50.) )
    # f_ = np.interp(log_l_,np.log10(l),np.log10(f))
    # l = 10**log_l_
    # f = 10**f_

    # print(l.shape)
    # stop
    # min_order = np.floor(np.log10(np.min(f)))
    # print(l.size)

    # Evaluate
    plot = True
    start1time = time.time()
    database_filters.evaluate_photometry( l,f,ind_filters=ind,plot=plot )
    end1time = time.time()
    print(f"... Time to plot photometric fluxes of spectrum: {end1time-start1time:<010.5} s")

    # photflux = database_filters.photflux
    # lamb_eff = database_filters.lamb_eff
    
    #print()
    #print( database_filters.mag_spec )
    #print( photflux )
    #print( lamb_eff )
    
    # # Evaluate from pyphot
    # start2time = time.time()
    # database_filters.evaluate_photometry_pyphot( l,f, plot=plot )
    # end2time = time.time()

    # plt.xscale('log')
    # plt.yscale('log')

    # plt.plot( l, f ,zorder=0 )
    # plt.scatter( lamb_eff, photflux, color='red',s=20.0,zorder=1 )
    # plt.scatter( database_filters.lamb_eff, database_filters.photflux , color='darkorange', s=7.5,zorder=2 )
    
    # plt.ylabel('$F_\lambda$')
    # plt.xlabel('log wavelength [$\AA$]')
    
    # print( "... Times: {0:} and {1:}".format(end1time-start1time,end2time-start2time) )
    # print( "... Ratio of times: {0:}".format( (end1time-start1time)/(end2time-start2time) ) )

    # import pyphot
    # # get the internal default library of passbands filters
    # lib = pyphot.get_library()
    
    
    # name_lib = [ 'GALEX_FUV',
    #                       'GALEX_NUV',
    #                       'SDSS_u',
    #                       'SDSS_g',
    #                       'SDSS_r',
    #                       'SDSS_i',
    #                       'SDSS_z',
    #                       '2MASS_H',
    #                       '2MASS_J',
    #                       '2MASS_Ks',
    #                       'WISE_RSR_W1',
    #                       'WISE_RSR_W2',
    #                       'WISE_RSR_W3',
    #                       'WISE_RSR_W4'     ]
    
    # # Compute photometric flux with pyphot   
    # n = len(name_lib) 
    # lambd = np.zeros( [n], dtype=float )
    # fluxes = np.zeros( [n], dtype=float )
   
    # for i in range( n ):
    #     filters = lib[ name_lib[ i ] ]
            
    #     #print( filters.info() )
        
    #     # compute the integrated flux through the filter f
    #     # note that it work on many spectra at once
    #     val = filters.get_flux( l, f, axis=1 )
    #     lbd = filters.leff.to("Angstrom")
        
    #     fluxes[ i ] = val.value
    #     lambd[ i ] = lbd.value
    
    # plt.scatter( np.log10(lambd), np.log10(fluxes) )
    
    # print(lambd)
    
    # convert to vega magnitudes
    #mags = -2.5 * np.log10(fluxes) - f.Vega_zero_mag
    # or similarly
    #mags = -2.5 * np.log10(fluxes / f.Vega_zero_flux)

# o = Filters()

# path = '../../data/'
# arq_fil1 = 'ListFilters.txt'

# o.plotfilter()
# o.ReadFilters( path,arq_fil1)
#o.plotfilter()

# o.evaluate_photometry()

#print(o.T_lambda.shape)
#print(repr(o.name_fil))

