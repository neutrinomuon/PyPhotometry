#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Revised interface on Sat Jan 30 12:07:21 2021

@author: Jean Gomes

RESUME :  Filters

Version: v0.0.6

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
from pylab import *
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
from PyPhotometry.flib import propfilters as prop
from PyPhotometry.flib import evalfilters as evalf

# import pyphot to check the results
import pyphot
from pyphot import Sun, Vega, unit

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
class create_database_filters(object):
    
    global Base
    Base = declarative_base()
    
    def __init__( self, database_path=None,database_filename='filters.db' ):
        self.store_filters = 0
        self.readfilters = 0
        
        # database_path by default is define as None, if so, look for installation of PyPhotometry
        # Specify the package name
        if database_path == None:
            package_name = 'PyPhotometry'

            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)

            # Get the base directory path of the package
            database_path = package_dist.location + '/' + package_name + '/data/'

            print("... Package directory:", database_path)
        
        self.database_path = database_path
        self.database_filename = database_filename
        self.database_filters = os.path.join(database_path, database_filename)
        
        # By using pkg_resources.resource_filename(__name__, 'filters.db'), the code retrieves the 
        # absolute file path of 'filters.db' within the package or module. This path is then assigned 
        # to self.database_filters for further use, such as creating an SQLAlchemy engine with the 
        # SQLite database located at that path. However, it does not work if where you're running it 
        # does not contain the filters.db file. So, a better approach is from above.
        #self.database_filters = pkg_resources.resource_filename(__name__, 'filters.db')
        
        self.Engine        = create_engine('sqlite:///' + self.database_filters, echo=False)
                
        self.Session       = sessionmaker(bind=self.Engine)
        self.session       = self.Session()

    # Create the table information object class
    class _filters_class( Base ):
        """Information to create the filters SQL table"""

        __tablename__ = 'filters'
        
        filterid = Column(Integer, primary_key=True)
        name_filter = Column(String, primary_key=True)
        detector =  Column(String)
        units =Column(String)
        
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
        """Information to create the filters SQL table"""

        __tablename__ = 'spectra'
        
        specid = Column(Integer, primary_key=True)
        name_spec = Column(String, primary_key=True)
        flux_units = Column(String)
        lbd_units = Column(String)
        
        wavelength = Column(PickleType)
        fluxes = Column(PickleType)

        N_lambda = Column(Integer)
                
    #def read_spectra( self, path_Vega='../../data/calibration_stars/VegaLR.dat',path_Sun='../../data/calibration_stars/Sun_LR.dat',path1Sun='../../data/calibration_stars/Sun.dat',path2Sun='/../../data/calibration_stars/sun_reference_stis_001.fits',path_BD='../../data/calibration_stars/BD+17d4708.dat', verbose=False ):
    def read_spectra( self, path_Vega=None,path_Sun=None,path1Sun=None,path2Sun=None,path_BD=None, verbose=False ):
        
        # Read Vega *******************************************************************
        if path_Vega == None:
            package_name = 'PyPhotometry'

            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)

            # Get the base directory path of the package
            path_Vega = package_dist.location + '/' + package_name + '/data/calibration_stars/VegaLR.dat'

            print("... path_Vega directory:", path_Vega)
        # Read Vega *******************************************************************

        # Read Sun ********************************************************************
        if path_Sun == None:
            package_name = 'PyPhotometry'

            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)

            # Get the base directory path of the package
            path_Sun = package_dist.location + '/' + package_name + '/data/calibration_stars/Sun_LR.dat'

            print("... path_Sun directory:", path_Sun)
        # Read Sun ********************************************************************
        
        # Read Sun_1 ******************************************************************
        if path1Sun == None:
            package_name = 'PyPhotometry'

            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)

            # Get the base directory path of the package
            path1Sun = package_dist.location + '/' + package_name + '/data/calibration_stars/Sun.dat'

            print("... path1Sun directory:", path1Sun)
        # Read Sun_1 ******************************************************************
        
        # Read Sun_2 ******************************************************************
        if path2Sun == None:
            package_name = 'PyPhotometry'

            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)

            # Get the base directory path of the package
            path2Sun = package_dist.location + '/' + package_name + '/data/calibration_stars/sun_reference_stis_001.fits'

            print("... path2Sun directory:", path2Sun)
        # Read Sun_2 ******************************************************************
        
        # Read BD *********************************************************************
        if path_BD == None:
            package_name = 'PyPhotometry'

            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)

            # Get the base directory path of the package
            path_BD = package_dist.location + '/' + package_name + '/data/calibration_stars/BD+17d4708.dat'

            print("... path_BD directory:", path2Sun)
        # Read BD *********************************************************************
        
        o = Filters( )
        o.ReadCalibrationStars(  path_Vega=path_Vega,path_Sun=path_Sun,path1Sun=path1Sun,path2Sun=path2Sun,path_BD=path_BD )
        
        self.lambVega = o.lambVega
        self.fluxVega =o.fluxVega
        self.lambvega = self.lambVega
        self.fluxvega = self.fluxVega
        
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

        if verbose:
            print("[read_spectra]")
            print("... Read spectra of calibration stars")
            print("[read_spectra]")
            
        return
        
    def read_filters( self, path_data=None,N_lambda=500,verbose=False ):
                
        if verbose:
            print("[read_filters]")
            print("... Reading filters")
        
        # database_path by default is define as None, if so, look for installation of PyPhotometry
        # Specify the package name
        if path_data == None:
            package_name = 'PyPhotometry'

            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)

            # Get the base directory path of the package
            path_data = package_dist.location + '/' + package_name + '/data/'

            print("... path_data directory:", path_data)
        
        # Verify if database contains models for Draine and Li (2007)
        try:
            read = self.session.query(self._filters_class).filter( np.size(self._filters_class.name_filter) > 0 )
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
            print("... Verified: {0:} entries in filters.db".format(np.count_nonzero(read.all())))
        except:
            print("... Failed")
                    
        if verbose:
            print("... Start reading filters database")
        
        try:
            if np.count_nonzero(read.all()) > 0:
                print("... Filters database already exists")
                self.store_filters = 1
                                
                meta = MetaData()
                
                filter_table = Table(
                                                'filters', meta, 
                                                Column('filterid', Integer, primary_key = True), 
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
                
                conn = self.Engine.connect( ) 
                s = filter_table.select()
                result = conn.execute(s)
                                
                filters_in_database = []
                for row in result:
                    filters_in_database.append(row[1])

                #filters_in_database = nparray( filters_in_database, dtype=object )
                N_filters_in_database =  np.size(filters_in_database)
                
                #print(filters_in_database)
            
                # Add filter if is not in database                
                path = path_data #'../../data/'
                arq_fil1 = 'ListFilters.txt' # Just an ascii file
                
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

                o = Filters( )
                path = path_data #'../../data/'
                arq_fil1 = 'ListFilters.txt'
                
                for j in name_list[ index_filters_in_database == 0 ]:
                    #print(j)
                    string_ = '%{ }%'.format(j)
                    s = filter_table.select().where( filter_table.c.name_filter.ilike(string_) )
                    result = conn.execute(s)
                    # print( string_ )
                    # print( result )

                    count = 0
                    for row in result:
                        # print (row[1])
                        count += 1
                        
                    # print( count )
                        
                    if count <= 0:
                        #ReadONEFilter( self,path,arq_fil1,N_lambda=500 )
                        path = path_data #'../../data/'
                        arq_fil1 = j + '.txt'
                        #print(path + arq_fil1)
                        o.ReadONEFilter( path=path, arq_fil1=arq_fil1)
                                                
                        d = o.onefilter[0]
                        filter_object = self._filters_class( name_filter=d[0], filterid=N_filters_in_database, detector=d[1], units=d[4]['units'], wavelength=np.array(d[2], dtype=float), transmission=d[3], N_lambda=int(d[4]['N_lambda']), t_l_area=d[4]['t_l_area'], t_n_area=d[4]['t_n_area'], vegaflux=d[4]['standard'], lamb_eff=d[4]['lamb_eff'], widtheff=d[4]['widtheff'], magabsys=d[4]['magabsys'], magtgsys=d[4]['magtgsys'] )
                        
                        #print(filter_object)
                        self.session.add( filter_object )
                        
                        try:
                            self.session.commit()
                            print("... The filter {0:} was included in the database filters.db.".format(d[0]))
                        except:
                            self.session.rollback()
                    
                self.store_filters = 1                                
                
                # Read database
                s = filter_table.select()
                result = conn.execute(s)
                
                rows = self.session.query( filter_table ).count()
                self.N_filters = rows
                self.Nfilters = self.N_filters

                filterid = np.zeros( [self.N_filters],dtype=int )
                t_l_area = np.zeros( [self.N_filters],dtype=float )
                t_n_area =  np.zeros( [self.N_filters],dtype=float )
                magabsys =  np.zeros( [self.N_filters],dtype=float )
                magtgsys =  np.zeros( [self.N_filters],dtype=float )
                standard =  np.zeros( [self.N_filters],dtype=float )
                numb_lbd =  np.zeros( [self.N_filters],dtype=int )
                lamb_eff =  np.zeros( [self.N_filters],dtype=float )
                widtheff =  np.zeros( [self.N_filters],dtype=float )
                units = np.zeros( [self.N_filters],dtype=object )
                name_filter = np.zeros( [self.N_filters],dtype=object )
                detector =  np.zeros( [self.N_filters],dtype=object )
                
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
                    units[ row[0] ] = row[3]
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
                    d['magtgsys'] = magtgsys[i]
                    d['standard'] = standard[ i ]
                    d['lamb_eff'] = lamb_eff[ i ]
                    d['widtheff'] = widtheff[ i ]
                    d['units'] = units[ i ]
                    
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
            print("... Need to construct filters database")
        
            metadata = Base.metadata.create_all(self.Engine)
        
            o = Filters()
            #print(o)
            path = path_data #'../../data/'
            arq_fil1 = 'ListFilters.txt'
            IsKeepOn = o.ReadFilters( path,arq_fil1,verbose=verbose )
        
            if IsKeepOn != 1:
                print("... Problem running filters database")
                return
        
            # 4 different entry for each filter
            N_entries = int( len(o.filters) / 4 )
            # print( len(o.filters) )
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
               
                #print(filter_object.name_filter)
                self.session.add( filter_object )
                
                # wave = np.array(d[2], dtype=float)
                # fluxes = d[3]
                wave = np.array(d[3], dtype=float)
                fluxes = d[4]
                j = wave.size
                self.T_lambda[ 0:j , i ] = wave[ 0:j ] 
                self.T_fluxes[ 0:j , i ] = fluxes[ 0:j ] 
                                
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
                self.t_l_area[ i ]  = d[5]['t_l_area']
                self.magabsys[ i ] = d[5]['magabsys']
                self.magtgsys[ i ] = d[5]['magtgsys']
                self.standard[ i ] = d[5]['standard']
                self.lamb_eff[ i ] = d[5]['lamb_eff']
                self.widtheff[ i ] =d[5]['widtheff']
            
                # print('AQUI',self.name_fil[ i ],i,N_entries)    
            
            self.numb_lbd = self.Numblbd
            # print( 'PASSOU',N_entries )
            
            #print(self.session.new)
            try:
                self.session.commit()
                print("... The filters were included in the database filters.db.")
                self.store_filters = 1
            
            except:
                self.session.rollback()
                print("... The filters are already included in the database filters.db.")
                
            self.filters = o.filters
            self.session.close()
            
            self.readfilters = 1
            return
    
    def FindFilters( self,word_search=None ):
    
        if self.readfilters == 0:
            print("... You need to first read the filters")
            return

        # Filters in the database
        if word_search == None:
            #N_filters = self.N_filters
            d = list( self.filters.keys() )
            # print(d)
            print("... List all filters - N_filters: {}".format(self.N_filters))
            for i in enumerate(d[::3]):
                print( "    {0:>04d} {1:}".format(i[0],self.filters[ i[1] ][0]) )
            i = d[::3]
            r = d[2::3]
        else:               
            r,i = self.search( self.filters, word_search )
            #print( r )
        
        self.ind_names = r
        self.ind_filters = i
        return r,i
    
    def search( self, values, searchFor ):
        
        isfilterinlist = 0
        keys = list( values.keys() )
        
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
            
            for k in enumerate(keys[2::3]):
                fullstring = str(k[1]).lower()
                #print( fullstring )
            
                if substring in fullstring:
                    #print(k)
                    indsearched.append( k[0] )
                    listsearched.append( k[1] )
                    isfilterinlist = 1
            
        if isfilterinlist:
            indsearched = np.array( indsearched,dtype=int )
            # listsearched = listsearched
            return listsearched,indsearched
        else:
            return None, None
                    
    def plotfilter( self, N_lambda=500 ):
        
        if self.readfilters == 0:
            print("... You need to first read the filters")
            return
                
        #self.T_lambda = self.filters
        #self.T_fluxes =
        
        Figure = plt.figure( figsize=(12,10),dpi=120,facecolor='w',edgecolor='w' )
        plt.subplots_adjust(bottom=.02, left=.06, right=.95, top=.98, wspace=0.0, hspace=0.0) 
        
        #ax = plt.subplot(111)

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

        #print("Continuar aqui tem que ler e estocar os espectros de referência ")
    
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

        # x_vec = np.arange(2,8.4,0.2)
        # for i in enumerate(self.T_lambda[0,0:self.Nfilters]):
        #     N_lambda = self.Numblbd[i[0]]
        #     print( self.T_lambda[N_lambda-1,i[0]] )
        positive_values = self.T_lambda[np.where(self.T_lambda > 0)]

        min_xvec = np.round( np.log10( np.min(positive_values) ), 1 )
        max_xvec  = np.round( np.log10( np.max(positive_values) ),1 ) 
        
        delta_xvec = np.round((max_xvec - min_xvec) / 12, 1)
        min_xvec -= delta_xvec
        max_xvec += delta_xvec + 1.0

        x_vec = np.arange( min_xvec,max_xvec,delta_xvec )
        
        #x_vec  =np.arange(1000,60000,3000) #10.0**x_vec
        #ax1_top.set_xticks( 10**x_vec ) #np.geomspace(1000, 20000 ,15).round() )
        ax1_top.xaxis.set_minor_formatter(FormatStrFormatter("%.2f"))

        ax1_top.set_xlim(x_vec[0],x_vec[-1])

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
        ax1_top.set_ylabel(r"log F$_\lambda$ [Normalized at {0:<3.0f} $\AA$]".format(lbd_norm))
        # Top plot ###########################################################

        # Bottom plot ########################################################
        ax1_bot = subplot2grid( (20,20), (11,0), colspan=20, rowspan=10 )                                           
        # Sets the position and size of the panel for Plot #01
        #ax1_top.axis('on')
        #ax1_top.axes.get_xaxis().set_visible(False)
        
        # x_vec = np.arange(np.min(),6.7,0.2)
        ax1_bot.set_xlim(x_vec[0],x_vec[-1])

        #x_vec = np.arange(3,4.8,0.1)
        #x_vec  =np.arange(1000,60000,3000) #10.0**x_vec
        #ax1_top.set_xticks( 10**x_vec ) #np.geomspace(1000, 20000 ,15).round() )
        ax1_top.xaxis.set_minor_formatter(FormatStrFormatter("%.2f"))

        Ny = 8
        Nx = 12
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
                        
            if self.detector[ i[0] ] == 'energy':
                ax1_bot.plot(np.log10(self.T_lambda[0:N_lambda,i[0]]),self.T_fluxes[0:N_lambda,i[0]]*self.T_lambda[0:N_lambda,i[0]],label=str(self.name_fil[i[0]]).split('.txt')[0])
            else:
                ax1_bot.plot(np.log10(self.T_lambda[0:N_lambda,i[0]]),self.T_fluxes[0:N_lambda,i[0]],label=str(self.name_fil[i[0]]).split('.txt')[0])
                              
            #print(N_lambda,str(self.name_fil[i[0]]).split('.txt')[0])
            #print(self.T_fluxes[0:N_lambda,i[0]])
            
        ax1_bot.legend(loc=1, prop={'size': 8}, title='Filters')
        ax1_bot.set_xlabel("log Wavelength")
        ax1_bot.set_ylabel("Transmission")
        # Bottom plot ########################################################

        plt.show()

        return
    
    def evaluate_photometry( self, l, f, ind_filters=None, plot=False  ):
        o_lambda = l #self.lamb_sun
        o_fluxes = f #self.flux_sun 
        
        print(o_lambda.shape,o_fluxes.shape)
        
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
        
        #print(numb_lbd)
        
        LSun = 3.839e33
        mag_spec,photflux,mcalibra,iskeepon = evalf( o_lambda,o_fluxes,t_lambda,t_fluxes,t_l_area,magabsys,magtgsys,standard,numb_lbd,int_type=int_type,verbosity=0 )
        
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

        self.mag_spec = mag_spec
        self.photflux    = photflux
        #self.iskeepon_photflux = iskeepon
        #self.lamb_eff =  self.lamb_eff

        #print(magabsys,mag_spec)
        
        if plot:
            # xlow = self.T_lambda.flatten()
            # xlow = xlow[ xlow>0.0 ].min() * 0.85
            # xupp = self.T_lambda.flatten()
            # xupp = xupp[ xupp>0.0 ].max() * 1.15
            
            #print(xlow,xupp)
            try:
                if ind_filters==None:
                    ind_filters = np.arange(0,self.N_filters,1)
            except:                
                if np.size(ind_filters) != 0:
                    ind_filters = ind_filters
                else: 
                    ind_filters = np.arange(0,self.N_filters,1)
            #print( ind_filters )
            
            zorder  = 0
            ind =  np.argsort(self.lamb_eff[ind_filters] )
            ind = ind_filters[ind]
            
            #print( self.lamb_eff[ ind ].size )
            
            xlow = self.T_lambda[:,ind].flatten()
            xlow = xlow[ xlow>0.0 ].min() * 0.85
            xupp = self.T_lambda[:,ind].flatten().flatten()
            xupp = xupp[ xupp>0.0 ].max() * 1.15
            
            n = ind.size
            colors  = pl.cm.jet(np.linspace(0,1,n))
            
            for i in enumerate(self.T_lambda[0,ind_filters]):
                #print(i)
                
                plt.plot( self.T_lambda[0:self.Numblbd[ind[i[0]]],ind[i[0]]],self.T_fluxes[0:self.Numblbd[ind[i[0]]],ind[i[0]]],alpha=0.6,zorder=zorder,label=self.name_fil[ind[i[0]]],color=colors[i[0]] )
                zorder += 1
                plt.fill_between(self.T_lambda[0:self.Numblbd[ind[i[0]]],ind[i[0]]],self.T_fluxes[0:self.Numblbd[ind[i[0]]],ind[i[0]]],color=colors[i[0]],alpha=0.4,zorder=zorder)
          
                # plt.plot( self.T_lambda[0:self.Numblbd[ind[i[0]]],ind[i[0]]],self.T_fluxes[0:self.Numblbd[ind[i[0]]],ind[i[0]]]/max(self.T_fluxes[0:self.Numblbd[ind[i[0]]],ind[i[0]]])*max(o_fluxes[  (o_lambda >= xlow) & (o_lambda <= xupp) ]),alpha=0.6,zorder=zorder,label=self.name_fil[ind[i[0]]],color=colors[i[0]] )
                # zorder += 1
                # plt.fill_between(self.T_lambda[0:self.Numblbd[ind[i[0]]],ind[i[0]]],self.T_fluxes[0:self.Numblbd[ind[i[0]]],ind[i[0]]]/max(self.T_fluxes[0:self.Numblbd[ind[i[0]]],ind[i[0]]])*max(o_fluxes[  (o_lambda >= xlow) & (o_lambda <= xupp) ]),color=colors[i[0]],alpha=0.4,zorder=zorder)
                        
            zorder += 1
        
            lbd =  self.lamb_eff[ ind ] # [ (self.lamb_eff >= xlow) & (self.lamb_eff <= xupp) ]
            jnd = np.argsort(lbd)
        
            #print(self.lamb_eff[ (self.lamb_eff >= xlow) & (self.lamb_eff <= xupp) ][ind]  )
        
            plt.plot(self.lamb_eff[ ind ][ jnd ],photflux[ ind ][ jnd ],color='black',zorder=zorder)
            
            #print(self.lamb_eff[ind],self.name_fil[ind])
            
            zorder += 1
            plt.plot(o_lambda[ (o_lambda >= xlow) & (o_lambda <= xupp) ],o_fluxes[  (o_lambda >= xlow) & (o_lambda <= xupp) ],color='blue',zorder=zorder)
            zorder += 1
            plt.scatter(self.lamb_eff[ ind ] ,photflux[ ind ],color='red',zorder=zorder)
              
            plt.xscale('log')
            plt.yscale('log')
    
            plt.xlim(xlow,xupp)
            
            legend = plt.legend(loc=1, prop={'size': 6}, facecolor='grey', fancybox=True, framealpha=1, edgecolor="black", title='Filters')
            legend.get_frame().set_alpha(None)
            legend.set_zorder( zorder + 102 ) 
            plt.xlabel("log Wavelength [$\AA$]")
            plt.ylabel("F$_\lambda$ [L$_\odot \AA^{-1}$]" )
            
            plt.show()
            
    def evaluate_photometry_pyphot( self, l, f, plot=False  ):
        o_lambda = l #self.lamb_sun
        o_fluxes = f #self.flux_sun 
        
        t_lambda = self.T_lambda
        t_fluxes = self.T_fluxes
        
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
            plt.plot(o_lambda[ (o_lambda >= xlow) & (o_lambda <= xupp) ],o_fluxes[  (o_lambda >= xlow) & (o_lambda <= xupp) ],color='blue',zorder=zorder)
            zorder += 1
            plt.scatter(self.lamb_eff[ (self.lamb_eff >= xlow) & (self.lamb_eff <= xupp) ] ,photflux[  (self.lamb_eff >= xlow) & (self.lamb_eff <= xupp) ],color='red',zorder=zorder)
              
            plt.xscale('log')
            plt.yscale('log')
    
            plt.xlim(xlow,xupp)
            
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
class Filters( object ):
    
    def __init__( self ):
        self.readfilters = 0
        self.filters = {}            

    def ReadFilters( self,path,arq_fil1,N_lambda=500,verbose=False ):
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
            
            for j in enumerate(r1filter):
                j_split = j[1].split()
                                
                if j_split[0][0] != '#':
                    k = self.Numblbd[i[0]]

                    self.T_lambda[k,i[0]] =j_split[0]
                    self.T_fluxes[k,i[0]] =j_split[1]
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
            
            #self.T_lambda[:,i[0]] = self.T_lambda[:,i[0]] * units.Angstrom
                        
            # If detector type is energy then
            #a = np.trapz( ifT[ind] * _sflux, _slamb[ind], axis=axis )
            #b = np.trapz( ifT[ind], _slamb[ind])
            if self.detector_type[i[0]] == 'energy':
                self.T_fluxes[0:self.Numblbd[i[0]],i[0]] /= self.T_lambda[0:self.Numblbd[i[0]],i[0]]
                
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
        t_l_area,t_n_area,magabsys,magtgsys,standard,lamb_eff,widtheff,iskeepon = prop( self.T_lambda,self.T_fluxes,self.Numblbd,self.lambvega,self.fluxvega,self.lamb_sun,self.flux_sun,self.lamb_fbd,self.flux_fbd,verbosity=0 )

        self.t_l_area = t_l_area
        self.t_n_area = t_n_area
        self.magabsys = magabsys
        self.magtgsys = magtgsys
        self.standard = standard
        self.numb_lbd = self.Numblbd
        self.lamb_eff = lamb_eff
        self.widtheff = widtheff
        
        for i in range(Nfilters):
            N_lambda = self.Numblbd[i]
            d = {}
            d['N_lambda'] = N_lambda
            d['t_l_area'] = t_l_area[ i ]
            d['t_n_area'] = t_n_area[ i ]
            d['magabsys'] = magabsys[ i ]
            d['magtgsys'] = magtgsys[i]
            d['standard'] = standard[ i ]
            d['lamb_eff'] = lamb_eff[ i ]
            d['widtheff'] = widtheff[ i ]
            d['units'] = self.units[ i ]
            
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
        print("[ReadFilters]")

        return 1
    
    def ReadONEFilter( self,path,arq_fil1,N_lambda=500 ):
        print("")
        print( "... Reading One Filter " )
        print( "... arq_fil1: {0:}".format(arq_fil1) )
                
        Nfilters =1
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
                        print("... Detector type: {0:}".format(j_split[2]))
                        self.detector1type[0] = j_split[2].lower()
                            
                    if det == 'units':
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
            d['magtgsys'] = magtgsys[i]
            d['standard'] = standard[ i ]
            d['lamb_eff'] = lamb_eff[ i ]
            d['widtheff'] = widtheff[ i ]
            d['units'] = self.units1[ i ]
            
            self.onefilter[ i ] = [ str(self.name1fil[i]).split('.txt')[0],self.detector1type[i],self.T1lambda[0:N_lambda,i], self.T1fluxes[0:N_lambda,i], d ]
            self.onefilter[ str(self.name1fil[i]) + '.txt' ] = self.onefilter[ i ]
            self.onefilter[ self.name1fil[i] ] = self.onefilter[ i ]
            
            #print( self.onefilter )
            return
            
    #def ReadCalibrationStars( self,path_Vega='../../data/calibration_stars/VegaLR.dat',path_Sun='../../data/calibration_stars/Sun_LR.dat',path1Sun='../../data/calibration_stars/Sun.dat',path2Sun='../../data/sun_reference_stis_001.fits',path_BD='../../data/calibration_stars/BD+17d4708.dat',verbose=0 ):
    def ReadCalibrationStars( self,path_Vega=None,path_Sun=None,path1Sun=None,path2Sun=None,path_BD=None,verbose=0 ):
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
            package_name = 'PyPhotometry'
        
            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)
        
            # Get the base directory path of the package
            path_Vega = package_dist.location + '/' + package_name + '/data/calibration_stars/VegaLR.dat'
        
            print("... path_Vega directory:", path_Vega)
# Read Vega *******************************************************************

# Read Sun ********************************************************************
        if path_Sun == None:
            package_name = 'PyPhotometry'
        
            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)
        
            # Get the base directory path of the package
            path_Sun = package_dist.location + '/' + package_name + '/data/calibration_stars/Sun_LR.dat'
        
            print("... path_Sun directory:", path_Sun)
# Read Sun ********************************************************************

# Read Sun_1 ******************************************************************
        if path1Sun == None:
            package_name = 'PyPhotometry'
        
            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)
        
            # Get the base directory path of the package
            path1Sun = package_dist.location + '/' + package_name + '/data/calibration_stars/Sun.dat'
        
            print("... path1Sun directory:", path1Sun)
# Read Sun_1 ******************************************************************

# Read Sun_2 ******************************************************************
        if path2Sun == None:
            package_name = 'PyPhotometry'
        
            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)
        
            # Get the base directory path of the package
            path2Sun = package_dist.location + '/' + package_name + '/data/calibration_stars/sun_reference_stis_001.fits'
        
            print("... path2Sun directory:", path2Sun)
# Read Sun_2 ******************************************************************

# Read BD *********************************************************************
        if path_BD == None:
            package_name = 'PyPhotometry'
        
            # Get the distribution object for the package
            package_dist = pkg_resources.get_distribution(package_name)
        
            # Get the base directory path of the package
            path_BD = package_dist.location + '/' + package_name + '/data/calibration_stars/BD+17d4708.dat'
        
            print("... path_BD directory:", path2Sun)
# Read BD *********************************************************************

#  *** VEGA spectrum ***************************************************************
#         Intrinsic Flux: erg/s/cm2/A                                                                                                                  !
######################################################################
        o = open(path_Vega,'r')
        r = o.readlines()
        r_split  = r[0].split()[1]
        
        self.NVegalbd = int(r_split)
        if verbose:
            print("... NVegalbd: {0:}".format(self.NVegalbd))
        self.lambVega = np.zeros([self.NVegalbd], dtype=float)
        self.fluxVega = np.zeros([self.NVegalbd], dtype=float)
        
        for i in enumerate(r[1:]):
            r_split = i[1].split()
            self.lambVega[i[0]] = r_split[0]
            self.fluxVega[i[0]] = r_split[1]
        o.close()
#  *** VEGA spectrum ***************************************************************

#  *** SUN spectrum ****************************************************************
#         Intrinsic Flux: erg/s/cm2/A                                                                                                                  !
######################################################################     
        o = open(path_Sun,'r')
        r = o.readlines()
        r_split  = r[0].split()[1]
        
        self.NSun_lbd = int(r_split)
        if verbose:
            print("... NSun_lbd: {0:}".format(self.NSun_lbd))
        self.lamb_Sun = np.zeros([self.NSun_lbd], dtype=float)
        self.flux_Sun = np.zeros([self.NSun_lbd], dtype=float)
        
        for i in enumerate(r[1:]):
            r_split = i[1].split()
            self.lamb_Sun[i[0]] = r_split[0]
            self.flux_Sun[i[0]] = r_split[1]
        o.close()
        
        o = open(path1Sun,'r')
        r = o.readlines()
        r_split  = r[0].split()[1]
        
        self.NSun1lbd = int(r_split)
        if verbose:
            print("... NSun1lbd: {0:}".format(self.NSun1lbd))
        self.lamb1Sun = np.zeros([self.NSun1lbd], dtype=float)
        self.flux1Sun = np.zeros([self.NSun1lbd], dtype=float)
        
        for i in enumerate(r[1:]):
            r_split = i[1].split()
            self.lamb1Sun[i[0]] = r_split[0]
            self.flux1Sun[i[0]] = r_split[1]
        o.close()
        
        #fits_open = fits.open(path2Sun)
        #header = fits_open[0].data
        #data = fits_open[1].data
        
        #self.lamb2Sun = data['WAVELENGTH']
        #self.flux2Sun = data['FLUX']
        #self.NSun2lbd = self.lamb2Sun.size
        
        #fits_open.close()
        
#  *** SUN spectrum ****************************************************************

# *** Reading of the spectrum of BD+17d4708 ******************************************
#        RESUME : F subdwarf used to calibrate the Thuan & Gunn system.                                            !
# Coordinates : 22:11:31.37 +18:05:34.1  ±   0.001                                                                                     !
######################################################################     
        o = open(path_BD,'r')
        r = o.readlines()
        r_split  = r[0].split()[1]
        
        self.NFBD_lbd = int(r_split)
        if verbose:
            print("... NFBD_lbd: {0:}".format(self.NFBD_lbd))
        self.lamb_FBD = np.zeros([self.NFBD_lbd], dtype=float)
        self.flux_FBD = np.zeros([self.NFBD_lbd], dtype=float)
        
        for i in enumerate(r[1:]):
            r_split = i[1].split()
            self.lamb_FBD[i[0]] = r_split[0]
            self.flux_FBD[i[0]] = r_split[1]
        o.close()
# *** Reading of the spectrum of BD+17o4708 ******************************************

        self.lambvega = self.lambVega
        self.fluxvega = self.fluxVega
        
        self.lamb_sun = self.lamb_Sun
        self.flux_sun  =self.flux_Sun
        self.lamb1sun = self.lamb1Sun
        self.flux1sun  =self.flux1Sun
        #self.lamb2sun = self.lamb2Sun
        #self.flux2sun  =self.flux2Sun
        
        self.lamb_fbd  = self.lamb_FBD
        self.flux_fbd = self.flux_FBD
        
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

        mag_spec,photflux,mcalibra,iskeepon = evalf( o_lambda,o_fluxes,t_lambda,t_fluxes,t_l_area,magabsys,magtgsys,standard,numb_lbd,int_type=int_type,verbosity=0 )
        self.mag_spec = mag_spec
        self.photflux = photflux
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

        x_vec = np.arange(2,6.7,0.2)
        #x_vec  =np.arange(1000,60000,3000) #10.0**x_vec
        #ax1_top.set_xticks( 10**x_vec ) #np.geomspace(1000, 20000 ,15).round() )
        ax1_top.xaxis.set_minor_formatter(FormatStrFormatter("%.2f"))

        ax1_top.set_xlim(x_vec[0],x_vec[-1])

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
        ax1_top.set_ylabel(r"log F$_\lambda$ [Normalized at {0:<3.0f} $\AA$]".format(lbd_norm))
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

        Ny = 8
        Nx = 12
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
    database_filters.read_spectra( verbose=verbose ) 
    
    # plot
    #database_filters.plotfilter()    

    # Read a given spectrum test
    l = []
    f = []
    file = '/home/jean/Codes/Pynoptic/DustyGlow/src/python/Reemission/SSPs/bc2003_hr_m62_chab_ssp_160.spec'
    o = open(file,'r')
    r = o.readlines()
    for i in r:
        i_split = i.split()
        if i_split[0][0] != '#':
            l.append(i_split[0])
            f.append(i_split[1])
            
    l = np.array(l, dtype=float)
    f = np.array(f, dtype=float)

    l_ = np.arange(500.,12001.,1., dtype=float)
    #print(l_[0],l_[-1])
    f_ = np.interp(l_,l,f)
    l = l_
    f = f_

    print(l.shape,f.shape)

    o.close()
    
    # Evaluate
    plot = True
    start1time = time.time()
    database_filters.evaluate_photometry( l,f, plot=plot )
    end1time = time.time()
    
    photflux = database_filters.photflux
    lamb_eff = database_filters.lamb_eff
    
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

