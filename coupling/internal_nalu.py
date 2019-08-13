#!/usr/bin/env python

'''
Tools for MMC via internal forcing
Based on SowfaInput (internal.py) by Dries Allaerts
'''

__author__ = "Colleen Kaul"
__date__   = "August 13, 2019"

import numpy as np
import pandas as pd
import os


class NaluInput(object):
    """
    Nalu input class for writing data to format that can be included in Nalu-Wind input files
    """
    def __init__(self,
                 dpath,
                 df,
                 dateref,
                 datefrom=None,
                 dateto=None):
        """
        Initialize Nalu-Wind input object

        Usage
        =====
        dpath : str
            Folder to write files to
        df : pandas.DataFrame
            Data (index should be called datetime)
        dateref : str
            Reference datetime, SOWFA uses seconds since this datetime
        datefrom : str
            Start date of the period that will be written out, if None
            start from the first timestamp in df
        dateto : str
            End date of the period that will be written out, if None end
            with the last timestamp in df
        """
        
        self.dpath = dpath
        # Create folder dpath if needed
        if not os.path.isdir(dpath):
            os.mkdir(dpath)

        # Use dataframe between datefrom and dateto
        if datefrom is None:
            datefrom = df.index[0]
        if dateto is None:
            dateto = df.index[-1]
        # Make copy to avoid SettingwithcopyWarning
        self.df = df.loc[(df.index>=datefrom) & (df.index<=dateto)].copy()
        assert(len(self.df.index.unique())>0), 'No data for requested period of time'
        
        # Store start date for ICs
        self.datefrom = datefrom

        # calculate time in seconds since reference date
        dateref = pd.to_datetime(dateref)
        tdelta = pd.Timedelta(1,unit='s')
        self.df.reset_index(inplace=True)
        self.df['t_index'] = (self.df['datetime'] - dateref) / tdelta
        self.df.set_index('datetime',inplace=True)

    '''
    this needs to be updated once Nalu-Wind code is finalized
    '''
    def write_BCs(self,
                  fname,
                  fieldname,
                  fact=1.0
                  ):
        """
        Write surface boundary conditions to Nalu-Wind readable input file
    
        Usage
        =====
        fname : str
            Filename
        fieldname : str
            Name of the field to be written out
        fact : float
            Scale factor for the field, e.g., to scale heat flux to follow
            OpenFOAM sign convention that boundary fluxes are positive if
            directed outward
        """
    
        # extract time and height array
        ts = self.df.t_index.values
        nt = ts.size
    
        # assert field exists and is complete
        assert(fieldname in self.df.columns), 'Field '+fieldname+' not in df'
        assert(~pd.isna(self.df[fieldname]).any()), 'Field '+fieldname+' is not complete (contains NaNs)'
    
        # scale field with factor,
        # e.g., scale heat flux with fact=-1 to follow OpenFOAM sign convention
        fieldvalues = fact * self.df[fieldname].values
    
        with open(os.path.join(self.dpath,fname),'w') as fid:
            fmt = ['    (%g', '%.12g)',]
            np.savetxt(fid,np.concatenate((ts.reshape((nt,1)),
                                          fieldvalues.reshape((nt,1))
                                          ),axis=1),fmt=fmt)
    
        return


    def write_ICs(self,
                  fname,
                  xmom = 'u',
                  ymom = 'v',
                  temp = 'theta',
                  ):
        """
        Write initial conditions to format that can be included
        in Nalu-Wind input files
    
        Usage
        =====
        fname : str
            Filename
        xmom : str
            Field name corresponding to the x-velocity
        ymom : str
            Field name corresponding to the y-velocity
        temp : str
            Field name corresponding to the potential temperature
        """
        
        # Make copy to avoid SettingwithcopyWarning
        df = self.df.loc[self.datefrom].copy()

        # set missing fields to zero
        fieldNames = [xmom, ymom, temp]
        for field in fieldNames:
            if not field in df.columns:
                df.loc[:,field] = 0.0
    
        # extract time and height array
        zs = df.height.values
        nz = zs.size
    
        # check data is complete
        for field in fieldNames:
            assert ~pd.isna(df[field]).any()
    
        # write data to Nalu-Wind format
        fname = 'ic_heights'
        with open(os.path.join(dpath,fname),'w') as fid:
            fmt = '- %4.2f '
            np.savetxt(fid,zs.reshape((nz,1)),fmt=fmt)
        fname = 'ic_theta'
        with open(os.path.join(dpath,fname),'w') as fid:
            fmt = '- %4.2f '
            np.savetxt(fid,df[temp].values.reshape((nz,1)),fmt=fmt)

        u = df[xmom].values.reshape((nz,))
        v = df[ymom].values.reshape((nz,))
        w = df[temp].values.reshape((nz,))

        fname = 'ic_velocities'
        with open(os.path.join(dpath,fname),'w') as fid:
            fmt = ' - [%4.2f, %4.2f, %4.2f] '
            np.savetxt(fid,np.column_stack((u,v,w)),fmt=fmt)


        return


    def write_timeheight(self,
                         fname,
                         xmom = 'u',
                         ymom = 'v',
                         zmom = 'w',
                         temp = 'theta',
                         ):
        """
        Write time-height data to format usable in Nalu-Wind input file
    
        Usage
        =====
        fname : str
            Filename
        xmom : str
            Field name corresponding to x momentum (field or tendency)
        ymom : str
            Field name corresponding to y momentum (field or tendency)
        zmom : str
            Field name corresponding to z momentum (field or tendency)
        temp : str
            Field name corresponding to potential temperature (field or tendency)
        """
    
        # extract time and height array
        zs = self.df.height.unique()
        ts = self.df.t_index.unique()
        nz = zs.size
        nt = ts.size
        
        ts2 = np.subtract(ts, ts[0])
        ts2 = np.reshape(ts2, (nt,1))
    
        # set missing fields to zero
        fieldNames = [xmom, ymom, zmom, temp]
        for field in fieldNames:
            if not field in self.df.columns:
                self.df.loc[:,field] = 0.0
    
        # pivot data to time-height arrays
        df_pivot = self.df.pivot(columns='height',values=fieldNames)
        # check data is complete
        for field in fieldNames:
            assert ~pd.isna(df_pivot[field]).any().any()


        srcMomX = self.df_pivot[xmom].values
        print(np.shape(srcMomX))

        with open(os.path.join(dpath,'abl_forcing_velocity_x'),'w') as fid:
            fmt = "- [ %g, " + "%.12g,"*(nz-1) + "%.12g]"
            np.savetxt(fid, np.block([ts2,srcMomX]),fmt=fmt)


        srcMomY = self.df_pivot[ymom].values
        print(np.shape(srcMomY))

        with open(os.path.join(dpath,'abl_forcing_velocity_y'),'w') as fid:
            fmt = "- [ %g, " + "%.12g,"*(nz-1) + "%.12g]"
            np.savetxt(fid, np.block([ts2,srcMomY]),fmt=fmt)


        srcMomZ = self.df_pivot[zmom].values
        print(np.shape(srcMomZ))

        with open(os.path.join(dpath,'abl_forcing_velocity_z'),'w') as fid:
            fmt = "- [ %g, " + "%.12g,"*(nz-1) + "%.12g]"
            np.savetxt(fid, np.block([ts2,srcMomZ]),fmt=fmt)


        srcTemp = self.df_pivot[temp].values
        print(np.shape(srcTemp))

        with open(os.path.join(dpath,'abl_forcing_temperature'),'w') as fid:
            fmt = "- [ %g, " + "%.12g,"*(nz-1) + "%.12g]"
            np.savetxt(fid, np.block([ts2,srcTemp]),fmt=fmt)

    
        return
