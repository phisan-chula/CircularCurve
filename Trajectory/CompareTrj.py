#
#
# CompareTrj : compare trajectories from an MMS resulting for various reference
#              base stations from differenct distances 0 km .. 50 km
#
import pandas as pd
import geopandas as gpd
import shapely
import pdal
import json
from pyproj import Geod
from shapely.geometry import LineString
from pyogrio import read_dataframe
from pathlib import Path

DIST_KM = { 'BASE'    : [ 'GNSS01','GNSS02','SBKK','PKKT','BPLE','OKRK'], 
            'dist_km' : [     0,       0,      8,    20,    30,    47  ]  }

COLUMNS = ['GPSTime', 'GPSTime.1', 'Latitude', 'Longitude', 'H-Ell', 'Roll',
       'Pitch', 'Heading', 'SDEast', 'SDNorth', 'SDHeight', 'RollSD',
       'PitchSD', 'HdngSD', 'AmbStatus', 'Q']

class CmpTrajectory:
    def __init__(self, SYS):
        self.SYS = SYS
        self.CACHE = Path('./CACHE')
        self.CACHE.mkdir( parents=True, exist_ok=True)
        self.PLOT = self.CACHE / f'COMPARE_{SYS.INSTRU}.gpkg'
        trj = list()
        for Speed in [30,60]: # kmh
            for Base in SYS.BASE:    
                FILE = self.TRJ.format( **{'Instru':SYS.INSTRU, 'Speed':Speed, 'Base':Base } )
                df = self.ReadTrj( FILE ) 
                #import pdb; pdb.set_trace()
                #df[['INSTRU','SPEED','BASE']]=[SYS.INSTRU,Speed,Base]
                df['INSTRU']=SYS.INSTRU  ; df['SPEED']=Speed;  df['BASE']=Base
                trj.append( df )
        self.dfTRJ = pd.concat( trj )

    def CalcAccuDist(self, df ):
        g = Geod( ellps='WGS84')
        dist_m = list() 
        for i in range( len(df) ):
            if i==0: 
                acc_dist = 0.0
            else:
                _,_,dist  = g.inv( df.iloc[i-1].Longitude, df.iloc[i-1].Latitude,
                                   df.iloc[i  ].Longitude, df.iloc[i  ].Latitude  ) 
                acc_dist += dist
            dist_m.append( acc_dist )
        df['dist_m'] = dist_m

    def MakeRefTrajectory( self, SPEED ):
        print( f'-----> MakeRefTrajectory( {SPEED}) ...')
        self.dfREF = self.dfTRJ[ (self.dfTRJ.SPEED==SPEED) & (self.dfTRJ.BASE==self.SYS.BASE[0])].copy()
        #self.CalcAccuDist( self.dfREF )
        self.dfREF.reset_index( drop=True, inplace=True )
        LS = LineString( self.dfREF[['Longitude','Latitude']].to_numpy() )
        self.gdfRefLS = gpd.GeoDataFrame( crs='EPSG:4326', geometry=[ LS, ] )
        print( f'Plotting reference trajectory speed={SPEED}kmh ...')
        self.gdfRefLS.to_file( self.PLOT , driver='GPKG' , layer=f'RefTraj_{SPEED}kmh' )

    def DoCompare(self):
        VC = self.dfTRJ[['SPEED','BASE']].value_counts()
        print( VC )
        diffs = list()
        for spd,spd_grp in self.dfTRJ.groupby('SPEED'):
            self.MakeRefTrajectory(spd)
            for base,base_grp in spd_grp.groupby('BASE'):
                gdf = gpd.GeoDataFrame( base_grp , crs='EPSG:4326', 
                          geometry=gpd.points_from_xy(base_grp.Longitude,base_grp.Latitude ) )
                #self.CalcAccuDist( gdf )
                gdf.reset_index( drop=True, inplace=True )
                self.MakeDiff( gdf )
                cols = ['SPEED','BASE','hor_mean', 'hor_std', 'hor_max', 'ver_mean','ver_std','ver_max']
                data = [spd,base] + gdf.HorDiff.describe()[['mean','std','max']].to_list() +\
                                    gdf.VerDiff.describe()[['mean','std','max']].to_list()
                df_diff = pd.DataFrame( [data], columns=cols )
                diffs.append( df_diff )
                print(f'Plotting {self.PLOT} v={spd}kmh base={base}...' )
                gdf.to_file( self.PLOT , driver='GPKG' , layer=f'v{spd}_{base}' ) 
        self.dfDIFF = pd.concat( diffs )
        self.dfDIFF = self.dfDIFF.merge( pd.DataFrame( DIST_KM ), on='BASE' )
        self.dfDIFF = self.dfDIFF.sort_values( by=['SPEED','dist_km'], ascending=True )
        #import pdb; pdb.set_trace()
    
    def MakeDiff(self, gdf ):
        def _MakeDiff(row, RefLS, RefTrj):
            lsdist = RefLS.project( row.geometry, normalized=False )
            pnt = RefLS.interpolate( lsdist, normalized=False )
            hor_diff = 111_000*shapely.distance( pnt, row.geometry ) 
            ver_diff = row['H-Ell'] - RefTrj.iloc[ row.name ]['H-Ell']
            #import pdb; pdb.set_trace()
            return hor_diff,ver_diff
        if 1:  
            LS = self.gdfRefLS.iloc[0].geometry 
            gdf[['HorDiff','VerDiff']] = gdf.apply( _MakeDiff, axis=1,
                                        result_type='expand', args=(LS,self.dfREF) )
        else:
            print( f'***DEBUG*** MakeDiff(self, gdf )')
            gdf[['HorDiff','VerDiff']] = 1.0,1.0  # debug !!!
            import pdb; pdb.set_trace()

#######################################################################################
class TrajectoryAU20( CmpTrajectory ):
    def __init__( self, SYS ):
        self.TRJ = './Trajectory_POC_MLS/{Instru}/Speed_{Speed:}/{Base:}/trajectory.csv'
        super().__init__( SYS )

    def ReadTrj( self, TRJFILE ):
        print( f'Reading {TRJFILE} ...' )
        df = pd.read_csv( TRJFILE, skiprows=[1,] )
        print( f'Epochs : {len(df):,} ...' )
        print( f'Columns : {len(df.columns)} ...' )
        return df

#######################################################################################
class TrajectoryMX9( CmpTrajectory ):
    def __init__( self, SYS ):
        self.TRJ = './Trajectory_POC_MLS/{Instru}/Speed_{Speed:}/{Base:}/Speed_{Speed:}_{Base:}.csv'
        super().__init__( SYS )

    def ReadTrj( self, TRJFILE ):
        HDR = " TIME, DISTANCE, EASTING, NORTHING, ELLIPSOID HEIGHT, LATITUDE, LONGITUDE, ELLIPSOID HEIGHT, ROLL, PITCH, HEADING, EAST VELOCITY, NORTH VELOCITY, UP VELOCITY, EAST SD, NORTH SD, HEIGHT SD, ROLL SD, PITCH SD, HEADING SD"
        HDR = [item.strip() for item in HDR.split(',')]
        print( f'Reading {TRJFILE} ...' )
        if TRJFILE[-17:]=='Speed_30_BPLE.csv':  # mistake !!!
            df = pd.read_csv( TRJFILE, skiprows=28, header=None )
        else:
            df = pd.read_fwf( TRJFILE, skiprows=28, header=None )
        df.columns = HDR
        #import pdb; pdb.set_trace()
        duplicate_columns = df.columns.duplicated()
        df = df.loc[:,~duplicate_columns]
        df = df.iloc[::32]  # sampling
        mapping = { 'TIME': 'GPSTime','LATITUDE':'Latitude','LONGITUDE':'Longitude',
                    'ELLIPSOID HEIGHT':  'H-Ell' }
        df.rename( columns=mapping , inplace=True )
        print( f'Epochs : {len(df):,} ...' )
        print( f'Columns : {len(df.columns)} ...' )
        return df

#######################################################################################
class TrajectoryM2X( CmpTrajectory):
    def __init__( self, SYS):
        self.TRJ = './Trajectory_POC_MLS/{Instru}/trjspeed{Speed:}_{Base:}/SBET_GTGV0M0.OUT'
        super().__init__( SYS )

    def ReadTrj( self, TRJFILE ):
        print( f'Reading {TRJFILE} ...' )
        pipeline = { "pipeline": [ { "type": "readers.sbet", "filename": TRJFILE } ] }
        pipeline = pdal.Pipeline(json.dumps(pipeline))
        try:
            pipeline.execute()
        except:
            import pdb; pdb.set_trace()
        arrays = pipeline.arrays[0]
        df = pd.DataFrame(arrays)
        print( f'Epochs : {len(df):,} ...' )
        print( f'Columns : {len(df.columns)} ...' )
        mapping = { 'GpsTime': 'GPSTime','Y':  'Latitude','X': 'Longitude','Z':'H-Ell' }
        df.rename( columns=mapping , inplace=True )
        #import pdb; pdb.set_trace()
        df = df.iloc[::16].copy()
        return df

#############################################################################
#############################################################################
#############################################################################
SYS_AU20 = pd.Series( { 'INSTRU' : 'AU20',   #  AU20
           'BASE' : ['GNSS01','SBKK','PKKT','BPLE','OKRK'] } ) # first BASE will referennced

SYS_MX9 = pd.Series( { 'INSTRU' : 'MX9',   #   MX9 
           'BASE' : ['GNSS02','SBKK','PKKT','BPLE','OKRK'] } ) 

SYS_M2X = pd.Series( { 'INSTRU' : 'ScoutM2X',   #  M2X
           'BASE' : ['GNSS01','SBKK', 'BPLE','OKRK'] } )   # PKKT error!

import argparse
parser = argparse.ArgumentParser( description=\
        'read trajectories in varios formats and compare with i'\
        'the slowest speed and nearest base reference station')
parser.add_argument("--au20",action="store_true", help="CHC AU20")
parser.add_argument("--mx9", action="store_true", help="Trimble MX9")
parser.add_argument("--m2x", action="store_true", help="Phoenix M2X")
args = parser.parse_args()
print( args )
if set(vars(args).values())==set([False,]):
    parser.print_help()
#import pdb ; pdb.set_trace()

if args.au20:
    cmp = TrajectoryAU20( SYS_AU20 )
    cmp.DoCompare()
    print(cmp.dfDIFF.to_markdown( floatfmt='.3f' ))
if args.mx9:
    cmp = TrajectoryMX9( SYS_MX9 )
    cmp.DoCompare()
    print(cmp.dfDIFF.to_markdown( floatfmt='.3f' ))
if args.m2x:  
    cmp = TrajectoryM2X( SYS_M2X )
    cmp.DoCompare( )
    print(cmp.dfDIFF.to_markdown( floatfmt='.3f' ))





