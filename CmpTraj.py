#
#
#
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString
from pyogrio import read_dataframe
from pathlib import Path

class CmpTrajectory:
    def __init__(self):
        self.INSTRU = 'AU20'  #  AU20, MX9 , M2X
        self.CACHE = Path('./CACHE')
        self.CACHE.mkdir( parents=True, exist_ok=True)
        self.PLOT = self.CACHE / f'COMPARE_AU20.gpkg'
        self.BASE = ['GNSS01','SBKK','PKKT','BPLE','OKRK']
        TRJ = './Trajectory_POC_MLS/{Instru}/Speed_{Speed:}/{Base:}/trajectory.csv'
        trj = list()
        for Speed in [30,60]: # kmh
            for Base in self.BASE:    
                FILE = TRJ.format( **{'Instru':self.INSTRU, 'Speed':Speed, 'Base':Base } )
                df = self.ReadTrj( FILE ) 
                df[['INSTRU','SPEED','BASE']]=[self.INSTRU,Speed,Base]
                trj.append( df )
        self.dfTRJ = pd.concat( trj )
        self.MakeReferenceTrajectory()

    def ReadTrj( self, TRJFILE ):
        print( f'Reading {TRJFILE} ...' )
        df = pd.read_csv( TRJFILE, skiprows=[1,] )
        print( f'Epocs : {len(df)} ...' )
        print( f'Columns : {len(df.columns)} ...' )
        return df

    def MakeReferenceTrajectory( self ):
        for spd in [30,60]:
            self.REF = self.dfTRJ[ (self.dfTRJ.SPEED==spd) & (self.dfTRJ.BASE=='GNSS01') ].copy()
            LS = LineString( self.REF[['Longitude','Latitude']].to_numpy() )
            gdf = gpd.GeoDataFrame( crs='EPSG:4326', geometry=[ LS, ] )
            print( f'Plotting reference trajectory speed={spd}kmh ...')
            gdf.to_file( self.PLOT , driver='GPKG' , layer=f'RefTraj_{spd}kmh' )
            #import pdb; pdb.set_trace()

    def DoCompare(self):
        VC = self.dfTRJ[['SPEED','BASE']].value_counts()
        print( VC )

        for [spd,bas],grp in self.dfTRJ.groupby( ['SPEED','BASE'] ):
            gdf = gpd.GeoDataFrame( grp , crs='EPSG:4326', geometry=gpd.points_from_xy( grp.Longitude, grp.Latitude ) )
            print(f'Plotting {self.PLOT} v={spd}kmh base={bas}...' )
            gdf.to_file( self.PLOT , driver='GPKG' , layer=f'v{spd}_{bas}' ) 
        #    import pdb; pdb.set_trace()
        return

#############################################################################
cmp = CmpTrajectory()
cmp.DoCompare()



