#
# MakeRoadSection.py : RoadSection reads alignment of MMS pathway and generate polygon of 
#               1-km (default) sections. The section will be buffered BUF meter to both-side. 
#               Polygon of section will be use for partitioning point-cloud and images into smaller chunks
#
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely import to_wkt
from shapely.geometry import LineString
from shapely.ops import substring
from pathlib import Path
from simplekml import Kml, Style 
from simplekml import Polygon as kmlPoly
import pyogrio
gpd.options.io_engine = "pyogrio"

class Section:
    def __init__(self, DATA, LS):
        self.LS = LS
        self.DATA = DATA
        self.dfSTA = self.MakeStation()
        self.dfTILE = self.MakeSection()
        #import pdb ; pdb.set_trace()

    def MakeStation( self ):
        LS_LEN = self.LS.length
        ndiv,rest = divmod( self.DATA.START_SECT + LS_LEN , self.DATA.DIV ) 
        dist_m = np.linspace( 0, ndiv*self.DATA.DIV, num=int(ndiv)+1, endpoint=True ) 
        if rest>10: dist_m = list(dist_m) + [ self.DATA.START_SECT + LS_LEN ]
        ndiv,rest = divmod( self.DATA.START_SECT, self.DATA.DIV)
        if (self.DATA.DIV-rest>10):  dist_m.append(  self.DATA.START_SECT )
        df = pd.DataFrame( dist_m , columns=['dist_m'] ).sort_values(by='dist_m')
        df = df[df.dist_m>=self.DATA.START_SECT].copy()
        df.reset_index( inplace=True, drop=True )
        df['ls_dist'] = df['dist_m']-self.DATA.START_SECT
        def mkpnt(row, DATA,LS):
            pnt = LS.interpolate( row.ls_dist ) 
            km,rest = divmod( row.dist_m,DATA.DIV)
            sta = f'{int(km):03.0f}+{rest:03.0f}'
            return sta,pnt
        df[['Name', 'geometry']] = df.apply( mkpnt, axis=1, result_type='expand',
                                   args=(self.DATA,self.LS ) )
        gdf = gpd.GeoDataFrame( df, crs=self.DATA.EPSG, geometry=df.geometry )
        return gdf

    def MakeSection( self ):
        tiles = list()
        for i in range( len(self.dfSTA)-1 ):
            fr  = self.dfSTA.iloc[i  ]
            to  = self.dfSTA.iloc[i+1]
            sect = substring( self.LS, fr.ls_dist, to.ls_dist , normalized=False)
            box  = sect.buffer( 30, cap_style='flat' )         
            tiles.append( [ fr.Name, to.Name, box  ] )
        df = pd.DataFrame( tiles, columns=[ 'BEG','END' ,'geometry' ] )      
        gdf = gpd.GeoDataFrame( df, crs=self.DATA.EPSG, geometry=df.geometry )
        return gdf

class MMS_Route:
    def __init__(self):
        self.CACHE = Path('CACHE')
        if not self.CACHE.exists():
            self.CACHE.mkdir( parents=True, exist_ok=True)
        DATA = pd.Series( { 'INSTRU':'M2X_1', 'ACQ_DATE':'20240621', 'DIV': 1000,  
            'ROUTE_NAME': 'DOH_351_PS_MANOO', 'CL_KML': 'DOH_351_CL.kml' , 'START_SECT' : 9_300 } )
        gdf = pyogrio.read_dataframe( DATA.CL_KML )
        route = list()
        for i in range(len(gdf)):
            if not isinstance( gdf.iloc[i].geometry, LineString ):
                print( gdf ); raise '***ERROR*** must be LineString'
            DATA['EPSG'] = gdf.estimate_utm_crs() # utm
            gdf = gdf.to_crs( DATA.EPSG )
            LS =  gdf.iloc[i].geometry   # current center line 
            LS = LineString([(x, y) for x, y, z in LS.coords])
            sect = Section( DATA,LS )  #  generate
            #import pdb ; pdb.set_trace()
            data = {'NAME': gdf.iloc[i]['Name'] ,'LS': LS, 'STA':sect.dfSTA, 'TILE': sect.dfTILE }
            route.append(data)
        self.dfROUTE = pd.DataFrame( route )
        self.DATA = DATA

    def PlotKML( self ):
        kml = Kml(name="Sections for MMS Routes", visibility=1 )
        STA = Style()
        STA.labelstyle.scale = 1  # Text twice as big
        STA.labelstyle.color = 'ff0000ff' 
        STA.iconstyle.color = 'ff0000ff'  
        STA.iconstyle.scale = 1  # Icon thrice as big
        STA.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/flag.png'
        POLY = Style()
        POLY.linestyle.color = 'ff0000ff'  # Red
        POLY.polystyle.color = '320000ff'  # Red
        POLY.polystyle.fill = 1  
        POLY.polystyle.outline = 1 

        for i,row in self.dfROUTE.iterrows():
            fold_name = kml.newfolder( name=row['NAME'])
            fold_sect = fold_name.newfolder( name='Section')
            tile_w84 = row.TILE.to_crs( crs='OGC:CRS84' ).copy()
            for i,tile in tile_w84.iterrows():
                coords =  [ [point[0],point[1]] for point in tile.geometry.exterior.coords]
                poly = fold_sect.newpolygon( name=f'{tile.BEG}-{tile.END}', outerboundaryis=coords)
                poly.style=POLY
            fol_sta = fold_name.newfolder( name='STATION', visibility=1 )
            sta_w84 = row.STA.to_crs( crs='OGC:CRS84' ).copy()
            for i,row in sta_w84.iterrows():
                pnt = fol_sta.newpoint( name=row['Name'], coords=[(row.geometry.x,row.geometry.y)] )
                pnt.style = STA
        KML =  f'{self.CACHE}/{self.DATA.INSTRU}_{self.DATA.ACQ_DATE}/{self.DATA.ROUTE_NAME}.kml'
        Path(KML).parent.mkdir(parents=True, exist_ok=True)
        print( f'Plotting {KML} ...' )
        kml.save( KML )

    def MakeFileStruct( self):
        REPORT =  f'{self.CACHE}/{self.DATA.INSTRU}_{self.DATA.ACQ_DATE}/'
        for rep in ['Trajectory.sbet', 'Trajectory.csv', 'MMS_System.pdf' ]:
            (Path(REPORT) / rep ).touch()

        PNC =  f'{self.CACHE}/{self.DATA.INSTRU}_{self.DATA.ACQ_DATE}/PntCloud'
        Path(PNC).mkdir(parents=True, exist_ok=True)
        for i,row_rt in self.dfROUTE.iterrows():
            for j,row_tile in row_rt.TILE.iterrows():
                TILE_NAME = f'{PNC}/{row_tile.BEG}_{row_tile.END}_{row_rt.NAME}'
                print( f'Writing PNTCLOUD {TILE_NAME} wkt|las ...' )
                with open( f'{TILE_NAME}.wkt', "w") as fd:
                    fd.write( to_wkt( row_tile.geometry, rounding_precision=2, trim=False ) )
                Path(f'{TILE_NAME}.las').touch()

        IMG =  f'{self.CACHE}/{self.DATA.INSTRU}_{self.DATA.ACQ_DATE}/Images'
        Path(IMG).mkdir(parents=True, exist_ok=True)
        for i,row_rt in self.dfROUTE.iterrows():
            for j,row_tile in row_rt.TILE.iterrows():
                IMG_DIR = f'{IMG}/{row_tile.BEG}_{row_tile.END}_{row_rt.NAME}'
                Path(IMG_DIR).mkdir(parents=True, exist_ok=True)
                print( f'Writing IMAGES {IMG_DIR}/xxxxx.jpg ..' )
                with open( f'{IMG_DIR}/{row_tile.BEG}_{row_tile.END}_{row_rt.NAME}.wkt',"w") as fd:
                    fd.write( to_wkt( row_tile.geometry, rounding_precision=2, trim=False ) )
                for k in list(range(1,5)):
                    JPG = f'{IMG_DIR}/img{k:05d}.jpg'
                    print( f'Toching .... {JPG}' )
                    Path(JPG).touch()

###############################################################################
mms =  MMS_Route()
mms.PlotKML()
mms.MakeFileStruct()
print('************ end of RoadTile *************' ) 
#import pdb ; pdb.set_trace()yy

