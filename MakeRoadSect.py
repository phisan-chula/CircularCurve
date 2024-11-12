#
# MakeRoadSection.py : RoadSection reads alignment of MMS pathway and generate polygon of 
#               1-km (default) sections. The section will be buffered BUF meter to both-side. 
#               Polygon of section will be use for partitioning point-cloud and images 
#               into smaller chunks
#
import argparse
import tomllib
import shutil
import pyogrio
import numpy as np
import pandas as pd
import geopandas as gpd
from sklearn.cluster import DBSCAN
from shapely import to_wkt
from shapely.geometry import LineString
from shapely.ops import substring
from pathlib import Path
from simplekml import Kml, Style 
from simplekml import Polygon as kmlPoly
from itertools import cycle
from pathlib import Path
gpd.options.io_engine = "pyogrio"

class Section:
    def __init__(self, DATA, LS):
        self.LS = LS
        self.DATA = DATA
        self.dfSTA    = self.MakeStation(self.DATA.DIV)  # station for data section
        self.dfSTA100 = self.MakeStation(100)  # for visualization
        self.dfTILE = self.MakeSection()

    def MakeStation( self, DIV ):
        TRUNCATE = 10  # meter
        ndiv,rest = divmod( self.DATA.START_SECT+self.LS.length ,DIV ) 
        dist_m = np.linspace( 0, ndiv*DIV, num=int(ndiv)+1, endpoint=True ) 
        if rest>TRUNCATE: 
            dist_m = list(dist_m) + [ self.DATA.START_SECT + self.LS.length ]
        else:
            dist_m = list(dist_m)
        ndiv,rest = divmod( self.DATA.START_SECT, DIV)
        if (DIV-rest)>TRUNCATE:
            dist_m.append( self.DATA.START_SECT )
        df = pd.DataFrame( dist_m , columns=['dist_m'] ).sort_values(by='dist_m')
        df = df[df.dist_m>=self.DATA.START_SECT].copy()
        df.drop_duplicates( 'dist_m', keep='first', inplace=True, ignore_index=True )
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
            buf = self.DATA.BUFFER
            sect = substring( self.LS, fr.ls_dist, to.ls_dist , normalized=False)
            geom  = sect.buffer( buf, cap_style='flat' )         
            tiles.append( [ fr.Name, to.Name, buf, geom  ] )
        df = pd.DataFrame( tiles, columns=[ 'BEG','END', 'buffer' ,'geometry' ] )      
        gdf = gpd.GeoDataFrame( df, crs=self.DATA.EPSG, geometry=df.geometry )
        return gdf

class MMS_Route:
    def __init__(self, DATA):
        self.DATA = DATA
        self.CACHE = Path('./CACHE')
        self.REP_DIR =  f'{self.CACHE}/{self.DATA.INSTRU}_{self.DATA.ACQ_DATE}'
        self.CACHE.mkdir( parents=True, exist_ok=True)
        Path( self.REP_DIR).parent.mkdir(parents=True, exist_ok=True)

        gdfCL = self.ReadKML_Valid()

        route = list()
        for i in range(len(gdfCL)):
            if not isinstance( gdfCL.iloc[i].geometry, LineString ):
                print( gdfCL ); raise '***ERROR*** must be LineString'
            LS =  gdfCL.iloc[i].geometry   # current center line 
            LS = LineString([(x, y) for x, y, z in LS.coords])
            sect = Section( DATA,LS )  #  generate
            data = {'NAME': gdfCL.iloc[i]['Name'] ,'LS': LS, 'STA':sect.dfSTA, 'STA100':sect.dfSTA100,
                     'TILE': sect.dfTILE }
            route.append(data)
        self.dfROUTE = pd.DataFrame( route )
        for rt,row in self.dfROUTE.iterrows():
            print( f'{30*"="} route : "{row.NAME}" {30*"="}')
            print( row.STA )
    
    def ReadKML_Valid(self):
        gdfCL = pyogrio.read_dataframe( self.DATA.CL_KML )
        self.DATA['EPSG'] = gdfCL.estimate_utm_crs() # auto utm
        gdfCL = gdfCL.to_crs( self.DATA.EPSG )
        gdfCL['Length'] = gdfCL.length
        print( '============ KML of MMS runs ============')
        print( gdfCL )
        print( 'KML: must have Name = { FWD | REV } and LINESTRING().....***')
        print( f'KML: lengths std = {gdfCL.Length.std():.1f} m')
        if 'LS_FLIP' in self.DATA.keys(): 
            for route in self.DATA.LS_FLIP:
                print( f'Found LS_FLIP , flipping "{route}" ...' )
                idx = gdfCL[gdfCL.Name==route].index[0]
                flipLS = gdfCL[gdfCL.Name==route].iloc[0].geometry.reverse()
                gdfCL.at[ idx,'geometry'] = flipLS
        begs=list();ends=list()
        for i in range( len(gdfCL) ):
            begs.append( gdfCL.iloc[i].geometry.coords[0] )
            ends.append( gdfCL.iloc[i].geometry.coords[-1] )
        clust_beg = DBSCAN( eps=200, min_samples=2).fit( np.array(begs) )
        clust_end = DBSCAN( eps=200, min_samples=2).fit( np.array(ends) )
        print(  'Check if begin/end nodes are clustered!...')
        print( f'Check begin nodes are clustered : {np.all(clust_beg.labels_==0)}...')
        print( f'Check end   nodes are clustered : {np.all(clust_end.labels_==0)}...')
        return gdfCL

    def PlotKML( self ):
        kml = Kml(name="Sections for MMS Routes", visibility=1 )
        STA = Style()
        STA.labelstyle.scale = 1  # Text twice as big
        STA.labelstyle.color = 'ff0000ff' 
        STA.iconstyle.color = 'ff0000ff'  
        STA.iconstyle.scale = 1  # Icon thrice as big
        STA.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/flag.png'
        STA100 = Style()
        STA100.labelstyle.scale = 0.5  # Text twice as big
        STA100.iconstyle.scale  = 0.5  # Icon thrice as big
        STA100.iconstyle.icon.href = 'https://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
        colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a']
        poly_cyc = list()
        for i in range( len(colors) ):
            POLY = Style()
            POLY.linestyle.color = colors[i]
            POLY.polystyle.color = colors[i].replace('#','B2')  # opacity 50%
            POLY.polystyle.fill = 1  
            POLY.polystyle.outline = 1 
            poly_cyc.append( POLY )
        #import pdb ; pdb.set_trace()
        POLY_CYC = cycle(poly_cyc)

        for i,row in self.dfROUTE.iterrows():
            fold_name = kml.newfolder( name=row['NAME'])
            fold_sect = fold_name.newfolder( name='Section')
            tile_w84 = row.TILE.to_crs( crs='OGC:CRS84' ).copy()
            for i,tile in tile_w84.iterrows():
                coords =  [ [point[0],point[1]] for point in tile.geometry.exterior.coords]
                poly = fold_sect.newpolygon( name=f'{tile.BEG}-{tile.END}', outerboundaryis=coords)
                poly.style = next(POLY_CYC)
            fol_sta = fold_name.newfolder( name='STATION', visibility=1 )
            sta = row.STA.to_crs( crs='OGC:CRS84' ).copy()
            for i,rw in sta.iterrows():
                pnt = fol_sta.newpoint( name=rw['Name'], coords=[(rw.geometry.x,rw.geometry.y)] )
                pnt.style = STA
            fol_sta100 = fold_name.newfolder( name='Sta100', visibility=0 )
            sta100 = row.STA100.to_crs( crs='OGC:CRS84' ).copy()
            for i,rw in sta100.iterrows():
                pnt = fol_sta100.newpoint( name=rw['Name'], coords=[(rw.geometry.x,rw.geometry.y)] )
                pnt.style = STA100
        MMS_Sections = Path( self.REP_DIR ) / 'ROUTE_SECTIONS.kml' 
        print( f'Plotting {MMS_Sections} ...' )
        kml.save( MMS_Sections )

    def MakeFileStruct(self):
        def epsg_fn():
            return f'EPSG_{self.DATA.EPSG.to_epsg()}'
        REP_DIR = self.REP_DIR
        Path(REP_DIR).mkdir(parents=True, exist_ok=True)
        for rep in ['Trajectory.out', 'Trajectory.csv', 'Trajectory.trj', 'MMS_MissionReport.pdf' ]:
            (Path(REP_DIR) / rep ).touch()
        #shutil.copy( Path(self.TILE_KML) , Path( f'{REP_DIR}/{Path(self.TILE_KML).name}' )  )
        shutil.copy( self.DATA.CL_KML , Path( f'{REP_DIR}/{self.DATA.CL_KML}' )  )

        PNC_DIR =  f'{self.CACHE}/{self.DATA.INSTRU}_{self.DATA.ACQ_DATE}/PntCloud'
        Path(PNC_DIR).mkdir(parents=True, exist_ok=True)
        for i,row_rt in self.dfROUTE.iterrows():
            with open( f'{PNC_DIR}/{epsg_fn()}.wkt',"w") as fd:
                fd.write( f'{self.DATA.EPSG.to_wkt(pretty=True)}' )
            for j,row_tile in row_rt.TILE.iterrows():
                TILE_NAME = f'{PNC_DIR}/{row_tile.BEG}_{row_tile.END}_{row_rt.NAME}'
                print( f'Writing PNTCLOUD {TILE_NAME} wkt|bnd|las ...' )
                with open( f'{TILE_NAME}.bnd', "w") as fd:
                    fd.write( to_wkt( row_tile.geometry, rounding_precision=2, trim=False ) )
                Path(f'{TILE_NAME}.las').touch()

        IMG_DIR =  f'{self.CACHE}/{self.DATA.INSTRU}_{self.DATA.ACQ_DATE}/Images'
        PANO_CNT = 1 
        Path(IMG_DIR).mkdir(parents=True, exist_ok=True)
        for i,row_rt in self.dfROUTE.iterrows():
            for j,row_tile in row_rt.TILE.iterrows():
                IMG_DIR_TILE = f'{IMG_DIR}/{row_tile.BEG}_{row_tile.END}_{row_rt.NAME}'
                Path(IMG_DIR_TILE).mkdir(parents=True, exist_ok=True)
                print( f'Writing IMAGES {IMG_DIR_TILE}/xxxxx.jpg ..' )
                #import pdb ; pdb.set_trace()
                with open( f'{IMG_DIR_TILE}/{epsg_fn()}.wkt',"w") as fd:
                    fd.write( f'{self.DATA.EPSG.to_wkt(pretty=True)}' )
                with open( f'{IMG_DIR_TILE}/{epsg_fn()}.bnd',"w") as fd:
                    fd.write( to_wkt( row_tile.geometry, rounding_precision=2, trim=False ) )
                for k in range(5):
                    JPG = f'{IMG_DIR_TILE}/img{PANO_CNT:05d}.jpg'
                    print( f'Toching .... {JPG}' )
                    Path(JPG).touch()
                    PANO_CNT += 1

###############################################################################
if __name__=="__main__":
    parser = argparse.ArgumentParser(
            description="Make MMS sections every 1-km, output CACHCE/* structure. ")
    parser.add_argument("toml", help="specified TOML file, KML file define fwd/rev MMS routes" )
    args = parser.parse_args()

    toml_file = Path( args.toml )
    with open( toml_file ,'rb' ) as f:
        #import pdb ; pdb.set_trace()
        TOML = pd.Series( tomllib.load(f) )
        TOML['CL_KML'] = toml_file.with_suffix('.kml')
    mms =  MMS_Route( TOML )
    mms.MakeFileStruct()
    mms.PlotKML()
    print('************ end of RoadTile *************' ) 
    #import pdb ; pdb.set_trace()
