#
# EstCureve.py : Estimate cicurlar from Linestring of road centerline
#                Input road centerline might incloud lead-in/lead-out straight
#                lines which come before and after the actural curve.
#
# Author : June,2024 Phisan Santitamnont ( phisan.s@cdg.co.th )
#
import numpy as np
import pandas as pd
import geopandas as gpd
import pyransac3d as pyrsc
from pathlib import Path
from skspatial.objects import Vector,Line
from shapely.wkt import dumps,loads
from shapely import ops,affinity
from shapely.geometry import LineString,Point
from pyogrio import read_dataframe
from CurvePnts import *


def drop_z(geometry):
    return ops.transform(lambda x, y, z=None: (x, y), geometry)

def AngLines(line1,line2):
    dir_vect1 = line1.direction
    dir_vect2 = line2.direction
    return dir_vect1.angle_signed( dir_vect2)

def Line2P_Len( p1, p2, Length ):
    dir_vect = Vector.from_points( p1,p2 ) 
    dir_mag = np.linalg.norm(dir_vect)
    unit_vect = dir_vect / dir_mag
    scaled_vector = unit_vect * Length 
    return Line( p1,scaled_vector)

########################################################################################
class EstimateCurve( CircularCurve ):
    def __init__(self, EPSG, ROAD_SECT, ANALY_DIV=0.5, THRES=0.2 ):
        if ROAD_SECT.has_z: ROAD_SECT=drop_z(ROAD_SECT)
        self.EPSG = EPSG
        self.ROAD_SECT = ROAD_SECT
        self.ANALY_DIV = ANALY_DIV
        self.THRES = THRES
        CURV_ALIGN,center,radius = self.CreateAlignment()
        super().__init__( self.EPSG ,CURV_ALIGN, radius, 2 ) 
        self.ROAD_SECT = gpd.GeoDataFrame( crs=self.EPSG, geometry=[ROAD_SECT,] )

    def CreateAlignment(self):
        center,axis,radius = self.FitCircRANSAC()
        PC = self.gdfInlier.iloc[0].geometry
        PT = self.gdfInlier.iloc[-1].geometry
        O_PC = Line2P_Len( center[:2], PC.coords[0], radius )
        O_PT = Line2P_Len( center[:2], PT.coords[0], radius )
        PC_PT_mid = [ (PC.x+PT.x)/2,(PC.y+PT.y)/2]
        delta = AngLines(O_PC,O_PT)
        print(f'AXIS = {axis}...,  DELTA = {np.degrees(delta)}')
        E = radius/np.cos(delta/2)-radius
        O_PI = Line2P_Len( center[:2], PC_PT_mid, radius+E )
        CURV_ALIGN = LineString([Point(O_PC.to_point()), Point(O_PI.to_point()), Point(O_PT.to_point()) ])
        CURV_ALIGN = affinity.scale( CURV_ALIGN, xfact=1.05, yfact=1.05, origin=CURV_ALIGN.coords[1] )
        #import pdb ; pdb.set_trace()
        return CURV_ALIGN,center,radius

    def FitCircRANSAC( self ):
        LEN = self.ROAD_SECT.length
        pnts = np.linspace( 0, self.ROAD_SECT.length, num=int(LEN/self.ANALY_DIV), endpoint=True )
        dfPnt = pd.DataFrame( pnts, columns=['dist_m'] )
        def MkPnt(row,ROAD_SECT):
            return ROAD_SECT.interpolate( row.dist_m, normalized=False )
        dfPnt['geometry'] = dfPnt.apply( MkPnt, axis=1, result_type='expand', 
                                                        args=(self.ROAD_SECT,) ) 
        #import pdb ; pdb.set_trace()
        gdfPnt = gpd.GeoDataFrame( dfPnt, crs=self.EPSG, geometry=dfPnt.geometry )
        pnts = np.array([[geom.x, geom.y] for geom in gdfPnt.geometry])
        pnts3d = np.insert(pnts, 2, 0, axis=1)
        circ = pyrsc.Circle()
        center,axis,radius,inliers = circ.fit( pnts3d, thresh=self.THRES, maxIteration=1000 )
        gdfPnt['INLIER'] = gdfPnt.index.isin( inliers )
        gdfInlier = gdfPnt[gdfPnt.INLIER==True].copy().reset_index()
        #import pdb ; pdb.set_trace()
        self.gdfInlier = gdfInlier
        return center,axis,radius

    def WriteGIS( self ):
        print( f'EstimateCureve:WriteGIS(): write {self.GPKG} layer Road/InlierPnt' )
        self.ROAD_SECT.to_file( self.GPKG, driver='GPKG', layer='Road' )
        self.gdfInlier.to_file( self.GPKG, driver='GPKG', layer='InlierPnt' )
        super().WriteGIS()

##############################################################
##############################################################
if __name__ == "__main__":
    if 0:
        Path('./CACHE_C').mkdir(parents=True, exist_ok=True)
        df = read_dataframe( 'CurvePrasert.kml' )
        EPSG = df.estimate_utm_crs().to_epsg()  # UTM
        df = df.to_crs( EPSG )
        for i in range( len(df)):
            #import pdb ; pdb.set_trace()
            EC = EstimateCurve( EPSG , df.iloc[i].geometry )
            print( dumps( drop_z(df.iloc[i].geometry), rounding_precision=7) )
            #EC.DoPlot()
            #EC.WriteGIS()
            #Path('./CACHE').rename( Path(f'./CACHE_C/curve_{i}')  )
    else:
        WKT = 'LINESTRING (681119.0450817 1527757.4696346, 681119.3968534 1527757.6504096, 681143.6579172 1527756.5851417, 681159.0136977 1527756.6929437, 681173.4631250 1527759.6319443, 681189.6597877 1527767.1997959, 681200.3085576 1527776.5055463, 681208.2744513 1527787.3780785, 681214.6224477 1527799.4765282, 681214.9790651 1527799.4806670, 681220.5513379 1527819.9152865, 681220.9922170 1527833.4072438, 681219.2690322 1527849.8834369, 681216.6632099 1527865.2979926, 681213.7011911 1527880.3589159, 681208.0542947 1527900.6247882)'
        EC = EstimateCurve( 32647 , loads(WKT) )

        import pdb ; pdb.set_trace()
        #EC.DoPlot()
