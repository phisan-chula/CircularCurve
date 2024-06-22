#
#
PROG="""\
 CurvePnts.py : Generate points on a designated curve from its 3-point alignment.
 Author : Phisan Santitamnont ( phisan.chula@gmail.com )
 History  25 Dec 2022 : version 0.1
          18 Jun 2024 : version 0.5
"""
import sys 
import numpy as np 
import pandas as pd
import geopandas as gpd
from skspatial.objects import Vector
from shapely.geometry import LineString,Point
from shapely.affinity import translate, rotate
from pathlib import Path
import matplotlib.pyplot as plt
from pygeodesy import dms

#####################################################################################
class CircularCurve:
    def __init__(self, EPSG, ALIGN, RADIUS, DIV, ROUND_ABOUT=False ):
        self.CACHE = Path( './CACHE' )
        self.CACHE.mkdir(parents=True, exist_ok=True)
        self.PLOT = self.CACHE.joinpath('Plot_Curve')
        PAR = pd.Series( {'EPSG':EPSG, 'ALIGN' : ALIGN, 'RADIUS' : RADIUS, 
                              'DIV': DIV, 'ROUND_ABOUT': ROUND_ABOUT } )
        assert( len(PAR.ALIGN.coords) ==3 ),'***ERROR*** limit 3 points on LS_ALIGN'
        pc,pi,pt = list(ALIGN.coords)
        vcPC = Vector.from_points( pc,pi )
        vcPI = Vector.from_points( pi,pt )
        sgDEFLEC = vcPI.angle_signed( vcPC )
        PAR['sgDEFL']   = sgDEFLEC
        PAR['DEFL']     = abs( sgDEFLEC )
        PAR['DEFLdms' ] = dms.toDMS( np.degrees(sgDEFLEC), prec=1 )
        print( f'Deflection angle : {np.degrees(PAR.sgDEFL)} ...' )
        print( f'Deflection angle : {PAR.DEFLdms} ...' )
        self.PAR = PAR
        #import pdb ;pdb.set_trace()
        self.GenNormArc()
        self.RotTransNormArc()
        LS_PNT = LineString( list(self.gdfPNT.geometry) )
        LS_PC = LineString( [PAR.ORIGIN,PAR.PC] ) 
        LS_PI = LineString( [PAR.ORIGIN,PAR.PI] ) 
        LS_PT = LineString( [PAR.ORIGIN,PAR.PT] ) 
        self.dfLS = gpd.GeoDataFrame( {'Type': ['Alignment', 'POC', 'O_PC', 'O_MO', 'O_PT'] },
                                   crs=EPSG, geometry=([PAR.ALIGN,LS_PNT,LS_PC,LS_PI,LS_PT ]) ) 

    def GenNormArc(self):
        ''' curvature at origin (0,0) draw arc clock-wise '''
        PAR = self.PAR
        PAR['TL']   = PAR.RADIUS*np.tan( PAR.DEFL/2 )   # tangent length
        PAR['LENCUR'] = PAR.RADIUS*PAR.DEFL          
        pc,pi,pt = list(PAR.ALIGN.coords)
        pi_pc = LineString( [pi,pc] ).length
        assert( pi_pc>=PAR.TL  ),'***ERROR** leadin PC too shore!'
        pi_pt = LineString( [pi,pt] ).length
        assert( pi_pt>=PAR.TL ),'***ERROR** leadout PT too shore!'
        #import pdb ;pdb.set_trace()
        if PAR.ROUND_ABOUT:
            PAR.LENCUR = 2*np.pi*PAR.RADIUS - PAR.LENCUR
        ndiv,rest = divmod(PAR.LENCUR, PAR.DIV)
        pnt_div = np.linspace(rest/2,PAR.LENCUR-rest/2,num=int(ndiv)+1,endpoint=True )
        pnts = np.concatenate( [np.array([0]), pnt_div, np.array([PAR.LENCUR]) ] )
        dfPNT = pd.DataFrame( pnts, columns=['cvDist'] )
        def DoCurve(row, PAR):
            theta = row.cvDist/PAR.RADIUS
            x,y = PAR.RADIUS*np.sin(theta),PAR.RADIUS*np.cos(theta)
            if PAR.sgDEFL<0.   : y=-y
            if PAR.ROUND_ABOUT : x=-x
            return [ f'{row.cvDist:03.0f}', f'{row.cvDist:.3f}', Point(x,y) ]
        dfPNT[['Name','cvDist', 'geometry']] = dfPNT.apply( DoCurve, 
                                  axis=1, result_type='expand', args=(PAR,) )
        gdfPNT = gpd.GeoDataFrame( dfPNT, crs=PAR.EPSG, geometry=dfPNT.geometry )
        self.gdfPNT=gdfPNT
        ############################################
        PC = LineString( [ pi,pc ]).interpolate( PAR.TL, normalized=False )
        PT = LineString( [ pi,pt ]).interpolate( PAR.TL, normalized=False )
        PAR['PC'] = PC
        PAR['PI'] = Point(pi) 
        PAR['PT'] = PT

    def RotTransNormArc(self):
        PAR = self.PAR
        p000  = self.gdfPNT.iloc[0].geometry
        dx,dy = PAR.PC.x-p000.x, PAR.PC.y-p000.y
        PC_PI = Vector.from_points(list(PAR.PC.coords[0]) ,list(PAR.PI.coords[0]) )
        sgRot = Vector([1,0]).angle_signed( PC_PI )
        def RotTr(row,dx,dy,sgRot,pntPC):
            x,y   = row.geometry.x, row.geometry.y
            x_,y_ = x+dx , y+dy
            return rotate( Point( x_,y_),sgRot, origin=(pntPC.x,pntPC.y), use_radians=True ) 
        self.gdfPNT['geometry'] = self.gdfPNT.apply( RotTr, axis=1, 
                            result_type='expand', args=( dx,dy,sgRot,PAR.PC ) )
        PAR['ORIGIN'] = rotate( Point( dx,dy),sgRot, origin=(PAR.PC.x,PAR.PC.y), use_radians=True ) 
        PAR['MO']  = LineString( [PAR.ORIGIN,PAR.PI] ).interpolate(PAR.RADIUS,normalized=False )

    def DoPlot(self, SUFFIX=None ):
        fig, ax = plt.subplots( figsize=(20,18))
        self.dfLS.plot( ax=ax ) 
        for i,row in self.gdfPNT.iterrows():
            x =  row.geometry.x ; y =  row.geometry.y
            ax.scatter( x, y, c='k', s=30, alpha=0.5 )
            ax.text( x,y, s=row['Name'], c='g', fontsize=15 )
        for pnt in ['PC','PI','PT','ORIGIN','MO']: 
            geom = self.PAR[pnt]
            ax.scatter( geom.x, geom.y, c='r', s=50 )
            ax.text( geom.x,geom.y, s=pnt, c='r', fontsize=20 )
        C_DATA = f'R = {self.PAR.RADIUS:.3f} m.\n\u03B4 = {self.PAR.DEFLdms}\n'\
                 f'LEN = {self.PAR.LENCUR:.3f} m.\nTangential (T) = {self.PAR.TL:.3f} m\n'\
                 f'Division: {self.PAR.DIV} m.'
        om = LineString([self.PAR.MO,self.PAR.ORIGIN]).centroid
        ax.text( om.x,om.y,s=C_DATA,c='r',fontsize=15, ha='center', va='center' )
        ax.tick_params(axis='x', rotation=90)
        ax.ticklabel_format( useOffset=False, style='plain' )
        plt.gca().set_aspect('equal')
        plt.grid()
        print(f'CircularCurve:DoPlot() Writing result "pdf|png" into ./{self.CACHE}/...')
        if SUFFIX is None: PLT = self.PLOT 
        else: PLT = f'{self.PLOT}_{SUFFIX}'
        plt.savefig( f'{PLT}.png' )
        plt.savefig( f'{PLT}.pdf' )

    def WriteGIS( self, SUFFIX=None ):
        print(f'CircularCurve:WriteGIS() "csv|gpkg" into ./{self.CACHE}/...')
        if SUFFIX is None: PLT = f'{self.PLOT}.gpkg' 
        else: PLT = f'{self.PLOT}_{SUFFIX}.gpkg' 
        self.gdfPNT.to_file( PLT , driver='GPKG', layer='CurveLoci' ) 
        self.dfLS.to_file( PLT, driver='GPKG', layer='Elements' ) 

###############################################################################
class CLI_CircCurve(CircularCurve):
    def __init__(self, ARGS ):
        self.ARGS = ARGS
        #import pdb; pdb.set_trace()
        if type(ARGS) is dict:  # self testing mode
            self.LS_ALIGN =  ARGS['CURVE']
            self.RADIUS, self.DIV = ARGS['RADIUS'], ARGS['DIV']
        else:
            self.LS_ALIGN = LineString(np.array( eval(args.align) ).tolist())
            self.RADIUS, self.DIV = float(ARGS.radius) , float(ARGS.division)
        super().__init__('EPSG:32647', self.LS_ALIGN,self.RADIUS,self.DIV)

###############################################################################
###############################################################################
###############################################################################
if __name__ == "__main__":
    USAGE = '''python3 CurvePnts.py -a [542939.592,1560557.148],[543219.123,1560612.552],[543408.493,1560534.688] -r 300 -d 20'''
    ###############################################################################
    ###############################################################################
    ###############################################################################
    if len(sys.argv)==1:
        xCURVE = LineString( [ [0,0],[0,1000],[-1000,2000] ] ), 1200 # PC-PI-PT, Radius
        xCURVE = LineString( [ [0,0],[0,1000],[1000,2000] ] ), 1200 # PC-PI-PT, Radius
        xCURVE = LineString( [ [0,0],[1000,0],[2000,-1000] ] ), 1200 # PC-PI-PT, Radius
        xCURVE = LineString( [ [0,2000],[0,0],[2000,0] ] ), 1200 # PC-PI-PT, Radius
        xCURVE = LineString( [ [2000,0],[0,0],[0,2000] ] ), 1200 # PC-PI-PT, Radius
        xCURVE = LineString( [ [2000,-300],[0,0],[0,2000] ] ), 1200 # PC-PI-PT, Radius
        xCURVE = LineString( [ [0,0],[1000,0],[1800,-1000] ] ), 500 # PC-PI-PT, Radius
        CURVE = LineString( [ [542939.592,1560557.148],[543219.123,1560612.552],
                              [543408.493,1560534.688] ] ) # PC-PI-PT, Radius
        args = { 'CURVE': CURVE , 'RADIUS':500 , 'DIV':10 }
    else:
        import argparse
        parser = argparse.ArgumentParser(description=PROG, usage=USAGE, 
                                        formatter_class=argparse.RawTextHelpFormatter )
        parser.add_argument( '-a','--align', action='store',             
                    help='3-point of alignment "[E,N],[E,N],[E,N]" for circular curve design' )
        parser.add_argument( '-r','--radius', action='store',type=float,
                    help='design value of the radius in meter' )
        parser.add_argument( '-d','--division', action='store',type=float,
                    help='desired division of the point-on-curve in meter' )
        args = parser.parse_args()
        print(args)
    cc = CLI_CircCurve( args )
    cc.DoPlot( )
    print( cc.gdfPNT )
    cc.WriteGIS()
    print('@@@@@@@@@@@@@@@@@ end of  CurvePnt @@@@@@@@@@@@@@@@@@@@')
    #import pdb; pdb.set_trace()
