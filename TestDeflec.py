#
#
#
import os
from CurvePnts import *
from shapely.geometry import LineString,Point
from shapely.affinity import rotate
from shapely.ops import linemerge


ls = LineString( [[542939.592,560557.148],[543219.123,1560612.552],[543408.493,1560534.688]] )
for def_ang in range(10, 360, 10 ):
    if def_ang == 180: continue
    LEN = 5000 # 
    PC = Point( [500_000      ,1_500_000] )
    PI = Point( [500_000+  LEN,1_500_000] )
    PT = Point( [500_000+2*LEN,1_500_000] )
    L1 = LineString( [PC,PI] )
    L2 = LineString( [PI,PT] )
    L2r = rotate( L2, def_ang, origin=PI, use_radians=False )
    ls = LineString( [ PC,PI, L2r.coords[1] ] )
    cc = CircularCurve( 32647, ls , 300, 20 )
    cc.DoPlot(PLOT_FILE=f'{def_ang:03d}' )
    #import pdb ; pdb.set_trace()
