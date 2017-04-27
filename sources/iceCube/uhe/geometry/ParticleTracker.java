package iceCube.uhe.geometry;

import iceCube.uhe.points.*;
import geometry.*;

import java.io.*;

/** 
    This class contains static methods concerning the particle tracking.
*/

public class ParticleTracker {

    public static void setInitialPoint(J3Line axis, IceCubeCoordinate iceCube){
        double REarth_ice = ParticlePoint.REarth + iceCube.getGlacierDepth();
        J3Utility.setJ3LineNegativeAxisLengthForGivenLength(axis, REarth_ice);
    }

    /** Check if the particle location represented by EarthCenter coordinate
	is inside the earth. 
        This 'earth' contains glacier of the the Antarctica.
	<pre>
	J3Vector r_center  :  Prarticle location defined by EarthCenterCoordinate
	iceCubeCoordinate iceCube  : The IceCube local coordinate
	EarthCenterCoordinate center : The Earth center coordinate
	Volume outVol : The outside volume - If the partile is outside this volume
	                AND the flag is 1, the particle is considered to be "outside"
                        the earth and may stop its tracking.
	int flag     :  control flag - when 1, then the volume outVol is involved.
	</pre>
    */
    public static boolean isInsideEarth(J3Vector r_center,
                                        IceCubeCoordinate iceCube,
                                        EarthCenterCoordinate center,
					Volume outVol, int flag){
        double REarth = ParticlePoint.REarth + iceCube.getGlacierDepth();
        J3Vector r_ice3 = iceCube.transformVectorToThisCoordinate(r_center,center);
        if(r_center.getLength() > REarth || ((flag ==1) && !outVol.isInsideVolume(r_ice3))){
            return false;
        }else{
            return true;
        }
    }
    public static boolean isInsideEarth(J3Vector r_center, J3Vector shift,
                                        IceCubeCoordinate iceCube,
                                        EarthCenterCoordinate center,
					Volume outVol, int flag){
        double REarth = ParticlePoint.REarth + iceCube.getGlacierDepth();
        J3Vector r_ice3 = iceCube.transformVectorToThisCoordinate(r_center,center);
        if(r_center.getLength() > REarth || ((flag ==1) && !outVol.isInsideVolume(r_ice3,shift))){
            return false;
        }else{
            return true;
        }
    }
}
