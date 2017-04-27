package geometry;

import java.util.*;

/**
<pre>
    This class defines the coordinate with its origin located at
     center of the earth. The z axis points to the North Pole
     and the y axis points to longitude 0 degree latitude 0 degree.
</pre>
*/

public class EarthCenterCoordinate extends Coordinate {

    /** Radius of the earth [cm] */
    public final static double REarth = 6.37814e8; 

    /** Constructor. */ 

    public EarthCenterCoordinate( ) {
        super();
    }


    /** Calculate longitude in degree.
	West : Positive
	East : Negative
    */
    public double getLongitude(J3Vector r){

	double alpha = getAzimuthInPolar(r);

	// corresction for the EAST
	if(alpha>= Math.PI) alpha -= 2.0*Math.PI;

	return Math.toDegrees(alpha);
    }

    /** Calculate latitude in degree.
        North : Positive
	South : Negative
    */

     public double getLatitude(J3Vector r){

	 double theta = getZenithInPolar(r);
	 double lat = 0.5*Math.PI-theta;

	 return Math.toDegrees(lat);
    }



    /** Get the point vector in this earth-center coordinate
	for a given longitude [deg] and latitude [deg].
	<pre>
	double height :  height [cm] from the earth center
        double longitude : longitude [deg]
	double latitude :  latitude [deg]
	</pre>
	Returns J3Vector defined in the EarthCenterCoordinate
    */
    public J3Vector getPointVectorFromLongitudeLatitude(double height, 
						       double longitude, 
						       double latitude){

	// to [rad] from [deg]
        double l = Math.toRadians(longitude);
        double b = Math.toRadians(latitude);

	// get zenith and azimuth angles
	double zenith = 0.5*Math.PI-b;
	double azimuth = l;

	// build a vector
        J3Vector a = getPointVectorFromPolarCoordinate(height,zenith,azimuth);
        return a;
    }

    /** Get the point vector on the earth surface.
	Longitude [deg] and latitude [deg] are given.
    */
    public J3Vector getSurfacePoint(double longitude, double latitude){
        return getPointVectorFromLongitudeLatitude(REarth,longitude,latitude);
    }

    /** Check if the point is inside the earth. If it's inside, returns true. */
    public static boolean IsInsideEarth(J3Vector r){
        if(r.getLength() <= REarth){
            return true;
        }else{
            return false;
        }
    }
         
}    
