package geometry;

import java.util.*;

/**
    This class defines local coordinates in/on the earth.
    Its z axis always points to direction from the earth center to the origin.
*/

public class EarthLocalCoordinate extends Coordinate {

    /** Constructor. Rotation angle in x'-y' plane should be given
	to define those axises directions.
	<pre>
	J3Vector origin :  Origin location vector. Should be represented in
                           EarthCenterCoordinate system.
	double alpha    :  the final rotation angle [rad] to determine
	                   the x and y axis directions. 0 correponds to
			   y axis pointing to local south, Math.PI to
			   y axis poing to local north.
	</pre>
    */ 
    public EarthLocalCoordinate(J3Vector origin, double alpha) {
	super(origin);
	J3UnitVector nz = new J3UnitVector(origin.getX(),
					   origin.getY(),origin.getZ());
	transformThisCoordinate(nz,alpha);
    }


    /** Constructor. 
	<pre>
	J3Vector origin :  Origin location vector. Should be represented in
                           EarthCenterCoordinate system.
	</pre>
	Y axis points to local south
    */ 
    public EarthLocalCoordinate(J3Vector origin) {
	this(origin, 0.0);
    }

    /** Constructor. This is a special constructor when the origin is (0, 0, 0).
	<pre>
	J3Vector nz   :  the unit vector to define the z-axis of this local coordinate.
	double alpha  :  the final rotation angle [rad] to determine
	                 the x and y axis directions. 0 correponds to
			 y axis pointing to local south, Math.PI to
			 y axis poing to local north.
	</pre>
    */ 
    public EarthLocalCoordinate(J3UnitVector nz, double alpha) {
	transformThisCoordinate(nz,alpha);
    }

    /** Constructor. Origin location described by longitude, latitude,
	hight from the earth center.
	<pre>
	double height :  height [cm] from the earth center
        double longitude : longitude [deg]
	double latitude :  latitude [deg]
	double alpha    :  the final rotation angle [rad] to determine
	                   the x and y axis directions. 0 correponds to
			   y axis pointing to local south, Math.PI to
			   y axis pointing to local north.
	</pre>
    */
    public EarthLocalCoordinate(double height, double longitude, double latitude,
				double alpha){
	super();

	// build a vector to the origin location. Represented 
	// in the EarthCenterCoordinate.
        EarthCenterCoordinate center = new EarthCenterCoordinate();
        origin = 
	    center.getPointVectorFromLongitudeLatitude(height,
						       longitude,
						       latitude);

	// setting this origin
	setOrigin(origin);

	// transform to the local coordinate planes.
	J3UnitVector nz = new J3UnitVector(origin.getX(),
					   origin.getY(),origin.getZ());
	transformThisCoordinate(nz,alpha);
    }

    /** Transform vector represented in this coordinate to EarthCenterCoordinate.
        Returns the transformed vector.
    */
    public J3Vector transformVectorToEarthCenter(J3Vector r_local){

	J3Vector x_ex = J3Vector.multipleFactor(r_local.getX(),ex);
        J3Vector y_ey = J3Vector.multipleFactor(r_local.getY(),ey);
        J3Vector z_ez = J3Vector.multipleFactor(r_local.getZ(),ez);

        J3Vector r = J3Vector.add((J3Vector.add(x_ex, y_ey)), z_ez);
        J3Vector r_center = J3Vector.add(r, origin);

        return r_center;
    }

    /** Transform unit vector (direction) J3UnitVector n represented 
        in this coordinate system to
        the one in EarthCenterCoordinate. Returns the transformed vector.
    */
    public J3UnitVector transformUnitVectorToEarthCenter(J3UnitVector n_local){

        J3Vector x_ex = J3Vector.multipleFactor(n_local.getX(),ex);
        J3Vector y_ey = J3Vector.multipleFactor(n_local.getY(),ey);
        J3Vector z_ez = J3Vector.multipleFactor(n_local.getZ(),ez);

        J3Vector direction_center = J3Vector.add((J3Vector.add(x_ex, y_ey)), z_ez);
        J3UnitVector n_center = new J3UnitVector(direction_center.getX(),direction_center.getY(),direction_center.getZ());

        return n_center;
    }

}
