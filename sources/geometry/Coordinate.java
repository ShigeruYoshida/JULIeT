package geometry;

import java.util.*;

/**
    This base class defines universal coordinates.
    The right-handed system is setup.
*/

public class Coordinate{

    J3Vector origin;
    J3Vector ex;
    J3Vector ey;
    J3Vector ez;
    
    private static double epsilon = 1.0e-8;

    /** Constructor. Default setting. */
    public Coordinate(){

        ex = new J3UnitVector(1.0, 0.0, 0.0);
        ey = new J3UnitVector(0.0, 1.0, 0.0);
        ez = new J3UnitVector(0.0, 0.0, 1.0);
        origin = new J3Vector(0.0, 0.0, 0.0);
    }

    /** Constructor. Create Cartesian coordinate system and set the origin. */
    public Coordinate(J3Vector origin){

        ex = new J3UnitVector(1.0, 0.0, 0.0);
        ey = new J3UnitVector(0.0, 1.0, 0.0);
        ez = new J3UnitVector(0.0, 0.0, 1.0);

        this.origin = origin;        
    }


    /** Constructor. Create Cartesian coordinate system and set the origin.
	x and z axis with location of the origin is given.
    */
    public Coordinate(J3Vector origin, J3UnitVector ex, J3UnitVector ez ){

        if(Math.abs(J3Vector.getDotProduct(ex,ez)) < epsilon){
            this.ex = ex;
            this.ez = ez;
            ey = J3Vector.getCrossProduct(ez,ex);

            this.origin = origin;
        }else{
            System.err.println("X axis and z axis must be perpendicular.");
            System.exit(0);
        }
    }        

    /** Convert Cartesian to Poler coordinate system.
        Returns r component.
    */
    public double getLengthInPolar(J3Vector a){

        double r;
        double x = a.getX();
        double y = a.getY();
        double z = a.getZ();

        r = Math.sqrt(x*x + y*y + z*z);
        return r;
    }
    
    /** Convert Cartesian to Poler coordinate system.
        Returns zenith angle [rad].
    */
    public double getZenithInPolar(J3Vector a){

        double r;
        double theta;
        double x = a.getX();
        double y = a.getY();
        double z = a.getZ();

        r = Math.sqrt(x*x + y*y + z*z);
        theta = Math.acos(z/r);
        return theta;
    }
    
    /** Convert Cartesian to Poler coordinate system.
        Returns azimuth angle [rad].
    */
    public double getAzimuthInPolar(J3Vector a){

        double alpha;
        double x = a.getX();
        double y = a.getY();

	if(y>0.0){
	    alpha = Math.atan(x/y);
	}else if(y<0.0){
	    alpha = Math.atan(x/y)+ Math.PI;
	}else{
	    if(x>=0.0) alpha = Math.PI/2.0;
	    else alpha = 3.0*Math.PI/2.0;
	}
	if(alpha<0.0) alpha += 2.0*Math.PI;

        return alpha;  
    }

    /** get a point vector for a given Poler system. Input r, and zenith angle
        (theta), azimuth angle (phi) in radian. Returns J3Vector.*/
    public J3Vector getPointVectorFromPolarCoordinate(double r, 
						      double theta, double phi){

        double x = r*Math.sin(theta)*Math.sin(phi);
        double y = r*Math.sin(theta)*Math.cos(phi);
        double z = r*Math.cos(theta);

        J3Vector a = new J3Vector(x,y,z);
        return a;
    }

    /** Rotate the coordinate system.
        It rotete the system theta [rad] around one axis.
        Axis number 1 -> x, 2 -> y, 3 -> z.
	<pre>
	Rotation angle "theta" sign + for clockwise
                                    - for anti-clockwise
	</pre>
    */
    public void rotate(int axis, double theta){

        double x, y, z;
        double sin = Math.sin(theta);
        double cos = Math.cos(theta);
        if(axis == 1){

	    J3Vector eyY = J3Vector.multipleFactor(cos,ey);
	    J3Vector eyZ = J3Vector.multipleFactor(-sin,ez);
	    J3Vector eyRotated = J3Vector.add(eyY,eyZ);

	    J3Vector ezY = J3Vector.multipleFactor(sin,ey);
	    J3Vector ezZ = J3Vector.multipleFactor(cos,ez);
	    J3Vector ezRotated = J3Vector.add(ezY,ezZ);

	    ey = eyRotated;
	    ez = ezRotated;

        }else if(axis == 2){

	    J3Vector ezZ = J3Vector.multipleFactor(cos,ez);
	    J3Vector ezX = J3Vector.multipleFactor(-sin,ex);
	    J3Vector ezRotated = J3Vector.add(ezZ,ezX);

	    J3Vector exZ = J3Vector.multipleFactor(sin,ez);
	    J3Vector exX = J3Vector.multipleFactor(cos,ex);
	    J3Vector exRotated = J3Vector.add(exZ,exX);

	    ez = ezRotated;
	    ex = exRotated;

        }else if(axis == 3){

	    J3Vector exX = J3Vector.multipleFactor(cos,ex);
	    J3Vector exY = J3Vector.multipleFactor(-sin,ey);
	    J3Vector exRotated = J3Vector.add(exX,exY);

	    J3Vector eyX = J3Vector.multipleFactor(sin,ex);
	    J3Vector eyY = J3Vector.multipleFactor(cos,ey);
	    J3Vector eyRotated = J3Vector.add(eyX,eyY);

	    ex = exRotated;
	    ey = eyRotated;

        }else{
            System.err.println("Axis number is 1 -> x, 2 -> y, 3 -> z.");
            System.exit(0);
        } 
    }

    /** transform this coordiante to one with the external z-axis vector
	nz. After the transformation, the new x-y plane is rotated
        by the given rotation angle alpha [rad] to redefine the directions
	of x and y axis.
	<pre>
	J3UnitVector nz    :  new z-axis direction. This vector has to be
                              described by unit base vectors of THIS coordinate
                              before applying this transformation.
	double alpha       :  the final rotation angle [rad] to determine
                              the x and y axis directions.
	</pre>
    */

    public void transformThisCoordinate(J3UnitVector nz, double alpha){

	// get the azimuth angle of nz 
	double phi = getAzimuthInPolar(nz);
	// get the zenith angle of nz
	double theta = getZenithInPolar(nz);

	// rotate around the z-axis
	rotate(3,phi);

	// rotate around the x'-axis
	rotate(1,theta);

	// rotate around the new z-axis by a given angle alpha
	rotate(3,alpha);

    }

    /** Set the origin of the coordinate system. */
    public void setOrigin(J3Vector origin){
        this.origin = origin;
    }

    /** Get x axis direction. */
    public J3Vector getEx(){
        return ex;
    }

    /** Get y axis direction. */
    public J3Vector getEy(){
        return ey;
    }

    /** Get z axis direction. */
    public J3Vector getEz(){
        return ez;
    }

    /** Get the origin of the coordinate system. */
    public J3Vector getOrigin(){
        return origin;
    }



    /** Transform vector point J3Vector r represented 
	in the external coordinate system to
	the one in this coordinate. Returns the transformed vector.
	<pre>
	J3Vector r :    point vector that will be transfomed to this corrdinate.
	Coordinate external: the cordinate system that describes the vector r.
	</pre>
    */ 
    public J3Vector transformVectorToThisCoordinate(J3Vector r, Coordinate external){

	// add the vector from the orgin in THIS coordinate
	J3Vector originToExternalOrigin = J3Vector.subtract(external.getOrigin(),
							    origin);
	J3Vector rFromThisOrigin = J3Vector.add(r,originToExternalOrigin);

	double x = J3Vector.getDotProduct(rFromThisOrigin,ex);
	double y = J3Vector.getDotProduct(rFromThisOrigin,ey);
	double z = J3Vector.getDotProduct(rFromThisOrigin,ez);

	J3Vector a = new J3Vector(x,y,z);

	return a;
    }

    /** Transform unit vector (direction) J3UnitVector n represented 
	in the different coordinate system to
	the one in this coordinate. Returns the transformed vector.
	<pre>
	J3UnitVector n :  direction vector that will be transfomed 
                          to this corrdinate.
	</pre>
    */ 
    public J3UnitVector transformUnitVectorToThisCoordinate(J3UnitVector n){

	double x = J3Vector.getDotProduct(n,ex);
	double y = J3Vector.getDotProduct(n,ey);
	double z = J3Vector.getDotProduct(n,ez);

	J3UnitVector m = new J3UnitVector(x,y,z);

	return m;
    }

}
