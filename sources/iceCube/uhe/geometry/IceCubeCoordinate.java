package iceCube.uhe.geometry;

import geometry.*;
import java.util.*;
import java.io.*;

/**
    This class defines IceCube local coordinate system.
    The y axis is "Grid North"(it is aligned with Prime Meridian,
    pointing toward Greeenwich), the x axis is "Grid East"(pointing
    90 degree clock-wise from Grid North), and the z axis is normal 
    to the Earth's surface, pointing "up".
*/

public class IceCubeCoordinate extends EarthLocalCoordinate {

    /** Nothing of the origin in the horizontal plane. [cm] 52200' East */
    public static double northing = 1.591056e6;

    /** Easting of the origin in the horizontal plane. [cm] 46500' East*/
    public static double easting = 1.41732e6;

    /** Elevation of the origin from earth ROCK surface (sea level). [cm] 
	9284.46ft */
    public static double elevation = 8.8392e4;

    /** The depth of glacier. [cm]*/
    double glacierDepth = 2.829903408e5;
    

    /** origin of the coordinate */
    public static double origin_x = -easting;
    public static double origin_y = northing;
    public static double origin_r = EarthCenterCoordinate.REarth+elevation;
    public static double origin_z = -Math.sqrt(origin_r*origin_r-origin_x*origin_x-origin_y*origin_y);

    
    /** constructor. */ 
    public IceCubeCoordinate(){
	// Point z axis to the earth surface from the IceCube coordinate origin.
	// No aditional rotation for the moment
	super(new J3Vector(origin_x,origin_y,origin_z));
	// Calculate the rotaion angle around the z-axis
	// rotate x-y plane so that the new y axis 
	// has no ex_earthCenter compoment. This operation results in
	// the y axis pararell to the Prime Meridian
	double xx = this.getEx().getX();
	double xy = this.getEx().getY();
	double yx = this.getEy().getX();
	double yy = this.getEy().getY();
        double rotationAngle = -Math.atan(yx/xx);
	if((xy*Math.sin(rotationAngle)+yy*Math.cos(rotationAngle))<0.0) 
	    rotationAngle += Math.PI;
	rotate(3,rotationAngle);
    }

    /** Set the depth of glacier [cm]. */
    public void setGlacierDepth(double depth){
	if(depth > 0.0) this.glacierDepth = depth;
    }


    /** Get the depth of glacier. [cm] */
    public double getGlacierDepth(){
        return glacierDepth;
    }
}
