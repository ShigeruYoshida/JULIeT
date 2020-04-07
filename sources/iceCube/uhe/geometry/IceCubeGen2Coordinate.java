package iceCube.uhe.geometry;

import geometry.*;
import java.util.*;
import java.io.*;

/**
    This class defines IceCube-Gen2 local coordinate system.
    The y axis is "Grid North"(it is aligned with Prime Meridian,
    pointing toward Greeenwich), the x axis is "Grid East"(pointing
    90 degree clock-wise from Grid North), and the z axis is normal 
    to the Earth's surface, pointing "up".
*/

public class IceCubeGen2Coordinate extends EarthLocalCoordinate {

    /** Nothing of the origin in the horizontal plane. [cm] 52200' East */
    public static double gen2northing = IceCubeCoordinate.northing + 1.84521e4;

    /** Easting of the origin in the horizontal plane. [cm] 46500' East*/
    public static double gen2easting = IceCubeCoordinate.easting -4.51064e4;

    /** Elevation of the origin from earth ROCK surface (sea level). [cm] 
	9284.46ft */
    public static double elevation = IceCubeCoordinate.elevation;

    /** Nothing of the origin in the horizontal plane. [cm] 52200' East */
    public static double ara_northing = IceCubeCoordinate.northing - 9.23e4;

    /** Easting of the origin in the horizontal plane. [cm] 46500' East*/
    public static double ara_easting = IceCubeCoordinate.easting -3.345e5;

    /** Elevation of the origin from earth ROCK surface (sea level). [cm] 
	9284.46ft */
    public static double ara_elevation = 1.4e5; 

    /** The depth of glacier. [cm]*/
    double glacierDepth = 2.829903408e5;


    /** origin of the coordinate */
    public static double origin_x = -gen2easting;
    public static double origin_y = gen2northing;
    public static double origin_r = EarthCenterCoordinate.REarth+elevation;
    public static double origin_z = -Math.sqrt(origin_r*origin_r-origin_x*origin_x-origin_y*origin_y);

    /** constructor. */ 
    public IceCubeGen2Coordinate(){
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

    /** Static method to swich the default parameters to those for ARA.
	You have to call the constructor after this method, 
	if you simulate the large volume including ARA 
    */
    public static void setARADimensionToDefaults(){
	origin_x = -ara_easting;
	origin_y = ara_northing;
	origin_r = EarthCenterCoordinate.REarth+ara_elevation;
	origin_z = -Math.sqrt(origin_r*origin_r-origin_x*origin_x-origin_y*origin_y);
    }

    /** Static method to swich the default parameters to those for gen2.
	You have to call the constructor after this method, 
	if you simulate the large volume including Gen2 
    */
    public static void setGen2DimensionToDefaults(){
	origin_x = -gen2easting;
	origin_y = gen2northing;
	origin_r = EarthCenterCoordinate.REarth+elevation;
	origin_z = -Math.sqrt(origin_r*origin_r-origin_x*origin_x-origin_y*origin_y);
    }

    /** Simple main method */

    public static void main(String args[]) {
	IceCubeGen2Coordinate gen2coordinate = new IceCubeGen2Coordinate();
	double x = gen2coordinate.origin_x;
	double y = gen2coordinate.origin_y;
	double r = gen2coordinate.origin_r;
	System.out.format("Gen2 coordinate origin x(%e) y(%e) r(%e)\n",x,y,r);
	IceCubeGen2Coordinate.setARADimensionToDefaults();
	x = gen2coordinate.origin_x;
	y = gen2coordinate.origin_y;
	r = gen2coordinate.origin_r;
	System.out.format("ARA coordinate origin x(%e) y(%e) r(%e)\n",x,y,r);
	IceCubeGen2Coordinate.setGen2DimensionToDefaults();
	x = gen2coordinate.origin_x;
	y = gen2coordinate.origin_y;
	r = gen2coordinate.origin_r;
	System.out.format("Gen2 coordinate origin x(%e) y(%e) r(%e)\n",x,y,r);
    }
 
}
