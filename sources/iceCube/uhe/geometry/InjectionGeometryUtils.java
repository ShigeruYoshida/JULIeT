package iceCube.uhe.geometry;

import geometry.*;

/**

   Collection of the static methods on calculation of the geometrical parameters
   for in-ice particle injection in JULIeT. This class is introduced
   for IceCube-Gen2 whose detector volume is no longer cubic.

   For (non-gen2) IceCube, the injection position can be simply calculated
   by the sphare with radis of  IceCubeCoordinate.elevation. This simple approach
   is not valid any more for Gen2. Insetad this class introduces the hypothetical
   of Cylinder with radius R and height Z

 */
public class InjectionGeometryUtils {

    /** default radis of the cylinder radius */
    public static double R_cylinder = 2.5e5; // 2,500 [m] = 2.5x10^5 [cm]. Default
    /** default radis of the cylinder height */
    public static double z_cylinder = 2.0*IceCubeGen2Coordinate.elevation;

    /**
       Calculate the radius of circle on which the juliet in-ice particles are injected
       in the IceCube gen2 simulation.
       <pre>
       double R  : radis of the hypothetical cylinder [cm]. Should be order of R_cylinder, 
                   the static member of this class.
       double z  : Height of the cylinder [cm]. 
       double theta : Zenith (or nadir) angle [Rad].
       </pre>
     */
    public static double getInjectionRadius(double R, double z, double theta){

	double cosZenith = Math.abs(Math.cos(theta));
	double sinZenith = Math.abs(Math.sin(theta));
	double r_injection = R*cosZenith + (z/2.0)*sinZenith;
	return r_injection;
    }

    /**
       Calculate the radius of circle on which the juliet in-ice particles are injected
       in the IceCube gen2 simulation. It considers the hypothetical cylinder
       with the default dimensions, R_cylinder and z_cylinder.
    */
    public static double getInjectionRadius(double theta){
	return getInjectionRadius(R_cylinder,z_cylinder,theta);
    }

    /**
       Calculate a rectangle in-ice injection area used in the JULIeT IceCube simulation.
       It is 4 * R_cylinder(R,z, 0.0 [rad]) * getInjectionRadius(R,z,theta) [cm2].
     */
    public static double getInIceRectangleInjectionArea(double R, double z, double theta){
	return 4.0*getInjectionRadius(R, z, theta)*getInjectionRadius(R, z, 0.0);
    }

    /**
       Calculate a rectangle in-ice injection area used in the JULIeT IceCube simulation.
       The defaul dimension is assumed.
       It is 4 * R_cylinder * getInjectionRadius(R_cylinder,z_cylinder,theta) [cm2].
     */
    public static double getInIceRectangleInjectionArea(double theta){
	return 4.0*getDefaultCylinderRadius()*getInjectionRadius(theta);
    }


    /**
       Calculate the distance of the injection location from the IceCube-gen2 center.
       An in-ince JULIeT particles starts its propagation from this location.
       <pre>
       double R  : radis of the hypothetical cylinder [cm]. Should be order of R_cylinder, 
                   the static member of this class.
       double z  : Height of the cylinder [cm]. 
       double theta : Zenith (or nadir) angle [Rad].
       </pre>
     */
    public static double getDistanceOfStartLocation(double R, double z, double theta){

	double cosZenith = Math.abs(Math.cos(theta));
	double tanZenith = Double.POSITIVE_INFINITY;
	if(cosZenith!=0.0){
	    tanZenith = Math.abs(Math.tan(theta));
	}

	double a = (2.0*R)/z;
	double d = 1.0;
	if(tanZenith>=a) {
	    double sinZenith = Math.abs(Math.sin(theta));
	    d = R/sinZenith;
	}else{
	    d = z/(2.0*cosZenith);
	}

	return d;
    }

    /**
       Calculate the distance of the injection location from the IceCube-gen2 center.
       An in-ince JULIeT particles starts its propagation from this location.
       It considers the hypothetical cylinder
       with the default dimensions, R_cylinder and z_cylinder.
    */
    public static double getDistanceOfStartLocation(double theta){
	return getDistanceOfStartLocation(R_cylinder, z_cylinder, theta);
    }

    public static double getDefaultCylinderRadius(){return R_cylinder;}
    public static double getDefaultCylinderHeight(){return z_cylinder;}


    /** Simple main method */

    public static void main(String args[]) {

	double theta = 0.0; // nadir angle
	double cm2m = 1.0e-2;
	while(theta<=180.0){
	    double l = InjectionGeometryUtils.getDistanceOfStartLocation(Math.toRadians(theta));
	    double R = InjectionGeometryUtils.getInjectionRadius(Math.toRadians(theta));
	    System.out.format("theta(%4.1f [deg]) R=%6.1f[m] l=%6.1f[m]\n",theta,R*cm2m,l*cm2m);
	    theta += 2.0; 
	}
    }

}
