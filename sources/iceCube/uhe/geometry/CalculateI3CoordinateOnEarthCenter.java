package iceCube.uhe.geometry;

import geometry.*;

import java.io.*;

/**
<pre>
    This is a demo program for calculating x, y, and z coordinate
    defined by the IceCobeCoordinate on the EarthCenterCoordinate
</pre>
*/


public class CalculateI3CoordinateOnEarthCenter {

    public static void main(String args[]) throws IOException {

        /** Components of point vector for fixed point in IceCube coordinate 
            system */
        double x_ice3 = 0.0;
        double y_ice3 = 0.0;
        double z_ice3 = 0.0;


	if(args.length == 3){
            x_ice3 = Double.valueOf(args[0]).doubleValue();
            y_ice3 = Double.valueOf(args[1]).doubleValue();
            z_ice3 = Double.valueOf(args[2]).doubleValue();
        }else{
            System.out.println("Usage: CalculateI3CoordinateOnEarthCenter x y z(IceCube coord.[cm])");
            System.exit(0);
        }

        /** Defines IceCube coordinate system. */
        IceCubeCoordinate iceCube = new IceCubeCoordinate();

        J3Vector origin = iceCube.getOrigin();
        J3Vector ex = iceCube.getEx();
        J3Vector ey = iceCube.getEy();
        J3Vector ez = iceCube.getEz(); 
        System.err.println("IceCube origin:(" + 
                           origin.getX() + ", " + 
                           origin.getY() + ", " + 
                           origin.getZ()+") ");
        System.err.println("IceCube ex:("+ex.getX()+", "+ex.getY()+", "+ex.getZ()+") ");
        System.err.println("IceCube ey:("+ey.getX()+", "+ey.getY()+", "+ey.getZ()+") ");
        System.err.println("IceCube ez:("+ez.getX()+", "+ez.getY()+", "+ez.getZ()+") ");
  

        // Calculation of the fixed point
        J3Vector fixedPoint_ice3 = new J3Vector(x_ice3,y_ice3,z_ice3);
        System.err.println("Fixed point on Ice3:(" + fixedPoint_ice3.getX()+ ", " + 
                           fixedPoint_ice3.getY()+", "+fixedPoint_ice3.getZ()+") ");

        /** Defines EarthCenterCoordinate system. */ 
        EarthCenterCoordinate c = new EarthCenterCoordinate();

        /** Transform to EarthCenterCoordinate of the fixed point and direction */
        J3Vector fixedPoint_center = iceCube.transformVectorToEarthCenter(fixedPoint_ice3);
        System.err.println("Fixed point on Earth center:(" + 
                           fixedPoint_center.getX() + ", " + 
                           fixedPoint_center.getY() + ", " + 
                           fixedPoint_center.getZ() + ")");

        /** Defines IceCubeGen2 coordinate system. */
        IceCubeGen2Coordinate gen2 = new IceCubeGen2Coordinate();

        /** Transform to IceCubeGen2 coordinate of the fixed point and direction */
	J3Vector fixedPoint_gen2 = gen2.transformVectorToThisCoordinate(fixedPoint_center,c);
        System.err.println("Fixed point on gen2:(" + 
                           fixedPoint_gen2.getX() + ", " + 
                           fixedPoint_gen2.getY() + ", " + 
                           fixedPoint_gen2.getZ() + ")");

        /** Transform BACK to IceCube coordinate of the fixed point and direction */
	J3Vector fixedPoint_ice3_back = iceCube.transformVectorToThisCoordinate(fixedPoint_center,c);
        System.err.println("Fixed point back on Ice3:(" + 
                           fixedPoint_ice3_back.getX() + ", " + 
                           fixedPoint_ice3_back.getY() + ", " + 
                           fixedPoint_ice3_back.getZ() + ")");

    }
}
