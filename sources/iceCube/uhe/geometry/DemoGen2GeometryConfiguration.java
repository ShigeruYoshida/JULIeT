package iceCube.uhe.geometry;

import geometry.*;
import iceCube.uhe.geometry.*;
import iceCube.uhe.points.*;

import java.io.*;
import java.util.*;

/**
    This is a demo program for configuration of
    the propagation geometry defined by IceCubeGen2Coordinate.
    This algorithm is implemented in 
    JulietEventGenerator3Gen2.configurePropagationGeometry().
*/


public class DemoGen2GeometryConfiguration {

    public static void main(String args[]) throws IOException {

        /** Components of point vector for fixed point in IceCubeGen2 coordinate 
            system */
        double x_gen2 = 0.0;
        double y_gen2 = 0.0;
        double z_gen2 = 0.0;
	double nadirAngle_gen2_in = 0.0;
	double azimuthAngle_gen2_in=0.0;


	if(args.length == 5){
            x_gen2 = Double.valueOf(args[0]).doubleValue();
            y_gen2 = Double.valueOf(args[1]).doubleValue();
            z_gen2 = Double.valueOf(args[2]).doubleValue();
            nadirAngle_gen2_in = Double.valueOf(args[3]).doubleValue();
            azimuthAngle_gen2_in = Double.valueOf(args[4]).doubleValue();
        }else{
            System.out.println("Usage: DemoGen2GeometryCoordinate x y z(IceCubeGen2 coord.[cm]) nadir az [deg]");
            System.exit(0);
        }

        /** Defines IceCube coordinate system. */
        IceCubeCoordinate ice3Coordinate = new IceCubeCoordinate();
        IceCubeGen2Coordinate gen2Coordinate = new IceCubeGen2Coordinate();
	EarthCenterCoordinate earthCoordinate =  new EarthCenterCoordinate();

        J3Vector passingPoint_J3Vector_gen2 = new J3Vector(x_gen2,y_gen2,z_gen2);
        J3Vector n_gen2 = gen2Coordinate.getPointVectorFromPolarCoordinate
            (1.0, Math.toRadians(nadirAngle_gen2_in), Math.toRadians(azimuthAngle_gen2_in));
        J3UnitVector direction_J3UnitVector_gen2 =
	    new J3UnitVector(n_gen2.getX(),n_gen2.getY(),n_gen2.getZ());
	J3Line particleAxis_J3Line_gen2 = new J3Line(passingPoint_J3Vector_gen2,
                                              direction_J3UnitVector_gen2);

	System.out.format("passing point in gen2 coordinate(%f %f %f)\n",
	    particleAxis_J3Line_gen2.getR0().getX(),
	    particleAxis_J3Line_gen2.getR0().getY(),
	    particleAxis_J3Line_gen2.getR0().getZ());
	System.out.format("direction in gen2 coordinate (%f %f %f)\n",
	    particleAxis_J3Line_gen2.getDirection().getX(),
	    particleAxis_J3Line_gen2.getDirection().getY(),
   	    particleAxis_J3Line_gen2.getDirection().getZ());

	J3Vector passingPoint_J3Vector_center =
	    gen2Coordinate.transformVectorToEarthCenter(particleAxis_J3Line_gen2.getR0());
	J3UnitVector direction_J3UnitVector_center =
            gen2Coordinate.transformUnitVectorToEarthCenter(particleAxis_J3Line_gen2.getDirection());
	J3Line particleAxis_J3Line_center = new J3Line(passingPoint_J3Vector_center,
						       direction_J3UnitVector_center);

	System.out.format("passing point in earth center coordinate(%e %e %e)\n",
	    particleAxis_J3Line_center.getR0().getX(),
	    particleAxis_J3Line_center.getR0().getY(),
	    particleAxis_J3Line_center.getR0().getZ());
	System.out.format("direction in earth center coordinate (%f %f %f)\n",
	    particleAxis_J3Line_center.getDirection().getX(),
	    particleAxis_J3Line_center.getDirection().getY(),
   	    particleAxis_J3Line_center.getDirection().getZ());

        ParticleTracker.setInitialPoint(particleAxis_J3Line_center,ice3Coordinate);
        J3Vector particleEntrance = new J3Vector(
                                                 particleAxis_J3Line_center.getX(),
                                                 particleAxis_J3Line_center.getY(),
                                                 particleAxis_J3Line_center.getZ());
        particleAxis_J3Line_center.setR0(particleEntrance);
        particleAxis_J3Line_center.setAxisLength(0.0); 

	System.out.format("entrance point in earth center coordinate(%e %e %e)\n",
	    particleEntrance.getX(),
	    particleEntrance.getY(),
	    particleEntrance.getZ());

	double r_earth = ParticlePoint.REarth + gen2Coordinate.getGlacierDepth();
	double ratio = particleEntrance.getLength()/r_earth;

	System.out.format("Distance from earth center to the entrance point/R-earth+glacir depth=%e\n",
			  ratio);

        double nadirAngleAtEntrance =
            Math.PI - J3Vector.getAngleInRadian(particleAxis_J3Line_center,
                                                direction_J3UnitVector_center);

	double nadirAngleAtCenter =
   	    Math.acos(particleAxis_J3Line_center.getDirection().getZ());

	double cosNadir =  direction_J3UnitVector_gen2.getZ();
	double nadirAngle_gen2 = Math.acos(cosNadir);

	System.err.format("nadir angle at the entrance %f [deg]\n",
			  Math.toDegrees(nadirAngleAtEntrance));
	System.err.format("nadir angle in the earth center coordinate %f [deg]\n",
			  Math.toDegrees(nadirAngleAtCenter));
	System.err.format("nadir angle in gen2 coordinate %f [deg]\n",
			  Math.toDegrees(nadirAngle_gen2));


        J3Vector particleEntrance_gen2 =
            gen2Coordinate.transformVectorToThisCoordinate(
							   particleEntrance,earthCoordinate);
        particleAxis_J3Line_gen2.setR0(particleEntrance_gen2);
	particleAxis_J3Line_gen2.setAxisLength(0.0);

	// now search for the starting point
	double trackLength = 0.0; 
	double trackLengthMax = 2.0*(ParticlePoint.REarth + gen2Coordinate.getGlacierDepth());
	double delta_length = 1.0e2; // 1 [m] step
	boolean foundTheInjectionPoint = false;
	while(trackLength < trackLengthMax){
	    trackLength += delta_length;
	    particleAxis_J3Line_gen2.setAxisLength(trackLength);
	    J3UnitVector n_primary = new J3UnitVector(particleAxis_J3Line_gen2.getX(),
						      particleAxis_J3Line_gen2.getY(),
						      particleAxis_J3Line_gen2.getZ());
	    double angle2primary  = Math.acos(n_primary.getZ());
	    double distance = InjectionGeometryUtils.getDistanceOfStartLocation(angle2primary);
	    double distance2primary = particleAxis_J3Line_gen2.getLength();
	    if(distance2primary<=distance){
		System.out.format("distance from the gen2 center=%e %e angle(%f)\n",
				  distance,distance2primary,Math.toDegrees(angle2primary));
		foundTheInjectionPoint = true;
		break;
	    }
	}
	if(!foundTheInjectionPoint){// outside earth+glacier!!
	    particleAxis_J3Line_gen2.setAxisLength(0.0);
	    // back to the earth entrance 
	}

        double startLocation = particleAxis_J3Line_gen2.getAxisLength();
	System.out.format("start location from the earth surface %e\n",startLocation);

	J3Vector startLocation_J3Vector_gen2 = new J3Vector(
						     particleAxis_J3Line_gen2.getX(),
                                                     particleAxis_J3Line_gen2.getY(),
                                                     particleAxis_J3Line_gen2.getZ());

	System.out.format("injection point in gen2 coordinate(%e %e %e)\n",
	    startLocation_J3Vector_gen2.getX(),
	    startLocation_J3Vector_gen2.getY(),
	    startLocation_J3Vector_gen2.getZ());

	particleAxis_J3Line_center.setAxisLength(startLocation);

	System.out.format("injection point in the earth center coordinate(%e %e %e)\n",
			  particleAxis_J3Line_center.getX(),
			  particleAxis_J3Line_center.getY(),
			  particleAxis_J3Line_center.getZ());

	J3Vector startLocation_J3Vector_ice3 =
	    ice3Coordinate.transformVectorToThisCoordinate(
							   particleAxis_J3Line_center,earthCoordinate);

	System.out.format("injection point in ice3 coordinate(%e %e %e)\n",
	    startLocation_J3Vector_ice3.getX(),
	    startLocation_J3Vector_ice3.getY(),
	    startLocation_J3Vector_ice3.getZ());

        J3Vector particleEntrance_ice3 =
            ice3Coordinate.transformVectorToThisCoordinate(
							   particleEntrance,earthCoordinate);
        J3UnitVector direction_J3UnitVector_ice3 =
            ice3Coordinate.transformUnitVectorToThisCoordinate(direction_J3UnitVector_center);
        J3Line particleAxis_J3Line_ice3 = new J3Line(particleEntrance_ice3,
                                              direction_J3UnitVector_ice3);
	particleAxis_J3Line_ice3.setAxisLength(startLocation);


	System.out.format("injection point in ice3 coordinate by J3Line_axis(%e %e %e)\n",
			  particleAxis_J3Line_ice3.getX(),
			  particleAxis_J3Line_ice3.getY(),
			  particleAxis_J3Line_ice3.getZ());

    }
}
