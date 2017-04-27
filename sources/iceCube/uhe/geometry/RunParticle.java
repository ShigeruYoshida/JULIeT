package iceCube.uhe.geometry;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import geometry.*;

import java.io.*;

/**
<pre>
    This is a demo program for running a particle in the earth.
    Mainly for debugging use.
    It runs particle and output the coordinate of particle's location.
</pre>
*/


public class RunParticle {

    public static void main(String args[]) throws IOException {

        /** Nadir angle [deg] and azimuth [deg] angle of particle in IceCube
            coordinate system */      
        double nadir = 0.0;        
        double azimuth = 0.0;

        /** Components of point vector for fixed point in IceCube coordinate 
            system */
        double x_ice3 = 0.0;
        double y_ice3 = 0.0;
        double z_ice3 = 0.0;

        if(args.length == 5){
	    nadir = Double.valueOf(args[0]).doubleValue();
            azimuth = Double.valueOf(args[1]).doubleValue();
            x_ice3 = Double.valueOf(args[2]).doubleValue();
            y_ice3 = Double.valueOf(args[3]).doubleValue();
            z_ice3 = Double.valueOf(args[4]).doubleValue();
        }else{
            System.out.println("Usage: RunParticle nadir[deg] azimuth[deg] x y z(IceCube coord.[cm])");
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

        /** Calculation of the direction */
        J3Vector direction_ice3 = iceCube.getPointVectorFromPolarCoordinate
	    (1.0, Math.toRadians(nadir), Math.toRadians(azimuth));
        J3UnitVector n_ice3 = new J3UnitVector(direction_ice3.getX(),
					       direction_ice3.getY(),
					       direction_ice3.getZ()); 
        System.err.println("Direction on Ice3:(" + 
			   n_ice3.getX() + ", " + 
			   n_ice3.getY() + ",  "+ 
			   n_ice3.getZ() + ") ");

        /** Defines EarthCenterCoordinate system. */ 
        EarthCenterCoordinate c = new EarthCenterCoordinate();

        J3Vector ex_center = c.getEx();
        J3Vector ey_center = c.getEy();
        J3Vector ez_center = c.getEz();
        System.err.println("Center ex:(" + 
			   ex_center.getX() + ", " + 
			   ex_center.getY() + ", " + 
			   ex_center.getZ() + ") ");
        System.err.println("Center ey:(" + 
			   ey_center.getX() + ", " + 
			   ey_center.getY() + ", " + 
			   ey_center.getZ() + ") ");
        System.err.println("Center ez:(" + 
			   ez_center.getX() + ", " + 
			   ez_center.getY() + ", " + 
			   ez_center.getZ() + ") ");

        /** Transform to EarthCenterCoordinate of the fixed point and direction */
        J3Vector fixedPoint_center = iceCube.transformVectorToEarthCenter(fixedPoint_ice3);
        System.err.println("Fixed point on center:(" + 
			   fixedPoint_center.getX() + ", " + 
			   fixedPoint_center.getY() + ", " + 
			   fixedPoint_center.getZ() + ")");

	J3UnitVector n_center = iceCube.transformUnitVectorToEarthCenter(n_ice3);
        System.err.println("Direction on center:(" + 
			   n_center.getX() + ", " + 
			   n_center.getY() + ", " + 
			   n_center.getZ() + ") ");
 
        fixedPoint_ice3 = iceCube.transformVectorToThisCoordinate(fixedPoint_center,c);
        System.err.println("Fixed point on Ice3 (again) :(" + 
			   fixedPoint_ice3.getX() + ", " + 
			   fixedPoint_ice3.getY() + ", " + 
			   fixedPoint_ice3.getZ() + ") ");

        n_ice3 = iceCube.transformUnitVectorToThisCoordinate(n_center);
        System.err.println("Direction on Ice3 (again) :(" + 
			   n_ice3.getX() + ", " + 
			   n_ice3.getY() + ", " +
			   n_ice3.getZ() + ") ");

        /** Generate J3Line object. */
        J3Line axis = new J3Line(fixedPoint_center,n_center);
 
        /** Set incident point.               
            The axisLength of J3Line is set to 0 at the incident point.*/        
        ParticleTracker.setInitialPoint(axis,iceCube);
        J3Vector incidentPoint = new J3Vector(axis.getX(),axis.getY(),axis.getZ());
        axis.setR0(incidentPoint);
        axis.setAxisLength(0.0);

        double nadir_center = Math.PI - J3Vector.getAngleInRadian(axis, n_center);

        System.err.println("Nadir angle: " + Math.toDegrees(nadir_center) + " [deg]");
        System.err.println("Start point on center:(" +
			   axis.getX() + ", " + 
			   axis.getY() + ", " + 
			   axis.getZ() + ") ");
	System.err.println("Distance from the earth center " +
			   axis.getLength() + " [cm]");
       

        /** Point vector of the particle */
        J3Vector r_center = axis;
        J3Vector r_ice3 = iceCube.transformVectorToThisCoordinate(r_center, c);
        System.err.println("Start point on ice3:(" + 
			   r_ice3.getX() + ", " + 
			   r_ice3.getY() + ", " +
			   r_ice3.getZ() +") ");
        Particle mu = new Particle(1,1);

        /** Generate ParticlePoint object. */
        ParticlePoint point = new ParticlePoint(0.0,nadir_center,1);
        System.err.println("Trajectory length: " + point.getAxisLength() + " [cm]");
        double pathLength = 1.0e4;

        IceCubeVolume iceVol = new IceCubeVolume();
        Volume outVol = new Volume(2.0e5);

        /** Set start point for tracking.  */
        double startDistance = 2.0e5;    //Distance from IceCube origin [cm]
        J3Line ice3axis = new J3Line(fixedPoint_ice3,n_ice3);
        J3Utility.setJ3LineNegativeAxisLengthForGivenLength(ice3axis, startDistance);
        J3Vector startPoint = iceCube.transformVectorToEarthCenter(ice3axis);
        double startLength = J3Vector.subtract(startPoint,incidentPoint).getLength();
        axis.setAxisLength(startLength);
        J3Vector axis_ice3 =
	    iceCube.transformVectorToThisCoordinate(axis, c);
	System.err.println("startLength = " + startLength + " [cm]");
	System.err.println("distance from IceCubeCenter " +
			   axis_ice3.getLength() +  " [cm]");


        System.out.println("zone 2 2");

        //X-Y
        System.out.println("titx X");
        System.out.println("tity Y");
        System.out.println("gwin 0.2 0.8 0.2 0.8");
        System.out.println("scal -1.0e5 1.0e5 -1.0e5 1.0e5");
        System.out.println("mksz 0.1");

        int i = 0;
        int flag = 0;
        double l = startLength;

        while(true) {

            l += pathLength;
            point.setParticleLocation(l);
            axis.setAxisLength(l);
  
	    if(!ParticleTracker.isInsideEarth(axis,iceCube,c,outVol,flag))    break; 
	    //Your particle has now emerged from underground or far away from detector.
	    //	     System.out.println("data "+axis.getX()+" 0.0 "+axis.getY()+" 0.0");
	    //       System.out.println("Density: " +point.getMediumDensity());

	    axis_ice3 = iceCube.transformVectorToThisCoordinate(axis,c);	    
            if(iceVol.isInsideVolume(axis_ice3)){
                i++;
                if(i == 1){
                    System.out.println("mkcl 5"); 
                }else{
                    System.out.println("mkcl 4");
                }
   
                System.out.println("data "+axis_ice3.getX()+" 0.0 "+axis_ice3.getY()+" 0.0");
                System.out.println("mkcl 0");
            }else if(outVol.isInsideVolume(axis_ice3)){
                flag = 1;
                System.out.println("data "+axis_ice3.getX()+" 0.0 "+axis_ice3.getY()+" 0.0");
            }
	    
	}
        System.out.println("plot");
        System.out.println("disp");
        System.out.println("cont");

	System.out.println("line -5.0e4 -5.0e4 -5.0e4 5.0e4");
        System.out.println("line -5.0e4 -5.0e4 5.0e4 -5.0e4");
	System.out.println("line 5.0e4 5.0e4 -5.0e4 5.0e4");
	System.out.println("line 5.0e4 5.0e4 5.0e4 -5.0e4");

        System.out.println("plot");
        System.out.println("disp");
        System.out.println("endg");

        //Y-Z
        System.out.println("titx Y");
        System.out.println("tity Z");
        System.out.println("gwin 0.2 0.8 0.2 0.8");
        System.out.println("scal -1.0e5 1.0e5 -1.0e5 1.0e5");
        System.out.println("mksz 0.1");

        i = 0;
        flag = 0;
        l = startLength;
        while(true) {

            l += pathLength;
            point.setParticleLocation(l);
            axis.setAxisLength(l);
            if(!ParticleTracker.isInsideEarth(axis,iceCube,c,outVol,flag))    break; 

	    //	     System.out.println("data "+axis.getY()+" 0.0 "+axis.getZ()+" 0.0");
	    //       System.out.println("Density: " +point.getMediumDensity());

	    axis_ice3 = iceCube.transformVectorToThisCoordinate(axis,c);	    
            if(iceVol.isInsideVolume(axis_ice3)){
                i++;
                if(i == 1){
                    System.out.println("mkcl 5"); 
                }else{
                    System.out.println("mkcl 4");
                }
   
                System.out.println("data "+axis_ice3.getY()+" 0.0 "+axis_ice3.getZ()+" 0.0");
                System.out.println("mkcl 0");
            }else if(outVol.isInsideVolume(axis_ice3)){
                flag = 1;
                System.out.println("data "+axis_ice3.getY()+" 0.0 "+axis_ice3.getZ()+" 0.0");
            }
	    
	}
        System.out.println("plot");
        System.out.println("disp");
        System.out.println("cont");

	System.out.println("line -5.0e4 -5.0e4 -5.0e4 5.0e4");
        System.out.println("line -5.0e4 -5.0e4 5.0e4 -5.0e4");
	System.out.println("line 5.0e4 5.0e4 -5.0e4 5.0e4");
	System.out.println("line 5.0e4 5.0e4 5.0e4 -5.0e4");

        System.out.println("plot");
        System.out.println("disp");
        System.out.println("endg");

        //X-Z
        System.out.println("titx X");
        System.out.println("tity Z");
        System.out.println("gwin 0.2 0.8 0.2 0.8");
        System.out.println("scal -1.0e5 1.0e5 -1.0e5 1.0e5");
        System.out.println("mksz 0.1");

        i = 0;
        flag = 0;
        l = startLength;
        while(true) {

            l += pathLength;
            point.setParticleLocation(l);
            axis.setAxisLength(l);
            if(!ParticleTracker.isInsideEarth(axis,iceCube,c,outVol,flag))    break;
         

	    //	     System.out.println("data "+axis.getX()+" 0.0 "+axis.getZ()+" 0.0");
	    //       System.out.println("Density: " +point.getMediumDensity());

            axis_ice3 = iceCube.transformVectorToThisCoordinate(axis,c);	    
            if(iceVol.isInsideVolume(axis_ice3)){
                i++;
                if(i == 1){
                    System.out.println("mkcl 5"); 
                }else{
                    System.out.println("mkcl 4");
                }
   
                System.out.println("data "+axis_ice3.getX()+" 0.0 "+axis_ice3.getZ()+" 0.0");
                System.out.println("mkcl 0");
            }else if(outVol.isInsideVolume(axis_ice3)){
                flag = 1;
                System.out.println("data "+axis_ice3.getX()+" 0.0 "+axis_ice3.getZ()+" 0.0");
            }
	    
	}
        System.out.println("plot");
        System.out.println("disp");
        System.out.println("cont");

	System.out.println("line -5.0e4 -5.0e4 -5.0e4 5.0e4");
        System.out.println("line -5.0e4 -5.0e4 5.0e4 -5.0e4");
	System.out.println("line 5.0e4 5.0e4 -5.0e4 5.0e4");
	System.out.println("line 5.0e4 5.0e4 5.0e4 -5.0e4");

        System.out.println("plot");
        System.out.println("disp");
        System.out.println("endg");


        //        System.out.println("data 0.0 "+ParticlePoint.REarth+" 0.0 0.0");
	//        System.out.println("mapc");
	//        System.out.println("disp");
	//        System.out.println("cont");

        //        System.out.println("mkcl 3");
	//        System.out.println("mksz 0.1");
	//        System.out.println("data 0.0 0.0 0.0 0.0");
        //        System.out.println("data "+(depth-ParticlePoint.REarth)+" 0.0 0.0 0.0");

        //Density
        System.out.println("titx l");
        System.out.println("tity Density");
        System.out.println("gwin 0.2 0.8 0.2 0.8");
        System.out.println("scal 0.0 "+point.getAxisLength()+" 0.0 15.0");
        System.out.println("mksz 0.1");

        flag = 0;
        l = startLength;
        pathLength = 1.0e5;

        while(true) {

            l += pathLength;
            point.setParticleLocation(l);
            axis.setAxisLength(l);
  
	    if(!ParticleTracker.isInsideEarth(axis,iceCube,c,outVol,flag))    break; 
	    //Your particle has now emerged from underground or far away from detector.

            axis_ice3 = iceCube.transformVectorToThisCoordinate(axis,c);   
	    if(outVol.isInsideVolume(axis_ice3)){
                flag = 1;
            }
             System.err.println("location; "+point.getParticleLocation());
             System.err.println("AxisLength: "+axis.getAxisLength());
                System.out.println("data "+axis.getAxisLength()+" 0.0 "+point.getMediumDensity()+" 0.0");
	}

        System.out.println("plot");
        System.out.println("disp");
        System.out.println("endg");

     }
        
}
  
