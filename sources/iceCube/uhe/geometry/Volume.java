package iceCube.uhe.geometry;

import geometry.*;

/**
    This class defines the cubic or cuboid volume.
    It is used for checking whether the particle is inside or outside of 
    the given volume.
*/

public class Volume{

    /** Size of the cubic volume [cm] */
    double sizeOfOneSide;      // case of cubic
    /** Size of the cuboid volume [cm] */
    private double xDimension;
    private double yDimension;
    private double zDimension;

    boolean isCubic = true;

    /** The vector from the origin to the center of a plane.
        It is nomal to the plane.
    */
    J3Vector p1,p2,p3,p4,p5,p6;

    /**
        Constructor. Set six planes.
        The origin of the IceCube coordinate (IceCubeCoordinate.java)
	is ASSUMED to be at the center of the volume.
    */
    public Volume(double sizeOfOneSide){

        double halfSize = sizeOfOneSide*0.5;
	this.sizeOfOneSide = sizeOfOneSide;

        p1 = new J3Vector(halfSize,0.0,0.0);
        p2 = new J3Vector(-halfSize,0.0,0.0);
        p3 = new J3Vector(0.0,halfSize,0.0);
        p4 = new J3Vector(0.0,-halfSize,0.0);
        p5 = new J3Vector(0.0,0.0,halfSize);
        p6 = new J3Vector(0.0,0.0,-halfSize);

    }

    /**
        Constructor. Set six planes.
        The origin of the IceCube coordinate (IceCubeCoordinate.java)
	is ASSUMED to be at the center of the volume.
	This constructor is for cuboid volume (e.g. for IceCube-gen2).
    */
    public Volume(double xD, double yD, double zD){

        p1 = new J3Vector(xD/2.0,0.0,0.0);
        p2 = new J3Vector(-xD/2.0,0.0,0.0);
        p3 = new J3Vector(0.0,yD/2.0,0.0);
        p4 = new J3Vector(0.0,-yD/2.0,0.0);
        p5 = new J3Vector(0.0,0.0,zD/2.0);
        p6 = new J3Vector(0.0,0.0,-zD/2.0);

	isCubic = false;
	xDimension = xD;
	yDimension = yD;
	zDimension = zD;

    }

    /**
       Tells if this volume is cubit. If not, it is a cuboid.
     */
    public boolean isCubic(){
	return isCubic;
    }

    /**
       return the X dimension [cm] of the cuboid volume
     */
    public double xDim(){
	if(isCubic) return xDimension;
	else{
	    System.err.println("This volume is cuboid!");
	    return 0.0;
	}
    }

    /**
       return the Y dimension [cm] of the cuboid volume
     */
    public double yDim(){
	if(isCubic) return yDimension;
	else{
	    System.err.println("This volume is cuboid!");
	    return 0.0;
	}
    }

    /**
       return the Z dimension [cm] of the cuboid volume
     */
    public double zDim(){
	if(isCubic) return zDimension;
	else{
	    System.err.println("This volume is cuboid!");
	    return 0.0;
	}
    }

    /**
        Check if the point described by IceCube coordinate is inside the volume.
        If it's inside, returns true.
    */
    public boolean isInsideVolume(J3Vector r){

        if(!isInsidePlane(r,p1))         return false;
        else if(!isInsidePlane(r,p2))    return false;
        else if(!isInsidePlane(r,p3))    return false;
        else if(!isInsidePlane(r,p4))    return false;
        else if(!isInsidePlane(r,p5))    return false;
        else if(!isInsidePlane(r,p6))    return false;
        else                             return true;     
    }

    /**
        Check if the point r - shift described by IceCube coordinate 
	is inside the volume.
        If it's inside, returns true.
    */
    public boolean isInsideVolume(J3Vector r, J3Vector shift){
	    J3Vector a = J3Vector.subtract(r,shift);
	    return isInsideVolume(a);
    }
    /** 
        Check if the point described by IceCube coordinate is inside the
        plane.
        If it's inside, returns true.
    */

    public static boolean isInsidePlane(J3Vector r, J3Vector plane){

        J3Vector rp = J3Vector.subtract(plane,r);
        if(J3Vector.getDotProduct(plane,rp) >= 0.0)   return true;
        else                                          return false;
    }

    /** Check if the given J3Line would pass inside the volome
	in the range of [axisLengthFrom, axisLengthTo] */
    public boolean isJ3LineInsideVolume(J3Line line, double axisLengthFrom,
					      double axisLengthTo) {
	if(axisLengthFrom>axisLengthTo){
	    System.err.println("axisLengthFrom must be smaller then axisLengthTo");
	    return false;
	}

	line.setAxisLength(axisLengthFrom);
	double deltaLength = 1.0e-3*(axisLengthTo-axisLengthFrom);
	double axisLength = axisLengthFrom;
	while(axisLength<axisLengthTo){
	    if(isInsideVolume(line)) return true;
	    axisLength += deltaLength;
	    line.setAxisLength(axisLength);
	}
	return false;
    }

    /** Check if the given J3Line - J3Vector would pass inside the volome
	in the range of [axisLengthFrom, axisLengthTo] */
    public boolean isJ3LineInsideVolume(J3Line line, J3Vector shift, 
					double axisLengthFrom, double axisLengthTo) {
	if(axisLengthFrom>axisLengthTo){
	    System.err.println("axisLengthFrom must be smaller then axisLengthTo");
	    return false;
	}

	line.setAxisLength(axisLengthFrom);
	double deltaLength = 1.0e-3*(axisLengthTo-axisLengthFrom);
	double axisLength = axisLengthFrom;
	while(axisLength<axisLengthTo){
	    if(isInsideVolume(line,shift)) return true;
	    axisLength += deltaLength;
	    line.setAxisLength(axisLength);
	}
	return false;
    }
}
