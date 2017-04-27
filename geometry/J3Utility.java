package geometry;

import java.util.*;
import java.io.*;

/**
<pre>

  This class provides some utility for common geometrical calculation.

</pre>

*/

public class J3Utility {

    /** returns distance between a line and a point */
    public static double getDistanceFromLineToPoint(J3Line line, J3Vector r){
        J3Line line_copy = new J3Line(line.getR0(),line.getDirection(),
                                      line.getAxisLength());
	J3UnitVector n = line_copy.getDirection();
	J3Vector r0ToR = J3Vector.subtract(r,line_copy.getR0());
	double axisAtPerpendicular = J3Vector.getDotProduct(n,r0ToR);
	line_copy.setAxisLength(axisAtPerpendicular);
	J3Vector axisToPoint = J3Vector.subtract(r,line_copy);
	double distance = axisToPoint.getLength();

	r0ToR = null;
	axisToPoint = null;
        line_copy = null;

	return distance;
    }

    /** Move the line vector to where it crosses a plane with
	a fixed point vector rPlane and plane direction nPlane */
    public static void moveLineVectorToPlane(J3Line line, 
					     J3UnitVector nPlane, 
					     J3Vector rPlane){
	J3UnitVector n = line.getDirection();
	J3Vector r0 = line.getR0();
	J3Vector d = J3Vector.subtract(r0, rPlane);
	double l = -J3Vector.getDotProduct(nPlane,d)/
	    J3Vector.getDotProduct(nPlane,n);

	line.setAxisLength(l);
	d = null;
    }

    /**
       Set J3Line's axis length so that its absolute length |J3Line|
       is equal to the given length.
       <pre>
       J3Line line        :   J3Line vector you want to set.
       lineLength         :   The value you have for |line|.
       </pre>
       Choose the positive value of axis length.
    */ 

    public static void setJ3LinePositiveAxisLengthForGivenLength(J3Line line,
								 double lineLength){

        J3Vector r0 = line.getR0();
        J3UnitVector n = line.getDirection();
        double r0Square = r0.getLength()*r0.getLength();
        double nSquare = n.getLength()*n.getLength(); // should be 1
        double dot = J3Vector.getDotProduct(r0,n);

        double l = (-dot+Math.sqrt(dot*dot-nSquare*(r0Square-lineLength*lineLength)))
                   /nSquare;
	line.setAxisLength(l);
    }

    /**
       Set J3Line's axis length so that its absolute length |J3Line|
       is equal to the given length.
       <pre>
       J3Line line        :   J3Line vector you want to set.
       lineLength         :   The value you have for |line|.
       </pre>
       Choose the negative value of axis length.
    */ 
    public static void setJ3LineNegativeAxisLengthForGivenLength(J3Line line,
								 double lineLength){

        J3Vector r0 = line.getR0();
        J3UnitVector n = line.getDirection();
        double r0Square = r0.getLength()*r0.getLength();
        double nSquare = n.getLength()*n.getLength(); // should be 1
        double dot = J3Vector.getDotProduct(r0,n);

        double l = -(dot+Math.sqrt(dot*dot-nSquare*(r0Square-lineLength*lineLength)))
                   /nSquare;
	line.setAxisLength(l);
    }

    /**
       Set J3Line's axis length so that its absolute length |J3Line-J3Vector|
       is equal to the given length.
       <pre>
       J3Line line        :   J3Line vector you want to set.
       J3Vector r         :   J3Vector from which the vector goes to J3Line above
       lineLength         :   The value you have for |line|.
       </pre>
       Choose the positive value of axis length.
    */ 

    public static void setJ3LinePositiveAxisLengthForGivenLength(J3Line line,
								 J3Vector r,
								 double lineLength){

        J3Vector r0 = J3Vector.subtract(line.getR0(),r);
        J3UnitVector n = line.getDirection();
        double r0Square = r0.getLength()*r0.getLength();
        double nSquare = n.getLength()*n.getLength(); // should be 1
        double dot = J3Vector.getDotProduct(r0,n);

        double l = (-dot+Math.sqrt(dot*dot-nSquare*(r0Square-lineLength*lineLength)))
                   /nSquare;
	line.setAxisLength(l);
    }

    /**
       Set J3Line's axis length so that its absolute length |J3Line-J3Vector|
       is equal to the given length.
       <pre>
       J3Line line        :   J3Line vector you want to set.
       J3Vector r         :   J3Vector from which the vector goes to J3Line above
       lineLength         :   The value you have for |line|.
       </pre>
       Choose the negative value of axis length.
    */ 
    public static void setJ3LineNegativeAxisLengthForGivenLength(J3Line line,
								 J3Vector r,
								 double lineLength){

        J3Vector r0 = J3Vector.subtract(line.getR0(),r);
        J3UnitVector n = line.getDirection();
        double r0Square = r0.getLength()*r0.getLength();
        double nSquare = n.getLength()*n.getLength(); // should be 1
        double dot = J3Vector.getDotProduct(r0,n);

        double l = -(dot+Math.sqrt(dot*dot-nSquare*(r0Square-lineLength*lineLength)))
                   /nSquare;
	line.setAxisLength(l);
    }
    
    /* Return the intersection of the line and its perpendicular which includes 
       the point "r"
       <pre>
       J3Line line :     J3Line vector for given line.
       J3Vector r  :     J3Vector for given point.
       </pre>
     */
      
    public static J3Vector getThePerpendicularOnLineFromPoint(J3Line line, J3Vector r){
        J3UnitVector n = line.getDirection();
        J3Vector r0 = line.getR0();
        J3Vector r0ToR = J3Vector.subtract(r,line.getR0());
    	double axisAtPerpendicular = J3Vector.getDotProduct(n,r0ToR);
    	J3Vector r0ToNearest = J3Vector.multipleFactor(axisAtPerpendicular,n);
    	J3Vector nearestPoint = J3Vector.add(r0,r0ToNearest);
    	return nearestPoint;
    }
}


