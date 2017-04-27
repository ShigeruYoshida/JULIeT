package geometry;

import java.util.*;
import java.io.*;

/**
<pre>

  This class describes line in three dimensional space.
  It provides a point vector along a given line
  represented as  r = r0 + l x n where l is the length
  long the line, n is the unit vector.
  and the relevant operations.

</pre>

*/

public class J3Line extends J3Vector {

    /** Length along the line */
    double axisLength;
    /** The unit vector to define the direction */
    J3UnitVector n;
    /** a fixed point where the line passes */
    J3Vector r0;


    /** Constructor */
    public J3Line(J3Vector fixedPoint, J3UnitVector n, double l){
	super();
	this.axisLength = l;
	this.r0 = fixedPoint;
	this.n = n;
	configure();
    }
    public J3Line(J3Vector fixedPoint, J3UnitVector n){
	this(fixedPoint,n,0.0);
    }

    /** Configuration. The point vector components defined by
	the gemetrical parameters are set.*/
    public void configure(){
	J3Vector d = multipleFactor(axisLength,n);
	J3Vector r = add(r0, d);
	putVector(r);
	r = null;
	d = null;
    }

    /** set a new axis length l */
    public void setAxisLength(double l){
	this.axisLength = l;
	configure();
    }

    /** set a new passing point */
    public void setR0(J3Vector r0){
	this.r0 = r0;
	configure();
    }

    /** set a new direction */
    public void setDirection(J3UnitVector n){
	this.n = n;
	configure();
    }

    /** return the current axis length l */
    public double getAxisLength(){
	return axisLength;
    }

    /** return the current direction */
    public J3UnitVector getDirection(){
	return n;
    }

    /** return the current passing point */
    public J3Vector getR0(){
	return r0;
    }
}



