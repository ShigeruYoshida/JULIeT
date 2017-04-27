package geometry;

import java.util.*;
import java.io.*;

/**
<pre>

  This class provides a base to represent/handle three vector.

</pre>

*/

public class J3Vector implements Serializable {

    /** The vector components in three deminsional space */
    private double x;
    private double y;
    private double z;

    /** Vector's length */
    private double length;

    /** Constructor. For XYZ coordinate system. */
    public J3Vector(double x, double y, double z){
	this.x = x;
	this.y = y;
	this.z = z;
	setLength();
    }

    /** Constructor. Default setting: x=y=0, z=1.0 */
    public J3Vector(){
	this(0.0,0.0,1.0);
    }


    /** set x component */
    public void setX(double x){
	this.x = x;
	setLength();
    }
    /** set y component */
    public void setY(double y){
	this.y = y;
	setLength();
    }
    /** set z component */
    public void setZ(double z){
	this.z = z;
	setLength();
    }
    /** set all xyz components */
    public void setAll(double x, double y, double z){
	this.x = x;
	this.y = y;
	this.z = z;
	setLength();
    }

    /** put a new vector a = b */
    public void putVector(J3Vector a){
	double xx = a.getX();
	double yy = a.getY();
	double zz = a.getZ();
	this.setX(xx);
	this.setY(yy);
	this.setZ(zz);
	setLength();
    }



    /** get x component */
    public double getX(){
	return x;
    }
    /** get y component */
    public double getY(){
	return y;
    }
    /** get z component */
    public double getZ(){
	return z;
    }

    /** calculate the vector length to set*/
    public void setLength(){
	double sqLength = x*x + y*y + z*z;
	if(sqLength>0.0) length = Math.sqrt(sqLength);
	else length = 0.0;
    }

    /** get the vector length */
    public double getLength(){
	return length;
    }

    /** calculate the dot-product */
    public static double getDotProduct(J3Vector a, J3Vector b){
	double dot = a.getX()*b.getX() + a.getY()*b.getY() + a.getZ()*b.getZ();
	return dot;
    }

    /** calculate the angle [rad] between two vectors.
        Return 0 if either of the two vectors is 0. */
    public static double getAngleInRadian(J3Vector a, J3Vector b){
	double dot = getDotProduct(a,b);
	double aLength = a.getLength();
	double bLength = b.getLength();
	double theta = 0.0;
	if(aLength>0.0 && bLength>0.0){
	    double cosTheta = dot/(aLength*bLength);
            if(cosTheta > 1.0){
                cosTheta = 1.0;
            }else if(cosTheta < -1.0){
                cosTheta = -1.0;
            }
	    theta = Math.acos(cosTheta);
	}
	return theta;
    }
    /** calculate the angle [deg] between two vectors.
        Return 0 if either of the two vectors is 0. */
    public static double getAngleInDegree(J3Vector a, J3Vector b){
	return Math.toDegrees(getAngleInRadian(a,b));
    }


    /** add vectors a+b */
    public static J3Vector add(J3Vector a, J3Vector b){
	double xx = a.getX() + b.getX();
	double yy = a.getY() + b.getY();
	double zz = a.getZ() + b.getZ();
	J3Vector c = new J3Vector(xx,yy,zz);
	return c;
    }

    /** sabtract vectors a-b */
    public static J3Vector subtract(J3Vector a, J3Vector b){
	double xx = a.getX() - b.getX();
	double yy = a.getY() - b.getY();
	double zz = a.getZ() - b.getZ();
	J3Vector c = new J3Vector(xx,yy,zz);
	return c;
    }

    /** increment vectors a+=b */
    public void increment(J3Vector a, J3Vector b){
	double xx = a.getX() + b.getX();
	double yy = a.getY() + b.getY();
	double zz = a.getZ() + b.getZ();
	a.setX(xx);
	a.setY(yy);
	a.setZ(zz);
    }

    /** decrement vectors a-=b */
    public void decrement(J3Vector a, J3Vector b){
	double xx = a.getX() - b.getX();
	double yy = a.getY() - b.getY();
	double zz = a.getZ() - b.getZ();
	a.setX(xx);
	a.setY(yy);
	a.setZ(zz);
    }

    /** calculate the cross Products axb */
    public static J3Vector getCrossProduct(J3Vector a, J3Vector b){
	double xx = a.getY()*b.getZ()-a.getZ()*b.getY();
	double yy = a.getZ()*b.getX()-a.getX()*b.getZ();
	double zz = a.getX()*b.getY()-a.getY()*b.getX();
	J3Vector c = new J3Vector(xx,yy,zz);
	return c;
    }

    /** calculate scalor x vector fxa */
    public static J3Vector multipleFactor(double f, J3Vector a){
	double xx = f*a.getX();
	double yy = f*a.getY();
	double zz = f*a.getZ();
	J3Vector c = new J3Vector(xx,yy,zz);
	return c;
    }

}
