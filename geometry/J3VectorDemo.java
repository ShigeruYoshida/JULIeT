package geometry;

import geometry.*;

import java.util.*;
import java.io.*;

/** A simple Demo program to show usage of the J3Vector class */

public class J3VectorDemo {

    public static void main(String[] args){

	double theta1 = 60.0; //[deg]
	double alpha1 = -30.0; // [deg]
	double theta2 = 45.0; //[deg]
	double alpha2 = -30.0; // [deg]

	double x1 = Math.sin(Math.toRadians(theta1))*Math.sin(Math.toRadians(alpha1));
	double y1 = Math.sin(Math.toRadians(theta1))*Math.cos(Math.toRadians(alpha1));
	double z1 = Math.cos(Math.toRadians(theta1));
	double x2 = Math.sin(Math.toRadians(theta2))*Math.sin(Math.toRadians(alpha2));
	double y2 = Math.sin(Math.toRadians(theta2))*Math.cos(Math.toRadians(alpha2));
	double z2 = Math.cos(Math.toRadians(theta2));

	J3Vector a = new J3Vector(x1,y1,z1);
	J3Vector b = new J3Vector(x2,y2,z2);

	System.out.println("Z(a) = " + a.getZ() );
	System.out.println("Y(b) = " + b.getY() );
	System.out.println("|a| = " + a.getLength() );
	System.out.println("|b| = " + b.getLength() );
	System.out.println("a*b = " + J3Vector.getDotProduct(a,b) );
	System.out.println("angle(a,b) = " + J3Vector.getAngleInDegree(a,b) + " [deg]");

	J3Vector c = J3Vector.multipleFactor(2.0,a);
	System.out.println("|2a| = " + c.getLength() );
	J3Vector d = J3Vector.getCrossProduct(a,b);
	System.out.println("|a x b| = " + d.getLength() );
	J3Vector e = J3Vector.add(a,b);
	System.out.println("|a + b| = " + e.getLength() );
	J3UnitVector f = new J3UnitVector(c.getX(),c.getY(),c.getZ());
	System.out.println("|a/|a|| = " + f.getLength() );

    }
}





