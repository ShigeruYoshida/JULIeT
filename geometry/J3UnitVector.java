package geometry;

import java.util.*;
import java.io.*;

/**
<pre>

  The class describes a unit vector.

</pre>

*/

public class J3UnitVector extends J3Vector {

    public J3UnitVector(double x, double y, double z){
	super(x,y,z);
	double l = this.getLength();
	if(l>0.0){
	    double norm = 1.0/l;
	    this.setAll(norm*x,norm*y,norm*z);
	}else{
	    System.err.println("This is a null vector!");
	    System.exit(0);
	}
    }
}
