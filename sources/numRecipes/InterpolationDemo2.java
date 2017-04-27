package numRecipes;

import numRecipes.*;
import java.io.*;

public class InterpolationDemo2 {


    public static void main(String[] args) throws IOException {

	char separator = ' ';
	double[ ] xx = new double[128];
	double[ ] yy = new double[128];
	double searchValue = 0;
	File fileName = null;

        if(args.length!=2){
            System.out.println("Usage: InterpolationDemo file-name x");
            System.exit(0);
        }else{
             fileName = new File(args[0]);
	     searchValue = Double.valueOf(args[1]).doubleValue( );
        }

	BufferedReader in =  
	    new BufferedReader(new FileReader(fileName));

	int n=0;
	String buffer;
	while ((buffer=in.readLine( ))!=null) {
	    int sep = buffer.indexOf(separator);
	    xx[n] = Double.valueOf(buffer.substring(0,sep)).doubleValue( );
	    sep = buffer.lastIndexOf(separator);
	    yy[n] = Double.valueOf(buffer.substring(sep+1)).doubleValue( );
	    n++;
	}



	int index = Interpolation.searchIndex(xx, searchValue, n);

	System.out.println("x= " + searchValue + 
			   " xx[" + index + "]= "+ xx[index]
			   + " yy[" + index + "]= "+ yy[index]
			   + " xx[" + (index+1) + "]= "+ xx[index+1]
			   + " yy[" + (index+1) + "]= "+ yy[index+1]);

	double yOut = 
	    Interpolation.mThPolynominalInterpolate(xx,yy,n,searchValue,6);

	System.out.println("Interpolated value = " + yOut);


    }

}
