package iceCube.uhe.analysis;

import java.io.*;

public class DrawEffAreaTable {

    public static void main(String[] args) throws IOException{

	int flavor = 1;
	int doublet = 0;
	double logE = 5.0;
	if(args.length!=3){
            System.out.println("Usage: DrawEffAreaTable flavor doublet logEnergy");
            System.exit(0);
	}else{
	    flavor = Integer.valueOf(args[0]).intValue();
	    doublet = Integer.valueOf(args[1]).intValue();
	    logE = Double.valueOf(args[2]).doubleValue();
	}

	EffAreaTable areaTable = new EffAreaTable(flavor,doublet);


        System.out.println("titx cosZenith ");
        System.out.println("tity Area [km^2!]");

	double cosZenith = -1.0;
	while(cosZenith <= 1.0){
	    double area = areaTable.getArea(logE,cosZenith);
	    if(area>0.0)
		System.out.println("data " + cosZenith + " 0.0 " +
				     area  + " 0.0");
	    cosZenith += 0.05;
	}

	System.out.println("logy");
	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");
    }


}
