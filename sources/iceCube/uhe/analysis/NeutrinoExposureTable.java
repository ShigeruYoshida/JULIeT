package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.analysis.*;
import numRecipes.*;

import java.io.*;
import java.util.*;

/**

Neutrino Exposure [cm^2 sec sr] is calculated from the numerical table
stored in the classes/iceCube/uhe/analysis.

Written by S. Yoshida 2012 May 30th

*/

public class NeutrinoExposureTable {

    private String TableFileName6yr = 
	"Exposure_6yAve_20150702.dat";
    private String pathToFile = 
	"iceCube/uhe/analysis/";
    private String[] headerFile= {"NuE","NuMu","NuTau"};

    private final int numberOfFlavor = 3;
    private int numLogEbin = 112;
    private int[] numLogEbinEffective = new int[numberOfFlavor];

    private double[][] areaArray = new double[numberOfFlavor][numLogEbin];
    private double[][] logEArray = new double[numberOfFlavor][numLogEbin];

    /** This is for nu-e : get rid of the glashow resonance */
    private boolean removeGlashowEnergyRange = false;
    private double logGlashowEnergyLower = 6.6;
    private double logGlashowEnergyUpper = 6.91;
    private double m2_to_cm2 = 1.0e4;
    private double km2_to_cm2 = 1.0e10;

    /**
       Constructor. Specify IC40-IC86III 6yr-based or not.
    */
    public NeutrinoExposureTable(boolean fullIceCube) throws IOException {

	if(!fullIceCube){
	    System.err.println(" must be fullICeCube=true");
	    System.exit(0);
	}
	
	// the flavor loop
	for(int flavor=0; flavor<numberOfFlavor;flavor++){
	    String tableFileName = pathToFile + headerFile[flavor] + TableFileName6yr;

	    // Reading out the table data
	    BufferedReader in =  
		new BufferedReader(new InputStreamReader(ClassLoader.getSystemResourceAsStream(tableFileName)));
	    char separator = ' ';int sepstart = 0; int sep = 0;
	    String buffer;
	    int iLogE=0;
	    while((buffer = in.readLine())!=null){
		boolean validRange = false;

		sepstart = 0;
		sep = buffer.indexOf(separator,sepstart+1);
		double logEnergy = Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double area = Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		if(!removeGlashowEnergyRange){
		    logEArray[flavor][iLogE] = logEnergy;
		    areaArray[flavor][iLogE] = area;
		    validRange = true;
		}else{
		    if(logGlashowEnergyLower > logEnergy ||logEnergy > logGlashowEnergyUpper){
			logEArray[flavor][iLogE] = logEnergy;
			areaArray[flavor][iLogE] = area;
			validRange = true;
		    }
		}

		//System.out.println(" logE= " + logEArray[flavor][iLogE] + 
		//		       " nu_area=" + areaArray[flavor][iLogE]);
		if(validRange) iLogE++;


	    }
	    numLogEbinEffective[flavor] = iLogE + 1;

	    // Place the edge point
	    logEArray[flavor][numLogEbinEffective[flavor]-1] = 12.0;
	    areaArray[flavor][numLogEbinEffective[flavor]-1]=
		areaArray[flavor][numLogEbinEffective[flavor]-2];
	    //System.err.format(" number of bins (flavor %d) = %d max logE=%f\n",
	    //		      flavor,numLogEbinEffective[flavor],
	    //		      logEArray[flavor][numLogEbinEffective[flavor]-2]);

	    in.close();
	}

    }


    /** Return the exposure of the primary inIce particle [cm^2 src sr].
	<pre>
	int flavor : 0 nu-e 1 nu-mu 2 nu-tau
	double logEnergy       : Log10(inIce Energy [GeV])
	</pre>
     */
    public double getExposure(int flavor, double logEnergy){

	if(logEnergy<5.4 ) return 0.0;


	if(logEnergy >= logEArray[flavor][numLogEbinEffective[flavor]-2]) logEnergy = logEArray[flavor][numLogEbinEffective[flavor]-2];
	if(logEnergy <= logEArray[flavor][0]) logEnergy = logEArray[flavor][0];

	int mTh = 7;
	if(flavor==0){
	    if(logEnergy >=logGlashowEnergyUpper) mTh=12;
	}

	double area = Interpolation.mThPolynominalInterpolate(logEArray[flavor],
		areaArray[flavor],numLogEbinEffective[flavor],logEnergy,mTh);
	if(area<0.0) area = 0.0;

	return area;

    }


    /** A simple main method for test */
    public static void main(String[] args) throws IOException{

	int flavor = 1;
	boolean isFullArray = true;
	if(args.length!=2){
	    System.out.println("Usage:NeutrinoEffAreaTable flavor isfullArray(1 for yes, 0 for IC40)");
	    System.exit(0);
	}else{
	    flavor = Integer.valueOf(args[0]).intValue();
	    if(Integer.valueOf(args[1]).intValue()==0) isFullArray=false;
	}

	NeutrinoExposureTable areaTable = new NeutrinoExposureTable(isFullArray);


	System.out.println("titx Energy [GeV]");
	System.out.println("tity Exposure [m^2! sec sr])");
	System.out.println("scal 1.0e5 1.0e11 1.0e7 5.0e+13");

	double logE = 5.0;
	while(logE <= 11.0){
	    double area = areaTable.getExposure(flavor,logE)/areaTable.m2_to_cm2;
	    double energy =  Math.pow(10.0,logE);
	    System.out.println("data " + energy + " 0.0 " +
				     area  + " 0.0");
	    logE += 0.02;
	}

	System.out.println("logx");
	System.out.println("logy");
	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");

	// plot the table data
	System.out.println("mksz 0.2");
	for(int i=0;i<areaTable.numLogEbinEffective[flavor];i++){
	    double energy = Math.pow(10.0,areaTable.logEArray[flavor][i]);
	    System.out.println("data " + energy + " 0.0 " +
			       areaTable.areaArray[flavor][i]/areaTable.m2_to_cm2 + " 0.0");
	}

	System.out.println("logx");
	System.out.println("logy");
	System.out.println("plot");
	System.out.println("disp");
	System.out.println("endg");
    }



}

