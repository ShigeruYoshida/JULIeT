package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.analysis.*;
import numRecipes.*;

import java.io.*;
import java.util.*;

/**

Neutrino Exposure by the HESE analysis [cm^2 sec sr] is calculated from the numerical table
stored in the classes/iceCube/uhe/analysis.

Written by S. Yoshida 2013 July 18th

*/

public class NeutrinoExposureTableHESE {

    private String TableFileNameHESE = 
	"iceCube/uhe/analysis/NeutrinoExposureHESE2012.dat";

    private final int numberOfFlavor = 3;
    private int numLogEbin = 51;

    private NeutrinoExposureTable eheExposure = null; // EHE exposure
    private boolean addEHE = false;

    private double[ ][ ] areaArray = new double[numberOfFlavor][numLogEbin];
    private double[ ] logEArray = new double[numLogEbin];

    private double m2_to_cm2 = 1.0e4;
    private double km2_to_cm2 = 1.0e10;


    /**
       Constructor.
    */
    public NeutrinoExposureTableHESE(boolean addEHE) throws IOException {

	this.addEHE = addEHE;

	String tableFileName = new String(TableFileNameHESE);

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
	    double nu_e_area = Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
	    sepstart = sep;

	    sep = buffer.indexOf(separator,sepstart+1);
	    double nu_mu_area = Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
	    sepstart = sep;

	    sep = buffer.indexOf(separator,sepstart+1);
	    double nu_tau_area = Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
	    sepstart = sep;

	    logEArray[iLogE] = logEnergy;
	    areaArray[0][iLogE] = nu_e_area;areaArray[1][iLogE] = nu_mu_area;areaArray[2][iLogE] = nu_tau_area;
	    //System.out.println(" logE= " + logEArray[iLogE] + 
	    //		       " nu_e_area=" + areaArray[0][iLogE] +
	    //		       " nu_mu_area=" + areaArray[1][iLogE] +
	    //		       " nu_tau_area=" + areaArray[2][iLogE]);


	    iLogE++;

	}
	numLogEbin = iLogE + 1;

	// Place the edge point
	logEArray[numLogEbin-1] = 7.0;
	for(int j=0;j<numberOfFlavor;j++) 
	    areaArray[j][numLogEbin-1]=areaArray[j][numLogEbin-2];

	in.close();

	// add EHE IC40+IC79+IC96 exposure
	if(addEHE) eheExposure = new NeutrinoExposureTable(true);

    }


    /** Return the exposure of the primary inIce particle [cm^2 src sr].
	<pre>
	int flavor : 0 nu-e 1 nu-mu 2 nu-tau
	double logEnergy       : Log10(inIce Energy [GeV])
	</pre>
     */
    public double getExposure(int flavor, double logEnergy){

	if(logEnergy<3.5 ) return 0.0;

	double area = 0.0;
	if(addEHE) area = eheExposure.getExposure(flavor,logEnergy);

	if(logEnergy >= logEArray[numLogEbin-2]) logEnergy = logEArray[numLogEbin-2];
	if(logEnergy <= logEArray[0]) logEnergy = logEArray[0];

	double areaHESE = Interpolation.mThPolynominalInterpolate(logEArray,
		areaArray[flavor],numLogEbin,logEnergy,7);
	if(areaHESE<0.0) areaHESE = 0.0;

	if(areaHESE>area) area = areaHESE;

	return area;

    }


    /** A simple main method for test */
    public static void main(String[] args) throws IOException{

	int flavor = 1;
	boolean addEHE = true;
	if(args.length!=2){
	    System.out.println("Usage:NeutrinoExposureTableHESE flavor addEHE(1 for yes, 0 no)");
	    System.exit(0);
	}else{
	    flavor = Integer.valueOf(args[0]).intValue();
	    if(Integer.valueOf(args[1]).intValue()==0) addEHE=false;
	}

	NeutrinoExposureTableHESE areaTable = new NeutrinoExposureTableHESE(addEHE);


	System.out.println("titx Energy [GeV]");
	System.out.println("tity Exposure [m^2! sec sr])");
	System.out.println("scal 1.0e3 1.0e11 1.0e8 3.0e+13");

	double logE = 3.0;
	while(logE <= 11.0){
	    double area = areaTable.getExposure(flavor,logE)/areaTable.m2_to_cm2;
	    double energy =  Math.pow(10.0,logE);
	    System.out.println("data " + energy + " 0.0 " +
				     area  + " 0.0");
	    logE += 0.1;
	}

	System.out.println("logx");
	System.out.println("logy");
	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");

	// plot the table data
	System.out.println("mksz 0.1");
	for(int i=0;i<areaTable.numLogEbin;i++){
	    double energy = Math.pow(10.0,areaTable.logEArray[i]);
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

