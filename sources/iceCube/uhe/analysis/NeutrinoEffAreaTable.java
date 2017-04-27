package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.analysis.*;
import numRecipes.*;

import java.io.*;
import java.util.*;

/**

Neutrino Effective Area is calculated from the numerical table
stored in the classes/iceCube/uhe/analysis.

The table data has been culculated by 
the Root Tree file by the weighting-module project
in the iceTray.

Written by S. Yoshida 2011 May 1st

*/

public class NeutrinoEffAreaTable {

    private String IC40TableFileName = 
	"iceCube/uhe/analysis/NeutrinoEffectiveArea_IC40.dat";
    private String TableFileName = 
	"iceCube/uhe/analysis/NeutrinoEffectiveArea_IC86.dat";

    private final int numberOfFlavor = 3;
    private int numLogEbin = 122;

    private double[ ][ ] areaArray = new double[numberOfFlavor][numLogEbin];
    private double[ ] logEArray = new double[numLogEbin];

    /** This is for nu-e : get rid of the glashow resonance */
    private boolean removeGlashowEnergyRange = false;
    private double logGlashowEnergyLower = 6.6;
    private double logGlashowEnergyUpper = 6.91;
    private double m2_to_cm2 = 1.0e4;

    /**
       Constructor. Specify IC40 or IC86;
       and the doublet variables in the argument. For example,
       flavor=1, doublet = 1, gives the effective area of inIce muons.
    */
    public NeutrinoEffAreaTable(boolean fullIceCube86) throws IOException {

	String tableFileName;
	if(fullIceCube86) tableFileName = new String(TableFileName);
	else  tableFileName = new String(IC40TableFileName);

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
	    double nu_e_area = Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( )*m2_to_cm2;
	    sepstart = sep;

	    sep = buffer.indexOf(separator,sepstart+1);
	    double nu_mu_area = Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( )*m2_to_cm2;
	    sepstart = sep;

	    sep = buffer.indexOf(separator,sepstart+1);
	    double nu_tau_area = Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( )*m2_to_cm2;
	    sepstart = sep;

	    if(!removeGlashowEnergyRange){
		logEArray[iLogE] = logEnergy;
		areaArray[0][iLogE] = nu_e_area;areaArray[1][iLogE] = nu_mu_area;areaArray[2][iLogE] = nu_tau_area;
		validRange = true;
	    }else{
		if(logGlashowEnergyLower > logEnergy ||logEnergy > logGlashowEnergyUpper){
		    logEArray[iLogE] = logEnergy;
		    areaArray[0][iLogE] = nu_e_area;areaArray[1][iLogE] = nu_mu_area;areaArray[2][iLogE] = nu_tau_area;
		    validRange = true;
		}
	    }

	    //System.out.println(" logE= " + logEArray[iLogE] + 
	    //		       " nu_e_area=" + areaArray[0][iLogE] +
	    //		       " nu_mu_area=" + areaArray[1][iLogE] +
	    //		       " nu_tau_area=" + areaArray[2][iLogE]);
	    if(validRange) iLogE++;


	}
	numLogEbin = iLogE + 1;

	// Place the edge point
	logEArray[numLogEbin-1] = 12.0;
	for(int j=0;j<numberOfFlavor;j++) 
	    areaArray[j][numLogEbin-1]=areaArray[j][numLogEbin-2];

	in.close();

    }


    /** Return the effective area of the primary inIce particle [cm^2].
	<pre>
	int flavor : 0 nu-e 1 nu-mu 2 nu-tau
	double logEnergy       : Log10(inIce Energy [GeV])
	</pre>
     */
    public double getArea(int flavor, double logEnergy){

	if(logEnergy<5.0 ) return 0.0;


	if(logEnergy >= logEArray[numLogEbin-2]) logEnergy = logEArray[numLogEbin-2];
	if(logEnergy <= logEArray[0]) logEnergy = logEArray[0];

	double area = Interpolation.mThPolynominalInterpolate(logEArray,
		areaArray[flavor],numLogEbin,logEnergy,7);
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

	NeutrinoEffAreaTable areaTable = new NeutrinoEffAreaTable(isFullArray);


	System.out.println("titx Energy [GeV]");
	System.out.println("tity Area [m^2!])");
	System.out.println("scal 1.0e5 1.0e11 1.0e-2 5.0e+4");

	double logE = 5.0;
	while(logE <= 11.0){
	    double area = areaTable.getArea(flavor,logE)/areaTable.m2_to_cm2;
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
			       areaTable.areaArray[flavor][i]/areaTable.m2_to_cm2  + " 0.0");
	}

	System.out.println("logx");
	System.out.println("logy");
	System.out.println("plot");
	System.out.println("disp");
	System.out.println("endg");
    }



}

