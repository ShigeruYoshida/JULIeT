package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.analysis.*;
import numRecipes.*;

import java.io.*;
import java.util.*;

/**

(In-ice) Effective Area is calculated from the numerical table
stored in the classes/iceCube/uhe/analysis.

The table data has been culculated by 
the Root Tree file by the weighting-module project
in the iceTray.

This class is used for calculating 
the inIce event rate, the detectoable neutrino flux
at the earth surface, and the neutrino yield [cm^2 sec sr]
in PropagationMatrixFlux.java in the analysis project.

You can alternativeldy calculate these stuff by
I3ParticleFlux that directly uses I3Particle events
filled with Propation Matrix by I3ParticlePropMatrixFiller.

Written by S. Yoshida 2007 April 14th

*/

public class EffAreaTable {

    private String muonTableFileName = 
	"iceCube/uhe/analysis/EffectiveAreaMuon.dat";

    private String tauTableFileName = 
	"iceCube/uhe/analysis/EffectiveAreaTau.dat";

    private String nuETableFileName = 
	"iceCube/uhe/analysis/EffectiveAreaNuE.dat";

    private String nuMuTableFileName = 
	"iceCube/uhe/analysis/EffectiveAreaNuMu.dat";

    private String nuTauTableFileName = 
	"iceCube/uhe/analysis/EffectiveAreaNuTau.dat";

    private int numCosZenithBin = 40;
    private int numLogEbin = 31; // 30 data points + one edge

    private double[ ][ ] areaArray = new double[numCosZenithBin][numLogEbin];
    private double[ ] logEArray = new double[numLogEbin];
    private double[ ] cosZenithArray = new double[numCosZenithBin];

    /** This is for nu-e : get rid of the glashow resonance */
    private boolean removeGlashowEnergyRange = false;
    private double logGlashowEnergyLower = 6.6;
    private double logGlashowEnergyUpper = 6.91;

    /**
       Constructor. Specify particle spiece by giving the flavor
       and the doublet variables in the argument. For example,
       flavor=1, doublet = 1, gives the effective area of inIce muons.
    */
    public EffAreaTable(int flavor, int doublet) 
	throws IOException {

	String tableFileName;
	switch (doublet) {

	case 1: // this is charged leption
	    if(flavor==1) tableFileName = new String(muonTableFileName);
	    else if(flavor==2) tableFileName = new String(tauTableFileName);
	    else tableFileName = null;
	    break;

	case 0: // this is neutrino
	    if(flavor==0){
		tableFileName = new String(nuETableFileName);
		removeGlashowEnergyRange = true;
	    }else if(flavor==1) tableFileName = new String(nuMuTableFileName);
	    else if(flavor==2) tableFileName = new String(nuTauTableFileName);
	    else tableFileName = null;
	    break;

	default: // WRONG!
	    tableFileName = null;
	}

	if(tableFileName == null){
	    System.err.println(Particle.particleName(flavor, doublet) +
			       " is not a inIce particle!");
	    System.exit(0);
	}


	// Reading out the table data
	BufferedReader in =  
	    new BufferedReader(new InputStreamReader(ClassLoader.getSystemResourceAsStream(tableFileName)));
        char separator = ' ';
	int sepstart = 0; int sep = 0;
        String buffer;
  
	int iLogE=0;
	for(int i=0;i<numLogEbin-1;i++){
	    boolean validRange = false;
	    for(int j=0;j<numCosZenithBin;j++){

		buffer=in.readLine( );

		sepstart = 0;
		//sep = buffer.indexOf(separator,sepstart+1);
		//sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		int ie =
  	        Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double area =
  	        Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		cosZenithArray[j] = -0.975+0.05*(double )j;
		double logEnergy = 5.1+0.2*(double )(ie-1);

		if(!removeGlashowEnergyRange){
		    logEArray[iLogE] = logEnergy;
		    areaArray[j][iLogE] = area;
		    validRange = true;
		}else{
		    if(logGlashowEnergyLower > logEnergy ||
		       logEnergy > logGlashowEnergyUpper){
			logEArray[iLogE] = logEnergy;
			areaArray[j][iLogE] = area;
			validRange = true;
		    }
		}

		//System.out.println(i + " " + j + 
		//		   " logE= " + logEArray[i] + " cosZ= " +
		//		   cosZenithArray[j]  + " area=" + areaArray[j][i]);

	    }
	    if(validRange) iLogE++;


	}
	numLogEbin = iLogE + 1;

	// Place the edge point
	logEArray[numLogEbin-1] = 12.0;
	for(int j=0;j<numCosZenithBin;j++) 
	    areaArray[j][numLogEbin-1]=areaArray[j][numLogEbin-2];

	in.close();

    }


    /** Return the effective area of the primary inIce particle [cm^2].
	Note: This method (and this entire class ) is concerned with
	the inIce particle. For connecting it with the earth surface particle,
	you shuold use PropagationMatrixFlux class that reads out
        propagation matrix data and convolute it with this class.
	<pre>
	double logEnergy       : Log10(inIce Energy [GeV])
	double cosZenith       : cosine of in-ice Zenith angle
	</pre>
     */
    public double getArea(double logEnergy, double cosZenith){

	if(logEnergy<6.0 ) return 0.0;

	if(cosZenith >= cosZenithArray[numCosZenithBin-1]) 
	    cosZenith = cosZenithArray[numCosZenithBin-1];
	if(cosZenith <= cosZenithArray[0]) 
	    cosZenith = cosZenithArray[0];

	int index = Interpolation.searchIndex(cosZenithArray,cosZenith,numCosZenithBin);

	if(logEnergy >= logEArray[numLogEbin-2]) logEnergy = logEArray[numLogEbin-2];
	if(logEnergy <= logEArray[0]) logEnergy = logEArray[0];

	// cosZenith loop
	int mData = 3;
	index -= (mData-1)/2; if(index<0) index = 0;
	double[] areaCosArray = new double[mData];
	double[] cosArray = new double[mData];
	int numberOfData = 0;
	for(int i= 0; i<mData && index<numCosZenithBin;i++){
	    cosArray[i] = cosZenithArray[index];
	    areaCosArray[i] = 
		Interpolation.mThPolynominalInterpolate(logEArray,
		areaArray[index],numLogEbin,logEnergy,7);
	    if(areaCosArray[i]<0.0) areaCosArray[i]=0.0;
	    numberOfData++;
	    index++;
	}
	double area = 0.0;
	if(numberOfData == mData){
	    area = Interpolation.mThPolynominalInterpolate(cosArray,
		   areaCosArray,numberOfData,cosZenith,numberOfData);
	}else{
	    area = areaCosArray[numberOfData-1];
	}

	//double area = 
	//    Interpolation.mThPolynominalInterpolate(logEArray,
	//					    areaArray[index],numLogEbin,
	//					    logEnergy,7);

	//System.out.println(" " + cosZenithLow + 
	//	   " areaLow= " + areaLow + " area= " + area);
	if(area<0.0) area = 0.0;

	return area;

    }


    /** A simple main method for test */
    public static void main(String[] args) throws IOException{

	int flavor = 1;
	int doublet = 0;
	double cosZenith = -1.0;
	if(args.length!=3){
            System.out.println("Usage: EffAreaTable flavor doublet cosTheta");
            System.exit(0);
	}else{
	    flavor = Integer.valueOf(args[0]).intValue();
	    doublet = Integer.valueOf(args[1]).intValue();
	    cosZenith = Double.valueOf(args[2]).doubleValue();
	}

	EffAreaTable areaTable = new EffAreaTable(flavor,doublet);


        System.out.println("titx Energy [GeV]");
        System.out.println("tity Area [km^2!])");
        if(doublet == 1) System.out.println("scal 1.0e5 1.0e11 1.0e-5 2.5");
	else  System.out.println("scal 1.0e5 1.0e11 1.0e-6 2.5e-2");

	double logE = 5.0;
	while(logE <= 11.0){
	    double area = areaTable.getArea(logE,cosZenith);
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
	int index = Interpolation.searchIndex(areaTable.cosZenithArray,cosZenith,
					      areaTable.numCosZenithBin);

	double cosZenithTable = areaTable.cosZenithArray[index];
	System.out.println("mssg cosZenithTable = " + cosZenithTable);
	for(int i=0;i<areaTable.numLogEbin;i++){
	    double energy = Math.pow(10.0,areaTable.logEArray[i]);
	    System.out.println("data " + energy + " 0.0 " +
			       areaTable.areaArray[index][i]  + " 0.0");
	}

	System.out.println("mksz 0.5");
	System.out.println("logx");
	System.out.println("logy");
	System.out.println("plot");
	System.out.println("disp");
	System.out.println("endg");
    }



}

