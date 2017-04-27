package iceCube.uhe.muonModel;

import iceCube.uhe.particles.*;

import java.io.*;

/**

   This class handles the table containing 
   the atmospheric muon flux calculated by PropagatingAtmMuonFlux
   class PropagatingAtmMuonFlux.getDFMuDLogE(logE,cosZenith)
   for combinations of alpha and Eth in AtmMuonBundleFlux.

   While direct use of PropagatingAtmMuonFlux must reads 
   propagation matrix data each time which takes I/O time,
   this class does not need to read the matrix data 
   that realizes faster processing. It is mainly used
   in AtmMuonBunfleFitter static class in the analysis package.


   The main method of MakeElbertFluxTable class generates the table file
   this class handles. In order to reduce the data size,
   fluxes with only the limited range of Eth in the Elbert model 
   that can be consistent  with the IceCube data are stored in the table.
   This range is stored in the protected vauluables 
   <pre>
    protected double[] muEThMin;       // [alpha]
    protected double[] muEThStepSize;  // [alpha]
    public int[] numberOfMuEThSteps;// [alpha]
   </pre>
   in each of alpha bins. The range of the alpha parameter of the Elbert formola
   is defined by the protected member valuables
   <pre>
    protected double alphaMin = 1.85;
    protected double alphaStepSize = 0.01;
    public final int numberOfAlphaSteps = 36;

   alpha =  alphaMin + alphaStepSize x alphabin
   </pre>

   Check out the simple main method for demonstration of this factory.

   Written by S. Yoshida 2008 June 2nd

 */


public class ElbertFluxTableFactory {

    private static final double ln10 = Math.log(10.0);
    private double roundoff = 1.0e-3;

    /** Path to the directory where the flux tables are stored.*/
    public static String tablePath = "../data/muonBundleFluxTable/";
    protected static final String[] tableFileName = 
    {
    "1e5_00cm.table", "1e5_05cm.table", "1e5_10cm.table", "1e5_15cm.table", "1e5_20cm.table", 
    "1e5_25cm.table", "1e5_30cm.table", "1e5_35cm.table", "1e5_40cm.table", "1e5_45cm.table", 
    "1e5_50cm.table", "1e5_55cm.table", "1e5_60cm.table", "1e5_65cm.table", "1e5_70cm.table", 
    "1e5_75cm.table", "1e5_80cm.table", "1e5_85cm.table", "1e5_90cm.table", "1e5_95cm.table", 
    "1e6_00cm.table", "1e6_10cm.table", "1e6_20cm.table", "1e6_30cm.table", "1e6_40cm.table", 
    "1e6_50cm.table", "1e6_60cm.table", "1e6_70cm.table", "1e6_80cm.table", "1e6_90cm.table", 
    "1e7_00cm.table", "1e7_10cm.table", "1e7_20cm.table", "1e7_30cm.table", "1e7_40cm.table", "1e7_50cm.table"};
    protected static final int numberOfDistance = tableFileName.length;

    protected static double[] icedist ={//distance in log10(distance[cm])
      5.00,  5.025, 5.075, 5.125, 5.175, 5.225, 5.275, 5.325, 5.375, 5.425, 
      5.475, 5.525, 5.575, 5.625, 5.675, 5.725, 5.775, 5.825, 5.875, 5.925, 
      5.975, 6.05,  6.15,  6.25,  6.35,  6.45,  6.55,  6.65,  6.75,  6.85, 
      6.95,  7.05,  7.15,  7.25,  7.35,  7.45, 7.50};

    // range of alpha
    /** Minimum bound of alpha in the tanle. any value of alpha larger than this
	works in this class. Determined by range of the flux table generated
        by the MakeElbertFluxTable in the muonModel package. */
    public static double alphaMin = 1.9;
    /** Step size of the alpha in the table */
    public static double alphaStepSize = 0.01;
    /** Number of bins of the alpha in the table */
    public static int numberOfAlphaSteps = 46;

    // range of muon Energy Threshold [GeV]
    protected double[] muEThMin;       // [alpha]
    protected double[] muEThStepSize;  // [alpha]
    public int[] numberOfMuEThSteps;// [alpha]
    /** Number of bins of the Muon Energy Threshold  in the table */
    public static int maxNumberOfMuEThSteps = 46;//

    // range of log10(muon energy [GeV])
    protected double logEmuonMin = Particle.getLogEnergyMinimum();
    protected double logEmuonStepSize = 5.0*Particle.getDeltaLogEnergy();
    protected int numberOfEmuon = (int )(Particle.getDimensionOfLogEnergyMatrix()/5.0+0.01);

    protected double[][][][] dFMuDLogEarray;//[alpha][muEth][distance][Emu]

    // bin number to specify the alpha value in the bundle flux. 
    protected int alphaBin = 0;
    // bin number to specify the muEth value in the bundle flux. 
    protected int muEThBin = 0;
    // index to tell whether those bin numbers are set by a user
    // with the method setElbertParameters
    private boolean elbertParametersHaveBeenSet = false;

    /** 
	Default constrctor : Reads the table data from all the files
        tablePath + tableFileName[i].
    */
    public ElbertFluxTableFactory() throws IOException{
	generateTableArrays();
	readTableData();
    }

    /** 
	This consrtuctor reads the single table data from the fine
	given in the argument.
	<pre>
	String fileName   : name of the table file name you read
	int distanceIndex : index of distance to contain the flux data in the array
	</pre>
    */
    public ElbertFluxTableFactory(String fileName, int distanceIndex) 
	throws IOException{
	generateTableArrays();
	readTableData(fileName,distanceIndex);
    }

    private void generateTableArrays(){
	muEThMin = new double[numberOfAlphaSteps];
	muEThStepSize = new double[numberOfAlphaSteps];
	numberOfMuEThSteps = new int[numberOfAlphaSteps];
	dFMuDLogEarray = new double[numberOfAlphaSteps][maxNumberOfMuEThSteps][numberOfDistance][numberOfEmuon];
    }


    /** Reading out the table data files from tableFileName[i] */
    protected void readTableData() throws IOException {

	for(int i =0 ; i< numberOfDistance ; i++){
	    String fileName = tablePath + tableFileName[i];
	    System.err.println("Reading the table from File " + fileName);
	    readTableData(fileName,i);
	}
    }

    /** Reading out the table data file with name fileName
	<pre>
	String fileName   : name of the table file name you read
	int distanceIndex : index of distance to contain the flux data in the array
	</pre>
    */
    protected void readTableData(String fileName, int distanceIndex) throws IOException {

	BufferedReader in = 
	    new BufferedReader(new FileReader(fileName));

        char separator = ' ';
	int sepstart = 0; int sep = 0;
        String buffer;

	buffer = in.readLine(); // read alpha range info : skip

	// alpha loop
	for(int ialpha =0;ialpha<numberOfAlphaSteps;ialpha++){
	    buffer = in.readLine();

	    sepstart = 0;
	    sep = buffer.indexOf(separator,sepstart+1);
	    muEThMin[ialpha]=
		Double.valueOf(buffer.substring(sepstart,sep)).doubleValue( );
	    //System.err.println("muEmin sepstart=" + sepstart + " sep="+ sep);

	    sepstart = sep;
	    sep = buffer.indexOf(separator,sepstart+1);
	    numberOfMuEThSteps[ialpha] =
		Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
	    //System.err.println("nStep sepstart=" + sepstart + " sep="+ sep);

	    sepstart = sep;
	    muEThStepSize[ialpha]=
		Double.valueOf(buffer.substring(sepstart+1)).doubleValue( );
	    //System.err.println("muEstepsize sepstart=" + sepstart + " sep="+ sep);

	    //System.err.println(muEThMin[ialpha] + " " + numberOfMuEThSteps[ialpha] + 
	    //	       " " + muEThStepSize[ialpha]);

	    // muE threshold loop
	    for(int ith = 0; ith< numberOfMuEThSteps[ialpha]; ith++){

		// logEmu loop
		for(int iLogE=0;iLogE<numberOfEmuon;iLogE++){
		    buffer = in.readLine();

		    sep =  buffer.lastIndexOf(separator);
		    dFMuDLogEarray[ialpha][ith][distanceIndex][iLogE]=
			Double.valueOf(buffer.substring(sep)).doubleValue( );

		    //System.err.println(" " + iLogE + 
		    //		       " " + 
		    //	       dFMuDLogEarray[ialpha][ith][distanceIndex][iLogE]);

		}// logEmu loop ends
	    }// muE threshold loop ends

	}//alpha loop ends




    }


    /** Returns index of [distance] for a given distance */
    private int indexOfDistance(double distance){
	double logDistance = Math.log(distance)/ln10;

	int i;
	for(i=0;i<icedist.length-1;i++){
	    if((logDistance>=icedist[i] )&&(logDistance<icedist[i+1])) break;
	}

	//if(i<tableFileName.length)
	//    System.err.println(" logDistance=" + logDistance + 
	//		       " " + tableFileName[i]);
	return i;
    }

    /** 
	Returns dFMu/dLogE for given indexes of alpha, muEth, distance, logEmu.
	<pre>
	ialpha = (alpha - alphaMin)/alphaStepSize
	iMuETh = (muETh - muEThMin)/muEThStepSize
	distance : Propagation distance [cm]
	logEmu : log(Muon Bundle Energy[GeV])
	</pre>
    */
    public double getDFMuDLogE(int ialpha, int iMuETh, double distance,
			       double logEmu){

	int idistance = indexOfDistance(distance);
	int iLogE = (int )((logEmu-logEmuonMin + roundoff)/logEmuonStepSize + 0.01);

	if(0<= idistance && idistance<numberOfDistance){
	    if(0<=iLogE && iLogE < numberOfEmuon){

		return dFMuDLogEarray[ialpha][iMuETh][idistance][iLogE];
	    }
	}

	//System.err.println("### logEmuon " + logEmu + " or distance "
	//		   + distance + " out of range ###");
	return 0.0;
    }

    /** 
	Returns dFMu/dLogE for given alpha, muEth, distance, logEmu.
	<pre>
	alpha : The Elbert formula's parameter - used in AtmMuonBundleFlux class
	muETh : The threshold energy of muon   - used in AtmMuonBundleFlux class
	distance : Propagation distance [cm]
	logEmu : log(Muon Bundle Energy[GeV])
	</pre>
    */
    public double getDFMuDLogE(double alpha, double muETh, double distance,
			       double logEmu){

	int ialpha = (int )((alpha - alphaMin + roundoff)/alphaStepSize + 0.01);
	int iMuETh = (int )((muETh - muEThMin[ialpha] + roundoff)/muEThStepSize[ialpha] + 0.01);

	if(0<= ialpha && ialpha<numberOfAlphaSteps){
	    if(0<=iMuETh && iMuETh < numberOfMuEThSteps[ialpha]){
		return getDFMuDLogE(ialpha,iMuETh,distance,logEmu);
	    }
	}

	System.err.println("### alpha " + alpha + " or muETh "
			   + muETh + " out of range ###");
	return 0.0;
    }

    /**
       Set the bin numbers of alpha and the threshold energy of Muon.
       Once you set these bins, calling getDFMuDLogE(distance, logEmu)
       gives the bundle flux. The relation between bin number
       and the values are :
	<pre>
	alphabin = (alpha - alphaMin)/alphaStepSize
	muEThbin = (muETh - muEThMin)/muEThStepSize
	</pre>

	They are must be within [0, numberOfAlphaSteps], 
	[0, numberOfMuEThSteps[alphabin]], which is determined
	by the calculated range of the table data.

    */
    public void setElbertParameters(int alphabin, int muEThbin){

	if(0<= alphabin && alphabin<numberOfAlphaSteps){
	    if(0<=muEThbin && muEThbin < numberOfMuEThSteps[alphabin]){
		this.alphaBin = alphabin;
		this.muEThBin = muEThbin;
		elbertParametersHaveBeenSet=true;
	    }
	}
	if(!elbertParametersHaveBeenSet){
	    System.err.println("### alpha " + alphabin + " or muETh "
			   + muEThbin + " out of range ###");
	}
    }


    /**
       Set alpha and the threshold energy of Muon.
       Once you set these bins, calling getDFMuDLogE(distance, logEmu)
       gives the bundle flux.
    */
    public void setElbertParameters(double alpha, double muETh){

	int alphabin = (int )((alpha - alphaMin + roundoff)/alphaStepSize + 0.01);
	int muEThbin = (int )((muETh - muEThMin[alphabin] + roundoff)/muEThStepSize[alphabin] + 0.01);
	setElbertParameters(alphabin,muEThbin);
    }

    /**
       Set the bin numbers of alpha. Eth is set to be the lowest value i.e. 
       that in the first bin.
       Once you set these bins, calling getDFMuDLogE(distance, logEmu)
       gives the bundle flux. The relation between bin number
       and the values are :
	<pre>
	alphabin = (alpha - alphaMin)/alphaStepSize
	</pre>

	They are must be within [0, numberOfAlphaSteps], 
	which is determined by the calculated range of the table data.

    */
    public void setElbertParameters(int alphabin){

	setElbertParameters(alphabin,0);
    }


    /**
       Set alpha.Eth is set to be the lowest value i.e. that in the first bin.
       Once you set these bins, calling getDFMuDLogE(distance, logEmu)
       gives the bundle flux.
    */
    public void setElbertParameters(double alpha){

	int alphabin = (int )((alpha - alphaMin + roundoff)/alphaStepSize + 0.01);
	setElbertParameters(alphabin,0);
    }


    /** Return the alpha value set by setElbertParameters */
    public double getAlpha(){
	if(!elbertParametersHaveBeenSet){
	    System.err.println("Has not set the alpha value by setElbertParameters!");
	    System.exit(0);
	}
	return alphaMin + alphaStepSize*(double )alphaBin;
    }


    /** Return the muETh value set by setElbertParameters */
    public double getMuETh(){
	if(!elbertParametersHaveBeenSet){
	    System.err.println("Has not set the alpha value by setElbertParameters!");
	    System.exit(0);
	}
	return muEThMin[alphaBin] + muEThStepSize[alphaBin]*(double )muEThBin;
    }

    /** 
	Returns dFMu/dLogE for given indexes of distance, logEmu.
	The parameters in the Elbert Formula, alpha and muETh must be
	set in advance by the method setElbertParameters.
	<pre>
	distance : Propagation distance [cm]
	logEmu : log(Muon Bundle Energy[GeV])
	</pre>
    */
    public double getDFMuDLogE(double distance, double logEmu){

	if(!elbertParametersHaveBeenSet){
	    System.err.println("Has not set the alpha value by setElbertParameters!");
	    System.exit(0);
	}
	return getDFMuDLogE(alphaBin,muEThBin,distance,logEmu);
    }


    /** Simple main function for debugging/demonstrating this class */
    public static void main(String[] args) throws IOException {

	ElbertFluxTableFactory atmMuonFluxTable= 
	    new ElbertFluxTableFactory(tablePath+tableFileName[0],0);
	// Reads only "1e5_00cm.table" for log10(distnce) = 5.0;

	double distance = 1.05e5; // distance muon propagated [cm]
	double alpha = 2.02;
	double muETh = 2500.0; 
	for(int iLogE=0;iLogE<10;iLogE++){
	    double logEmu = Particle.getLogEnergyMinimum( ) + 
		50.0*Particle.getDeltaLogEnergy( )*(double )iLogE;
	    System.out.println("logEmuon=" + logEmu + " flux=" +
	    	       atmMuonFluxTable.getDFMuDLogE(alpha,muETh,distance,logEmu));
	}

	// Or different way to set Alpha and MuETh that are within the range
	// of the table.
	atmMuonFluxTable.setElbertParameters(10,12); // alphaBin muEThbin
	System.out.println("alpha=" + atmMuonFluxTable.getAlpha() + 
			   " muETh=" + atmMuonFluxTable.getMuETh());
	for(int iLogE=0;iLogE<10;iLogE++){
	    double logEmu = Particle.getLogEnergyMinimum( ) + 
		50.0*Particle.getDeltaLogEnergy( )*(double )iLogE;
	    System.out.println("logEmuon=" + logEmu + " flux=" +
	    	       atmMuonFluxTable.getDFMuDLogE(distance,logEmu));
	}

    }


}


