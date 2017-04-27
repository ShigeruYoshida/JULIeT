package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import numRecipes.*;

import hep.aida.*;
import hep.aida.ext.*;
import hep.aida.util.*;
import hep.aida.util.comparison.*;

import java.io.*;
import java.util.*;

/** 
    Optimize NPE threshold to minimize either the model rejection potential
    or the signal discovery potential. Rely on numRecipes.FeldmanCousins.java

    Written originally by S. Yoshida for the IceCube EHE analysis.
    2010/3/29
*/
public class NPEthresholdOptimizer {

    // Jaida FreeHep objects
    IAnalysisFactory jaidaFactory = null;
    ITree jaidaTree = null;
    IHistogramFactory jaidaHistoFactory = null;

    /** dimension of the parameter y. usually cosZenith, but depends on analysis */
    protected int dimensionY = 20;
    /** maximum range of the parameter y (ex. cos(Zenith)) */ 
    protected double yMax = 1.0;
    /** mininum range of the parameter y (ex. cos(Zenith)) */ 
    protected double yMin = -1.0;


    /** dimension of the logNPE. bin */
    protected int dimensionNPE = 300;
    /** maximum range of logNPE */ 
    protected double logNPEMax = 4.0; // 10^4 < NPE
    /** mininum range of logNPE */ 
    protected double logNPEMin = 7.0; // 10^7 < NPE

    /** histogram to store number of signal events */
    protected IHistogram2D signalRate2D = null;
    protected IHistogram1D signalRate1D = null;
    /** histogram to store number of background events */
    protected IHistogram2D bgRate2D = null;
    protected IHistogram1D bgRate1D = null;

    /** Histogram to store the model rejection factor */
    protected IHistogram1D modelRejectHisto = null;
    /** Histogram to store the signal discovery potential */
    protected IHistogram1D signalDiscoverPotentialHisto = null;

    /** name of the output Aida Tree file */
    protected String aidaFileName = "NPEoptimizer.aida";

    private boolean isSignalFilled = false;
    private boolean isBGFilled = false;
    private double epsilon = 1.0e-8;

    /** Constructor. Creaate the Jaida Histogram [logNPE][y] to store
        number of background and signal events. 
	<pre>
	int dimensionY   dimension of the parameter y (ex. cos(zenith) in histograms
	double yMin      mininum rang of the parameter y (ex. cos(zenith) in histograms
	double yMax      maximum rang of the parameter y (ex. cos(zenith) in histograms
	int dimensionNPE dimension of log(NPE)  in histograms
	double logNPEMin    mininum rang of log(NPE) in histograms
	double logNPEMax    maximum rang of log(NPE) in histograms
	</pre>
    */
    public NPEthresholdOptimizer(int dimensionY, double yMin, double yMax,
				 int dimensionNPE, double logNPEMin,double logNPEMax) 
	throws IOException{

	this.dimensionY = dimensionY;
	this.yMin = yMin;
	this.yMax = yMax;
	this.dimensionNPE = dimensionNPE;
	this.logNPEMin = logNPEMin;
	this.logNPEMax = logNPEMax;

	createHistograms();
    }

    /** Constructor. Creaate the Jaida Histogram [logNPE][y] to store
        number of background and signal events. Range of bins and dimensions 
	are set by default values */
    public NPEthresholdOptimizer() throws IOException{
	createHistograms();
    }

    private void createHistograms() throws IOException{

        // Jaida FreeHep objects
        jaidaFactory = IAnalysisFactory.create();
        jaidaTree = jaidaFactory.createTreeFactory().
            createTree(aidaFileName,"xml",ITreeFactory.RECREATE);
	jaidaHistoFactory = 
	    jaidaFactory.createHistogramFactory(jaidaTree);

	signalRate2D = 
	    jaidaHistoFactory.createHistogram2D("Signal Rate",dimensionNPE,logNPEMin,logNPEMax,
						dimensionY,yMin,yMax);
	bgRate2D = 
	    jaidaHistoFactory.createHistogram2D("Background Rate",dimensionNPE,logNPEMin,logNPEMax,
						dimensionY,yMin,yMax);

	signalRate1D = 
	    jaidaHistoFactory.createHistogram1D("Signal Rate 1D",
						dimensionNPE,logNPEMin,logNPEMax);

	bgRate1D = 
	    jaidaHistoFactory.createHistogram1D("Background Rate 1D",
						dimensionNPE,logNPEMin,logNPEMax);

	modelRejectHisto = 
	    jaidaHistoFactory.createHistogram1D("Model Rejection Factor",
						dimensionNPE,logNPEMin,logNPEMax);

	signalDiscoverPotentialHisto = 
	    jaidaHistoFactory.createHistogram1D("Signal Discovery Potential",
						dimensionNPE,logNPEMin,logNPEMax);
    }

    /**
       Reading the number of signal and background events binned in logNPE and y(such as cos(Zenith)
    */
    public void readAndFillEventRate(DataInputStream in) throws IOException{

	double y = 0.0;

	// Reading data
	BufferedReader  d = new BufferedReader(new InputStreamReader(in));
	String buffer; int sep = 0; int sepstart = 0;
        char separator = ' ';
	while((buffer = d.readLine())!=null){
	    // line -- logNpe, signal, bg
	    sepstart = 0;

	    sep = buffer.indexOf(separator,sepstart+1);
	    double logNPE =
		Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
	    sepstart = sep;

	    sep = buffer.indexOf(separator,sepstart+1);
	    double sigRate =
		Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
	    sepstart = sep;

	    sep = buffer.indexOf(separator,sepstart+1);
	    double bgRate =
		Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
	    sepstart = sep;

	    signalRate2D.fill(logNPE-epsilon,y-epsilon,sigRate);
	    bgRate2D.fill(logNPE-epsilon,y-epsilon,bgRate);
	    //System.err.format(" logNPE=%5.2f sig(%f) BG(%f)\n",logNPE,sigRate,bgRate);
	}

	isSignalFilled = true; isBGFilled = true;

    }

    /**
       Reading the number of signal and background events binned in logNPE and y(such as cos(Zenith).
       Data format is different from the method above.
    */
    public void readAndFillEventRate(DataInputStream in, boolean wideFormat) throws IOException{

	if(!wideFormat){
	    readAndFillEventRate(in);
	    return;
	}

	double y = 0.0;

	// Reading data
	BufferedReader  d = new BufferedReader(new InputStreamReader(in));
	String buffer; int sep = 0; int sepstart = 0;
        char separator = ' ';
	while((buffer = d.readLine())!=null){
	    // line -- logNpe bg: rate1, rate2,..., rateN sig: rate1,....,rateN
	    sepstart = 0;

	    sep = buffer.indexOf(separator,sepstart+1);
	    double logNPE =
		Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
	    sepstart = sep;

	    sep = buffer.indexOf(separator,sepstart+1);
	    String bg = buffer.substring(sepstart+1,sep);
	    sepstart = sep;

	    for(int i=0;i<dimensionY;i++){ // loop over y-bin
		y = (double )i;
		sep = buffer.indexOf(separator,sepstart+1);
		double bgRate =
		    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		bgRate2D.fill(logNPE-epsilon,y+epsilon,bgRate);
		//System.err.format(" logNPE=%5.3f y=%5.2f BG(%f)\n",logNPE,y,bgRate);
	    }

	    sep = buffer.indexOf(separator,sepstart+1);
	    String sig = buffer.substring(sepstart+1,sep);
	    sepstart = sep;

	    for(int i=0;i<dimensionY;i++){ // loop over y-bin
		y = (double )i;
		sep = buffer.indexOf(separator,sepstart+1);
		double sigRate =
		Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		signalRate2D.fill(logNPE-epsilon,y+epsilon,sigRate);
		//System.err.format(" logNPE=%5.3f y=%5.2f sig(%f)\n",logNPE,y,sigRate);
	    }

	}

	isSignalFilled = true; isBGFilled = true;

    }

    /** Calculate the model rejection factor for each NPE threshold
	and fill the 1D histogram(ogNPE) with the factors.
	You have to fill the background and signal number historam first
	by FillSignalRate() and FillBGRate()
	<pre>
	double y  : presumably cosZenith, a parameter on y-axis. logNPE is on x-axis.
	boolean isIntegral  : if the filled BG/Signal rate is integral rate, set true.
	</pre>
    */
    public void calculateAndFillModelRejectionFactor(double y, boolean isIntegral){
	if(!isSignalFilled || !isBGFilled){
	    System.err.println(" You have to fill the signal/background rate histogram");
	    System.exit(0);

	}else{
	    modelRejectHisto.reset();  // initialize
	    signalRate1D.reset();
	    bgRate1D.reset();

	    // NPE loop
	    int iy = signalRate2D.coordToIndexY(y);
	    int nbin = signalRate2D.xAxis().bins();
	    for(int iNPE = 0; iNPE < nbin; iNPE++){

		// Extract number of signal and background events
		double nSignal = 0.0;
		double nBG = 0.0;
		if(isIntegral){
		    nSignal = signalRate2D.binHeight(iNPE,iy);
		    nBG = bgRate2D.binHeight(iNPE,iy);
		}else{
		    for(int ix = iNPE; ix< nbin; ix++){
			nSignal += signalRate2D.binHeight(ix,iy);
			nBG += bgRate2D.binHeight(ix,iy);
		    }
		}

		// average Feldman-Cousins upper limit
		double mu = FeldmanCousins.getAverageUpperLimit(nBG);
		if(nSignal<=1.0e-12) nSignal = 1.0e-12;
		double modelRejectionFactor = mu/nSignal;

		// fill the model rejection factor
		double logNPE =  signalRate2D.xAxis().binCenter(iNPE);
		modelRejectHisto.fill(logNPE,modelRejectionFactor);
		signalRate1D.fill(logNPE,nSignal);
		bgRate1D.fill(logNPE,nBG);
		System.err.format(" log(NPE)=%5.2f Sig(%f) BG(%f) MRF(%f)\n",
				  logNPE,nSignal,nBG,modelRejectionFactor);

	    }// NPE loop ends
	}

    }

    /** Calculate the signal discovery potential for each NPE threshold
	and fill the 1D histogram(ogNPE) with the factors.
	You have to fill the background and signal number historam first
	by FillSignalRate() and FillBGRate()
	<pre>
	double y  : presumably cosZenith, a parameter on y-axis. logNPE is on x-axis.
	double significanceOfDiscovery : significance to define "discovery" of signals. 5.0 e.x., 5sigma
	boolean isIntegral  : if the filled BG/Signal rate is integral rate, set true.
	</pre>
    */
    public void calculateAndFillSignalDiscoveryPotential(double y, double significanceOfDiscovery,
							 boolean isIntegral) throws IOException{
	if(!isSignalFilled || !isBGFilled){
	    System.err.println(" You have to fill the signal/background rate histogram");
	    System.exit(0);

	}else{
	    signalDiscoverPotentialHisto.reset();  // initialize

	    // NPE loop
	    int iy = signalRate2D.coordToIndexY(y);
	    int nbin = signalRate2D.xAxis().bins();
	    for(int iNPE = 0; iNPE < nbin; iNPE++){

		// Extract number of signal and background events
		double nSignal = 0.0;
		double nBG = 0.0;
		if(isIntegral){
		    nSignal = signalRate2D.binHeight(iNPE,iy);
		    nBG = bgRate2D.binHeight(iNPE,iy);
		}else{
		    for(int ix = iNPE; ix< nbin; ix++){
			nSignal += signalRate2D.binHeight(ix,iy);
			nBG += bgRate2D.binHeight(ix,iy);
		    }
		}

		// least signal for discovery - 5sigma significance
		double mu = FeldmanCousins.getLeastSignalForDiscovery(nBG,significanceOfDiscovery);
		if(nSignal<=1.0e-12) nSignal = 1.0e-12;
		double signalDiscoveryPotential = mu/nSignal;

		// fill the model rejection factor
		double logNPE =  signalRate2D.xAxis().binCenter(iNPE);
		signalDiscoverPotentialHisto.fill(logNPE,signalDiscoveryPotential);
		System.err.format(" log(NPE)=%5.2f Sig(%f) BG(%f) SDP(%f)\n",
				  logNPE,nSignal,nBG,signalDiscoveryPotential);

	    }// NPE loop ends
	}
        jaidaTree.commit();

    }


    /** Plot the model rejection factor and the signal discovery poteintial
     as a function of logNPE */
    public void plot(){
        IPlotterFactory plotterFactory = jaidaFactory.createPlotterFactory();
        IPlotter plotter = plotterFactory.create("MRF and SDP distribution");
        plotter.destroyRegions();
        plotter.createRegion(0,0,0.66,1);
        plotter.createRegions(2,1);

	IPlotterStyle mrfStyle = plotter.region(0).style();
        JaidaPlotStyleSetter.setPlotStyle(mrfStyle,
                                          "Log(Npe)","MRF(SDP)");
        plotter.region(0).plot(modelRejectHisto);
        plotter.region(0).plot(signalDiscoverPotentialHisto);

	IPlotterStyle neventStyle = plotter.region(1).style();
        JaidaPlotStyleSetter.setPlotStyle(neventStyle,
                                          "Log(Npe)","Number of Events");
        plotter.region(1).plot(signalRate1D);
        plotter.region(1).plot(bgRate1D);

	plotter.show();

	int nbin = modelRejectHisto.axis().bins();
	double mrfMin = 1.0e21; double npeMRF = modelRejectHisto.axis().binCenter(0);
	double sdpMin = 1.0e21; double npeSDP = signalDiscoverPotentialHisto.axis().binCenter(0);
	for(int iNPE = 0; iNPE < nbin; iNPE++){
	    double mrf = modelRejectHisto.binHeight(iNPE);
	    double sdp = signalDiscoverPotentialHisto.binHeight(iNPE);
	    double logNPE =  modelRejectHisto.axis().binCenter(iNPE);
	    if(mrf<mrfMin){
		mrfMin = mrf;
		npeMRF = logNPE;
	    }
	    if(sdp<sdpMin){
		sdpMin = sdp;
		npeSDP = logNPE;
	    }
	}
	System.err.format(" log(NPE)=%5.2f SDP(%f)\n",npeSDP,sdpMin);
	System.err.format(" log(NPE)=%5.2f MRF(%f)\n",npeMRF,mrfMin);
    }


    /** A simple main method for test */
    public static void main(String[] args) throws IOException{
 
	String eventRateFileName = null;
	double logNPEmin = 3.5;
	double logNPEmax = 7.0;
	double y = 0.0;
	double yMin = 0.0;
	double yMax = 7.0;

        if(args.length!=5){
            System.out.println("Usage: NPEthresholdOptimizer event-ratefile-name y yMax logNPEmin logNPEmax");
            System.exit(0);
        }else{
	    eventRateFileName = args[0];
	    y = Double.valueOf(args[1]).doubleValue();
	    yMax = Double.valueOf(args[2]).doubleValue();
	    logNPEmin = Double.valueOf(args[3]).doubleValue();
	    logNPEmax = Double.valueOf(args[4]).doubleValue();
        }

	int nYbin = (int)yMax;
	int nLogNPEbin = (int )((logNPEmax-logNPEmin)/0.1 + 0.001);

	//NPEthresholdOptimizer optimizer = new NPEthresholdOptimizer(5,0.0,5.0,25,4.0,6.5);
	NPEthresholdOptimizer optimizer = new NPEthresholdOptimizer(nYbin,yMin,yMax,nLogNPEbin,logNPEmin,logNPEmax);

        DataInputStream in = new DataInputStream (new FileInputStream(eventRateFileName));
	optimizer.readAndFillEventRate(in,true);
	optimizer.calculateAndFillModelRejectionFactor(y,true);
	optimizer.calculateAndFillSignalDiscoveryPotential(y,5.0,true); // 5sigma
	optimizer.plot();
    }
}
