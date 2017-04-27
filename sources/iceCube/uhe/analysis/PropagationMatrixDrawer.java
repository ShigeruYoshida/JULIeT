package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.propagation.*;
import iceCube.uhe.muonModel.*;
import iceCube.uhe.analysis.*;

import hep.aida.*; 

import java.io.*;
import java.util.*;

/** 
    This class makes a 2D/1D histogram of Energy correlations 
    of In-ice Energy Vs Surface Energy
    (using PropagationMatrixFactory.java in the propahgation package)

    Written by S. Yoshida 2007 June 24th
*/

public class PropagationMatrixDrawer {

    // FreeHep objects
    protected IAnalysisFactory jaidaFactory;
    private ITree jaidaTree;
    private IHistogramFactory jaidaHistoFactory;
    private IDataPointSetFactory jaidaDpsFactory;

    // Histogram Parameters
    protected double minLogE = Particle.getLogEnergyMinimum(); // min JULIeT energy
    protected double maxLogE = Particle.getLogEnergyMinimum() +
	Particle.getDeltaLogEnergy()*
	(double )(Particle.getDimensionOfLogEnergyMatrix()); // max JULIeT energy
    protected int dimensionLogE = Particle.getDimensionOfLogEnergyMatrix()/10;

    private final static double ln10 = Math.log(10.0); 

    /** Default Constructor */
    public PropagationMatrixDrawer(){
	// Jaida FreeHep objects
        jaidaFactory = IAnalysisFactory.create();
        jaidaTree = jaidaFactory.createTreeFactory().create();
        jaidaHistoFactory = 
            jaidaFactory.createHistogramFactory(jaidaTree);
        jaidaDpsFactory = 
            jaidaFactory.createDataPointSetFactory(jaidaTree);

	jaidaTree.mkdir("/energy2D"); 
	jaidaTree.cd("/energy2D"); 

    }

    /** Build Jaida 2D Histogram by propMatrixFactory */
    private IHistogram2D build2DHistogramByPropMatrix(Particle inIceParticle,
						      Particle surfaceParticle) 
	throws IOException{

	// Propagation Matrix 
	PropagationMatrixFactory propMtx = new PropagationMatrixFactory();
	IHistogram2D h2 = 
	    jaidaHistoFactory.createHistogram2D("IceCube depth Vs Earth surface",
						dimensionLogE,minLogE,maxLogE,
						dimensionLogE,minLogE,maxLogE);


	for(int itheta=0;itheta<I3ParticleWeightFiller.matrixFileName[2].length;
	    itheta++){
	    // Read the serialized object of the Neutrino Charged Interaction Matrix
	    String fileName = 
		I3ParticleWeightFiller.pathname[0] + 
		I3ParticleWeightFiller.matrixFileName[2][itheta];
	    DataInputStream in = new DataInputStream(new FileInputStream(fileName));
	    propMtx.readMatrix(in);
	    in.close( );
	    System.err.println("Reading the matrix from " + fileName + " done.");

	    // Solid angle calculation
	    double radiansUp = 0.0; double radiansDown = 0.0;
	    double cosZenith = 0.0;
	    radiansUp = Math.toRadians(I3ParticleWeightFiller.icerange[itheta]);
	    radiansDown = 
		Math.toRadians(I3ParticleWeightFiller.icerange[itheta+1]);
	    cosZenith = Math.cos(radiansUp);
	    double solidAngle = 
	    2.0*Math.PI*Math.abs(Math.cos(radiansDown)-Math.cos(radiansUp));
	    System.err.println(" Solid angle " + solidAngle);

	    // Earth-surface Neutrino Energy loop
	    for(int iLogE = 0;iLogE<Particle.getDimensionOfLogEnergyMatrix();
		iLogE++){

		double logSurfaceEnergy = Particle.getLogEnergyMinimum()
		    + Particle.getDeltaLogEnergy()*(double )iLogE;

		surfaceParticle.putLogEnergy(logSurfaceEnergy);

		// InIce Lepton Energy loop
		for(int jLogE= 0;jLogE<Particle.getDimensionOfLogEnergyMatrix();
		    jLogE ++) {

		    double logInIceEnergy = Particle.getLogEnergyMinimum()
			+ Particle.getDeltaLogEnergy()*(double )jLogE;
		    inIceParticle.putLogEnergy(logInIceEnergy);

		    double flux = 
			propMtx.getDF(surfaceParticle,inIceParticle);

		    h2.fill(logSurfaceEnergy,logInIceEnergy,flux*solidAngle);
		}
	    }
	}

	return h2;
    }

    /** Build Jaida 1D Histogram by propMatrixFactory */
    private IHistogram1D build1DHistogramByPropMatrix(Particle inIceParticle,
						      Particle surfaceParticle)

	throws IOException{

	// Propagation Matrix 
	PropagationMatrixFactory propMtx = new PropagationMatrixFactory();
	IHistogram1D h1 = 
	    jaidaHistoFactory.createHistogram1D("Energy Distribution",
						dimensionLogE,minLogE,maxLogE);



	for(int itheta=0;itheta<I3ParticleWeightFiller.matrixFileName[2].length;
	    itheta++){
	    // Read the serialized object of the Neutrino Charged Interaction Matrix
	    String fileName = 
		I3ParticleWeightFiller.pathname[0] + 
		I3ParticleWeightFiller.matrixFileName[2][itheta];
	    DataInputStream in = new DataInputStream(new FileInputStream(fileName));
	    propMtx.readMatrix(in);
	    in.close( );
	    System.err.println("Reading the matrix from " + fileName + " done.");

	    // Solid angle calculation
	    double radiansUp = 0.0; double radiansDown = 0.0;
	    double cosZenith = 0.0;
	    radiansUp = Math.toRadians(I3ParticleWeightFiller.icerange[itheta]);
	    radiansDown = 
		Math.toRadians(I3ParticleWeightFiller.icerange[itheta+1]);
	    cosZenith = Math.cos(radiansUp);
	    double solidAngle = 
	    2.0*Math.PI*Math.abs(Math.cos(radiansDown)-Math.cos(radiansUp));
	    System.err.println(" Solid angle " + solidAngle);

	    // Earth-surface Neutrino Energy loop
	    for(int iLogE = 0;iLogE<Particle.getDimensionOfLogEnergyMatrix();
		iLogE++){

		double logSurfaceEnergy = Particle.getLogEnergyMinimum()
		    + Particle.getDeltaLogEnergy()*(double )iLogE;

		surfaceParticle.putLogEnergy(logSurfaceEnergy);

		inIceParticle.putLogEnergy(logSurfaceEnergy);

		double flux = 
		    propMtx.getDF(surfaceParticle,inIceParticle);

		h1.fill(logSurfaceEnergy,flux*solidAngle);
	    }
	}

	return h1;
    }


    /** Main method -- Draw the 2D histogram to show the colleration of energies */
    public static void main(String[] args) throws IOException{

	PropagationMatrixDrawer histMaker = new PropagationMatrixDrawer();

	// Interactive session to set the inIceParticle flavor and doublets
        DataInputStream input = new DataInputStream(System.in); 
        BufferedReader  d     = new BufferedReader(new InputStreamReader(input)); 
        String buffer; 

	System.err.print("inIceParticle flavor ->"); 
	buffer   = d.readLine(); 
	int flavor = Integer.valueOf(buffer).intValue();

	System.err.print("inIceParticle doublet ->"); 
	buffer   = d.readLine(); 
	int doublet = Integer.valueOf(buffer).intValue();
	Particle inIceParticle = new Particle(flavor,doublet);

	System.err.print("surfaceParticle flavor ->"); 
	buffer   = d.readLine(); 
	flavor = Integer.valueOf(buffer).intValue();

	System.err.print("surfaceParticle doublet ->"); 
	buffer   = d.readLine(); 
	doublet = Integer.valueOf(buffer).intValue();

	Particle surfaceParticle = new Particle(flavor,doublet);
	System.err.println("From " + 
			   surfaceParticle.particleName(surfaceParticle.getFlavor(),
						surfaceParticle.getDoublet()) +
			   " to " +
			   inIceParticle.particleName(inIceParticle.getFlavor(),
						inIceParticle.getDoublet()));

	System.err.print("2D Histogram?  [yes(1)/no(0)] ->"); 
	buffer   = d.readLine(); 
	boolean draw2D = false;
	IHistogram2D h2_energy = null;
	IHistogram1D h1_energy = null;
	if(Integer.valueOf(buffer).intValue()==1){
	    System.err.println(" Drawing 2D Histogram..");
	    h2_energy =	histMaker.build2DHistogramByPropMatrix(inIceParticle,surfaceParticle);
	    draw2D = true;
	}else{
	    System.err.println(" Drawing 1D Histogram..");
	    h1_energy =	histMaker.build1DHistogramByPropMatrix(inIceParticle,surfaceParticle);
	}

	//
	// plotting 
	// 
	IPlotterFactory plotterFactory = 
	    histMaker.jaidaFactory.createPlotterFactory(); 
	IPlotter plotter = plotterFactory.create("Energy Correlation"); 
 
	IPlotterStyle distStyle = plotter.region(0).style(); 

	if(draw2D){ // drawing the 2D Histogram
	    JaidaPlotStyleSetter.setPlotStyle(distStyle, 
				  "Log(Energy at Surface [GeV])",
				  "Log(Energy at IceCube depth[GeV])",
				  "hist2DStyle","colorMap");

	    plotter.region(0).plot(h2_energy); 
	}else{ // 1D Historgram
	    JaidaPlotStyleSetter.setPlotStyle(distStyle, 
			  "Log(Energy [GeV])",
  		          "Probability");
	    plotter.region(0).plot(h1_energy); 
	}

 
	plotter.show(); 

	if(args.length>1){ // f2k out
	    if(draw2D){
		int nbinx = h2_energy.xAxis().bins();
		int nbiny = h2_energy.yAxis().bins();
		System.out.println(nbinx + " " + histMaker.minLogE + " " + 
				   histMaker.maxLogE);
		System.out.println(nbiny + " " + histMaker.minLogE + " " +
				   histMaker.maxLogE);
		for(int iy = 0; iy< nbiny ; iy++){
		    for(int ix = 0; ix< nbinx; ix++){
			System.out.println(h2_energy.binHeight(ix,iy));
		    }
		}
	    }else{
		int nbinx = h1_energy.axis().bins();
		System.out.println(nbinx + " " + histMaker.minLogE + " " + 
				   histMaker.maxLogE);
		for(int ix = 0; ix< nbinx ; ix++){
		    System.out.println(h1_energy.binHeight(ix));
		}
	    }
	}
 

    }
}
