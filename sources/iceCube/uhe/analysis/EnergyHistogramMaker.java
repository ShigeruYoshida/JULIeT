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
    This class makes a 2D histogram of Energy correlations 
    like Muon In-ice Energy Vs Primary Cosmic Ray enegry 
    (obtained by AtmMuonBundleFlux.java in the MuonModel package)
    or Muon In-ice Energy Vs Muon Surface Energy
    (using PropagationMatrixFactory.java in the propahgation package)

    Written by S. Yoshida 2007 June 24th
*/

public class EnergyHistogramMaker {

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

    // parameter of the bundle flux
    protected double alpha = 2.04;
    protected double muEth = 3730.0;

    /** Default Constructor */
    public EnergyHistogramMaker(){
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



    /** Make 2D Histogram 
	<pre>
	string option : "inice-surface"  - inice muon Vs surface muon
                        "surface-primary" - surface muon Vs primary CR
                        "inice-primary"  - inice muon Vs primary CR
	</pre>
    */
    public IHistogram2D makeJaidaHistogram(String option) throws IOException{

	// Interactive session to set the inIceParticle flavor and doublets
        DataInputStream input = new DataInputStream(System.in); 
        BufferedReader  d     = new BufferedReader(new InputStreamReader(input)); 
        String buffer; 

	if(option.startsWith("inice-surf")){ // inice particle Vs surface particle
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

	    boolean isWeighted = false;
	    System.err.print("Flux Weighted? [yes(1)/no(0)]->"); 
	    buffer   = d.readLine(); 
	    if(Integer.valueOf(buffer).intValue()==1) 
		isWeighted= true;
	    if(!isWeighted) System.err.print(" no weighted"); 

	    Particle surfaceParticle = new Particle(flavor,doublet);
	    System.err.println("From " + 
			     surfaceParticle.particleName(surfaceParticle.getFlavor(),
							  surfaceParticle.getDoublet()) +
			     " to " +
			     inIceParticle.particleName(inIceParticle.getFlavor(),
							inIceParticle.getDoublet()));

	    IHistogram2D h2 = 
		build2DHistogramByPropMatrix(inIceParticle,surfaceParticle,isWeighted);
	    return h2;

	}else{

	    System.err.print("Elbert alpha ->"); // Elbert alpha parameter
	    buffer   = d.readLine(); 
	    alpha = Double.valueOf(buffer).doubleValue();

	    System.err.print("Muon E threshold [GeV]->"); // Elbert muE threshold
	    buffer   = d.readLine(); 
	    muEth = Double.valueOf(buffer).doubleValue();

	    boolean includeFluctuationEffects = false;
	    System.err.print("Include the fluctuation effect? [yes(1)/no(0)]->"); 
	    buffer   = d.readLine(); 
	    if(Integer.valueOf(buffer).intValue()==1) 
		includeFluctuationEffects = true;

	    boolean isWeighted = false;
	    System.err.print("Flux Weighted? [yes(1)/no(0)]->"); 
	    buffer   = d.readLine(); 
	    if(Integer.valueOf(buffer).intValue()==1) 
		isWeighted= true;
	    if(!isWeighted) System.err.print(" no weighted"); 

	    // Muon Bundle Flux model
	    AtmMuonBundleFlux atmFlux = new AtmMuonBundleFlux(alpha,1.0,muEth);
	    atmFlux.setFluxCalculationMode(false,false);

	    IHistogram2D h2 = 
		build2DHistogramByAtmMuonBundleFlux(option,atmFlux,
						    includeFluctuationEffects,
						    isWeighted);
	    return h2;
	}
    }

    /** Build Jaida 2D Histogram by propMatrixFactory */
    private IHistogram2D build2DHistogramByPropMatrix(Particle inIceParticle,
						      Particle surfaceParticle,
						      boolean isWeighted)
	throws IOException{

	// Propagation Matrix 
	PropagationMatrixFactory propMtx = new PropagationMatrixFactory();
	IHistogram2D h2 = 
	    jaidaHistoFactory.createHistogram2D("IceCube depth Vs Earth surface",
						dimensionLogE,minLogE,maxLogE,
						dimensionLogE,minLogE,maxLogE);

	AtmMuonBundleFlux atmFlux = null;
	double muonThresholdEnergy = 3750.0;
	if(isWeighted){
	    // Muon Bundle Flux model
	    atmFlux = new AtmMuonBundleFlux();
	    atmFlux.setFluxCalculationMode(true,false);
	    muonThresholdEnergy = atmFlux.getMuonThresholdEnergy();
	}

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

	    if(isWeighted){
		double logReferenceEnergy = 7.0; // 10^7 GeV
		double thresholdE = (muonThresholdEnergy+atmFlux.criticalEnergy)/
		    propMtx.getAverageMuonEnergyLossAfterPropagation(logReferenceEnergy)-atmFlux.criticalEnergy;
		atmFlux.setMuonThresholdEnergy(thresholdE);
		System.err.println(" JULIET: Muon E threshold = " + 
				   atmFlux.getMuonThresholdEnergy());
	    }


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
		double muonBundleFlux = 1.0;
		if(isWeighted){
		    muonBundleFlux = atmFlux.getDFDLogE(logSurfaceEnergy,cosZenith)*1.0e15;
		    
		}

		// InIce Lepton Energy loop
		for(int jLogE= 0;jLogE<Particle.getDimensionOfLogEnergyMatrix();
		    jLogE ++) {

		    double logInIceEnergy = Particle.getLogEnergyMinimum()
			+ Particle.getDeltaLogEnergy()*(double )jLogE;
		    inIceParticle.putLogEnergy(logInIceEnergy);

		    double flux = 
			propMtx.getDF(surfaceParticle,inIceParticle)*muonBundleFlux;

		    h2.fill(logSurfaceEnergy,logInIceEnergy,flux*solidAngle);
		}
	    }
	}

	return h2;
    }


    /** Build Jaida 2D Histogram by AtmMuonBundleFlux */
    private IHistogram2D build2DHistogramByAtmMuonBundleFlux(String option,
					     AtmMuonBundleFlux atmFlux,
					     boolean includeFluctuationEffects,
					     boolean isWeighted){
	IHistogram2D h2 = 
	    jaidaHistoFactory.createHistogram2D(option,
						dimensionLogE,minLogE,maxLogE,
						dimensionLogE,minLogE,maxLogE);

	CosmicRayFlux crFlux = new CosmicRayFlux();
	crFlux.setCutOffFeature(true);

	//Zenith Angle loop
	double cosZenith = 1.0; double deltaCos = 0.01;
	double solidAngle =  2.0*Math.PI*deltaCos;
	double detectorDepth = 1.4e5;  // Detector Depth.. 1400m = 1.4e5 cm below sea level
	ParticlePoint s = new ParticlePoint(0.0,0.0,0); // ice
	while(cosZenith>=0.17){
	    // InIce Muon Energy loop
	    for(int jLogE= 0;jLogE<Particle.getDimensionOfLogEnergyMatrix();
		jLogE ++) {

		double logInIceEnergy = Particle.getLogEnergyMinimum()
		    + Particle.getDeltaLogEnergy()*(double )jLogE;
		double inIceEnergy = Math.pow(10.0,logInIceEnergy);

		double beta_loss = 4.4619776009127244E-6;
		if(logInIceEnergy >= 7.0) beta_loss = CELbeta.getBeta(logInIceEnergy);

		double cosmicRayEnergy;
		double atmMuonDFDLogE;

		if(option.startsWith("surface-pri")){ // surface muon Vs primary CR
		    cosmicRayEnergy = 
			atmFlux.getEffectiveEnergyOfCRs(inIceEnergy,cosZenith);
		    atmMuonDFDLogE = atmFlux.getDFDLogE(inIceEnergy,cosZenith);
		}else{ // inice muon Vs primary CR

		    double sq_term = 
			Math.sqrt((ParticlePoint.REarth-detectorDepth)*
				  (ParticlePoint.REarth-detectorDepth)
				  *cosZenith*cosZenith + 
				  2.0*ParticlePoint.REarth*detectorDepth-
				  detectorDepth*detectorDepth);
		    double trajectoryLength = 
			sq_term - (ParticlePoint.REarth-detectorDepth)*cosZenith;

		    double slantDepth = trajectoryLength*s.getMediumDensity();

		    cosmicRayEnergy = 
			atmFlux.getEffectiveEnergyOfCRs(logInIceEnergy,cosZenith,
						       beta_loss,slantDepth);
		    atmMuonDFDLogE = atmFlux.getDFDLogE(inIceEnergy,cosZenith,
						       beta_loss,slantDepth);
		}

		double logCREnergy = Math.log(cosmicRayEnergy)/ln10;
		double flux = solidAngle;
		if(isWeighted) flux = atmMuonDFDLogE*1.0e10*1.0e2;

		if(!includeFluctuationEffects){
		    h2.fill(logCREnergy,logInIceEnergy,flux);
		}else{
		    boolean asInIce = true;
		    if(option.startsWith("surface-pri")) asInIce = false;
		    double logR = atmFlux.cascade.getLogMuOverCREnergyMin(logCREnergy,asInIce);
		    while(logR<=atmFlux.cascade.getLogMuOverCREnergyMax(logCREnergy,asInIce)){
			double prob = 
			    atmFlux.cascade.getProbability(logR,logCREnergy,asInIce);
			double logEnergy = logCREnergy-logR;
			if(isWeighted) 
			    flux = (1.0/(alpha-1.0))*crFlux.getDFDLogE(logCREnergy - logR)*1.0e10*1.0e2;;
			h2.fill(logEnergy,logInIceEnergy,
				flux*prob*Particle.getDeltaLogEnergy());
			logR += Particle.getDeltaLogEnergy();
		    }
		}


	    } // inIce muon loop ends

	    cosZenith -= deltaCos;

	} // cos(Zenith) loop ends

	return h2;
    }


    /** Main method -- Draw the 2D histogram to show the colleration of energies */
    public static void main(String[] args) throws IOException{

	if(args.length<1){
            System.out.println("Usage: EnergyHistogramMaker  option (f2kout)");
	    System.out.println("inice-surface  - inice muon Vs surface muon");
            System.out.println("surface-primary - surface muon Vs primary CR");
            System.out.println("inice-primary  - inice muon Vs primary CR");
            System.exit(0);
	}

	EnergyHistogramMaker histMaker = new EnergyHistogramMaker();

        IHistogram2D h2_energy = histMaker.makeJaidaHistogram(args[0]);

	//
	// plotting 
	// 
	IPlotterFactory plotterFactory = histMaker.jaidaFactory.createPlotterFactory(); 
	IPlotter plotter = plotterFactory.create("Energy Correlation"); 
 
	IPlotterStyle distStyle = plotter.region(0).style(); 

	if(args[0].startsWith("inice-surf")){ // inice particle Vs surface particle
	    JaidaPlotStyleSetter.setPlotStyle(distStyle, 
				  "Log(Muon Energy at Surface [GeV])",
				  "Log(Muon Energy at IceCube depth[GeV])",
				  "hist2DStyle","colorMap");
	}else{
	    JaidaPlotStyleSetter.setPlotStyle(distStyle, 
					      "Log(Primary CR Energy [GeV])",
					      "Log(Muon Bundle Energy [GeV])",
					      "hist2DStyle","colorMap");
	}
	plotter.region(0).plot(h2_energy); 
 
	plotter.show(); 

	if(args.length>1){ // f2k out
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
	}
 

    }
}
