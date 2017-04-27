package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.geometry.*;
import geometry.*;

import java.io.*;
import java.util.*;

/** 
    Build IceCube event(s) reading from the f2k data via standard input.
    By this way, You can leave the IceTray/i3 file framwork into
   the pure java world.

   The built I3Particle objects are outputed via outputStream.

   Orignally written by S. Yoshida for the IceCube EHE analysis
*/

public class I3ParticleBuilder {

    /** flag to read out MC events with MCTruth or real events */
    private boolean isMCTruth;
    private boolean fillMCWeight;
    private final double ln10 = Math.log(10.0);
    private double mcWeightOfThisData;
    private double log_EnergyMinimum;   // log[GeV]
    private double log_EnergyMaximum;   // log[GeV]
    private double powerLaw;

    private boolean isNeutrino = false;  // Case of neutrio primary
    private boolean isCorsika = false;   // Case of Corsika mmc primary
    private boolean isNuE   = false;  // if case of nu-e for the Glashow Resonane
    private boolean hasGR = false;
    private InteractionsMatrix nuCCmtx = null;   // CC cross section for weighting
    private InteractionsMatrix nuNCmtx = null;  // NC cross section for weighting
    private GlashowResonanceLeptonic grLepton = null;// GR lepton cross section for weighting
    private GlashowResonanceHadronic grHadron = null;// GR Hadron cross section for weighting
    private double enhancementFactor = 1.0;// cross section enhancement factor
    private double trackLength = 0.0;      // neutrino track length for weighting
    ParticlePoint s = null;


    /** Constructor. You are forced to choose MC or real data mode. */
    public I3ParticleBuilder(boolean isMCTruth){
	this.isMCTruth = isMCTruth;
	fillMCWeight = false;
	mcWeightOfThisData = 1.0;
	log_EnergyMinimum = 5.0;   // 10^5 GeV
	log_EnergyMaximum = 11.0;  // 10^11 GeV
	powerLaw = 1.0;
    }

    /**
       Set range of MC primary spectrum dN/dE 
       <pre>
       double logEnergyMinimum :  log(Energy Minimum [GeV]) in the spectral range
       double logEnergyMaximum :  log(Energy Maximum [GeV]) in the spectral range
       </pre>
    */
    public void setRangeOfMCSpectrum(double logEnergyMinimum, 
				     double logEnergyMaximum){
	log_EnergyMinimum = logEnergyMinimum;
	log_EnergyMaximum = logEnergyMaximum;
    }

    /** decide if you fills the MC spectrum weight with
        the dN/dLogE of the MC data. If true, then the method
        process() calls I3ParticleAnalysisFactor.getDNDLogE()
	and fills the MC weight with the retuned value 
	<pre>
	double powerLaw    :   dN/dE = E**(-powerLaw)
	</pre>
    */
    public void fillMCPrimarySpectrumWeight(double powerLaw){
	if(isMCTruth){
	    this.powerLaw = powerLaw;
	    fillMCWeight = true;
	}else{
	    System.err.println("You try to fill MC weight to the real data!");
	    System.exit(0);
	}
    }

    /** 
	You must call this method if the primary is a (interaction-weighted)
	neutrino. Then MCPrimarySpectrum weight is convoluted with the interaction
        weight. You have to call fillMCPrimarySpectrumWeight() a priori.
	InteractionsMatrix of NeutrinoCharge and NeutrinoNeutral objects 
	in the interaction package
	(iceCube.uhe.interactions) are set for calculating the relevant cross section.
    */
    public void fillNeutrinoWeight(InteractionsMatrix nuCCmtx, 
				   InteractionsMatrix nuNCmtx){
	if(fillMCWeight){
	    this.nuCCmtx = nuCCmtx; this.nuNCmtx = nuNCmtx;
	    isNeutrino = true; s = new ParticlePoint(0.0,0.0,0); //ice
	}else{
	    System.err.println("You must call fillMCPrimarySpectrumWeight()!");
	    System.exit(0);
	}
    }

    /**
       Call this method before executing process() when you read out
       the Corsika-mmc track primary.
     */
    public void corsikaMMC(){
	isCorsika = true;
	isNeutrino = false;
    }
      

    /** the method to build I3Particle objects with data from the DataInputStream. 
	Generated I3Particle objects are serialized and written out to OutputStream.
     */
    public void process(DataInputStream in, OutputStream out, boolean isFullData) 
	throws IOException {

	// Reading data
	BufferedReader  d     = new BufferedReader(new InputStreamReader(in));
	String buffer; int sep = 0; int sepstart = 0;
        char separator = ' ';
	while((buffer = d.readLine())!=null){
	    try{

		// 1st line -- eventNumber, flavor, doublet
		sepstart = 0;
		sep = buffer.indexOf(separator,sepstart+1);
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		int eventNumber =
		    Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		int flavor =
		    Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		int doublet =
		    Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		sepstart = sep;

		if(flavor == 0 && doublet == 0){
		    isNuE = true;
		    if(fillMCWeight && !hasGR){
			grLepton = new GlashowResonanceLeptonic(s,1);
			grHadron = new GlashowResonanceHadronic(s); 
			hasGR = true; // Has read the Glashow Resonance Matrix
		    }
		}
		else{
		    isNuE = false;
		}

		// 2nd line energy distance (energy reco)
		buffer = d.readLine();
		sepstart = 0;
		sep = buffer.indexOf(separator,sepstart+1);
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double energy =
		    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double distance =
		    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double energyReco = -1.0;
		if(sep!=-1){
		    energyReco =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		}
		sepstart = sep;

		// 3rd line -- direction (nx,ny,nz) (fitQuality)
		buffer = d.readLine();
		sepstart = 0;
		sep = buffer.indexOf(separator,sepstart+1);
		sepstart = sep;

		double[] n = new double[3];
		for(int jj=0;jj<3;jj++){
		    sep = buffer.indexOf(separator,sepstart+1);
		    n[jj] =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    sepstart = sep;
		}
		J3UnitVector direction = new J3UnitVector(n[0],n[1],n[2]);

		sep = buffer.indexOf(separator,sepstart+1);
		double fitQuality = -1.0;
		if(sep!=-1){
		    fitQuality =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		}
		sepstart = sep;

		// 4th line -- vertex (rx,ry,rz)
		buffer = d.readLine();
		sepstart = 0;
		sep = buffer.indexOf(separator,sepstart+1);
		sepstart = sep;

		double[] r = new double[3];
		for(int jj=0;jj<3;jj++){
		    sep = buffer.indexOf(separator,sepstart+1);
		    r[jj] =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    sepstart = sep;
		}
		J3Vector vertex = new J3Vector(r[0],r[1],r[2]);
		J3Line axis = new J3Line(vertex,direction);

		// Neutrino Line for rewrighting if this is a netruno
		if(isNeutrino){ // neutrino!
		    buffer = d.readLine();
		    sepstart = 0;
		    sep = buffer.indexOf(separator,sepstart+1);
		    sepstart = sep;
		    
		    sep = buffer.indexOf(separator,sepstart+1);
		    enhancementFactor =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    sepstart = sep;
		    
		    sep = buffer.indexOf(separator,sepstart+1);
		    trackLength =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    sepstart = sep;
		}

		// Corsika-mmc Line (if exists)
		if(isCorsika){
		    // CRenergy bundleEnergyAtEarthSurface
		    buffer = d.readLine();
		    sepstart = 0;
		    sep = buffer.indexOf(separator,sepstart+1);
		    sepstart = sep;

		    sep = buffer.indexOf(separator,sepstart+1);
		    energyReco =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    sepstart = sep;

		    sep = buffer.indexOf(separator,sepstart+1);
		    double energyAtsurface =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    sepstart = sep;

		    // bundle's (nx ny nz) -- skipped
		    buffer = d.readLine();
		}
		// Additional line (if exists)
		J3Line axisReco = null;
		if(isFullData){
		    // 5th line -- direction (nx,ny,nz) (fitQuality)
		    buffer = d.readLine();
		    sepstart = 0;
		    sep = buffer.indexOf(separator,sepstart+1);
		    sepstart = sep;

		    for(int jj=0;jj<3;jj++){
			sep = buffer.indexOf(separator,sepstart+1);
			n[jj] =
			    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
			sepstart = sep;
		    }
		    direction = new J3UnitVector(n[0],n[1],n[2]);

		    sep = buffer.indexOf(separator,sepstart+1);
		    if(sep!=-1){
			fitQuality =
			Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		    }
		    sepstart = sep;


		    // 6th line -- vertex (rx,ry,rz)
		    buffer = d.readLine();
		    sepstart = 0;
		    sep = buffer.indexOf(separator,sepstart+1);
		    sepstart = sep;

		    for(int jj=0;jj<3;jj++){
			sep = buffer.indexOf(separator,sepstart+1);
			r[jj] =
			    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
			sepstart = sep;
		    }
		    vertex = new J3Vector(r[0],r[1],r[2]);
		    axisReco = new J3Line(vertex,direction);
		}

		// 7th line IceCube data - ATWDNpe ATWDNDOMs FADCNPE FADCNDOMs BestNPE BestNDOMs
		buffer = d.readLine();
		sepstart = 0;
		sep = buffer.indexOf(separator,sepstart+1);
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double npeATWD =
		    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		int nDOMsATWD =
		    Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double npeFADC =
		    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		int nDOMsFADC =
		    Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		double npeBest =
		    Double.valueOf(buffer.substring(sepstart+1,sep)).doubleValue( );
		sepstart = sep;

		sep = buffer.indexOf(separator,sepstart+1);
		int nDOMsLaunch =
		    Integer.valueOf(buffer.substring(sepstart+1,sep)).intValue( );
		sepstart = sep;

		// Fills the MC Weight
		if(fillMCWeight){
		    mcWeightOfThisData =
			I3ParticleAnalysisFactory.getDNDLogE(powerLaw,energy,
							     log_EnergyMinimum,
							     log_EnergyMaximum);
		    if(isNeutrino){ // this is a neutrino. Mutiply the interaction weight.
			
			double sigma = 0.0; // cross section
			if(energy>0.0){
			    double logEnergy = Math.log(energy)/ln10;
			    int iLogE = (int)((logEnergy-Particle.getLogEnergyMinimum())/
					      Particle.getDeltaLogEnergy());
			    sigma = nuCCmtx.getSigmaMatrix(iLogE) +
				nuNCmtx.getSigmaMatrix(iLogE);
			    if(hasGR && isNuE){
				grLepton.setIncidentParticleEnergy(energy);
				grHadron.setIncidentParticleEnergy(energy);
				sigma +=0.5*(3.0*grLepton.getSigma()+grHadron.getSigma());
			    }
			}
			double density = s.getMediumDensity();
			double enhancement = (double )enhancementFactor;

			double weight = 
			    Math.exp(s.NA*density*trackLength*sigma*(enhancement-1.0))/
			    enhancement;
			mcWeightOfThisData = mcWeightOfThisData/weight;
			//System.err.println("sigma=" + sigma + " E=" + energy +
			//	   " eFac(" + enhancementFactor +
			//	   ") w=" + weight);
		    }

		}    
		//
		// Now generate an I3Particle object
		//
		I3Particle iceParticle = null;
		if(!isFullData){
		    iceParticle = generateI3Particle(flavor,doublet,energy,distance,axis,
						     eventNumber,
						     npeFADC,npeATWD,npeBest,
						     nDOMsFADC,nDOMsATWD,
						     nDOMsLaunch);
		    if(fillMCWeight)
			iceParticle.setMCPrimarySpectrumWeight(mcWeightOfThisData);
		}else{
		    iceParticle = generateI3Particle(flavor,doublet,energy,energyReco,
						     distance,
						     mcWeightOfThisData,
						     axis,axisReco,
						     eventNumber,
						     npeFADC,npeATWD,npeBest,
						     nDOMsFADC,nDOMsATWD,
						     nDOMsLaunch);
		}
								       
		// First guess results exists
		if(fitQuality!= -1.0) iceParticle.setFirstGuessQuality(fitQuality);
		// Write out I3Particle object
		I3ParticleOutputStream.outputI3Particle(iceParticle, out);

	    }catch (EOFException e){
		buffer = null;
		break;
	    }
	}

    }

    /** Generate I3Particle from a set of the given valuables. */
    protected I3Particle generateI3Particle(int flavor, int doublet,
					    double energy, double distance,
					    J3Line axisInIce3,
					    int eventNumber,
					    double npeFADC,double npeATWD, double npeBest,
					    int nDOMsFADC, int nDOMsATWD, 
					    int nDOMsLaunch){
	// Generate an I3Particle object
	I3Particle iceParticle = new I3Particle(flavor,doublet);

	if(isMCTruth) iceParticle.switchToMCTruth();
	else iceParticle.switchToReco();

	// Sets the Geometry
	iceParticle.setParticleAxisInIceCubeCoordinate(axisInIce3);
	iceParticle.transformParticleAxisToEarthCenterCoordinate();

	// Sets the energy
	if(isMCTruth){
	    iceParticle.putEnergy(energy);
	    if(energy>0.0){
		double logEnergy = Math.log(energy)/ln10;
		iceParticle.putLogEnergy(logEnergy);
	    }else{
		iceParticle.putLogEnergy(Double.NEGATIVE_INFINITY);
	    }
	}else{
	    iceParticle.putRecoEnergy(energy);
	}

	// Sets the distance from the Earth surface
	iceParticle.setDistanceFromEarthSurfaceToIceCube(distance);

	// Fills the IceCube data
	iceParticle.getIceCubeData().setEventNumber(eventNumber);
	iceParticle.getIceCubeData().setNpeATWD(npeATWD);
	iceParticle.getIceCubeData().setNpeFADC(npeFADC);
	iceParticle.getIceCubeData().setBestNpe(npeBest);
	iceParticle.getIceCubeData().setNDOMsATWD(nDOMsATWD);
	iceParticle.getIceCubeData().setNDOMsFADC(nDOMsFADC);
	iceParticle.getIceCubeData().setNDOMsLaunch(nDOMsLaunch);

	return iceParticle;

    }

    /** Generate I3Particle from a set of the given valuables.
	This method will be used when you build a MC event.
     */
    protected I3Particle generateI3Particle(int flavor, int doublet,
					    double energyMCTruth, double energyReco,
					    double distance, double mcWeight,
					    J3Line axisInIce3MCTruth,
					    J3Line axisInIce3Reco,
					    int eventNumber,
					    double npeFADC,double npeATWD, double npeBest,
					    int nDOMsFADC, int nDOMsATWD, 
					    int nDOMsLaunch){

	// Generate an I3Particle object
	isMCTruth = true;
	I3Particle iceParticle = 
	    generateI3Particle(flavor, doublet, energyMCTruth, distance,
			       axisInIce3MCTruth,
			       eventNumber, npeFADC, npeATWD, npeBest,
			       nDOMsFADC, nDOMsATWD, nDOMsLaunch);
	// Set the reconsturcted geometry
	iceParticle.switchToReco();
	iceParticle.setParticleAxisInIceCubeCoordinate(axisInIce3Reco);
	iceParticle.transformParticleAxisToEarthCenterCoordinate();

	// set the (reco) energy
	iceParticle.putRecoEnergy(energyReco);

	// set the MC primary spectrum weight
	iceParticle.switchToMCTruth();
	iceParticle.setMCPrimarySpectrumWeight(mcWeight);

	return iceParticle;

    }	

}
