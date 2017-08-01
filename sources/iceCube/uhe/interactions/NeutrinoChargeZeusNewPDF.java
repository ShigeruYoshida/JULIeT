package iceCube.uhe.interactions;

import java.io.*;
import java.util.*;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import numRecipes.*;

/** The Neutrino-isoscalar nucleon charged current interactions
taking place in undergound rock and ice are calculated in this subclass.

The class variables  Particle, and Point are necessary for the methods
described here because the cross section depends on the particle propaty
and the medium like Z and A(atomic number).

*/
public class NeutrinoChargeZeusNewPDF extends Interactions implements Function{
    private static final long serialVersionUID = 5427600879174779867L;
    private String neutrinoFile = "iceCube/uhe/interactions/diffy_nu_CC_HERAPDF15NLO.dat";
    private static final int numberOfEnergyBin = 91;
    private static final int numberOfYBin = 102;
    private double[] logEArray = new double[numberOfEnergyBin];
    private double[][] logDyArray = new double[numberOfEnergyBin][numberOfYBin];
    private double[][] yDsigmaArray = new double[numberOfEnergyBin][numberOfYBin];
    protected static double yLowerLimit = 1.0e-13;
    protected static double dsLowerLimit = 1.0e-20; // [mb]
    protected static final double mb2cm2 = 1.0e-27;  //[cm2/mb]
	
	/** Constructor: Register the Particle and ParticlePoint classes.
	It also reads the pre-calculated y/E * dsigma/dy
        from the data file. */
	public NeutrinoChargeZeusNewPDF(Particle p,ParticlePoint s) throws IOException{
		super(p,s,3);
		
		File fileName = new File(neutrinoFile);
		   BufferedReader in = 
			   new BufferedReader(new FileReader(fileName));
		   
		char separator = ' ';
		int n = 0;
		int sep = 0;  int sepstart = 0;
		String buffer;
		
		for(n = 0; n < numberOfEnergyBin ; n++){
		    for(int i = 0; i < numberOfYBin; i++){
			buffer = in.readLine();
			sepstart = 0;

			sep = buffer.indexOf(separator,sepstart);
			double energyData = Double.valueOf(buffer.substring(sepstart,sep)).doubleValue();
			if(i == 0){
			    logEArray[n] = Math.log(energyData);
			}
			sepstart = sep;
			
			sep = buffer.indexOf(separator,sepstart+1);
			double dy = 
			    Double.valueOf(buffer.substring(sepstart,sep)).doubleValue();
			sepstart = sep;
			if(dy == 0.0) dy= yLowerLimit;
			logDyArray[n][i] = Math.log(dy);

			sep = buffer.indexOf(separator,sepstart+1);
			double ds = 
			    Double.valueOf(buffer.substring(sepstart,sep)).doubleValue(); // [mb]
			sepstart = sep;
			if(ds == 0.0) ds = dsLowerLimit;
			yDsigmaArray[n][i] = Math.log(ds*dy/energyData);

			//if(n%10==0 && i%50==0) System.err.format("energy(%e) dy (%e) ds (%e)\n",energyData,dy,ds);
			
		    }

		}
		in.close();
		
	}
	
	 /** Differential cross section dsigma/dy [cm^2]
	y = 1 - Erecoiling/Eincoming -- inelasticity parameter  */
	public double getDSigmaDy(double y){
		double logEnergy = Math.log(energy);
		int index = Interpolation.searchIndex(logEArray,logEnergy,numberOfEnergyBin);
		if(index>=numberOfEnergyBin-2) index = numberOfEnergyBin-2;
		double logELow = logEArray[index];
		double logEUp = logEArray[index+1];
		if(!isValidInelasticity(y)) return 0.0;
		
		double logY = Math.log(y);
		double dSigmaDyLow = 
			Interpolation.mThPolynominalInterpolate(logDyArray[index],
			yDsigmaArray[index],numberOfYBin,logY,4);
		double dSigmaDyUp = 
			Interpolation.mThPolynominalInterpolate(logDyArray[index+1],
			yDsigmaArray[index+1],numberOfYBin,logY,4);
		double yEdSigmaDy = dSigmaDyLow + (dSigmaDyUp-dSigmaDyLow)/
		    (logEUp-logELow)*(logEnergy-logELow);
		double dSigmaDy = Math.exp(yEdSigmaDy)*energy/y*mb2cm2;
		return dSigmaDy;
	}
	
	
	/** Checking the range of the given inelasticity y 
	that is determined in an individual interaction channel.
    */
	public boolean isValidInelasticity(double y){
		if((getYmin()-roundOffError) <= y && 
				(y <= getYmax()+roundOffError)) return true;
		else return false;
	}
	
	
	/** Getting the range of allowed y for a given interaction */
	public double getYmin(){
		double logEnergy = Math.log(energy);
		int index = Interpolation.searchIndex(logEArray,logEnergy,numberOfEnergyBin);
		if(index>=numberOfEnergyBin-2) index = numberOfEnergyBin-2;
	        double yMin = Math.exp(logDyArray[index][1]);
		return yMin;
	}
	public double getYmax(){
		double logEnergy = Math.log(energy);
		int index = Interpolation.searchIndex(logEArray,logEnergy,numberOfEnergyBin);
		if(index>=numberOfEnergyBin-2) index = numberOfEnergyBin-2;
		double yMax = Math.exp(logDyArray[index][numberOfYBin-1]);
		return yMax;
	}
	
	
	/** Checking the particle kind involved with
	a given interaction. Only neutrinos
	are allowed to be involved with this interaction
	in the medium.
    */
	public boolean isValidParticle(Particle p){
		if(p.getDoublet()==0 && p.getFlavor()>=0 && 3>p.getFlavor())
			return true;
		else return false;
	}
	
	
	/** Show Name of the Interaction */
	public String interactionName(){
		String channel = "Neutrino-Nuclen Charged Interaction ";
		String incidentParticle = 
			p.particleName(p.getFlavor(),p.getDoublet());
		String name = channel.concat("from").concat(incidentParticle);
		return name;
	}
	
	
}
