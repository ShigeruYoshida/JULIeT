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
public class NeutrinoChargeZeusPDF extends Interactions implements Function{
    private static final long serialVersionUID = 5427600879174779867L;
	private String neutrinoFile = "iceCube/uhe/interactions/nuCC_ZeusPDF.dat";
	private static final int numberOfEnergyBin = 22;
	private static final int numberOfYBin = 321;
	private double[] logEArray = new double[numberOfEnergyBin];
	private double[] logDyArray = new double[numberOfYBin];
	private double[][] yDsigmaArray = new double[numberOfEnergyBin][numberOfYBin];
	
	/** Constructor: Register the Particle and ParticlePoint classes.
	It also reads the pre-calculated y/E * dsigma/dy
        from the data file. */
	public NeutrinoChargeZeusPDF(Particle p,ParticlePoint s) throws IOException{
		super(p,s,3);
		
		File fileName = new File(neutrinoFile);
		   BufferedReader in = 
			   new BufferedReader(new FileReader(fileName));
		   
		char separator = ' ';
		int n = 0;
		int sep = 0;  int sepstart = 0;
		String buffer;
		
		for(n = 0; n < numberOfYBin ; n++){
			buffer = in.readLine();
			sepstart = 0;
			sep = buffer.indexOf(separator,sepstart);
			double dy = 
				Double.valueOf(buffer.substring(sepstart,sep)).doubleValue();
			logDyArray[n] = Math.log(dy);
			for(int i = 0; i < numberOfEnergyBin; i++){
				sepstart = sep;
				sep = buffer.indexOf(separator,sepstart+1);
			    double energy =
					  Double.valueOf(buffer.substring(sepstart,sep)).doubleValue();
			    if(n == 0){
				  logEArray[i] = Math.log(energy);
				}
				sepstart = sep;
				sep = buffer.indexOf(separator,sepstart+1);
				double ds = 
					Double.valueOf(buffer.substring(sepstart,sep)).doubleValue();
				yDsigmaArray[i][n] = Math.log(ds*dy/energy);
			}
			
		}
		in.close();
		
	}
	
	 /** Differential cross section dsigma/dy [cm^2]
	y = 1 - Erecoiling/Eincoming -- inelasticity parameter  */
	public double getDSigmaDy(double y){
		if(!isValidInelasticity(y)) return 0.0;
		double logEnergy = Math.log(energy);
		int index = Interpolation.searchIndex(logEArray,logEnergy,numberOfEnergyBin);
		double logELow = logEArray[index];
		double logEUp = logEArray[index+1];
		
		double logY = Math.log(y);
		double dSigmaDyLow = 
			Interpolation.mThPolynominalInterpolate(logDyArray,
			yDsigmaArray[index],numberOfYBin,logY,4);
		double dSigmaDyUp = 
			Interpolation.mThPolynominalInterpolate(logDyArray,
			yDsigmaArray[index+1],numberOfYBin,logY,4);
		double yEdSigmaDy = dSigmaDyLow + (dSigmaDyUp-dSigmaDyLow)/
		    (logEUp-logELow)*(logEnergy-logELow);
		double dSigmaDy = Math.exp(yEdSigmaDy)*energy/y*1.0e-36;
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
	        double yMin = Math.exp(logDyArray[0]);
		return yMin;
	}
	public double getYmax(){
		double yMax = Math.exp(logDyArray[numberOfYBin-1]);
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
