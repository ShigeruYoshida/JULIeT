package iceCube.uhe.interactions;

import java.io.*;
import java.util.*;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import numRecipes.*;

/** 
<pre>
    The Interection Abstract class to provide several methods
    concerning particle interactions. 

    The Pair Creation, Bremsstrahlung, Photo nuclear interactions, and 
    the Week interactions involving UHE leptons and neutrinos 
    propagating undergound rock and ice are calculated in this subclass.

    The classvariables  Particle, and Point are necessary for the methods
    described here because the cross section depends on the particle propaty
    and the medium like Z and A(atomic number).

    The flavor of the produced particle, producedFlavor, is also set
    in the constructor in order to distingush mu+mu- pair creation 
    from e+e- pair for example.

</pre>
*/
abstract class Interactions implements Function, Serializable{

    /** Fine Stracture Const.*/
    static final double Alpha = 1.0D/137.035989561D;
    /** Classical Electron Radius [cm].*/
    static final double Re = 2.81794e-13;
    /** Compton wavelength [cm].*/
    static final double LambdaE = 3.86159e-11;

    /** the base of natural logarithm.*/
    static final double E = Math.E;

    /** RoundOff Error.*/
    static final double roundOffError = 1.0e-9;

    /** The Particle class involved with this interaction.*/
    Particle p;
    /** The ParticlePoint class involced with this interaction.*/
    ParticlePoint s;
    /**
       <pre>
      Flavor of the produced particle.
      For example 0 for e+e- creation
      For example 1 for mu+mu- creation
      </pre>
    */
    int producedFlavor;

    /** The incident particle energy [GeV].*/
    double energy;
    /** Mass of produced particle [GeV]. */
    double producedMass;
    /** Mass of Incoming particle [GeV]. */
    double mass;
    /** energy to define the integral range of y. */
    double energyCut = 0.0;

    /** for interface getFunction( ). */
    double[] parameters;
    double a;


    /** Constructor: Register the Particle and ParticlePoint classes.
	Also check if a given particle can be involved in this intereaction.
	<pre>

	Particle p       : Particle class involved in this interaction.
	ParticlePoint s  : ParticlePoint class to provide the propaty of medium.
	int flavor       : flavor of the produced particle. 0 for e+e- production
                           1 for mu+mu- production for example.
			   This is specified in the individual interaction
			   constructor whenever necessary. For Bremsstrahlung
			   for example, it calls just super(p,s,0) because
			   the produced particle's mass (in this case 0 -Photon)
			   does not appear in the cross section form.
        </pre>
    */
    public Interactions(Particle p, ParticlePoint s, int flavor){
	// flavor specifies produced lepton
	if(isValidParticle(p)){
	    this.p = p;
	    this.s = s;
	    this.producedFlavor = flavor;
	    this.energy = p.getEnergy();
  	    producedMass = Particle.particleMasses[flavor][1];
	    mass = p.getMass( );
	    a = energy;
	    parameters = new double[1];
	}else{
	    System.err.println("This particle " + 
		   p.particleName(p.getFlavor(), p.getDoublet()) + 
		   " should not be involved in this interaction");
	    System.exit(0);
	}
    }

    /** 
	<pre>
	Set The incident particle energy [GeV].
	The default value has been given by the constructor
	Interactions( ) with p.getEnergy( ) where p is 
	the Particle class. You might want to set, however,
	a different value such as 
	logE = logEnergyMinimum + deltaLogE*ilogE
	where ilogE is i'th index of the logEnergyMatrix
	in the Particle class p. This method provides you
	with a way to put the incident particle energy.
	</pre>
    */
    public void setIncidentParticleEnergy(double energy){
	this.energy = energy;
	p.putEnergy(energy);
    }
    public void setIncidentParticleEnergy(int iLogE){
        double logEnergy = p.getLogEnergyMinimum()+
	    p.getDeltaLogEnergy()*(double )iLogE;
        double energy = Math.pow(10.0,logEnergy);
	this.energy = energy;
	p.putEnergy(energy);
	p.putLogEnergy(logEnergy);
    }
    public double getIncidentParticleEnergy( ){
	return energy;
    }



    /** 
	<pre>
	Differential cross section dsigma/dy [cm^2]
	y = 1 - Erecoiling/Eincoming -- inelasticity parameter
	</pre>
    */
    abstract double getDSigmaDy(double y);



    /** 
	<pre>
	Differential cross section dsigma/dz [cm^2]
	z = Erecoiling/Eincoming -- inelasticity parameter
	</pre>
    */
    public double getDSigmaDz(double z){
	return getDSigmaDy(1.0-z);
    }


    /** Checking the range of the given inelasticity y 
	that is determined in an individual interaction channel.
    */
    abstract boolean isValidInelasticity(double y);

    /** Getting the range of allowed y for a given interaction */
    abstract double getYmin( );
    abstract double getYmax( );

    /** Energy Cut Parameter in integration to obtain the total cross section.*/
    public void setEnergyCut(double cutEnergy){
	energyCut = cutEnergy;
    }
    public double getEnergyCut(){
	return energyCut;
    }


    /** Checking the particle kind involved with
	a given interaction.
    */
    abstract boolean isValidParticle(Particle p);


    /** Total cross section [cm^2] */
    public double getSigma( ){
        int functionIndex =1; // Specifies dSigma/dy
	Interactions interactions = this;
	double yMin = getYmin( )+2.0*roundOffError;
	double yMax = getYmax( )-2.0*roundOffError;
	double sum = Integration.RombergIntegral(interactions, 
                                 functionIndex, parameters, yMin, yMax);
	return sum;
    }


    /** Integral dSigma/dy over a given range to obtain a partial
	cross section */
    public double integralDSigmaDy(double lowerY, double upperY){
	int functionIndex = 1; // Specifies dSigma/dy
	Interactions interactions = this;
	double sum = 0.0;
	if(isValidInelasticity(lowerY) && isValidInelasticity(upperY)){
	    sum = Integration.RombergIntegral(interactions, 
			      functionIndex, parameters, lowerY, upperY);
	}
	else{
	    showIntegralErrorMessage(lowerY,upperY);
	}
	return sum;
    }
    /** Integral dSigma/dz over a given range to obtain a partial
	cross section */
    public double integralDSigmaDz(double lowerZ, double upperZ){
	int functionIndex = 2; // Specifies dSigma/dz
	Interactions interactions = this;
	double sum = 0.0;
	if(isValidInelasticity(1.0-lowerZ) && isValidInelasticity(1.0-upperZ)){
	    sum = Integration.RombergIntegral(interactions, 
			      functionIndex, parameters, lowerZ, upperZ);
	}
	else{
	    showIntegralErrorMessage(1.0-upperZ,1.0-lowerZ);
	}
	return sum;
    }
    /** Integral y*dSigma/dy over a given range to obtain the inelasticity
	distribution. The energy transfer probability would also require
        this value. */
    public double getYDSigmaDy(double lowerY, double upperY){
	int functionIndex = 3; // Specifies y x dSigma/dy
	Interactions interactions = this;
	double sum = 0.0;
	if(isValidInelasticity(lowerY) && isValidInelasticity(upperY)){
	    sum = Integration.RombergIntegral(interactions, 
			      functionIndex, parameters, lowerY, upperY);
	}
	else{
	    showIntegralErrorMessage(lowerY,upperY);
	}
	return sum;
    }
    /** Integral z*dSigma/dz over a given range to obtain the inelasticity
	distribution. The energy transfer probability would also require
        this value.  z = 1-y */
    public double getZDSigmaDZ(double lowerZ, double upperZ){
	int functionIndex = 4; // Specifies z x dSigma/dz
	Interactions interactions = this;
	double sum = 0.0;
	if(isValidInelasticity(1.0-lowerZ) && isValidInelasticity(1.0-upperZ)){
	    sum = Integration.RombergIntegral(interactions, 
			      functionIndex, parameters, lowerZ, upperZ);
	}
	else{
	    showIntegralErrorMessage(1.0-upperZ,1.0-lowerZ);
	}
	return sum;
    }


    /** 
	<pre>
	Method for interface <Function>. 
        Interface the differntial cross sections given here
        to the utility methods such as the Romberg
        Integration code that is desinged for a genereal
        function in form of Func(x).

	functionIndex     1     dsigma/dy
	functionIndex     2     dsigma/dz  z = y-1
	functionIndex     3     y x dsigma/dy 
	functionIndex     4     z x dsigma/dz 
	</pre>
    */
    public double getFunction(int functionIndex, double[] parameters, 
	double x){ 
	double dSigma;
	double y;

        switch (functionIndex) {
        case 1 : 
	    dSigma = getDSigmaDy(x);
	    break;
	case 2 :
	    dSigma = getDSigmaDz(x);
	    break;
	case 3:
	    dSigma = getDSigmaDy(x)*x;
	    break;
	case 4:
	    dSigma = getDSigmaDz(x)*x;
	    break;
        default:
            dSigma = 0.0;
            System.err.println("Illegal parameters! Index" + functionIndex);
            System.exit(0);
        }

        return dSigma;

    }


    


    /** Show Name of the Interaction */
    abstract String interactionName();


    /** Error message utility **/
    public void showIntegralErrorMessage(double lowerY, double upperY){

	System.err.println(interactionName());
	System.err.println("The integral range out of the allowed y");
	System.err.println("Ymin " + getYmin() + " Ymax " + getYmax());
	System.err.println("lowerY " + lowerY + " upperY " + upperY);
	System.err.println("The Incident Particle " +
			   p.particleName(p.getFlavor(),p.getDoublet()));
	System.err.println("The Incident Energy " + energy + " [GeV]");
    }

}


