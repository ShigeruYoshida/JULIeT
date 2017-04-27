package iceCube.uhe.interactions;

import java.io.*;
import java.util.*;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import numRecipes.*;

/** The Pair Creation involving UHE leptons
    propagating undergound rock and ice are calculated in this subclass.

    The class variables  Particle, and Point are necessary for the methods
    described here because the cross section depends on the particle propaty
    and the medium like Z and A(atomic number).

*/

public class PairCreation extends Interactions implements Function{

    private double massRatio;
    private double producedScale;
    private double[] para = new double[2];



    /** Constructor: Register the Particle and ParticlePoint classes */
    public PairCreation(Particle p, ParticlePoint s, int flavor){
	super(p, s, flavor);
	massRatio = mass/producedMass;
	producedScale = Particle.particleMasses[0][1]/producedMass;
	// Ratio of electron mass (e+e-) to produced leptons (e+e- mu+mu- etc)
	System.err.println("producedScale " + producedScale);
    }



    /** Differential cross section dsigma/dy/drho [cm^2]
	y = 1 - Erecoiling/Eincoming -- inelasticity parameter
	The assymetry factor is defined by rho = (E^+ - E^-)/(E^+ + E^-) 
    */
    public double getDSigmaDyDrho(double y, double rho){
	if(!isValidInelasticity(y)) return 0.0;
	double factor = Alpha*Alpha*Re*Re*producedScale*producedScale*
	    2.0/(3.0*Math.PI)*(1.0-y)/y;
	double chargeFactor = 0.0;
	for(int i=0;i<s.NumberOfSpecies[s.getMaterialNumber( )];i++){
	    if(y<getYmaxCharge(i))
  	    chargeFactor += s.getCharge(i)*(s.getCharge(i)+getScreenFactor(i))
              *(double )(s.getNumberOfAtoms(i))*getAsymmetryTerm(rho,y,i);
	}
	double dSigmaDy = factor*chargeFactor;
	return dSigmaDy;
    }




    /** Differential cross section dsigma/dy [cm^2]
	y = 1 - Erecoiling/Eincoming -- inelasticity parameter    */
    public double getDSigmaDy(double y){
	if(!isValidInelasticity(y)) return 0.0;
	double factor = Alpha*Alpha*Re*Re*producedScale*producedScale*
	    2.0/(3.0*Math.PI)*(1.0-y)/y;
	double chargeFactor = 0.0;
	double rhoBound = (1.0-6.0*mass*mass/(energy*energy*(1.0-y)))*
	    Math.sqrt(1.0-4.0*producedMass/(energy*y));

	for(int i=0;i<s.NumberOfSpecies[s.getMaterialNumber( )];i++){
	    PairCreation paircreation = this;
	    // The asymmetry term integration
	    para[0]=y;para[1]=(double )i; 
	    double sum;
	    if(y<getYmaxCharge(i)){
		sum = Integration.RombergIntegral(paircreation, 
                               5, para, -rhoBound, rhoBound);
	    }else{
		sum = 0.0;
	    }
	    chargeFactor += s.getCharge(i)*(s.getCharge(i)+getScreenFactor(i))
		            *(double )(s.getNumberOfAtoms(i))*sum;
	}
	double dSigmaDy = factor*chargeFactor;
	//System.err.println("y= " + y + " dSigma/dy = " + dSigmaDy);
	return dSigmaDy;
    }




    /** Differential cross section dsigma/dy [cm^2]
	yPlus = Epositron/Eincoming  */
    public double getDSigmaDyPlus(double yPlus){
	para[0] = yPlus;
        double rhoMax = 
	(1.0-6.0*mass*mass/(energy*energy*(1.0-getYmax())))*
	Math.sqrt(1.0-4.0*producedMass/(energy*getYmax()));
	double rhoUpperBound;
	double rhoLowerBound;
	double yMin = getYmin( );
	double yMax = getYmax( );
	if((yMin*(1.0+rhoMax)/2.0)<= yPlus &&
	    yPlus <= (yMax*(1.0+rhoMax)/2.0)){

	    rhoUpperBound = rhoMax;
	    rhoLowerBound = 2.0*yPlus/yMax-1.0;

	}else if((yMax*(1.0-rhoMax)/2.0) <= yPlus && 
	    yPlus< (yMin*(1.0+rhoMax)/2.0)){

	    rhoUpperBound = 2.0*yPlus/yMin-1.0;
	    rhoLowerBound = 2.0*yPlus/yMax-1.0;

	}
	else if((yMin*(1.0-rhoMax)/2.0) <= yPlus && 
	    yPlus< (yMax*(1.0-rhoMax)/2.0)){

	    rhoUpperBound = 2.0*yPlus/yMin-1.0;
	    rhoLowerBound = -rhoMax;

	}
	else{
	    rhoUpperBound = 0.5;
	    rhoLowerBound = 0.5;
	}
        PairCreation paircreation = this;
        double sum = 
	Integration.RombergIntegral(paircreation, 6, para, 
				    rhoLowerBound, rhoUpperBound);
	return sum;
    }
        



    /** Integral dSigma/dyPlus over a given range to obtain a partial
	cross section */
    public double integralDSigmaDyPlus(double lowerY, double upperY){
	int functionIndex = 7; // Specifies dSigma/dy+
        PairCreation paircreation = this;
	double sum = 0.0;
	if(isValidInelasticityPlus(lowerY) && isValidInelasticityPlus(upperY)){
	    sum = Integration.RombergIntegral(paircreation, 
			      functionIndex, parameters, lowerY, upperY);
	}
	return sum;
    }



    /** Integral yPlus*dSigma/dyPlus over a given range to obtain 
	the inelasticity
	distribution. The energy transfer probability would also require
        this value. */
    public double getYPlusDSigmaDyPlus(double lowerY, double upperY){
	int functionIndex = 8; // Specifies yPlus x dSigma/dyPlus
        PairCreation paircreation = this;
	double sum = 0.0;
	if(isValidInelasticityPlus(lowerY) && isValidInelasticityPlus(upperY)){
	    sum = Integration.RombergIntegral(paircreation, 
			      functionIndex, parameters, lowerY, upperY);
	}
	return sum;
    }



    /** 
	<pre>
	Override Method for interface <Function>. 
        Interface the differntial cross sections given here
        to the utility routiner such as the Romberg
        Integration code that is desinged for a genereal
        function in form of Func(x).

	functionIndex     1     dsigma/dy
	functionIndex     2     dsigma/dz  z = y-1
	functionIndex     3     y x dsigma/dy 
	functionIndex     4     z x dsigma/dz 
	functionIndex     5     dn/drho (The assymetry term)
	functionIndex     6     dsigma/dy-drho x 2/(1+rho)
	functionIndex     7     dsigma/dyPlus
	functionIndex     8     dsigma/dyPlus x yPlus
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
        case 5 : 
	    y = para[0];
	    int ithSpecies = (int )para[1];
	    dSigma = getAsymmetryTerm(x,y,ithSpecies);
	    break;
	case 6 :
	    double yPlus = para[0];
	    y = 2.0*yPlus/(1.0+x);
	    dSigma = getDSigmaDyDrho(y,x)*2.0/(1.0+x);
	    break;
	case 7 :
	    dSigma = getDSigmaDyPlus(x);
	    break;
	case 8 :
	    dSigma = getDSigmaDyPlus(x)*x;
	    break;
        default:
            dSigma = 0.0;
            System.err.println("Illegal parameters! Index" + functionIndex);
            System.exit(0);
        }

        return dSigma;

    }




    /** Calculate the term on the asymmetry factor of the pair creation.
	The assymetry factor is defined by rho = (E^+ - E^-)/(E^+ + E^-) */
    public double getAsymmetryTerm(double rho, double y, int ithSpecies){
	double beta = y*y/(2.0*(1.0-y));
	double qsi = (y*massRatio/2.0)*(y*massRatio/2.0)*(1.0-rho*rho)/(1.0-y);
	double Yrec = (5.0-rho*rho+4.0*beta*(1.0+rho*rho))/
          	     (2.0*(1.0+3.0*beta)*Math.log(3.0+1.0/qsi)-rho*rho-
                      2.0*beta*(2.0-rho*rho));
        double Yin = (4.0+rho*rho+3.0*beta*(1.0+rho*rho))/
	    ((1.0+rho*rho)*(1.5+2.0*beta)*Math.log(3.0+qsi)+1.0-1.5*rho*rho);
	double chargeTerm = Math.pow(s.getCharge(ithSpecies),-1.0/3.0);
	double Lrec = Math.log(s.getRadiation(ithSpecies)*
        chargeTerm*Math.sqrt((1.0+qsi)*(1.0+Yrec))/
        (1.0+2.0*producedMass*Math.sqrt(E)*s.getRadiation(ithSpecies)*chargeTerm*
	 (1.0+qsi)*(1.0+Yrec)/(energy*y*(1.0-rho*rho)))) -
         0.5*Math.log(1.0+(1.5*1.5/(massRatio*massRatio*chargeTerm*chargeTerm))*
		      (1.0+qsi)*(1.0+Yrec));
	double Lin = Math.log(s.getRadiation(ithSpecies)*chargeTerm*chargeTerm*
        2.0/3.0*massRatio/(1.0+2.0*producedMass*Math.sqrt(E)*
        s.getRadiation(ithSpecies)*chargeTerm*(1.0+qsi)*(1.0+Yin)/
        (energy*y*(1.0-rho*rho))));

	double recoileTerm = 
          (((2.0+rho*rho)*(1.0+beta)+qsi*(3.0+rho*rho))*Math.log(1.0+1.0/qsi)+
	  (1.0-rho*rho-beta)/(1.0+qsi)-(3.0+rho*rho))*Lrec;
	double incomingTerm = 
          (((1.0+rho*rho)*(1.0+1.5*beta)-(1.0+2.0*beta)*(1.0-rho*rho)/qsi)*
          Math.log(1.0+qsi)+
	  qsi*(1.0-rho*rho-beta)/(1.0+qsi)+(1.0+2.0*beta)*(1.0-rho*rho))*Lin;
	double asymmetryTerm = recoileTerm+incomingTerm/massRatio/massRatio;
	return asymmetryTerm;
    }




    /** Calculate the factor due to the screening effect
	on an atomic ellectron */
    public double getScreenFactor(int ithSpecies){
	double screenFac;
	if(energy<= 35.0*mass){
	    screenFac = 0.0;
	}else{
	    double gamma1;
	    double gamma2;
	    if(s.getCharge(ithSpecies)==1.0){ // H
		gamma1 = 4.4e-5;
		gamma2 = 4.8e-5;
	    }else{                            // otherwise
		gamma1 = 1.95e-5;
		gamma2 = 5.30e-5;
	    }
	    screenFac = 
	    (0.073*Math.log((energy/mass)/
            (1.0+gamma1*Math.pow(s.getCharge(ithSpecies),2.0/3.0)*
            energy/mass))-0.26)/
	    (0.058*Math.log((energy/mass)/
            (1.0+gamma2*Math.pow(s.getCharge(ithSpecies),1.0/3.0)*
             energy/mass))-0.14);
	    if(screenFac<0.0) screenFac = 0.0;
	}

	return screenFac;
    }





    /** Checking the range of the given inelasticity y 
	that is determined in an individual interaction channel.
    */
    public boolean isValidInelasticity(double y){ 
	if((getYmin()+roundOffError)<= y && 
	   y <= (getYmax()-roundOffError)) return true;
	else return false;
    }
    public boolean isValidInelasticityPlus(double yPlus){ 
	double rhoUpperBound = 
	(1.0-6.0*mass*mass/(energy*energy*(1.0-getYmax())))*
	Math.sqrt(1.0-4.0*producedMass/(energy*getYmax()));
	if((producedMass/energy+roundOffError)<= yPlus && 
	   yPlus <= (getYmax()*(1.0+rhoUpperBound)/2.0 - roundOffError)) 
	    return true;
	else return false;
    }




    /** Getting the range of allowed y for a given interaction */
    public double getYmin( ){
	if(energyCut == 0.0){
	    double yMin = 4.0*producedMass/energy;
	    return yMin;
	}
	else {
	    return energyCut/energy;
	}
    }
    public double getYmax( ){
	double chargeTerm =1.0; // Pure Hydrogen as an upperbound
        double yMax = 1.0-3.0/4.0*Math.sqrt(E)*mass/energy*chargeTerm;
	return yMax;
    }
    public double getYmaxCharge(int ithSpecies){
	double chargeTerm = Math.pow(s.getCharge(ithSpecies),1.0/3.0);
	double yMax = 1.0-3.0/4.0*Math.sqrt(E)*mass/energy*chargeTerm;
	return yMax;
    }

    public double getYPlusMin( ){
	double rhoUpperBound = 
	(1.0-6.0*mass*mass/(energy*energy*(1.0-getYmax())))*
	Math.sqrt(1.0-4.0*producedMass/(energy*getYmax()));
	double yMin = (getYmin()*(1.0-rhoUpperBound)/2.0);
	return yMin;
    }
    public double getYPlusMax( ){
	double rhoUpperBound = 
	(1.0-6.0*mass*mass/(energy*energy*(1.0-getYmax())))*
	Math.sqrt(1.0-4.0*producedMass/(energy*getYmax()));
	double yMax = (getYmax()*(1.0+rhoUpperBound)/2.0);
	return yMax;
    }



    /** Checking the particle kind involved with
	a given interaction. Only muon and tauon
	but not neutrinos and pions
	are allowed to be involved with the pair creation
	in the medium.
    */
    public boolean isValidParticle(Particle p){
	if(p.getDoublet()==1 && p.getFlavor()>=0 && 3>p.getFlavor())
	return true;
	else return false;
    }




    /** Show Name of the Interaction */
    public String interactionName(){
	String channel = "Pair Creation ";
	String producedParticle = Particle.particleName(producedFlavor,1);
	String incidentParticle = 
	    p.particleName(p.getFlavor(), p.getDoublet());
	String name = channel.concat("from ").
	    concat(incidentParticle).concat(" to ").concat(producedParticle);
	return name;
    }
}


