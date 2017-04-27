package iceCube.uhe.muonModel;

import iceCube.uhe.muonModel.*;
import iceCube.uhe.points.*;
import iceCube.uhe.particles.*;
import iceCube.uhe.geometry.*;
import geometry.*;

import java.io.*;
import java.util.*;

/** 

    The muon bundle flux table relative to those given by
    the original parameters (alpha=2.04 Eth=3730=GeV - settable by the constructor
    in this class) to a given alpha with the fixed Eth.
    The ELbertFluxFactory object takes care of the original flux calculation
    and those with given alpha, respectively.


 */
public class RelativeElbertFluxTableMaker {

    protected double alphaReference= 2.04;
    protected double muEThresholdReference =  3730.0;
    protected double[][] dFMuDLogEReferece; // array to store the reference flux

    protected static String table1DPath = "../data/muonBundleFluxTable2TeVthresholdNoCutoffIron/";
    double alpha = 1.95;
    ElbertFluxTableFactory muonFluxTable = null;
    private boolean hasReadTable = false;

    protected int numberOfCosZenith = 100;
    protected double stepCosZenith = 0.01;

    IceCubeCoordinate iceCube = new IceCubeCoordinate();

    /** default constructor */
    public RelativeElbertFluxTableMaker()throws IOException{
	muonFluxTable = new ElbertFluxTableFactory();
	muonFluxTable.setElbertParameters(alphaReference,muEThresholdReference);
	fillReferenceFluxArray();
    }

    /** constructor with the reference parameters which gives the original fluxe*/
    public RelativeElbertFluxTableMaker(double alphaReference, double muEThresholdReference) throws IOException{
	this.alphaReference= alphaReference;
	this.muEThresholdReference = muEThresholdReference;
	muonFluxTable = new ElbertFluxTableFactory();
	muonFluxTable.setElbertParameters(alphaReference,muEThresholdReference);
	fillReferenceFluxArray();
    }

    private void fillReferenceFluxArray(){
	dFMuDLogEReferece = new double[numberOfCosZenith][muonFluxTable.numberOfEmuon];

	double cosZenith = 1.0;
	for(int icos=0;icos<numberOfCosZenith;icos++){
	    double distance = distanceFromCosZenith(cosZenith);

	    for(int iLogE=0;iLogE<muonFluxTable.numberOfEmuon;iLogE++){
		double logEmuon = muonFluxTable.logEmuonMin + 
		    muonFluxTable.logEmuonStepSize*(double )iLogE;
		dFMuDLogEReferece[icos][iLogE]= 
		    muonFluxTable.getDFMuDLogE(distance,logEmuon);
	    }
	    //System.err.println("cosZenith=" + cosZenith + " log(distance)=" + 
	    //	       Math.log(distance)/Math.log(10.0));
	    cosZenith -= stepCosZenith;

	}
    }

    private double distanceFromCosZenith(double cosZenith){
	double rEarth_ice = ParticlePoint.REarth + iceCube.getGlacierDepth();
	double detectorDepth = iceCube.getGlacierDepth()-iceCube.elevation;
	//880m from the icecube center
	double distance = -(rEarth_ice-detectorDepth)*cosZenith +
		Math.sqrt((rEarth_ice-detectorDepth)*(rEarth_ice-detectorDepth)*(cosZenith*cosZenith-1.0) + rEarth_ice*rEarth_ice);
	return distance;
    }


    public void setAlpha(double alpha){
	this.alpha = alpha;
	hasReadTable = false;
    }

    public double getAlpha(double alpha){
	return alpha;
    }

    public double getRelativeFlux(double logEmuon, double cosZenith) throws IOException{
	if(!hasReadTable){
	    ElbertFluxTableFactory.tablePath= table1DPath;
	    ElbertFluxTableFactory.alphaMin= 1.8;
	    ElbertFluxTableFactory.numberOfAlphaSteps = 56;
	    ElbertFluxTableFactory.maxNumberOfMuEThSteps = 1;
	    muonFluxTable = new ElbertFluxTableFactory();
	    muonFluxTable.setElbertParameters(alpha);
	    hasReadTable = true;
	}

	double distance = distanceFromCosZenith(cosZenith);
	int icos = (int )((1.0-cosZenith)/stepCosZenith + 0.01);
	int iLogE = (int )((logEmuon - muonFluxTable.logEmuonMin)/
			   muonFluxTable.logEmuonStepSize + 0.01);
	if(0<=iLogE && iLogE < muonFluxTable.numberOfEmuon && 
	   icos <numberOfCosZenith){
	    double referenceFlux = dFMuDLogEReferece[icos][iLogE];
	    if(referenceFlux>0.0)
	    return muonFluxTable.getDFMuDLogE(distance, logEmuon)/referenceFlux;
	    else return 0.0;
	}else{
	    return 0.0;
	}

    }



    public static void main(String[] args) throws IOException{

	double cosZenith = 1.0;
	double alpha = 1.0;
	double alphaReference= 1.97;
	double muEThresholdReference =  1505.5;
	if(args.length==0){
	    System.out.println("Usage: RelativeElbertFluxMaker alpha (alphaRef muEThRef)");
	    System.exit(0);
        }else {
	    alpha = Double.valueOf(args[0]).doubleValue();
	    if(args.length==3){ // use the given parameter values for the reference flux
		alphaReference = Double.valueOf(args[1]).doubleValue();
		muEThresholdReference = Double.valueOf(args[2]).doubleValue();
	    }
	}

	RelativeElbertFluxTableMaker tableMaker = null;
	if(args.length!=3) {
	    tableMaker = new RelativeElbertFluxTableMaker();
	}else {
	    System.err.println("alphaRef(" + alphaReference + ") EthRef(" + muEThresholdReference +
			       ") ");
            RelativeElbertFluxTableMaker.table1DPath = "../data/muonBundleFluxTable1-5TeVthresholdNoCutoff/";
	    tableMaker = new RelativeElbertFluxTableMaker(alphaReference,muEThresholdReference);
	}
	tableMaker.setAlpha(alpha);

	while(cosZenith>0.0){
	    for(int iLogE=0;iLogE<120;iLogE++){
		double logEmuon =Particle.getLogEnergyMinimum()+
		    5.0*Particle.getDeltaLogEnergy()*(double )iLogE;
		double relativeFlux = tableMaker.getRelativeFlux(logEmuon,cosZenith);
		System.out.println(cosZenith + " " + logEmuon + " " + relativeFlux);
	    }
	    cosZenith -= tableMaker.stepCosZenith;
	}
    }



}
