package iceCube.uhe.particles;

import iceCube.uhe.particles.*;
import geometry.*;

import java.io.*;
import java.util.*;

/** I3Particle.class Demo program. */
public class I3ParticleStreamDemo {

    public static void main(String[] args) throws IOException{

	int doublet =1;
	int flavor = 1;
	String fileName = null;
	I3Particle iceTrack;

        if(args.length!=1){
            System.out.println("Usage: I3ParticleStreamDemo file-name");
	    System.exit(0);
        }else{
             fileName = args[0];
        }

	// Generate the particle class.
	iceTrack = 
	    new I3Particle(flavor, doublet, 1.0e8);

	iceTrack.generateLogEnergyMatrix();
	for(int i=0;i<I3Particle.getDimensionOfLogEnergyMatrix();i++){
	    double logEnergy = I3Particle.getLogEnergyMinimum()+
		I3Particle.getDeltaLogEnergy()*(double )i;
	    double energy = Math.pow(10.0,logEnergy);
	    double sigma=6.9542e-34*Math.pow(energy/1.0e6,0.402);
	    iceTrack.putLogEnergyMatrix(i,sigma);
	}

	// Sets the track geometry
	J3Vector cob = new J3Vector(0.0,0.0,0.0);
	J3UnitVector direction = new J3UnitVector(1.0,1.0,1.0);
	J3Line axis = new J3Line(cob,direction);
	iceTrack.setParticleAxisInIceCubeCoordinate(axis);
	iceTrack.transformParticleAxisToEarthCenterCoordinate();

	System.out.println("The I3Particle Name is " + 
			   iceTrack.particleName(iceTrack.getFlavor(), 
						 iceTrack.getDoublet()));

	System.out.println("The I3Particle Mass (" + 
			   iceTrack.getMass( )+ ") [GeV]");
	System.out.println("The I3Particle Energy (" + 
			   iceTrack.getEnergy( ) + ") [GeV] " + 
			   iceTrack.getLogEnergy( ));
	System.out.println("The I3Particle lifetime (" + 
			   iceTrack.getLifeTime( ) + ") [sec]");

	System.out.println("Matrix(" + 100 +")=" +
			   iceTrack.getLogEnergyMatrix(100));

	J3UnitVector trackDirection = iceTrack.getDirectionInIceCubeCoordinate();
	System.out.println("The I3Particle direction in ice3(" +
			   trackDirection.getX() + ", " + trackDirection.getY() +
			   ", " + trackDirection.getZ());

	trackDirection = iceTrack.getDirectionInEarthCenterCoordinate();
	System.out.println("The I3Particle direction in EarthCenter(" +
			   trackDirection.getX() + ", " + trackDirection.getY() +
			   ", " + trackDirection.getZ());

	J3Vector cob_center = iceTrack.getR0InEarthCenterCoordinate();
	System.out.println("The I3Particle location in EarthCenter(" +
			   cob_center.getX() + ", " + cob_center.getY() +
			   ", " + cob_center.getZ());

	/** Fills the IceCube data */
	iceTrack.getIceCubeData().setNpeATWD(1.0e2);
	iceTrack.getIceCubeData().setNpeFADC(1.1e2);
	iceTrack.getIceCubeData().setNDOMsATWD(300);
	iceTrack.getIceCubeData().setNDOMsFADC(270);

	/** Fills the GZK neutrino flux */
	String modelName = new String("GZK m=4 Zmax=4");
	iceTrack.setGZKNeutrinoFlux(5.0e-16,modelName);
	iceTrack.setGZKNeutrinoFlux(4.5e-16,"GZK m=2 Zmax=4");



	/** Generate 2nd I3Particle */
	I3Particle iceCascade = new I3Particle(0,1,1.0e7);
	cob = new J3Vector(0.0,0.0,0.0);
	direction = new J3UnitVector(0.0,0.0,1.0);
	iceCascade.generateLogEnergyMatrix();
	iceCascade.setParticleAxisInIceCubeCoordinate(new J3Line(cob,direction));
	iceCascade.transformParticleAxisToEarthCenterCoordinate();
	iceCascade.getIceCubeData().setNpeATWD(2.0e2);
	iceCascade.getIceCubeData().setNpeFADC(3.1e2);
	iceCascade.getIceCubeData().setNDOMsATWD(400);
	iceCascade.getIceCubeData().setNDOMsFADC(385);
	iceCascade.setAtmosphericMuonFlux(3.0e-14,"Corsika Model");
	iceCascade.setAtmosphericMuonFlux(3.5e-14,"Elbert Model");


	/** Output the serialized I3Particles */

	FileOutputStream out = new FileOutputStream(fileName);
	I3ParticleOutputStream.outputI3Particle(iceTrack, out);
	I3ParticleOutputStream.outputI3Particle(iceCascade, out);
	out.close( );

	System.out.println("The I3Particle output done.");



	/** Input the serialized I3Particle class */

	FileInputStream in = new FileInputStream(fileName);
	I3Particle iceParticle = null;
	while((iceParticle = I3ParticleInputStream.inputI3Particle(in))
	      !=null){


	    System.out.println("- The I3Particle Name is " + 
			       iceParticle.particleName(iceParticle.getFlavor(), 
							iceParticle.getDoublet()));

	    System.out.println("- The I3Particle Mass (" + 
			       iceParticle.getMass( )+ ") [GeV]");
	    System.out.println("- The I3Particle Energy (" + 
			       iceParticle.getEnergy( ) + ") [GeV] " + 
			       iceParticle.getLogEnergy( ));
	    System.out.println("- The I3Particle lifetime (" + 
			       iceParticle.getLifeTime( ) + ") [sec]");

	    System.out.println("- Matrix(" + 100 +")=" +
			       iceParticle.getLogEnergyMatrix(100));


	    J3UnitVector n = iceParticle.getDirectionInIceCubeCoordinate();
	    System.out.println("- The I3Particle direction in ice3(" + 
			       n.getX() + ", " + n.getY() + ", " + n.getZ());

	    n = iceParticle.getDirectionInEarthCenterCoordinate();
	    System.out.println("- The I3Particle direction in EarthCenter(" + 
			       n.getX() + ", " + n.getY() + ", " + n.getZ());

	    J3Vector r = iceParticle.getR0InEarthCenterCoordinate();
	    System.out.println("- The I3Particle location in EarthCenter(" + 
			       r.getX() + ", " + r.getY() + ", " + r.getZ());

	    /** Loop over all the GZK flux weights */
	    Iterator gzkIterator = iceParticle.iteratorOfGZKNeutrinoFlux();
	    while(gzkIterator.hasNext()){
		Map.Entry entry = (Map.Entry )(gzkIterator.next());
		Double flux = (Double )(entry.getValue());
		String name = (String )(entry.getKey());
		double dFdLogE = flux.doubleValue();
		System.out.println("- GZK flux = " + dFdLogE + 
				   " [/cm^2 sec sr] (" + name + ")");
	    }

	    /** Loop over all the Atmospheric Muon weights */
	    Iterator atmMuonIterator = iceParticle.iteratorOfAtmosphericMuonFlux();
	    while(atmMuonIterator.hasNext()){
		Map.Entry entry = (Map.Entry )(atmMuonIterator.next());
		Double flux = (Double )(entry.getValue());
		String name = (String )(entry.getKey());
		double dFdLogE = flux.doubleValue();
		System.out.println("- Atm Muon flux = " + dFdLogE + 
				   " [/cm^2 sec sr] (" + name + ")");
	    }

	    String thisModel = new String("GZK m=2 Zmax=4");
	    System.out.println("- this GZK flux = " + 
			       iceParticle.getGZKNeutrinoFlux(thisModel) + 
			       " [/cm^2 sec sr] (" + thisModel + ")");
	    iceParticle.removeGZKNeutrinoFlux(thisModel);
	    /** Loop over all the GZK flux weights again */
	    gzkIterator = iceParticle.iteratorOfGZKNeutrinoFlux();
	    while(gzkIterator.hasNext()){
		Map.Entry entry = (Map.Entry )(gzkIterator.next());
		Double flux = (Double )(entry.getValue());
		String name = (String )(entry.getKey());
		double dFdLogE = flux.doubleValue();
		System.out.println(" -- GZK flux = " + dFdLogE + 
			       " [/cm^2 sec sr] (" + name + ")");
	    }

	    String thatModel = new String("Elbert Model");
	    System.out.println("- this Atm Muon flux = " + 
			       iceParticle.getAtmosphericMuonFlux(thatModel) + 
			       " [/cm^2 sec sr] (" + thatModel + ")");

	    System.out.println("- Npe by ATWD = " + 
			       iceParticle.getIceCubeData().getNpeATWD());
	    System.out.println("- Log Npe by ATWD = " + 
			       iceParticle.getIceCubeData().getLogNpeATWD());
	    System.out.println("- NDOMs with  ATWD = " + 
			       iceParticle.getIceCubeData().getNDOMsATWD());
	    System.out.println("- NDOMs with FADC = " + 
			       iceParticle.getIceCubeData().getNDOMsFADC());
	}

	in.close();

	System.out.println("Demo program ended --");
    }
}
