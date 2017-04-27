package iceCube.uhe.particles;

import iceCube.uhe.particles.*;

/** ParticleArray.class Demo program. */
public class ParticleArrayDemo {

    public static void main(String[] args){

	ParticleArray uheParticleArray = new ParticleArray( );
	Particle uheParticle;

	for(int flavor=0;flavor<Particle.NumberOfFlavor;flavor++){
	    for(int doublet=0;doublet<Particle.NumberOfDoublet;doublet++){
		for(int iLogE=0;iLogE<Particle.getDimensionOfLogEnergyMatrix();
		    iLogE++){

		    uheParticle = 
			uheParticleArray.getParticle(flavor,doublet,iLogE);
		    for(int jLogE=0;jLogE<=iLogE;jLogE++){
			double logEnergy = Particle.getLogEnergyMinimum()+
			Particle.getDeltaLogEnergy()*(double )jLogE;
			double energy = Math.pow(10.0,logEnergy);
			double sigma=6.9542e-34*Math.pow(energy/1.0e6,0.402);
			uheParticle.putLogEnergyMatrix(jLogE,sigma);
		    }
		}

	        uheParticle = 
		uheParticleArray.getParticle(flavor,doublet,0);
		System.out.println("The Particle flavor (" + 
				   uheParticle.getFlavor( ) + ")");
		System.out.println("The Particle doublet (" + 
				   uheParticle.getDoublet( ) + ")");
		System.out.println("The Particle Mass (" + 
				   uheParticle.getMass( ) + ") [GeV]");
		System.out.println("The Particle Energy (" + 
				   uheParticle.getEnergy( ) + ") [GeV] " 
				   + uheParticle.getLogEnergy( ));
		System.out.println("The Particle lifetime (" + 
				   uheParticle.getLifeTime( ) + ") [sec]");
		System.out.println("The Particle Name is " + 
	        uheParticle.particleName(uheParticle.getFlavor(), 
					 uheParticle.getDoublet()));

		for(int iLogE=0;iLogE<Particle.getDimensionOfLogEnergyMatrix();
		    iLogE+=100){
		    uheParticle = 
			uheParticleArray.getParticle(flavor,doublet,iLogE);
		    double logEnergy = uheParticle.getLogEnergy( );
		    System.out.println("logE("+logEnergy+")Matrix("+iLogE+")=" +
				       uheParticle.getLogEnergyMatrix(iLogE));
		}
	    }
	}

    }

}

