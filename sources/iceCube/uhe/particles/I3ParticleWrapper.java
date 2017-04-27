package iceCube.uhe.particles;

import iceCube.uhe.particles.*;
import geometry.*;

public class I3ParticleWrapper {

    static public J3Line getParticleAxisInIceCubeCoordinate(I3Particle particle){

	if(particle.isMCTruth()){
	    return  particle.particleAxis_J3Line_ice3_MC;
	}else{
	    return  particle.particleAxis_J3Line_ice3_reco;
	}
    }
}