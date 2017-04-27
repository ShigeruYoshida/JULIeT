package iceCube.uhe.muonModel;

/** 

Abstract class for all the particle flux subclasses
in the package muonModel. Implementations of the flux models
will be made in each of subclasses.

*/


import java.util.*;

abstract class ParticleFlux {

    /** Calculate dF/dLogE [/cm^2 sec sr] as a function
        of log E [GeV] and cosine of zenith angle */
    abstract double getDFDLogE(double logEnergy, double cosTheta);

    /** Calculate energy flux E^2 dF/dE [GeV/cm^2 sec sr] */
    abstract double  getEFlux(double logEnergy, double cosTheta);
}
