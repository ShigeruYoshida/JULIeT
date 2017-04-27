package iceCube.uhe.geometry;

import geometry.*;

/**
    This class defines geometry of IceCube detecter.
    To be simple, we define cubic volume.
    The coordinate is defined by IceCube local coordinate IceCubeCoordinate.java
*/

public class IceCubeVolume extends Volume{


    /**  Constructor. */
    public IceCubeVolume(){

        super(1.0e5);
    }
}        
