package numRecipes;

import numRecipes.*;

import java.io.*;
import java.util.*;

public class PoissonBinnedLHExtracter {


    public static PoissonBinnedLikelihoodCalculator getPoisonBinnedLikelihoodCalculator(Map plbMap){
	Iterator dataIterator = plbMap.entrySet().iterator();

	Map.Entry entryData = (Map.Entry )(dataIterator.next());
	PoissonBinnedLikelihoodCalculator cal = 
	    (PoissonBinnedLikelihoodCalculator)(entryData.getKey());
	return cal;
    }
}