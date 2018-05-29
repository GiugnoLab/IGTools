package igtools.analyses.recurrences.carpena;

import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;



/**
 * Carpena2009 methods factory
 * 
 * @author vbonnici
 *
 */
public class Carpena2009Factory {

	
	public static Carpena2009 factory(NELSA nelsa, boolean byParent, boolean noSigmaNor){
		
		if(byParent){
			if(noSigmaNor){
				return new Carpena2009ByParentC(nelsa){
					@Override
					public double C(IELSAIterator it){
						return (sigma(it) - sigma_nor_Poisson(it.multiplicity())) / sd_sigma_nor_Poisson(it.multiplicity());
					}
				};
			}
			else{
				return new Carpena2009ByParentC(nelsa);
			}
		}
		else{
			if(noSigmaNor){
				return new Carpena2009(nelsa){
					@Override
					public double C(IELSAIterator it){
						return (sigma(it) - sigma_nor_Poisson(it.multiplicity())) / sd_sigma_nor_Poisson(it.multiplicity());
					}
				};
			}
			else{
				return new Carpena2009(nelsa);
			}
		}
	}
}
