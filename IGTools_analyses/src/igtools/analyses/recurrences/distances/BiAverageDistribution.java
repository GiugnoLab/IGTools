package igtools.analyses.recurrences.distances;

import java.util.Map;
import java.util.TreeMap;

import igtools.common.distributions.UpdatableDistribution;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.INELSA;

public class BiAverageDistribution {

	
	private INELSA nelsa;
	private Map<Integer, UpdatableDistribution> kMap;
	private BiRecurrenceDistancesExtractor extractor;
	private boolean keepDiagonal;
	private boolean onlySup;
	
	public BiAverageDistribution(INELSA nelsa, BiRecurrenceDistancesExtractor extractor, boolean keepDiagonal, boolean onlySup){
		this.nelsa = nelsa;
		this.extractor = extractor;
		this.keepDiagonal = keepDiagonal;
		this.onlySup = onlySup;
		kMap = new TreeMap<Integer, UpdatableDistribution>();
	}
	
	public Map<Double,Double> getCachedAvgDistribution(int k){
		UpdatableDistribution udistr = kMap.get(k);
		if(udistr == null){
			udistr = getAvgUpdatableDistribution(nelsa, k, extractor, keepDiagonal, onlySup);
			kMap.put(k, udistr);
		}
		return udistr.getDistribution();
	}
	
	public static UpdatableDistribution getAvgUpdatableDistribution(INELSA nelsa, int k, BiRecurrenceDistancesExtractor extr, boolean keepDiagonal, boolean onlySup){
		UpdatableDistribution udistr = new UpdatableDistribution();
		IELSAIterator it_a = nelsa.begin(k);
		while(it_a.next()){
			IELSAIterator it_b;
			
			if(onlySup)
				it_b = it_a.clone();
			else
				it_b = nelsa.begin(k);
			
			if(keepDiagonal){
				if(!onlySup){
					udistr.updateLazy( extr.recurrenceDistanceDistributionMap(it_a,it_b));
				}
				while(it_b.next()){
					udistr.updateLazy( extr.recurrenceDistanceDistributionMap(it_a,it_b));
				}
			}
			else{
				if(onlySup){
					while(it_b.next()){
						udistr.updateLazy( extr.recurrenceDistanceDistributionMap(it_a,it_b));
					}
				}
				else{
					while(it_b.next()){
						if(it_a.compare(it_b) != 0){
							udistr.updateLazy( extr.recurrenceDistanceDistributionMap(it_a,it_b));
						}
					}
				}
			}
			
		}
		udistr.flushLazy();
		return udistr;
	}
	
	public static Map<Double,Double> getAvgDistribution(INELSA nelsa, int k, BiRecurrenceDistancesExtractor extr, boolean keepDiagonal, boolean onlySup){
		return getAvgUpdatableDistribution(nelsa, k, extr,keepDiagonal, onlySup).getDistribution();
	}
}
