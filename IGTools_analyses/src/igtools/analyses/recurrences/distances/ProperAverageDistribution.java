package igtools.analyses.recurrences.distances;

import java.util.Map;
import java.util.TreeMap;

import igtools.common.distributions.UpdatableDistribution;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.INELSA;

public class ProperAverageDistribution {

	
	private INELSA nelsa;
	private Map<Integer, UpdatableDistribution> kMap;
	ProperRecurrenceDistancesExtractor extractor;
	
	public ProperAverageDistribution(INELSA nelsa, ProperRecurrenceDistancesExtractor extractor){
		this.nelsa = nelsa;
		this.extractor = extractor;
		kMap = new TreeMap<Integer, UpdatableDistribution>();
	}
	
	public Map<Double,Double> getCachedAvgDistribution(int k){
		UpdatableDistribution udistr = kMap.get(k);
		if(udistr == null){
			udistr = getAvgUpdatableDistribution(nelsa, k, extractor);
			kMap.put(k, udistr);
		}
		return udistr.getDistribution();
	}
	
	public static UpdatableDistribution getAvgUpdatableDistribution(INELSA nelsa, int k, ProperRecurrenceDistancesExtractor extr){
		UpdatableDistribution udistr = new UpdatableDistribution();
		IELSAIterator it = nelsa.begin(k);
		while(it.next()){
			udistr.updateLazy( extr.recurrenceDistanceDistributionMap(it));
		}
		udistr.flushLazy();
		return udistr;
	}
	
	public static Map<Double,Double> getAvgDistribution(INELSA nelsa, int k, ProperRecurrenceDistancesExtractor extr){
		return getAvgUpdatableDistribution(nelsa, k, extr).getDistribution();
	}
}
