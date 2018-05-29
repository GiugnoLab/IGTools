package igtools.analyses.recurrences.distances;

import igtools.dictionaries.elsa.IELSAIterator;

import java.util.Map;
import java.util.TreeMap;

public abstract class BiRecurrenceDistancesExtractor {

	
	/**
	 * 
	 * @param it
	 * @return the complete list of recurrence distances
	 */
	public abstract int[] recurrencesDistances(IELSAIterator it_a, IELSAIterator it_b);
	
	/**
	 * Recurrence distance distribution
	 * 
	 * ret[i][0] = recurrence distance
	 * ret[i][1] = count
	 * 
	 * @param it
	 * @return
	 */
	public double[][] recurrenceDistanceDistribution(IELSAIterator it_a, IELSAIterator it_b){
		int[] dists = recurrencesDistances(it_a, it_b);
		
//		System.out.println("nof recurrences: "+dists.length);
		
		Map<Double, Double> distMap = new TreeMap<Double, Double>();
		Double v;
		for(int i=0; i<dists.length; i++){
			v = distMap.get((double)dists[i]);
			if(v == null) distMap.put((double)dists[i], 1.0);
			else distMap.put((double)dists[i], v+1.0); 
			//distMap.put(new Double(dists[i]), v + 1.0);
		}
		
		int j = 0;
		double[][] distr = new double[distMap.size()][2];
		for(Map.Entry<Double, Double> entry : distMap.entrySet()){
			//System.out.println("("+entry.getKey()+","+entry.getValue()+")");
			distr[j][0] = entry.getKey();
			distr[j][1] = entry.getValue();
			j++;
		}
		
		return distr;
	}
	
	/**
	 * Recurrence distance distribution
	 * 
	 * @param it
	 * @return
	 */
	public Map<Double,Double> recurrenceDistanceDistributionMap(IELSAIterator it_a, IELSAIterator it_b){
		int[] dists = recurrencesDistances(it_a, it_b);
		
//		System.out.println("nof recurrences: "+dists.length);
		
		Map<Double, Double> distMap = new TreeMap<Double, Double>();
		Double v;
		for(int i=0; i<dists.length; i++){
			v = distMap.get((double)dists[i]);
			if(v == null) distMap.put((double)dists[i], 1.0);
			else distMap.put((double)dists[i], v+1.0); 
			//distMap.put(new Double(dists[i]), v + 1.0);
		}
		return distMap;
	}
	
	
	
	
	
	
	
	public static BiRecurrenceDistancesExtractor factory(boolean minimal, boolean overlapping, boolean skipNgaps){
//		System.out.println("bifactory "+overlapping+" "+skipNgaps);
		if(minimal){
			return BiMinimalRecurrenceDistancesExtractor.factory(overlapping, skipNgaps);
		}
		else{
			return LimitedBiGlobalRecurrenceDistancesExtractor.factory(overlapping, skipNgaps);
		}
	}
}
