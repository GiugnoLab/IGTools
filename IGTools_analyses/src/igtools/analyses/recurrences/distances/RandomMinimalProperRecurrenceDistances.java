package igtools.analyses.recurrences.distances;

import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;

import igtools.dictionaries.elsa.IELSAIterator;


/**
 * Random recurrence distance generator
 * starting from a real recurrence distance distribution of a word. 
 * 
 * 
 * @author vbonnici
 *
 */
public class RandomMinimalProperRecurrenceDistances {

	
	
	/**
	 * randomize a complete list of recurrence distance of  a real word 
	 * 
	 * @param it
	 * @return
	 */
	private static int[] randomRecurrencesDistances(IELSAIterator it){
		int[] poss = it.sortedPositions();
		
		int min = poss[0];
		int max = poss[poss.length -1];
		double mm = (double)(max -min);
		
		System.out.println("mm "+mm);
		
		for(int i=0; i<poss.length; i+=2){
			int npos = min + (int)Math.floor((Math.random() * mm));
			//System.out.println(npos);
			poss[i] = npos;
		}
		
		
		for(int i=1; i<poss.length; i+=2){
			int npos = min + (int)Math.floor((Math.random() * mm));
			//System.out.println(npos);
			poss[i] = npos;
		}
		
		Arrays.sort(poss);
//		for(int i=0; i<poss.length; i++){
//			System.out.print(poss[i]+" ");
//		}
//		System.out.println();
		
		
		int k = it.k();
		
		int count = 0;

		
		for(int i=1; i<poss.length; i++){
			if(poss[i] - poss[i-1] >=k){
				count++;
			}
		}
		
		System.out.println("count "+count);
		
		int[] recs = new int[count];
		if(count > 0){
			int j = 0;
			for(int i=1; i<poss.length; i++){
				if(poss[i] - poss[i-1] >=k){
					recs[j] = poss[i] - poss[i-1];
					if(recs[j] == 0) recs[j] = 1;
					j++;
				}
			}
		}
		
		return recs;
	}
	
	
	/**
	 * Randomized ecurrence distance distribution
	 * 
	 * ret[i][0] = recurrence distance
	 * ret[i][1] = count
	 * 
	 * @param it
	 * @return
	 */
	public static double[][] extractDistribution(IELSAIterator it){
		int[] distances = randomRecurrencesDistances(it); 
		System.out.println("ds: "+distances.length);
		
		
		Map<Double, Double> distMap = new TreeMap<Double, Double>();
		Double v;
		for(int i=0; i<distances.length; i++){
			v = distMap.get((double)distances[i]);
			if(v == null) distMap.put((double)distances[i], 1.0);
			else distMap.put((double)distances[i], v+1.0); 
		}
		
		System.out.println("dm: "+distMap.size());
		
		
		int j = 0;
		double[][] distr = new double[distMap.size()][2];
		for(Map.Entry<Double, Double> entry : distMap.entrySet()){
			distr[j][0] = entry.getKey();
			distr[j][1] = entry.getValue();
			j++;
		}
		
		return distr;
	}
	
	
	/**
	 * Randomized recurrence distance distribution
	 *  
	 * @param it
	 * @return
	 */
	public static TreeMap<Double,Double> extractDistributionMap(IELSAIterator it){
		TreeMap<Double,Double> dmap = new TreeMap<Double,Double>();
		
		double[][] distr = extractDistribution(it);
		System.out.println(distr.length);
		
		for(int i=0; i<distr.length; i++){
			dmap.put(distr[i][0], distr[i][1]);
		}
		
		return dmap;
	}
}
