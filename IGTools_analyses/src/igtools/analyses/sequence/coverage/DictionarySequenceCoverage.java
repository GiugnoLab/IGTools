package igtools.analyses.sequence.coverage;

import java.util.Map;

public class DictionarySequenceCoverage {
	public static void coverageDistribution(double[] coverage, Map<Double,Double> dist){
		Double v;
		for(int i=0; i<coverage.length; i++){
			v = dist.get(coverage[i]);
			if(v == null){
				dist.put(coverage[i], 1.0);
			}
			else{
				dist.put(coverage[i], v + 1.0);
			}
		}
	}
}
