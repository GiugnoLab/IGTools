package igtools.common.distributions.distance;


import java.util.Map;

public class KLDistance implements DistributionDistance{

	
	/**
	 * d1 and d2 do not need to have the same domain
	 */
	@Override
	public Double distance(Map<Double,Double> d1, Map<Double,Double> d2){
		Double dist = 0.0;
		
		//Double laste2 = Double.MIN_VALUE;
		
		//Double ee;
		Double e1,e2,v;
		for(Map.Entry<Double, Double> en : d1.entrySet()){
			e1 = en.getValue();
			e2 = d2.get(en.getKey());
			//ee=e2;
			
			if(e1 > 0.0){
				if(e2 == null || e2 == 0.0){
					//e2 = laste2;
					//break;
				}
				else{
					//laste2 = e2;
					
					//System.out.println(en.getKey() +"\t"+ e1 +"\t"+ ee +"\t"+ e2 + "\t"+dist);
					//System.out.println(en.getKey() +"\t"+ e1  +"\t"+ e2 + "\t"+dist);
					
					v = e1 * Math.log(e1 / e2);
					if((!v.isInfinite()) && (!v.isNaN()))
						dist += e1 * Math.log(e1 / e2);
				}
				
				//System.out.println(en.getKey() +"\t"+ e1 +"\t"+ ee +"\t"+ e2 + "\t"+dist);
				
				//dist += e1 * Math.log(e1 / e2);
			}
		}
		
		//return Math.abs(dist);
                return dist;
	}
	
	
	public static double pDistances(Map<Double,Double> d1, Map<Double,Double> d2, Map<Double,Double> dists){
		Double dist = 0.0;
		Double e1,e2,v;
		for(Map.Entry<Double, Double> en : d1.entrySet()){
			e1 = en.getValue();
			e2 = d2.get(en.getKey());
			
			if(e1 > 0.0){
				if(e2 == null || e2 == 0.0){
					dists.put(en.getKey(), 0.0);
				}
				else{
					v = e1 * Math.log(e1 / e2);
					if((!v.isInfinite()) && (!v.isNaN())){
						dist += v;
						dists.put(en.getKey(), v);
					}
					else{
						dists.put(en.getKey(), 0.0);
					}
				}
			}
			else{
				dists.put(en.getKey(), 0.0);
			}
		}
		//return Math.abs(dist);
                return dist;
	}
	
	
	public static class MaxKLDistance extends KLDistance{
		@Override
		public Double distance(Map<Double,Double> d1, Map<Double,Double> d2){
			Double a2b = super.distance(d1, d2);
			Double b2a = super.distance(d2, d1);
			
			//System.out.println("KLD: "+a2b+" "+b2a);
			
			if(a2b > b2a)	return a2b;
			return b2a;
		}
	}
	
}
