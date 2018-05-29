package igtools.common.distributions.distance;


import java.util.Map;

public class KSDistance implements DistributionDistance{

	
	/**
	 * d1 and d2 do not need to have the same domain
	 */
	@Override
	public Double distance(Map<Double, Double> d1, Map<Double, Double> d2) {
		Double ret = 0.0;
		
		Double d1_cum = 0.0;
		Double d2_cum = 0.0;
		
		Double[] x1 = new Double[d1.size()];
		d1.keySet().toArray(x1);
		
		Double[] x2 = new Double[d2.size()];
		d2.keySet().toArray(x2);
		
		
		int i1=0, i2=0;
		
		while(i1<x1.length && i2<x2.length){
			//System.out.println("["+i1+","+x1[i1]+"] ["+i2+","+x2[i2]+"] v("+d1.get(x1[i1])+","+d2.get(x2[i2])+")  cum("+d1_cum+","+d2_cum+")");
			
			d1_cum += d1.get(x1[i1]);
			
			if(x1[i1] == x2[i2]){
				d2_cum += d2.get(x2[i2]);
				i2++;
			}
			else if(x1[i1] > x2[i2]){
				while(i2<x2.length && x1[i1] >= x2[i2]){
					d2_cum += d2.get(x2[i2]);
					i2++;
				}
			}
			
			i1++;
			
			if(Math.abs(d1_cum - d2_cum) > ret)
				ret += Math.abs(d1_cum - d2_cum);
		}
		
		while(i1<x1.length){
			d1_cum += d1.get(x1[i1]);
			i1++;
		}
		while(i2<x2.length){
			d2_cum += d2.get(x2[i2]);
			i2++;
		}
		
		if(Math.abs(d1_cum - d2_cum) > ret)
			ret += Math.abs(d1_cum - d2_cum);
		
		
		return ret;
	}

}
