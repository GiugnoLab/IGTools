package igtools.analyses;

import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.INELSA;



/**
 * Dictionary trends: size, max/min multiplicity, hapax count, repeats (repeatenes avg/sd)
 * 
 * 
 * 
 * @author vbonnici
 *
 */
public class Trends {

	
	public static int[] DictionarySizes(INELSA nelsa){
		return DictionarySizes(nelsa, Inform.Mrl(nelsa));
	}
	public static int[] DictionarySizes(INELSA nelsa, int toK){
		int[] trends = new int[toK+1];
		trends[0] = 0;
		for(int k=1; k<=toK; k++){
			int n = 0;
			IELSAIterator it = nelsa.begin(k);
			while(it.next()){
				n += it.multiplicity();
			}
			trends[k] = n;
		}
		return trends;
	}
	
	
	
	
	public static class TrendInfoRet{
		public int d_size = 0;
		public int m_size = 0;
		public int max_multiplicity = 0;
		public int min_multiplicity = Integer.MAX_VALUE;
		public int nof_hapaxes = 0;
		public int nof_repeates = 0;
		public double avg_repeatenes = 0;
		public double sd_repeatenes = 0;
		
		public void reset(){
			d_size = 0;
			m_size = 0;
			max_multiplicity = 0;
			min_multiplicity = Integer.MAX_VALUE;
			nof_hapaxes = 0;
			nof_repeates = 0;
			avg_repeatenes = 0;
			sd_repeatenes = 0;
		}
	}
	public static void trendInfo(INELSA nelsa, int k, TrendInfoRet retInfo){
		//retInfo.reset();
		IELSAIterator it = nelsa.begin(k);
		
		int m;
		while(it.next()){
			retInfo.d_size++;
			m = it.multiplicity();
			retInfo.m_size += m;
			if(m==1){
				retInfo.nof_hapaxes++;
			}
			else{
				retInfo.nof_repeates++;
				retInfo.avg_repeatenes += m;
			}
			
			if(m > retInfo.max_multiplicity)
				retInfo.max_multiplicity = m;
			if(m <retInfo.min_multiplicity)
				retInfo.min_multiplicity = m;
		}
		
		retInfo.avg_repeatenes /= (double)retInfo.nof_repeates;
		
		double dm;
		it = nelsa.begin(k);
		while(it.next()){
			dm = (double)it.multiplicity();
			if(dm!=1){
				retInfo.sd_repeatenes += (dm - retInfo.avg_repeatenes) * (dm - retInfo.avg_repeatenes); 
			}
		}
		retInfo.sd_repeatenes = Math.sqrt( retInfo.sd_repeatenes / (double)retInfo.nof_repeates );
		
	}
	
	
	
}
