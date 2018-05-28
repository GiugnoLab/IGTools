package igtools.analyses;

import java.util.Map;
import java.util.TreeMap;

import igtools.dictionaries.elsa.IELSA;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.INELSA;
import igtools.dictionaries.elsa.NELSA;


/**
 * Classic InfoGenomic measurements 
 * 
 * @author vbonnici
 *
 */
public class Inform {

	
	/**
	 * 
	 * @param nelsa
	 * @return Maximal repeat length
	 */
	public static int Mrl(INELSA nelsa){
		int mrl = 0; 
		int[] lcp = nelsa.lcp();
		int[] ns = nelsa.ns();
		for(int i=1; i<lcp.length; i++){
			if(ns[i] > lcp[i]){
				if(lcp[i] > mrl)
					mrl = lcp[i];
			}
			else{
				if(ns[i] > mrl){
					mrl = ns[i] - 1;
				}
			}
//			if(lcp[i] > mrl && ns[i] > lcp[i]) //new nelsa version
//				mrl = lcp[i];
		}
		return mrl;
	}
	
	
	/**
	 * 
	 * @param nelsa
	 * @return minimal hapax length
	 */
	public static int mhl(INELSA nelsa){
		int mhl = 1;
		while(true){
			IELSAIterator it = nelsa.begin(mhl);
			if(it.next()){
				if(it.multiplicity() == 1)
					return mhl;
				while(it.next()){
					if(it.multiplicity() == 1)
						return mhl;
				}
			}
			else{
				return -1;
			}
			mhl++;
		}
	}
	
	
	
	/**
	 * 
	 * @param nelsa
	 * @return the last k for which the sequence is k-complete
	 */
	public static int kCompleteness(IELSA nelsa){
		int k = 0;
		double k4 = 1;
		
		while(true){
			k++;
			k4 *= 4;
			int nof = nelsa.nof(k);
			if(nof != k4)
				break;
		}
		return k-1;
	}	
	
	
	
	public static int FofanofM(int length){
		int k = 0;
		double dlength = (double)length;
		double eq;
		while(true){
			eq = (dlength * dlength) + ((1-k)*(dlength + Math.pow(4,k)));
			if(eq < 0)
				break;
			else
				k++;
		}
		return k;
	}
	
	
	/**
	 * 
	 * @param nelsa
	 * @param k
	 * @return  empirical k-entropy
	 */
	public static double entropy(IELSA nelsa, int k){
		double ff = 0;
		IELSAIterator it = nelsa.begin(k);
		while(it.next()){
			ff += (double)it.multiplicity();
		}
		double entr = 0;
		double p;
		it = nelsa.begin(k);
		while(it.next()){
			p = ((double)it.multiplicity()) / ff;
			entr += p * Math.log(p); 
		}
		entr *= -1;
		entr /= Math.log(2);
		return entr;
	}
	
	
	
	public static int[] DSizes(INELSA nelsa){
		return DSizes(nelsa, Inform.Mrl(nelsa));
	}
	/**
	 * 
	 * @param nelsa
	 * @param toK
	 * @return the dictionary sizes for k in ]0,toK[
	 */
	public static int[] DSizes(INELSA nelsa, int toK){
		int[] trends = new int[toK+1];
		trends[0] = 0;
		for(int k=1; k<=toK; k++){
			trends[k] = nelsa.nof(k);
		}
		return trends;
	}
	
	
	
	public static int[] MSizes(INELSA nelsa){
		return MSizes(nelsa, Inform.Mrl(nelsa));
	}
	/**
	 * 
	 * @param nelsa
	 * @param toK
	 * @return the total multiplicity for k in ]0,toK[
	 */
	public static int[] MSizes(INELSA nelsa, int toK){
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
	
	
	
	public static TreeMap<Integer,Integer> multiplicityDistribution(INELSA nelsa){
		int Mrl = Inform.Mrl(nelsa);
		TreeMap<Integer,Integer> distr = new TreeMap<Integer,Integer>();
		for(int k=1; k<=Mrl; k++){
			multiplicityDistribution(nelsa, k, distr);
		}
		return distr;
	}
	public static TreeMap<Integer,Integer> multiplicityDistribution(INELSA nelsa, int k){
		TreeMap<Integer,Integer> distr = new TreeMap<Integer,Integer>();
		multiplicityDistribution(nelsa, k,distr);
		return distr;
	}
	public static void multiplicityDistribution(INELSA nelsa, int k, TreeMap<Integer,Integer> distr){
		Integer cm;
		int mult;
		IELSAIterator it = nelsa.begin(k);
		while(it.next()){
			mult = it.multiplicity();
			cm = distr.get(mult);
			if(cm == null){
				distr.put(mult, 1);
			}
			else{
				distr.put(mult, cm+1);
			}
		}
	}
	
	
}
