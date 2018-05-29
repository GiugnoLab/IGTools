package igtools.analyses.recurrences.carpena;

import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;

import igtools.common.ds.StaticB2Tree;
import igtools.dictionaries.elsa.ExtensionNELSAIterator;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;




/**
 * Original Carpena2009 algorithm and measurements
 * 
 * @author vbonnici
 *
 */
public class Carpena2009 {

	protected NELSA nelsa;
	protected int[] sa;
	
	
	public Carpena2009(NELSA nelsa){
		this.nelsa = nelsa;
		this.sa = nelsa.sa();
	}

	
	public double[] get_Cs(int k){
		int nof = 0;
		IELSAIterator bit = nelsa.begin(k);
		while(bit.next()){
			if(bit.multiplicity() > 1)
				nof++;
		}
		
		double[] cs = new double[nof];
		bit = nelsa.begin(k);
		int i=0;
		while(bit.next()){
			if(bit.multiplicity() > 1){
				cs[i] = C(bit);
				i++;
			}
		}
		return cs;
	}
	
	
	public double[] get_sigma_nors(int k){
		int nof = 0;
		IELSAIterator bit = nelsa.begin(k);
		while(bit.next()){
			if(bit.multiplicity() > 1)
				nof++;
		}
		
		double[] cs = new double[nof];
		bit = nelsa.begin(k);
		int i=0;
		while(bit.next()){
			if(bit.multiplicity() > 1){
				cs[i] = sigma_nor(bit);
				i++;
			}
		}
		return cs;
	}
	
	
	public double[] get_sigmas(int k){
		int nof = 0;
		IELSAIterator bit = nelsa.begin(k);
		while(bit.next()){
			if(bit.multiplicity() > 1)
				nof++;
		}
		
		double[] cs = new double[nof];
		bit = nelsa.begin(k);
		int i=0;
		while(bit.next()){
			if(bit.multiplicity() > 1){
				cs[i] = sigma(bit);
				i++;
			}
		}
		return cs;
	}
	
	
	public double get_C0_percetile(double percentile, int k){
		double[] cs = get_Cs(k);
		Arrays.sort(cs);
		
		int n = cs.length;
		for(int i=cs.length-1; i>=0; i--){
			if(!Double.isNaN(cs[i]))
				break;
			n--;
		}
		
		int p = (int)((percentile * ((double)n)) + (0.5));
		
//		for(int i=(cs.length-p-1); i<cs.length; i++)
//			System.out.print(cs[i]+" ");
//		System.out.println();
//		
//		System.out.println("nof "+cs.length+"; p "+p);
		
		//return cs[cs.length -p -1];
		return cs[n -p -1];
	}
	
	
	
	
	
	
	public int select(int initial_k, double C0, boolean usePercentile, double C_delta, StaticB2Tree words){
		int count = 0;
		IELSAIterator  bit = nelsa.begin(initial_k);
		while(bit.next()){
			double C = C(bit);
			if((!usePercentile) || (C >= C0 - (C0 * C_delta))){
				count += recursive_select(initial_k, bit, C0, 0.0, words);
			}
		}
		return count;
	}
	
	protected int recursive_select(int initial_k, IELSAIterator it, double C0, double C_delta, StaticB2Tree words){
		int count = 0;
		double[] cs = {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
		ExtensionNELSAIterator[] cs_its = {null, null, null, null};
		int i;
		double max = Double.NEGATIVE_INFINITY;
		
		ExtensionNELSAIterator eit = new ExtensionNELSAIterator(nelsa, it);
		i=0;
		while(eit.next()){
			cs_its[i] = (ExtensionNELSAIterator)eit.clone();
			cs[i] = C(cs_its[i]);
			if(cs[i] > max){
				max = cs[i];
			}
			i++;
		}
		
		
//		System.out.print("C: ");
//		for(i=0; i<cs.length; i++){
//			if(cs[i] != Double.MIN_VALUE){
//				System.out.print(B3Nucleotide.toString(cs_its[i].kmer())+" ");
//			}
//		}
//		System.out.println();
//		System.out.print("C: ");
//		for(i=0; i<cs.length; i++){
//			if(cs[i] != Double.MIN_VALUE){
//				System.out.print(cs[i]+" ");
//			}
//		}
//		System.out.println();
//		System.out.println();
		
		//i=0   there are no extensions
		//max<=C0  do not extend words having C<=C0
		if(i>0 && max>=C0){
			for(i=0; i<cs.length; i++){
				if(cs[i] >= max){
					if(cs_its[i].multiplicity()==1){
						words.add(nelsa.b3seq(), sa[cs_its[i].istart()], sa[cs_its[i].istart()] + it.k()+1);
						count++;
					}
					else{
						count += recursive_select(initial_k, cs_its[i], C0, 0.0, words);
					}
				}
			}
		}
		else if(it.k() != initial_k){
			//set selected
			//System.out.println("S: "+ B3Nucleotide.toString(it.kmer()) +"\n" );
			words.add(nelsa.b3seq(), sa[it.istart()], sa[it.istart()] + it.k());
			count++;
		}
		return count;
	}	
	
	
	
	private Map<Integer, Integer> nof = new TreeMap<Integer, Integer>();
	private int nof(int k){
		Integer n = nof.get(k);
		if(n == null){
			n = nelsa.nof_mults(k);
			nof.put(k, n);
		}
		return n;
	}
	
	
	
//	private int[] ordered_positions(IELSAIterator it){
//		int[] pos = Arrays.copyOfRange(nelsa.sa(), it.istart(), it.iend());
//		Arrays.sort(pos);
//		return pos;
//	}
//	private int[] proper_codistances(IELSAIterator it){
//		//int[] pos = ordered_positions(it);
//		int[] pos = it.sortedPositions();
//		
//		int k = it.k();
//		int count = 0;
//		int last_i = 0;
//		for(int i=1; i<pos.length; i++){
//			if(pos[i] - pos[last_i] >= k){
//				count++;
//				last_i = i;
//			}
//		}
//		
//		int[] codists = new int[count];
//		last_i = 0;
//		int j = 0;
//		for(int i=1; i<pos.length; i++){
//			if(pos[i] - pos[last_i] >= k){
//				codists[j] = pos[i] - pos[last_i];
//				j++;
//				last_i = i;
//			}
//		}
//		
//		return codists;
//	}
	
	
	private double avg(int[] a){
		double avg = 0;
		for(int i=0; i<a.length; i++){
			avg += a[i];
		}
		avg /= (double) a.length;
		return avg;
	}
	
	private double sd(int[] a, double avg){
		double sd = 0;
		for(int i=0; i<a.length; i++){
			sd += (((double)a[i]) - avg) * (((double)a[i]) - avg); 
		}
		sd /= (double)a.length;
		sd = Math.sqrt(sd);
		return sd;
	}
	
	
	
	
	public double sigma(int[] codists){
//		System.out.println(codists.length);
		double avg = avg(codists);
		//double sigma = sd(codists, avg) / avg;
		return sd(codists, avg) / avg;
	}
	
	
	public double sigma(IELSAIterator it){
		//int[] codists = proper_codistances(it);
		int[] codists = DistalRecurrence.proper_codistances(it);
		return sigma(codists);
	}
	
	
	
	public double prob(IELSAIterator it){
		return ((double)it.multiplicity()) / ((double)nof(it.k()));
	}
	public double sigma_nor(IELSAIterator it){
//		System.out.println(sigma(it) +" / "+ ( Math.sqrt(1 - prob(it)))          );
		return sigma(it) / ( Math.sqrt(1 - prob(it)));
	}
	
//	private double C(double sigma_nor, int n){
//		return (sigma_nor - sigma_nor_Poisson(n)) / sd_sigma_nor_Poisson(n);
//	}
	public double C(IELSAIterator it){
//		System.out.println(sigma_nor(it) +" - "+ sigma_nor_Poisson(it.multiplicity()) +" / "+ sd_sigma_nor_Poisson(it.multiplicity()));
		return (sigma_nor(it) - sigma_nor_Poisson(it.multiplicity())) / sd_sigma_nor_Poisson(it.multiplicity());
	}
	
	
	private double sigma_slow(int[] a){
		//sigma is the standard deviation of the normalized (by average) distances
		double i_avg = avg(a);
		
		double avg = 0;
		for(int i=0; i<a.length; i++){
			avg += ((double)a[i]) / i_avg;
		}
		avg /= (double) a.length;
		
		double sigma = 0;
		for(int i=0; i<a.length; i++){
			sigma += ((((double)a[i]) / i_avg) - avg) * ((((double)a[i]) / i_avg) - avg); 
		}
		sigma /= (double)a.length;
		sigma = Math.sqrt(sigma);
		
		return sigma;
	}
	
	
	private double sigma_slow(IELSAIterator it){
		//int[] codists = proper_codistances(it);
		int[] codists = DistalRecurrence.proper_codistances(it);
		return sigma_slow(codists);
	}
	
	
	
	
	
	public static double sigma_nor_Poisson(double n){
		return (2*n -1) / (2*n +2);
	}
	
	protected static double sd_sigma_nor_Poisson(double n){
		return 1 / (Math.sqrt(n) * (1+ 2.8 * Math.pow(n, -0.865)));
	}
	
}
