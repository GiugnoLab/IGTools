package igtools.cli.anagrams;

import java.util.Map;
import java.util.TreeMap;

//import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.common.util.Arrays;
import igtools.common.util.AvgSD;
import igtools.common.util.Maths;
import igtools.common.util.Vector4Int;
import igtools.dictionaries.elsa.NELSA;
import igtools.analyses.parikh.*;




/**
 * Calculate the avg and sd of the ratio  theoretical_possible_anagrams / own_anagrams of D_k(G).
 * Namely, for each word in D_k, it calculates the ratio between the theoretical number of anagrams and the real number of anagrams owned by G. 
 *  
 * @author vbonnici
 *
 */
public class KAvgSD {

	
	public static void usage(){
		System.out.println("Usage: cmd b3seq nelsa k");
	}
	
	public static void main(String[] args){
		String a_b3seq = "";
		String a_nelsa = "";
		int k = 0;
		
		try{
			a_b3seq = args[0];
			a_nelsa = args[1];
			
			k = Integer.parseInt(args[2]);
			if(k<1)
				throw new Exception();
			
		}catch(Exception e){
			usage();
			System.exit(1);
		}
		
		try{
			
			System.out.println("Loading sequence...");
			B3LLSequence b3seq = B3LLSequence.load(a_b3seq);
			System.out.println("done");
			System.out.println("Loading NELSA...");
			NELSA nelsa = new NELSA();
			nelsa.load(a_nelsa);
			System.out.println("done");
			nelsa.setSequence(b3seq);
			
			
			Map<Vector4Int, Integer> parikhClasses = new TreeMap<Vector4Int, Integer>();
			Map<Vector4Int, Integer> parikhMults = new TreeMap<Vector4Int, Integer>();
			
			ParikhClasses.getClassesNMultiplicities(nelsa, k, parikhClasses, parikhMults);
			
			
			
//			Vector4Int s_vect = new Vector4Int();
//			Vector4Int c_vect;
//			B3Nucleotide[] kmer = new B3Nucleotide[k];
//			Integer g, m;
			long t;
			
			double[] tmp = new double[parikhClasses.size()];
			int tt = 0;
			
			for(Map.Entry<Vector4Int, Integer> entry : parikhClasses.entrySet()){
				t = Maths.nofAnagrams(entry.getKey().values, k);
				
				System.out.println(entry.getKey()+" "+ entry.getValue()+" "+t+" "+ (((double)entry.getValue())/((double)t)) );
				
				tmp[tt] = (((double)entry.getValue())/((double)t));
				tt++;
			}
			
			int tot_k = nelsa.nof(k);
			
			System.out.println("nof k "+tot_k);
			System.out.println("nof p "+parikhClasses.size());
			
			double avg = AvgSD.avg(tmp);
			double sd = AvgSD.sd(tmp,  avg);
			System.out.println("# "+k+" "+tot_k+" "+(tot_k / Math.pow(4.0, k))+" "+parikhClasses.size()+" "+Arrays.min(tmp)+" "+Arrays.max(tmp)+" "+avg +" "+sd);
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
}
