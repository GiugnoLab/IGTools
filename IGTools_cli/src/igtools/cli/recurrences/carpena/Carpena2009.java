package igtools.cli.recurrences.carpena;

import java.util.Map;
import java.util.TreeMap;

import igtools.analyses.recurrences.carpena.Carpena2009Factory;
import igtools.common.ds.StaticB2Tree;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.NELSA;


/**
 * Run the Carpena2009 algorithm for clustered words extraction
 * 
 * @author vbonnici
 *
 */
public class Carpena2009 {
	private static void usage(){
		System.out.println("Usage: cmd sequence.3bit nelsa.nelsa [original|parent] [sigma|sigma_nor] C_delta percetile k... ");
	}
	
	
	public static void main(String[] args){
		String a_iseq = "";
		String a_inelsa = "";
		String a_type = "";
		String s_type = "";
		double c_delta = 0.0;
		double a_percentile = 0.0;
		int[] ks = null;
		
		
		boolean original = false;
		boolean sigma_nor = false;
		
		
		try{
			a_iseq = args[0];
			a_inelsa = args[1];
			
			a_type = args[2];
			if(a_type.compareTo("original")==0)
				original = true;
			else if(a_type.compareTo("parent")==0)
				original = false;
			else
				throw new Exception();
			
			s_type = args[3];
			if(s_type.compareTo("sigma")==0)
				sigma_nor = false;
			else if(s_type.compareTo("sigma_nor")==0)
				sigma_nor = true;
			else
				throw new Exception();
			
			c_delta = Double.parseDouble(args[4]);
			
			a_percentile = Double.parseDouble(args[5]);
			
			if(c_delta < 0 || c_delta>1)
				throw new Exception("c_delta must be in ]0,1[");
			
			if(a_percentile < 0 || a_percentile>1)
				throw new Exception("a_percentile must be in ]0,1[");
			
			if(args.length <= 6)
				throw new Exception();
			
			
			System.out.println(args.length +" "+ (args.length - 6));
			
			ks = new int[args.length - 6];
			for(int i=6; i<args.length; i++){
				ks[i-6] = Integer.parseInt(args[i]);
			}
		}catch(Exception e){
			e.printStackTrace();
			usage();
			System.exit(1);
		}
		
		try{
			
			System.out.println("Loading sequence...");
			B3LLSequence b3seq = B3LLSequence.load(a_iseq);
			System.out.println("done");
			System.out.println("Loading NELSA...");
			NELSA nelsa = new NELSA();
			nelsa.load(a_inelsa);
			System.out.println("done");
			nelsa.setSequence(b3seq);
			
			
			StaticB2Tree words = new StaticB2Tree();
			//igtools.analyses.words.Carpena2009 c2009 = new igtools.analyses.words.Carpena2009(nelsa);
			igtools.analyses.recurrences.carpena.Carpena2009 c2009 = Carpena2009Factory.factory(nelsa, !original, !sigma_nor);
			
			for(int i=0; i<ks.length; i++){
				int k = ks[i];
				int count = 0;
				
				System.out.println("k "+k);
				
				double C0 = 0.0;
				//if(original){
					System.out.println("getting C0 with "+a_percentile+" percentile...");
					C0 = c2009.get_C0_percetile(a_percentile, k);
					System.out.println("C0 = "+C0);
					System.out.println("C0 - delta = "+ (C0 - (C0 * c_delta)));
				//}
				
				System.out.println("running...");
				count = c2009.select(k, C0, true, c_delta, words);
				System.out.println("count "+count);
			}
			System.out.println("printing...");
			words.printWords("S: ");
			
			/*System.out.println("comults...");
			TreeMap<Integer,Integer> comults = new TreeMap<Integer,Integer>();
			words.stats(comults);
			for(Map.Entry<Integer, Integer> entry : comults.entrySet()){
				System.out.println(entry.getKey()+" "+entry.getValue());
			}*/
			
			
		}catch(Exception ex){
			ex.printStackTrace();
		}
	}
}
