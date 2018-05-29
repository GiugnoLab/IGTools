package igtools.cli.recurrences.carpena;

import igtools.analyses.recurrences.carpena.Carpena2009Factory;
import igtools.common.ds.StaticB2Tree;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.NELSA;

import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;

public class CarpenaMeasures {
	private static void usage(){
		System.out.println("Usage: cmd sequence.3bit nelsa.nelsa [original|parent] k ");
	}
	
	
	public static void main(String[] args){
		String a_iseq = "";
		String a_inelsa = "";
		String a_type = "";
		int k = 0;
		
		boolean original = false;		
		
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
			
			k = Integer.parseInt(args[3]);
			if(args.length <= 3)
				throw new Exception();
			
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
			igtools.analyses.recurrences.carpena.Carpena2009 c2009 = Carpena2009Factory.factory(nelsa, !original, false);
			
			System.out.println("sigma sigma_nor C");
			
			
			double[] sigmas = c2009.get_sigmas(k);
			double[] sigma_nors = c2009.get_sigma_nors(k);
			double[] cs = c2009.get_Cs(k);
			
			for(int i=0; i<sigmas.length; i++){
				System.out.println(sigmas[i] +" "+ sigma_nors[i] +" "+ cs[i]);
			}
			
			System.out.println("\nSORTED\n");
			System.out.println("sigma sigma_nor C");
			Arrays.sort(sigmas);
			Arrays.sort(sigma_nors);
			Arrays.sort(cs);
			for(int i=0; i<sigmas.length; i++){
				System.out.println(sigmas[i] +" "+ sigma_nors[i] +" "+ cs[i]);
			}
			
		}catch(Exception ex){
			ex.printStackTrace();
		}
	}
}