package igtools.cli.distributions;

import java.util.Map;
import java.util.TreeMap;

import igtools.analyses.Inform;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

public class RepeatLengthDistribution {

	public static void usage(){
		System.out.println("Usage: cmd m iseq.3bit iseq.nelsa");
		System.out.println("output the repeat length distribution that, given a multiplicity m,  assigns to each length the number of words having such length and multiplicity euqal to m");
		System.out.println("words are taken from D_1 to D_{MRL}");
		System.out.println("distribution etries are in numerical order in the space separated format");
		System.out.println("\t# length number_of_words");
	}
	
	public static void main(String[] args){
		
		int m = 0;
		String a_iseq = "";
		String a_inelsa = "";
		
		try{
			m = Integer.parseInt(args[0]);
			if(m < 0) throw new Exception();
			a_iseq = args[1];
			a_inelsa = args[2];
		}catch(Exception  e){
			usage();
			System.exit(0);
		}
		
		
		try{

			System.out.println("seq");
			B3LLSequence b3seq = B3LLSequence.load(a_iseq);
			System.out.println("nelsa");
			NELSA nelsa = new NELSA();
			nelsa.load(a_inelsa);
			nelsa.setSequence(b3seq);
			
			int mrl = Inform.Mrl(nelsa);
			
			TreeMap<Integer,Integer> distr = new TreeMap<Integer,Integer>();
			for(int k=1; k<=mrl; k++){
				IELSAIterator it = nelsa.begin(k);
				while(it.next()){
					if(it.multiplicity() == m){
						if(distr.containsKey(k)){
							distr.put(k, distr.get(k) + 1);
						}
						else{
							distr.put(k, 1);
						}
					}
				}
			}
			
			for(Map.Entry<Integer,Integer> entry : distr.entrySet()){
				System.out.println("# "+entry.getKey()+" "+entry.getValue());
			}
			
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
		
	}
	
	
}
