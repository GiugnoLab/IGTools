package igtools.cli.distributions;

import java.util.Map;
import java.util.TreeMap;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

public class MultiplicityDistribution {

	public static void usage(){
		System.out.println("Usage: cmd k iseq.3bit iseq.nelsa");
		System.out.println("output the multiplicity distribution that assigns to each multiplicity value the number of words having such multiplicity");
		System.out.println("distribution etries are in numerical order in the space separated format");
		System.out.println("\t# multiplicity number_of_words");
	}
	
	public static void main(String[] args){
		
		int k = 0;
		String a_iseq = "";
		String a_inelsa = "";
		
		try{
			k = Integer.parseInt(args[0]);
			if(k < 0) throw new Exception();
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
		
			
			
			TreeMap<Integer,Integer> distr = new TreeMap<Integer,Integer>();
			IELSAIterator it = nelsa.begin(k);
			int m;
			while(it.next()){
				m = it.multiplicity();
				if(distr.containsKey(m)){
					distr.put(m, distr.get(m) + 1);
				}
				else{
					distr.put(m,1);
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
