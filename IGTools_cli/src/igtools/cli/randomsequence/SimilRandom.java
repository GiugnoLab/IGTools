package igtools.cli.randomsequence;

import igtools.analyses.random.Random;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.common.sequence.B3Sequence;

public class SimilRandom {
	
	
	
	public static void usage(){
		System.out.println("Usage: cmd iseq.3bit oseq.3bit");
	}
	
	
	public static void main(String[] args){
		String a_iseq = "";
		String a_oseq = "";
		
		try{
			
			a_iseq = args[0];
			a_oseq = args[1];
			
		}catch(Exception e){
			usage();
			System.exit(1);
		}
		
		
		try{
			
			
			System.out.println("Loading sequence...");
			B3LLSequence iseq = B3LLSequence.load(a_iseq);
			System.out.println("done");
			
			B3LLSequence oseq = (B3LLSequence)Random.generateUniform(iseq.length());
			
			for(int i=0; i<iseq.length(); i++){
				if(iseq.getB3(i) == B3Nucleotide.N_CODE)
					oseq.setB3(i, B3Nucleotide.N_CODE);
			}
			
			oseq.save(a_oseq);
			
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
}
