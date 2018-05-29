package igtools.cli.util;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;

public class GenomeStats {

	
	private static void usage(){
		System.out.println("Usage: i_seq.3bit");
	}
	
	public static void main(String[] args){
		
		String a_i_seq = null;
				
		try{
			a_i_seq = args[0];
		}catch(Exception e){
			usage();
			System.exit(1);
		}
		
		try{
			
			B3LLSequence b3seq = B3LLSequence.load(a_i_seq);
			
			int nsInslandCount = 0;
			int max = 0;
			int l = 0;
			for(int i=0; i<b3seq.length(); i++){
				if(b3seq.getB3(i) == B3Nucleotide.N_CODE){
					l++;
				}
				else{
					if(l!=0)
						nsInslandCount++;
					if(l>max)
						max = l;
					l = 0;
				}
			}
			
			System.out.println(a_i_seq +"\t"+ b3seq.length() + "\t"+  b3seq.countBads() +"\t"+nsInslandCount+"\t"+ max);
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
	
}
