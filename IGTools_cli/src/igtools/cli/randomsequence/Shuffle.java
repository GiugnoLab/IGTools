package igtools.cli.randomsequence;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.common.util.Timer;

public class Shuffle {

	
	public static void usage(){
		System.out.println("Usage: cmd iseq.3bit oseq.3bit perc k");
	}
	
	
	public static void main(String[] args){
		String a_iseq = "";
		String a_oseq = "";
		double a_perc = 0;
		int a_k = 0;
		
		try{
			a_iseq = args[0];
			a_oseq = args[1];
			a_perc = Double.parseDouble(args[2]);
			a_k = Integer.parseInt(args[3]);
		}catch(Exception e){
			usage();
			System.exit(1);
		}
		
		Timer timer = new Timer();
		
		try {
			System.out.println("Loading sequence...");
			B3LLSequence b3seq = B3LLSequence.load(a_iseq);
			System.out.println("done "+timer.getElapsedSecs() +"sec.\n");
			
			double seq_length = b3seq.length();
			double to_shuffle = seq_length * a_perc;
			
			
			if(a_k == 1){
				for(int i=0; i<to_shuffle; i++){
					int p1 = (int) (Math.random() * seq_length);
					
					int c1 = b3seq.getB3(p1);
					while(c1 == B3Nucleotide.N_CODE){
						p1 = (int) (Math.random() * seq_length);
						c1 = b3seq.getB3(p1);
					}
					int p2 = (int) (Math.random() * seq_length);
					int c2 = b3seq.getB3(p2);
					while(c2 == B3Nucleotide.N_CODE){
						p2 = (int) (Math.random() * seq_length);
						c2 = b3seq.getB3(p2);
					}
					b3seq.setB3(p1, c2);
					b3seq.setB3(p2, c1);
				}
			}
			else{
				int[] c1 = new int[a_k];
				int[] c2 = new int[a_k];
				
				int p1, p2, j;
				
				for(int i=0; i<to_shuffle; i++){
					
					while(true){
						p1 = (int) (Math.random() * (seq_length-a_k));
						for(j=0; j<a_k; j++){
							c1[j] = b3seq.getB3(p1+j);
							if(c1[j] == B3Nucleotide.N_CODE){
								break;
							}
						}
						if(j == a_k)
							break;
					}
					
					
					while(true){
						p2 = (int) (Math.random() * (seq_length-a_k));
						for(j=0; j<a_k; j++){
							c2[j] = b3seq.getB3(p2+j);
							if(c2[j] == B3Nucleotide.N_CODE){
								break;
							}
						}
						if(j == a_k)
							break;
					}
					
					for(j=0; j<a_k; j++){
						b3seq.setB3(p1+j, c2[j]);
						b3seq.setB3(p2+j, c1[j]);
					}
				}
			}
			
			b3seq.save(a_oseq);
			
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println(e);
			System.exit(1);
		}
		
		
		
		
	}
	
	
}
