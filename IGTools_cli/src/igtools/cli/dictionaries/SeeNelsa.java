package igtools.cli.dictionaries;

import igtools.common.sequence.B3LLSequence;
import igtools.common.util.Timer;
import igtools.dictionaries.elsa.NELSA;

public class SeeNelsa {

	public static void usage(){
		System.out.println("Usage: cmd inputSequence.3bit k");
	}
	
	
	public static void main(String[] args){
		
		String i_file = "";
		int k  = 0;
		
		try{
			i_file = args[0];
			k = Integer.parseInt(args[1]);
		}catch(Exception e){
			usage();
			System.exit(1);
		}
		
		try{
			Timer timer = new Timer();
			
			timer.reset();
			System.out.println("Loading sequence...");
			B3LLSequence b3seq = B3LLSequence.load(i_file);
			System.out.println("done "+timer.getElapsedSecs()+" sec.");
			
			timer.reset();
			System.out.println("Building NELSA...");
			NELSA nelsa = new NELSA(b3seq);
			System.out.println("done "+timer.getElapsedSecs()+" sec.");
			
			
			
			int[] sa = nelsa.sa();
			int[] lcp = nelsa.lcp();
			int[] ns = nelsa.ns();
			
			for(int i=0; i<sa.length; i++){
				//System.out.println(i+"\t"+sa[i]+"\t"+lcp[i]+"\t"+ns[i]+"\t"+ b3seq.subSequence(sa[i], b3seq.length()) );
				System.out.println(i+"\t"+sa[i]+"\t"+lcp[i]+"\t"+ns[i]+"\t"+ b3seq.subSequence(sa[i], sa[i]+k) );
			}
			
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
		
	}
}
