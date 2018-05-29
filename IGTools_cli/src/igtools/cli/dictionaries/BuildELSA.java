package igtools.cli.dictionaries;

import igtools.common.sequence.B3LLSequence;
import igtools.common.sequence.B3Sequence;
import igtools.common.util.Timer;
import igtools.dictionaries.elsa.ELSA;

public class BuildELSA {

	
	
	public static void usage(){
		System.out.println("Usage: cmd inputSequence.3bit outputFile.elsa");
	}
	
	
	public static void main(String[] args){
		
		String i_file = "";
		String o_file = "";
		
		try{
			i_file = args[0];
			o_file = args[1];
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
			System.out.println("Building ELSA...");
			ELSA elsa = new ELSA(b3seq);
			System.out.println("done "+timer.getElapsedSecs()+" sec.");
			
			
			timer.reset();
			System.out.println("Saving ELSA...");
			elsa.save(o_file);
			System.out.println("done "+timer.getElapsedSecs()+" sec.");
			
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
		
	}
}
