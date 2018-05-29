package igtools.cli.util;

import java.io.File;
import java.io.PrintWriter;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;


/**
 * Convert a sequence in the 3bit format to the FASTA format.
 * 
 * @author vbonnici
 *
 */
public class B3ToFASTA {

	
	
	public static void usage(){
		System.out.println("Usage: cmd  iseq.3bit oseq.fasta");
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
			B3LLSequence b3seq = B3LLSequence.load(a_iseq);
			System.out.println("done");
			
			PrintWriter writer = new PrintWriter(new File(a_oseq));
			
			int col = 0;
			for(int i=0; i<b3seq.length(); i++){
				if(col!=0 && col%80==0)
					writer.println("");
				writer.print(B3Nucleotide.charFor(b3seq.getB3(i)));
				col++;
			}
			
			writer.flush();
			writer.close();
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
}
