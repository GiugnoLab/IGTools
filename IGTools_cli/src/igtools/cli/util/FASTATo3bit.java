package igtools.cli.util;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.common.streams.BufferedFASTANucleotideStream;
import igtools.common.streams.INucleotideStream;


public class FASTATo3bit {
	
	private static void usage(){
		System.out.println("Usage: cmd isequence.fa osequence.3bit");
	}
	
	public static void main(String[] args){
		
		String iseq=null, oseq=null;
		
		try{
			iseq = args[0];
			oseq = args[1];
		}catch(Exception e){
			usage();
			System.exit(1);
		}
		
		
		try{
			INucleotideStream<B3Nucleotide> fastaStream = new BufferedFASTANucleotideStream(iseq);
			//while(fastaStream.next()!=null);
			while(fastaStream.nextCode() != INucleotideStream.NULL_CODE);
			System.out.println("sequence length " + (fastaStream.position()));
			B3LLSequence bitgenome = new B3LLSequence(fastaStream.position() + 1);
			fastaStream.close();
			fastaStream = new BufferedFASTANucleotideStream(iseq);
			bitgenome.fillBy(fastaStream);
			fastaStream.close();
			bitgenome.save(oseq);
		}catch(Exception e){
			System.out.println(e);
			e.printStackTrace();
			System.exit(2);
		}
		System.exit(0);
	}
}
