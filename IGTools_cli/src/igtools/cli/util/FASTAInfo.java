package igtools.cli.util;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.streams.BufferedFASTANucleotideStream;
import igtools.common.streams.INucleotideStream;

public class FASTAInfo {

	
	private static void usage(){
		System.out.println("Usage: cmd isequence.fasta");
	}
	
	
	public static void main(String[] args){
		String a_ifile = "";
		
		try{
			a_ifile = args[0];
		}catch(Exception e){
			usage();
		}
		
		
		try{
			System.out.println(a_ifile);
			
			INucleotideStream<B3Nucleotide> fastaStream = new BufferedFASTANucleotideStream(a_ifile);
			while(fastaStream.nextCode() != INucleotideStream.NULL_CODE);
			System.out.println("sequence length " + (fastaStream.position()));
			
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
}
