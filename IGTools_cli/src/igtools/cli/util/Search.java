package igtools.cli.util;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.streams.BufferedFASTANucleotideStream;
import igtools.common.streams.INucleotideStream;

public class Search {

	
	private static void usage(){
		System.out.println("Usage: cmd fastaSequence  k-mer");
		System.out.println("example: chr1.fa  AATGGCCGTA");
	}
	
	private static void print(B3Nucleotide[] a){
		for(int i=0; i<a.length; i++)
			System.out.print(a[i]);
	}
	
	
	private static boolean equals(B3Nucleotide[] pattern, B3Nucleotide[] target, int s){
//		System.out.print("equals(");
//		print(pattern);
//		System.out.print(",");
//		print(target);
//		System.out.println(", "+s+")");
		
		int j=s;
		for(int i=0;  i<pattern.length; i++){
//			System.out.println("\tp["+i+"] t["+(j%pattern.length)+"]");
			if(pattern[i] != target[j%pattern.length])
				return false;
			j++;
		}
		
//		for(int i=0; i+s<pattern.length; i++){
//			System.out.println("\tp["+i+"] t["+(i+s)+"]");
//			if(pattern[i] != target[i+s])
//				return false;
//		}
//		
//		for(int i=pattern.length - s; i<pattern.length; i++){
//			System.out.println("\t.p["+i+"] t["+(i-s)+"]");
//			if(pattern[i] != target[i-s])
//				return false;
//		}
		
		
		return true;
	}
	
	
	public static void main(String[] args){
		
		String ifile = null;
		String iword = null;
		
		
		try{
			ifile = args[0];
			iword = args[1];
		}catch(Exception e){
			usage();
			System.exit(1);
		}
		
		System.out.println("Searching in "+ifile);
		System.out.println("searching for "+iword);
		
		
		try{
			
			B3Nucleotide[] pattern = new B3Nucleotide[iword.length()];
			for(int i=0; i<iword.length(); i++){
				pattern[i] = B3Nucleotide.by(iword.charAt(i));
				if(pattern[i] == null){
					throw new Exception("Unknown symbol "+iword.charAt(i));
				}
			}
			
			INucleotideStream<B3Nucleotide> fastaStream = new BufferedFASTANucleotideStream(ifile);
			B3Nucleotide[] window = new B3Nucleotide[iword.length()];
			B3Nucleotide r;
			int position = 0;
			
			for(int i=0; i<pattern.length; i++){
				r = fastaStream.next();
				if(r == null){
					System.out.println("Input sequence is smaller than searched word");
					System.exit(1);
				}
				window[i] = r;
			}
			if(equals(pattern, window, 0))
				System.out.println(position);
			
			int shift = 0;
			int s = 1;
			while((r = fastaStream.next()) != null){
				position++;
				window[shift] = r;
//				System.out.println("pos " + position);
				if(equals(pattern, window, s))
					System.out.println(position);
				
				shift++;
				if(shift == pattern.length)
					shift = 0;
				s++;
				if(s == pattern.length)
					s = 0;
			}
			
			
			System.out.println("final position "+position);
			
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
}
