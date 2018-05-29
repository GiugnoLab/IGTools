package igtools.cli.wordset;

import igtools.common.ds.CountingB3Tree;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;

import java.io.File;
import java.util.Scanner;

public class DominantRC {

	public static void usage(){
		System.out.println("Usage: cmd seqs.fa");
		System.out.println("OUTPUT <W: word>");
	}
	
	public static void main(String[] args){
		String a_file = "";
		try{
			a_file = args[0];
		}catch(Exception e){
			usage();
			System.exit(0);
		}
		
		try{
			
			CountingB3Tree tree = new CountingB3Tree();
			
			B3Nucleotide[] ns;
			
			String line;
			
			Scanner scanner = new Scanner(new File(a_file));
			while(scanner.hasNextLine()){
				line = scanner.nextLine().trim();
				if(line.length() > 0){
					ns = B3Nucleotide.toB3(line);
					if(B3Nucleotide.compareToRC(ns) > 0){
						B3Nucleotide.reverseComplement(ns);
					}
					tree.add(ns);
				}
			}
			scanner.close();
			
			CountingB3Tree.WordListener list = new CountingB3Tree.WordListener() {
				@Override
				public void word(String w, int w_count, int c_count) {
					System.out.println("W: "+w);
				}
			};
			tree.listWords(list);
			
		}catch(Exception e){
			e.printStackTrace(System.err);
			System.err.println(e);
		}
	}
	
}
