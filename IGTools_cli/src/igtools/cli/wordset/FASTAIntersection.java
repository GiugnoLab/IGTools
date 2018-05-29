package igtools.cli.wordset;

import igtools.common.ds.StaticB2Tree;
import igtools.common.sequence.B3LLSequence;

import java.io.File;
import java.util.Scanner;

/**
 * Get intersection of exact words.
 * 
 * @author vbonnici
 *
 */
public class FASTAIntersection {
	public static class WCounter implements StaticB2Tree.WordListener{
		public int wordCount = 0;
		public int totalLength = 0;
		public WCounter(){}
		@Override
		public void word(String s) {
			wordCount++;
			totalLength += s.length() + 1;
		}
	}
	
	
	public static void usage(){
		System.out.println("Usage: cmd {file list}");
	}
	
	public static void main(String[] args){
		
		
		/*
		 * remove duplicate and sub-included fasta polymers
		 */
		
		
		try{
			
			StaticB2Tree b2tree = new StaticB2Tree();
			int nofWords = 0;
			
			System.out.println("concatenating...");
			for(int i=0; i<args.length; i++){
				System.out.println(args[i]);
				Scanner scanner = new Scanner(new File(args[i]));
				String line;
				while(scanner.hasNextLine()){
					line = scanner.nextLine().trim();
					if(line.length() > 0){
						nofWords++;
						b2tree.add(new B3LLSequence(line), 0, line.length());
					}
				}
				scanner.close();
			}
			
			System.out.println("nof input words: "+ nofWords);
			
			
			int isize = 0;
			int nof = args.length;
			Scanner scanner = new Scanner(new File(args[0]));
			String line;
			while(scanner.hasNextLine()){
				line = scanner.nextLine().trim();
				if(line.length() > 0){
					if(b2tree.nof(new B3LLSequence(line)) == nof){
						isize++;
						System.out.println("W: "+line);
					}
				}
			}
			scanner.close();
			
			System.out.println("isize "+isize);
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
}
