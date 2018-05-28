import java.util.Arrays;
import java.util.Stack;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.ExtensionNELSAIterator;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

public class Es {

	
	public static void main(String[] args){
		
		//this is a single nucletodie in object representation
		//It can be A,C,G,T,N,NULL
		B3Nucleotide cc = B3Nucleotide.C;
		cc.code();
		
		byte cc_code = B3Nucleotide.C_CODE; 
		
		//this is the kmer ACG defined as an array of nucleotides
		B3Nucleotide[] kmer = {B3Nucleotide.A, B3Nucleotide.C, B3Nucleotide.G}; 
		
		//the static method B3Nucleotide.toString can be used to print B3Nucleotide arrays, namely kmers
		System.out.println(  B3Nucleotide.toString(kmer)  );
		
		//hot to create a kmer froma string
		B3Nucleotide[] mer = B3Nucleotide.toB3("ACGNTT");
		System.out.println(  B3Nucleotide.toString(mer)  );
		
		//another way to create kmer in object representation.
		B3Nucleotide[] mer2 = new B3Nucleotide[2];
		B3Nucleotide.toB3("AC", mer2);
		
		//---------------------------------------------------------------------------------------
		System.out.println("================================================================================");
		
		String  g = "ACGNTTACGATTGGCCAGTTCAANCGTTCCAGTTTGGGAAAACACGTGCAAGTTCAG";
		System.out.println(g);
		
		//This is a genomic sequence stored in a succint 3-bit representation.
		 B3LLSequence b3seq = new B3LLSequence(g);
		 
		 //a way to create a 3-bit sequence from an nucleotide array.
		 //b3seq = new B3LLSequence( kmer );
		 
		 //how to create a NELSA from a 3-bit nucleotide sequence
		 NELSA nelsa = new NELSA(b3seq);
		 
		 int k = 1;
		 
		 B3Nucleotide[] kk = new B3Nucleotide[k];
		 
		 //get the number of k-mers in our NELSA
		 System.out.println("|D_"+k+"| = "+ nelsa.nof(k));
		 //get the multiplicity table
		 System.out.println("|T_"+k+"| = "+ nelsa.nof_mults(k));
		 
		 //iterate over k-mers occurring in our sequence
		 IELSAIterator it = nelsa.begin(k);
		 while(it.next()){
			//B3Nucleotide[] kk = it.kmer();
			 it.kmer(kk);
			
			//print the kmer and its multiplicity
			System.out.print(B3Nucleotide.toString(kk)+"\t");
			System.out.println(it.multiplicity());
		 }
		
		 System.out.println("================================================================================");
			
			
		//search for a specific kmer by mean of the suffix-array
		 IELSAIterator itf = nelsa.find( kmer );
		 //if null, the kmer ib absent in our sequence
		 if(itf != null){
			System.out.print(B3Nucleotide.toString( itf.kmer() )+"\t");
			System.out.println(itf.multiplicity());
			
			//get the positions where the kmer occurrs, sorted
			int[] pos = itf.sortedPositions();
			//sorting the positions requires an extra time, but one can also get the unsorted positions
			pos = itf.positions();
			
			for(int i=0; i<pos.length; i++){
				System.out.print(pos[i]+" ");
			}
			System.out.println();
			
			
			//int[] apos = Arrays.copyOf(pos, pos.length);
			
			//both array of positions are a copy of the corresponding array slice made by Arrays,copyOf
			//however, one can retrieve the start and end (not inclusive) indexes on the suffix array
			System.out.println("["+itf.istart()+","+itf.iend()+"] = ["+nelsa.first(kmer)+"]");
			
			//a way to search for the reverce complement of the given kmer
			nelsa.find_rc(kmer);
		 }
		 
		System.out.println("================================================================================");

		//printing the suffix array by just showing onyl the first 12 nucleotides of every suffix
		k = 12;
		kk = new B3Nucleotide[k];
		int SA[] = nelsa.sa();
		int LCP[] = nelsa.lcp();
		int NS[] = nelsa.ns();
		for(int i=0; i<SA.length; i++){
			//B3LLSequence.get(position, array) retrieves the subsequece at a given position and as long as the given nucleotide array
			b3seq.getB3(SA[i], kk);
			System.out.println(i +"\t"+ SA[i]+"\t"+ LCP[i]+"\t"+ NS[i] +"\t"+ B3Nucleotide.toString(kk) );
		}
		

		System.out.println("================================================================================");
	
	}
	

}
