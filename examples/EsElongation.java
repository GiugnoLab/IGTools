import java.util.Arrays;
import java.util.Stack;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.ExtensionNELSAIterator;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

public class EsElongation {

	
	public static void main(String[] args){
		//The following example show how to elongate root kmers by their elongations occuring in the input sequence
		//The goal is obtained by making use of the ExtensionNELSAIterator
		//The first 3 solution are recursive methods, thus their cannot elongate for long runs due to practical limits
		//the last example simulates recursion stack by means of a stack data structure, in this way it has no limit on extension length
		
		System.out.println("================================================================================");
		
		String  g = "ACGNTTACGATTGGCCAGTTCAANCGTTCCAGTTTGGGAAAACACGTGCAAGTTCAG";
		System.out.println(g);
		
		B3LLSequence b3seq = new B3LLSequence(g);
		NELSA nelsa = new NELSA(b3seq);
		
		int SA[] = nelsa.sa();
		int LCP[] = nelsa.lcp();
		int NS[] = nelsa.ns();
		for(int i=0; i<SA.length; i++){
			//print the suffixes as they are sorted in the suffix array
			System.out.println(i +"\t"+ SA[i]+"\t"+ LCP[i]+"\t"+ NS[i] +"\t"+ b3seq.subSequence(SA[i], b3seq.length) );
		}
		

		System.out.println("================================================================================");
		
		//The ExtensionNELSAIterator iterates over the elongations of a given kmer
		//If the input iterator is representing the word w,
		//the the ExtensionNELSAIterator lists of the extensions in w{A,C,G,T} occurring in the NELSA, thus in the sequence.
		
		int k = 5;
		IELSAIterator it = nelsa.begin(k);
		while(it.next()){
			System.out.println(B3Nucleotide.toString(it.kmer()));
			ExtensionNELSAIterator nit = new ExtensionNELSAIterator(nelsa, it);
			while(nit.next()){
				System.out.println("\t"+ B3Nucleotide.toString(nit.kmer()));
			}
		}
		
		
		//elongate from 1-mers to the leafs of the suffix tree
		k = 1;
		it = nelsa.begin(k);
		elongate(nelsa, it.clone());

		System.out.println("================================================================================");

		 
		 //elongate until the kmer is not an hapax
		 k = 1;
		 it = nelsa.begin(k);
		 while(it.next()){
			 elongate_tohapax(nelsa, it.clone());
		 }

		System.out.println("================================================================================");
		 
		 
		//elongat until kmers are repeat and their elongability is greater than 1
		 k = 1;
		 it = nelsa.begin(k);
		 while(it.next()){
			 elongate_LR(nelsa, it);
		 }
		 
		 
		System.out.println("================================================================================");

		//elongate according to multiplicity and right elongability by simulatin recursion on a stack
		 k = 1;
		 it = nelsa.begin(k);
		 while(it.next()){
			 iterative_elongate_seed(nelsa, it);
		 }
		 
		 
		System.out.println("================================================================================");
	}
	

	private static void elongate(NELSA nelsa,IELSAIterator eit){
		while(eit.next()){
			System.out.println(eit.istart() +"\t"+ eit.k() +"\t"+ B3Nucleotide.toString(eit.kmer()) );
			elongate(nelsa, new ExtensionNELSAIterator(nelsa, eit.clone()));
		}
	}


	private static void elongate_tohapax(NELSA nelsa, IELSAIterator it){
		if(it.multiplicity() == 1){
			System.out.println(B3Nucleotide.toString(it.kmer()));
		}
		else{
			ExtensionNELSAIterator eit = new ExtensionNELSAIterator(nelsa, it);
			while(eit.next()){
				elongate_tohapax(nelsa, eit);
			}
		}
	}	


	private static void elongate_LR(NELSA nelsa, IELSAIterator it){
		int count = 0;
		ExtensionNELSAIterator eit = new ExtensionNELSAIterator(nelsa, it);
		while(eit.next()){
			count++;
		}
		
		if(count < 2){
			if(it.multiplicity() > 1){
				System.out.println(B3Nucleotide.toString(it.kmer()));
			}
		}
		else{
		
			eit = new ExtensionNELSAIterator(nelsa, it);
			
			while(eit.next()){
				if(eit.multiplicity() > 1){
					elongate_LR(nelsa, eit.clone());
				}
//				else{
//					System.out.println(B3Nucleotide.toString(eit.kmer()));
//				}
			}
		}
	}





	public static int elongability_R(NELSA nelsa, IELSAIterator it){
		int count = 0;
		ExtensionNELSAIterator eit = new ExtensionNELSAIterator(nelsa, it);
		while(eit.next()) count++;
		return count;
	}
	
	
	public static class Pair <T1, T2>{
		T1 first;
		T2 second;
		public Pair(T1 first, T2 second){
			this.first = first;
			this.second = second;
		}
	}



	public static void iterative_elongate_seed(NELSA nelsa, IELSAIterator _it){
		Stack< Pair<IELSAIterator,ExtensionNELSAIterator> > stack = new Stack< Pair<IELSAIterator,ExtensionNELSAIterator> >();
		stack.push(new Pair<IELSAIterator, ExtensionNELSAIterator>(_it, new ExtensionNELSAIterator(nelsa, _it)) );
		
		IELSAIterator it;
		ExtensionNELSAIterator eit;
		while(! stack.isEmpty()){
			it = stack.lastElement().first;
			if((it.multiplicity() > 1)){
				eit = stack.lastElement().second;
				if(eit.next()){
					if(elongability_R(nelsa, eit) <= 1){
						//stop recursion
						System.out.println("HEAD\t"+B3Nucleotide.toString(eit.kmer()));
					}
					else{
						//simulate recursive call
						stack.push( new Pair<IELSAIterator, ExtensionNELSAIterator>(eit.clone(), new ExtensionNELSAIterator(nelsa, eit))  );
					}
				}
				else{
					stack.pop();
				}
			}
			else{
				stack.pop();
			}
		}
	}

}
