package igtools.cli.teaching;


import java.util.Arrays;
import java.util.Stack;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.ExtensionNELSAIterator;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

public class Es {

	
	public static void main(String[] args){
		System.out.println("hello");
		
		
		B3Nucleotide cc = B3Nucleotide.C;
		cc.code();
		
		byte cc_code = B3Nucleotide.C_CODE; 
		
		B3Nucleotide[] kmer = {B3Nucleotide.A, B3Nucleotide.C, B3Nucleotide.G}; 
		
		System.out.println(  B3Nucleotide.toString(kmer)  );
		
		B3Nucleotide[] mer = B3Nucleotide.toB3("ACGNTT");
		System.out.println(  B3Nucleotide.toString(mer)  );
		
		
		B3Nucleotide[] mer2 = new B3Nucleotide[2];
		B3Nucleotide.toB3("AC", mer2);
		
		//---------------------------------------------------------------------------------------
		System.out.println("================================================================================");
		
		String  g = "ACGNTTACGATTGGCCAGTTCAANCGTTCCAGTTTGGGAAAACACGTGCAAGTTCAG";
		System.out.println(g);
		
		 B3LLSequence b3seq = new B3LLSequence(g);
		 //b3seq = new B3LLSequence( kmer );
		 
		 NELSA nelsa = new NELSA(b3seq);
		 
		 int k = 1;
		 
		 B3Nucleotide[] kk = new B3Nucleotide[k];
		 
		 System.out.println("|D_"+k+"| = "+ nelsa.nof(k));
		 System.out.println("|T_"+k+"| = "+ nelsa.nof_mults(k));
		 
		 IELSAIterator it = nelsa.begin(k);
		 while(it.next()){
			//B3Nucleotide[] kk = it.kmer();
			 it.kmer(kk);
			
			System.out.print(B3Nucleotide.toString(kk)+"\t");
			System.out.println(it.multiplicity());
		 }
		
		 System.out.println("================================================================================");
			 
		 IELSAIterator itf = nelsa.find( kmer );
		 if(itf != null){
			 System.out.print(B3Nucleotide.toString( itf.kmer() )+"\t");
				System.out.println(itf.multiplicity());
				
			int[] pos = itf.sortedPositions();
			pos = itf.positions();
			
			for(int i=0; i<pos.length; i++){
				System.out.print(pos[i]+" ");
			}
			System.out.println();
			
			//int[] apos = Arrays.copyOf(pos, pos.length);
			
			System.out.println("["+itf.istart()+","+itf.iend()+"] = ["+nelsa.first(kmer)+"]");
			
			
			nelsa.find_rc(kmer);
		 }
		 
 System.out.println("================================================================================");

		k = 12;
		kk = new B3Nucleotide[k];
		int SA[] = nelsa.sa();
		for(int i=0; i<SA.length; i++){
			b3seq.getB3(SA[i], kk);
			System.out.println(i +"\t"+ SA[i] +"\t"+ B3Nucleotide.toString(kk) );
		}
		

 System.out.println("================================================================================");
		 
		 k = 1;
		 it = nelsa.begin(k);
		 elongate(nelsa, it.clone());

		 System.out.println("================================================================================");

		 
		 k = 1;
		 it = nelsa.begin(k);
		 while(it.next()){
			 elongate_tohapax(nelsa, it.clone());
		 }

		 System.out.println("================================================================================");
		 
		 k = 1;
		 it = nelsa.begin(k);
		 while(it.next()){
			 elongate_LR(nelsa, it);
		 }
		 
		 
		  System.out.println("================================================================================");


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
