import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.ExtensionNELSAIterator;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;



public class EsCreods {
	
	public static void main(String[] args){

		String iseq = "";
		String inelsa = "";
		int k_seed = 0;
	
		try{
			iseq = args[0];
			inelsa = args[1];
			k_seed = Integer.parseInt(args[2]);
		}catch(Exception e){
			System.out.println("Usage: cmd iseq.b3seq inelsa.nelsa k");
			System.exit(1);
		}

		try{
			
			//int max_length = Integer.parseInt(args[3]);
			int max_length = Integer.MAX_VALUE; //only for debug 
		
			
			
			System.out.println("################################################################################");
			System.out.println("# LOADING");
			System.out.println("################################################################################");
			
			B3LLSequence seq = B3LLSequence.load(iseq);
			NELSA nelsa = new NELSA();
			nelsa.load(inelsa);
			nelsa.setSequence(seq);
			
			
			
//			System.out.println("################################################################################");
//			System.out.println("# RECURSIVE SEARCH");
//			System.out.println("################################################################################");
//			
//			final Set<String> rec_found_heads = new TreeSet<String>(); //only for debug purpose
//			final Set<String> rec_found_creods = new TreeSet<String>(); //only for debug purpose
//			
//			ElongationListener listener = new ElongationListener(){
//				public void head(IELSAIterator it) {
//					System.out.println("REC HEAD "+B3Nucleotide.toString(it.kmer()));
//					
//					rec_found_heads.add(B3Nucleotide.toString(it.kmer()));
//				}
//				public void creod(IELSAIterator it, int head_length) {
//					System.out.println("REC CREOD "+head_length+" "+B3Nucleotide.toString(it.kmer()));
//					
//					rec_found_creods.add(B3Nucleotide.toString(it.kmer()));
//				}
//			};
//			
//			
//			IELSAIterator it = nelsa.begin(k_seed);
//			while(it.next()){
//				elongate_seed(nelsa, it, max_length, listener);
//			}
			
			
			
			System.out.println("################################################################################");
			System.out.println("# ITERATIVE SEARCH");
			System.out.println("################################################################################");
			
//			final Set<String> iter_found_heads = new TreeSet<String>(); //only for debug purpose
//			final Set<String> iter_found_creods = new TreeSet<String>(); //only for debug purpose
			
			ElongationListener iter_listener = new ElongationListener(){
				public void head(IELSAIterator it) {
					System.out.println("ITE HEAD "+B3Nucleotide.toString(it.kmer()));
					
//					String s = B3Nucleotide.toString(it.kmer());
//					if(s.length() < max_length){
//						iter_found_heads.add(s);
//					}
				}
				public void creod(IELSAIterator it, int head_length) {
					System.out.println("ITE CREOD "+head_length+" "+B3Nucleotide.toString(it.kmer()));
					
//					String s = B3Nucleotide.toString(it.kmer());
//					if(s.length() < max_length){
//						iter_found_creods.add(s);
//					}
				}
				
				public void n_head(IELSAIterator it) {
					System.out.println("ITE NHEAD "+B3Nucleotide.toString(it.kmer()));
				}
				public void n_creod(IELSAIterator it, int head_length) {
					System.out.println("ITE NCREOD "+head_length+" "+B3Nucleotide.toString(it.kmer()));
				}
			};
			
			IELSAIterator it;
			it = nelsa.begin(k_seed);
			while(it.next()){
				iterative_elongate_seed(nelsa, it, iter_listener);
			}
			
			
			
//			System.out.println("################################################################################");
//			System.out.println("# COMPARISON: RECURSIVE VS ITERATIVE VERSION");
//			System.out.println("################################################################################");
//			
//			
//			System.out.println("Nof recursive heads "+rec_found_heads.size());
//			System.out.println("Nof recursive creods "+rec_found_creods.size());
//			System.out.println("Nof iterative heads "+iter_found_heads.size());
//			System.out.println("Nof iterative creods "+iter_found_creods.size());
//			
//			
//			//not lets to compare the two resulting sets
//			int difference_count = 0;
//			int intersection_count = 0;
//			for(String s : rec_found_heads){
//				if(iter_found_heads.contains(s)){
//					intersection_count++;
//				}
//				else{
//					difference_count++;
//				}
//			}
//			System.out.println("Nof recursive heads not in iterative heads "+difference_count);
//			System.out.println("Nof recursive head in interative heads "+intersection_count);
//			
//			
//			difference_count = 0;
//			intersection_count = 0;
//			for(String s : rec_found_creods){
//				if(iter_found_creods.contains(s)){
//					intersection_count++;
//				}
//				else{
//					difference_count++;
//				}
//			}
//			System.out.println("Nof recursive creods not in iterative creods "+difference_count);
//			System.out.println("Nof recursive creods in interative creods "+intersection_count);
			
			
//			System.out.println("################################################################################");
//			System.out.println("# CREODIC LENGTHS");
//			System.out.println("################################################################################");
//			double mean = 0.0;
//			Map<Integer,Integer> colength = new TreeMap<Integer,Integer>();
//			for(String s : iter_found_creods){
//				mean += s.length();
//				colength.put(s.length(), colength.getOrDefault(s.length(), 0) + 1);
//			}
//			mean /= (double)iter_found_creods.size();
//			System.out.println("Average creodic length "+mean);
//			System.out.println("Creodic co-length: [length , nof creods]");
//			for(Map.Entry<Integer, Integer> entry : colength.entrySet()){
//				System.out.println(entry.getKey()+" "+entry.getValue());
//			}
//			
//			System.out.println("################################################################################");
//			System.out.println("# CREODIC HEAD LENGTHS");
//			System.out.println("################################################################################");
//			mean = 0.0;
//			//Map<Integer,Integer> colength = new TreeMap<Integer,Integer>();
//			colength.clear();
//			for(String s : iter_found_heads){
//				mean += s.length();
//				colength.put(s.length(), colength.getOrDefault(s.length(), 0) + 1);
//			}
//			mean /= (double)iter_found_heads.size();
//			System.out.println("Average creodic head length "+mean);
//			System.out.println("Creodic head co-length: [length , nof heads]");
//			for(Map.Entry<Integer, Integer> entry : colength.entrySet()){
//				System.out.println(entry.getKey()+" "+entry.getValue());
//			}
//			
//			System.out.println("################################################################################");
//			System.out.println("# CREODIC COVERAGE");
//			System.out.println("################################################################################");
//			
//			int[] coverage = new int[seq.length];
//			IELSAIterator fit;
//			int[] positions;
//			int creod_length;
//			for(String s : iter_found_creods){
//				creod_length = s.length();
//				fit = nelsa.find(new B3LLSequence(s));
//				positions = fit.positions();
//				for(int i=0; i<positions.length; i++){
//					for(int j=0; j<creod_length; j++){
//						coverage[ positions[i] + j  ]++;
//					}
//				}
//			}
//			
//			int nof_ns = seq.countBads();
//			int covered_positions = 0;
//			Map<Integer, Integer> cocoverage = new TreeMap<Integer, Integer>();
//			for(int i=0; i<coverage.length; i++){
//				if(coverage[i] > 0){
//					covered_positions++;
//					cocoverage.put(coverage[i], cocoverage.getOrDefault(coverage[i], 0) + 1 );
//				}
//			}
//			System.out.println("sequence length "+seq.length);
//			System.out.println("Nof Ns "+nof_ns);
//			System.out.println("Creodic coverage "+covered_positions);
//			System.out.println("Creodic coverage ratio "+(  ((double)covered_positions) / ((double)(seq.length - nof_ns))  ) );
//			System.out.println("Creodic head co-coverage: [positional coverage value , nof covered positions]");
//			for(Map.Entry<Integer, Integer> entry : cocoverage.entrySet()){
//				System.out.println(entry.getKey()+" "+entry.getValue());
//			}
//			
//			
//			System.out.println("################################################################################");
//			System.out.println("# CREODIC HEAD COVERAGE");
//			System.out.println("################################################################################");
//			
//			//coverage = new int[seq.length];
//			for(int i=0; i<coverage.length; i++){
//				coverage[i] = 0;
//			}
//			for(String s : iter_found_heads){
//				creod_length = s.length();
//				fit = nelsa.find(new B3LLSequence(s));
//				positions = fit.positions();
//				for(int i=0; i<positions.length; i++){
//					for(int j=0; j<creod_length; j++){
//						coverage[ positions[i] + j  ]++;
//					}
//				}
//			}
//			
//			covered_positions = 0;
//			//cocoverage = new TreeMap<Integer, Integer>();
//			cocoverage.clear();
//			for(int i=0; i<coverage.length; i++){
//				if(coverage[i] > 0){
//					covered_positions++;
//					cocoverage.put(coverage[i], cocoverage.getOrDefault(coverage[i], 0) + 1 );
//				}
//			}
//			System.out.println("sequence length "+seq.length);
//			System.out.println("Nof Ns "+nof_ns);
//			System.out.println("Creodic head coverage "+covered_positions);
//			System.out.println("Creodic head coverage ratio "+(  ((double)covered_positions) / ((double)(seq.length - nof_ns))  ) );
//			System.out.println("Creodic head co-coverage: [positional coverage value , nof covered positions]");
//			for(Map.Entry<Integer, Integer> entry : cocoverage.entrySet()){
//				System.out.println(entry.getKey()+" "+entry.getValue());
//			}
			
			
		}catch(Exception e){
			e.printStackTrace(System.err);
			System.err.println(e);
		}
	}
	
	
	
	
	public static interface ElongationListener{
		public void head(IELSAIterator it);
		public void creod(IELSAIterator it, int head_length);
		
		public void n_head(IELSAIterator it);
		public void n_creod(IELSAIterator it, int head_length);
	}
	
	
	
	/*
	 * We define the elongability of a k-mer it as the number of elongations of it, from it to it.x such that x in Gamma,
	 * that are  present within the genome G.
	 * 
	 * A k-mer is defined as creodic if and only if the elongability of it is equal to 1 and it is a repeat.
	 * 
	 * Let it to be a k-mer and eit to be an elongation of it,
	 * eit is the head of a creod if and only if:
	 *  - it is a repeat
	 *  - it is not creodic
	 *  - eit is creodic
	 *  
	 * A seed k-mer, it, is elongated, in order to discover one elongation of it that is a creodic head, until:
	 *  - it is no longer than a fixed max length
	 *  - it is not an hapax
	 *  - its elongability is greater than 1
	 *  
	 *  After a creodic head have been discovered, it is elongated, in order to identify the entire creod, until:
	 *   - it is no longer than a fixed max length
	 *   - it is not an hapax
	 *   - its elongability is equal to 1 
	 */
	
	
//	/*
//	 * This method elongates seed k-mers until they became a creodic head or until some constraints fail.
//	 */
//	public static void elongate_seed(NELSA nelsa, IELSAIterator it, int max_length, ElongationListener listener){
//		if((it.multiplicity() > 1) && (it.k() < max_length)){
//			ExtensionNELSAIterator eit = new ExtensionNELSAIterator(nelsa, it);
//			while(eit.next()){
//				if(elongability_R(nelsa, eit) <= 1){
//					//we found a head
//					//System.out.println("HEAD "+B3Nucleotide.toString(eit.kmer()));
//					listener.head(eit);
//					elongate_head(nelsa, eit, max_length, eit.k(), listener);
//				}
//				else{
//					elongate_seed(nelsa, eit.clone(), max_length, listener);
//				}
//			}
//		}
//	}
//	
//	/*
//	 * This method elongates creodic heads until their elongability is equal to 1.
//	 * It is supposed that at the first time it is called, elongability(it) == 1
//	 */
//	public static void elongate_head(NELSA nelsa, IELSAIterator it, int max_length, int head_length, ElongationListener listener){
//		ExtensionNELSAIterator eit = new ExtensionNELSAIterator(nelsa, it);
//		while(eit.next()){
//			if(	(eit.multiplicity() == 1) || (eit.k() == max_length) || (elongability_R(nelsa, eit) > 1)){
//				//stop elongation and output the creod
//				//System.out.println("CREOD "+head_length+" "+B3Nucleotide.toString(eit.kmer()));
//				listener.creod(eit, head_length);
//			}
//			else{
//				//there is still something to elongate
//				elongate_head(nelsa, eit.clone(), max_length, head_length, listener);
//			}
//		}
//	}
	
	
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
	
	/*
	 * This method elongates seed k-mers until they became a creodic head or until some constraints fail.
	 */
	public static void iterative_elongate_seed(NELSA nelsa, IELSAIterator _it, ElongationListener listener){
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
						//we found a head
						//System.out.println("HEAD "+B3Nucleotide.toString(eit.kmer()));
						listener.head(eit);
						iterative_elongate_head(nelsa, eit, eit.k(), listener);
					}
					else{
						//elongate_seed(nelsa, eit.clone(), max_length);
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
	
	public static void iterative_elongate_head(NELSA nelsa, IELSAIterator it, int head_length, ElongationListener listener){
		Stack<ExtensionNELSAIterator> stack = new Stack<ExtensionNELSAIterator>();
		stack.push(new ExtensionNELSAIterator(nelsa, it));
		
		ExtensionNELSAIterator eit;
		while(! stack.isEmpty()){
			eit = stack.lastElement();
			if(eit.next()){
				if(	(eit.multiplicity() == 1) || (elongability_R(nelsa, eit) > 1)){
					//stop elongation and output the creod
					//System.out.println("CREOD "+head_length+" "+B3Nucleotide.toString(eit.kmer()));
					listener.creod(eit, head_length);
				}
				else{
					//there is still something to elongate
					//elongate_head(nelsa, eit.clone(), max_length, head_length);
					stack.push(new ExtensionNELSAIterator(nelsa, eit));
				}
			}
			else{
				stack.pop();
			}
		}
	}
	
	
	
}
