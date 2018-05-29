package igtools.cli.coverage;

import java.util.BitSet;
import java.util.Set;
import java.util.TreeSet;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.common.util.Timer;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;



/**
 * 
 * Get the coverage of words having the same multiplicity
 * 
 * @author vbonnici
 *
 */
public class Coverage {
	
	public static void usage(){
		System.out.println("Usage:  cmd  isequence.3bit isequence.nelsa k...");
	}
	
	private static void print(int[] series){
		for(int i=0; i<series.length; i++)
			System.out.print(series[i]+"\t");
	}
	private static void println(int[] series){
		print(series);
		System.out.println();
	}
	
	private static void print(B3Nucleotide[] kmer){
		for(int i=0; i<kmer.length; i++)
			System.out.print(kmer[i]);
	}
	private static void println(B3Nucleotide[] kmer){
		print(kmer);
		System.out.println();
	}
	
	
	public static void main(String[] args){
		
//		DecimalFormat df = new DecimalFormat("#.##");
		
		String a_iseq = null;
		String a_inelsa = null;
		
		Set<Integer> kk = new TreeSet<Integer>();
		
		try{
			a_iseq = args[0];
			a_inelsa = args[1];
			
			if(args.length == 2){
				usage();
				System.exit(1);
			}
			
			for(int i=2; i<args.length; i++){
				kk.add(Integer.parseInt(args[i]));
			}
			
		}catch(Exception e){
			usage();
			System.exit(1);
		}
		
		
		try{
			
			Timer timer = new Timer();
			
			System.out.println(a_iseq);
			System.out.println("Loading sequence...");
			B3LLSequence b3seq = B3LLSequence.load(a_iseq);
			System.out.println("done "+timer.getElapsedSecs() +"sec.\n");
			
			int s_length = b3seq.length();
			System.out.println("Sequence length: " + s_length);
			
			timer.reset();
			System.out.println("Loading NELSA...");
			NELSA nelsa = new NELSA();
			nelsa.load(a_inelsa);
			System.out.println("done "+timer.getElapsedSecs()+" sec.");
			nelsa.setSequence(b3seq);
			
			
			BitSet coverage = new BitSet(b3seq.length());
			int cardinality;
			Set<Integer> multiplicities = new TreeSet<Integer>();
			int[] positions;
			
			for(int k : kk){
				
				//B3Nucleotide[] kmer = new B3Nucleotide[k];
				
				IELSAIterator it = nelsa.begin(k);
				while(it.next()){
					//positions = it.sortedPositions();
					//kmer = it.kmer();
					
					multiplicities.add(it.multiplicity());
					
					if(it.isMinimalHapax()){
						multiplicities.add(-1);
					}
					else if(it.isGlobalMaximalRepeat()){
						multiplicities.add(-2);
					}
					
//					if(it.multiplicity() == 1){
//						for(int i=0; i<positions.length; i++){
//							coverage.set(positions[i], positions[i]+k);
//						}
//					}
				}
				
				
				for(int m : multiplicities){
					coverage.clear();
					
					if(m==-1){
						it = nelsa.begin(k);
						while(it.next()){
							if(it.isMinimalHapax()){
								positions = it.sortedPositions();
								for(int i=0; i<positions.length; i++){
									coverage.set(positions[i], positions[i]+k);
								}
							}
						}
						cardinality = coverage.cardinality();
						System.out.println("chr\t"+k+"\t"+m+"\t"+cardinality+"\t"+ ( ((double)cardinality) / ((double)s_length)  ));
					}
					else if(m==-2){
						it = nelsa.begin(k);
						while(it.next()){
							if(it.isGlobalMaximalRepeat()){
								positions = it.sortedPositions();
								for(int i=0; i<positions.length; i++){
									coverage.set(positions[i], positions[i]+k);
								}
							}
						}
						cardinality = coverage.cardinality();
						System.out.println("chr\t"+k+"\t"+m+"\t"+cardinality+"\t"+ ( ((double)cardinality) / ((double)s_length)  ));
					}
					else{
						it = nelsa.begin(k);
						while(it.next()){
							if(it.multiplicity() == m){
								positions = it.sortedPositions();
								for(int i=0; i<positions.length; i++){
									coverage.set(positions[i], positions[i]+k);
								}
							}
						}
						cardinality = coverage.cardinality();
						System.out.println("chr\t"+k+"\t"+m+"\t"+cardinality+"\t"+ ( ((double)cardinality) / ((double)s_length)  ));
					}
				}
			}
			
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}

}
