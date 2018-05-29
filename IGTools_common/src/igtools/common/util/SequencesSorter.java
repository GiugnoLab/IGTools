package igtools.common.util;

import java.util.Arrays;
//import java.util.Set;
//import java.util.TreeSet;

public class SequencesSorter {
	
	
	public static void globalBtSort(String[] seqs, double[][] scores){
		//TODO
		int n = seqs.length;
		
		
		double[][] max_scores = new double[n][n];
		double max_score = Double.MIN_VALUE;
		for(int i=0; i<n; i++){
			for(int j=0; j<n; j++){
				if(scores[i][j] > max_score){
					max_score = scores[i][j];
				}
			}
		}
		for(int i=0; i<n; i++){
			for(int j=0; j<n; j++){
				max_scores[i][j] = max_score / scores[i][j];
			}
		}
		
		
		
		
		int[] candidates = new int[n];
		for(int i=0; i<n; i++){
			candidates[i] = -1;
		}
		boolean[] matched = new boolean[n];
		int[] best_solution = new int[n];
		double min_score = Double.MAX_VALUE;
		double partial_score = 0;
		
		
		int si = 0;
		int ci = n;
		while(si > 0){
			
			//next_cand;
			if(si == 0){
				ci = candidates[si] + 1;
			}
			else{
				ci = candidates[si];
				while(ci<n){
					ci++;
					if(	!matched[ci] && 
						(partial_score + max_scores[ci][si-1])<min_score ){
						break;
					}
				}
			}
			
			if(ci == n){
				//backtrack
				if(si != 0)
					partial_score -= max_scores[si][si-1];
				matched[candidates[si]] = false;
				candidates[si] = -1;
				si--;
			}
			else{
				if(si == n){
					//found new min global solution
					//no tracking operations
					min_score = partial_score;
					for(int i=0; i<n; i++){
						best_solution[i] = candidates[i];
					}
				}
				else{
					//set partial
					//forward
					partial_score += max_scores[ci][si-1];
					matched[ci] = true;
					candidates[si] = ci;
					si++;
				}
			}
		}
		
		
		String[] ordered = new String[seqs.length];
//		String temp;
		for(int i=0; i<n; i++){
//			if(best_solution[i] != i){
//				temp = seqs[i];
//				seqs[i] = seqs[best_solution[i]];
//				seqs[best_solution[i]] = temp;
//			}
			ordered[i] = seqs[best_solution[i]];
		}
		for(int i=0; i<seqs.length; i++){
			seqs[i] = ordered[i];
		}
	}
	
	public static void localSort(String[] seqs, double[][] scores){
		for(int i=0; i<seqs.length; i++){
			System.out.println("["+i+"]"+seqs[i]);
		}
		
		int n = seqs.length;
		
		boolean[] choosen = new boolean[n];
		
		int max_i=0, max_j=0; 
		double max_score = Double.MIN_VALUE;
		
		for(int i=0; i<n; i++){
			for(int j=0; j<n; j++){
				if(i != j){
					if(scores[i][j] > max_score){
						max_score = scores[i][j];
						max_i = i;
						max_j = j;
					}
				}
			}
		}
		
		System.out.println("First max score i["+max_i+"] j["+max_j+"] score["+scores[max_i][max_j]+"]");
		
		String[] ordered = new String[seqs.length];
		
//		String temp;
//		temp = seqs[0];
//		seqs[0] = seqs[max_i];
//		seqs[max_i] = temp;
		choosen[max_i] = true;
		ordered[0] = seqs[max_i];
		
		for(int i=1; i<n; i++){
			max_score = Double.MIN_VALUE;
			max_j = 0;
			for(int j=0; j<n; j++){
				if(!choosen[j] && i!=j){
					if(scores[max_i][j] > max_score){
						max_score = scores[max_i][j];
						max_j = j;
					}
				}
			}
			max_i = max_j;
			
//			temp = seqs[i];
//			seqs[i] = seqs[max_i];
//			seqs[max_i] = temp;
			choosen[max_i] = true;
			ordered[i] = seqs[max_i];
			
		}
		
		for(int i=0; i<seqs.length; i++){
			seqs[i] = ordered[i];
		}
		
		for(int i=0; i<seqs.length; i++){
			System.out.println("["+i+"]"+seqs[i]);
		}
	}
	
	
	public static void prefixSort(String[] seqs){
		//Set<String> prefix = new TreeSet<String>();
		Arrays.sort(seqs);
	}
}
