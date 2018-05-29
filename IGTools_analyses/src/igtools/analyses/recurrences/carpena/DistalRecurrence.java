package igtools.analyses.recurrences.carpena;

import java.util.Arrays;
import java.util.TreeMap;

import igtools.dictionaries.elsa.IELSAIterator;

public class DistalRecurrence {

	public static int[] proper_codistances(IELSAIterator it){
		return proper_codistances(it.sortedPositions(), it.k(), true);
	}
	public static int[] proper_codistances(int[] pos, int k, boolean already_sorted){
		if(!already_sorted)
			Arrays.sort(pos);
		//pos must be sorted
		int count = 0;
		
		int last_i = 0;
		for(int i=1; i<pos.length; i++){
			if(pos[i] - pos[last_i] >= k){
				count++;
				last_i = i;
			}
		}
		
		int[] codists = new int[count];
		last_i = 0;
		int j = 0;
		for(int i=1; i<pos.length; i++){
			if(pos[i] - pos[last_i] >= k){
				codists[j] = pos[i] - pos[last_i];
				j++;
				last_i = i;
			}
			
		}
		return codists;
	};		
	
	
	public static
	double[] double_proper_codistances(IELSAIterator it)
	{
		return double_proper_codistances(it.sortedPositions(), it.k(), true);
	}
	
	public static 
	double[] double_proper_codistances(int[] pos, int k, boolean already_sorted)
	{
		if(!already_sorted)
			Arrays.sort(pos);
		//pos must be sorted
		int count = 0;
		
		int last_i = 0;
		for(int i=1; i<pos.length; i++){
			if(pos[i] - pos[last_i] >= k){
				count++;
				last_i = i;
			}
		}
		
		double[] codists = new double[count];
		last_i = 0;
		int j = 0;
		for(int i=1; i<pos.length; i++){
			if(pos[i] - pos[last_i] >= k){
				codists[j] = (double)(pos[i] - pos[last_i]);
				j++;
				last_i = i;
			}
			
		}
		return codists;
	};	
	
	
	
	
	
	
	
	
	
	
	
	
	
	public static int[] proper_overlapping_codistances(IELSAIterator it){
		return proper_overlapping_codistances(it.sortedPositions(),true);
	}
	public static int[] proper_overlapping_codistances(int[] pos, boolean already_sorted){
		//pos must be sorted
		if(!already_sorted)
			Arrays.sort(pos);
		int[] dists = new int[pos.length -1];
		for(int i=1; i<pos.length; i++){
			dists[i-1] = pos[i] - pos[i-1];
		}
		return dists;
	};
	
	
	
	public static TreeMap<Integer,Integer> proper_global_codistances(IELSAIterator it, int maxDist){
		TreeMap<Integer,Integer> codists = new TreeMap<Integer,Integer>();
		
		int[] pos = it.sortedPositions();
		
		Integer cc;
		int d;
		for(int i=0; i<pos.length-1; i++){
			for(int j=i+1; j<pos.length; j++){
				d = pos[j] - pos[i];
				if(d>maxDist)
					break;
				cc = codists.get( d );
				if(cc == null)
					codists.put(d,1);
				else
					codists.put(d, cc+1);
			}
		}
		
		return codists;
	}
	
	
	
	public static int[] codistances(IELSAIterator ita, IELSAIterator itb){
		return codistances(ita.sortedPositions(), itb.sortedPositions(), ita.k(), true);
	};
	
	public static int[] codistances(int[] pos_a, int[] pos_b, int k, boolean already_sorted){
		//pos_a and pos_b must be sorted
		int[][] poss = new int[2][];
				poss[0] = pos_a;
				poss[1] = pos_b;
		if(!already_sorted){
			Arrays.sort(poss[0]);
			Arrays.sort(poss[1]);
		}
		
		if(poss[0].length > 0 && poss[1].length > 0){
			
			//TODO
			
			return null;
		}
		else{
			int[] dists = {};
			return dists;
		}
	}
	
	
	public static int[] overlapping_codistances(IELSAIterator ita, IELSAIterator itb){
		return overlapping_codistances(ita.sortedPositions(), itb.sortedPositions(), true);
	}
	public static int[] overlapping_codistances(int[] pos_a, int[] pos_b, boolean already_sorted){
		//pos_a and pos_b must be sorted
		int[][] poss = new int[2][];
				poss[0] = pos_a;
				poss[1] = pos_b;
		if(!already_sorted){
			Arrays.sort(poss[0]);
			Arrays.sort(poss[1]);
		}
		
		if(poss[0].length > 0 && poss[1].length > 0){
			//pointers in poss[i]
			int[] poss_i = {0,0};
			
			int a_i = poss[0][0] <= poss[1][0] ? 0 : 1;
	        int a_j = poss[1][0] < poss[0][0] ? 0 : 1; 
			
	        int count = 0;
	        
			while(poss_i[0] < poss[0].length  &&  poss_i[1] < poss[1].length ){
				while(poss_i[a_i] < poss[a_i].length &&  
	                    poss[a_i][poss_i[a_i]] <= poss[a_j][poss_i[a_j]]
	                    ){
	                poss_i[a_i]++;
	            }
				
				count++;
				//dist = (poss[a_j][poss_i[a_j]]) - (poss[a_i][poss_i[a_i] - 1]);
				a_i = a_i == 0 ? 1 : 0;
	            a_j = a_j == 0 ? 1 : 0;
			}
			
			
			int[] dists = new int[count];
	        
	        poss_i[0] = 0;
			poss_i[1] = 0;
			a_i = poss[0][0] <= poss[1][0] ? 0 : 1;
	        a_j = poss[1][0] < poss[0][0] ? 0 : 1;
	        int i = 0;
			while(poss_i[0] < poss[0].length  &&  poss_i[1] < poss[1].length ){
				while(poss_i[a_i] < poss[a_i].length &&  
	                    poss[a_i][poss_i[a_i]] <= poss[a_j][poss_i[a_j]]
	                    ){
	                poss_i[a_i]++;
	            }
				dists[i] = (poss[a_j][poss_i[a_j]]) - (poss[a_i][poss_i[a_i] - 1]);
				i++;
				a_i = a_i == 0 ? 1 : 0;
	            a_j = a_j == 0 ? 1 : 0;
			}
		
			
			return dists;
		}
		else{
			int[] dists = {};
			return dists;
		}
		
	};
	
	
	
	
	
	/*
	public static TreeMap<Integer,Integer> co_recurrences_map(int[] recs){
		TreeMap<Integer,Integer> corecs = new TreeMap<Integer,Integer>();
		co_recurrences_map(recs,corecs);
		return corecs;
	}
	public static void co_recurrences_map(int[] recs, TreeMap<Integer,Integer> corecs){
		Integer corec;
		for(int i=0; i<recs.length; i++){
			corec = corecs.get(recs[i]);
			if(corec == null){
				corecs.put(recs[i], 1);
			}
			else{
				corecs.put(recs[i], corec+1);
			}
		}
	}
	
	public static int[][] co_recurrences_array(int[] recs){
		TreeMap<Integer,Integer> corecsMap = co_recurrences_map(recs);
		int[][] corecs = new int[2][corecsMap.size()];
		int i=0;
		for(Map.Entry<Integer,Integer> entry : corecsMap.entrySet()){
			corecs[0][i] = entry.getKey();
			corecs[1][i] = entry.getValue();
			i++;
		}
		return corecs;
	}
	*/
	
	
	
	public static int[][] co_recurrences_array(int[] recurrences, boolean already_sorted){
		if(!already_sorted)
			Arrays.sort(recurrences);
		
		int count = 0;
		int prev_value = recurrences[0];
		for(int i=0; i<recurrences.length; i++){
			if(recurrences[i] != prev_value){
				count++;
				prev_value = recurrences[i];
			}
		}
		count++;
		
		int[][] corecs = new int[2][count];
	
		
		count = 1;
		int p = 0;
		prev_value = recurrences[0];
		
		for(int i=0; i<recurrences.length; i++){
			if(recurrences[i] != prev_value){
				corecs[0][p] = prev_value;
				corecs[1][p] = count;
				p++;
				prev_value = recurrences[i];
				count = 1;
			}
			else{
				count++;
			}
		}
		//if(count>0){
		corecs[0][p] = recurrences[recurrences.length-1];
		corecs[1][p] = count;
		
		return corecs;
	}
	
	
}
