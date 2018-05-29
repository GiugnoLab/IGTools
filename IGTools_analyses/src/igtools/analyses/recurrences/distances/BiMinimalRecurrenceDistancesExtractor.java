package igtools.analyses.recurrences.distances;


import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.INELSA;


/**
 * TODO i think it's a dynamic programming problem, not just linear
 * TODO it does not work when it_a and it_b are the same word or it_a is a factor of it_b (and vice versa) 
 * 
 * @author vbonnici
 *
 */
public abstract class BiMinimalRecurrenceDistancesExtractor extends BiRecurrenceDistancesExtractor {

	
	public static BiMinimalRecurrenceDistancesExtractor factory(boolean overlapping, boolean skipNgaps){
		//TODO
//		System.out.println("BIMinimalFactory "+overlapping+" "+skipNgaps);
		if(overlapping){
			if(skipNgaps){
				//overlap skip
				
				return new BiMinimalRecurrenceDistancesExtractor() {
					@Override
					public int[] recurrencesDistances(IELSAIterator it_a, IELSAIterator it_b) {
						
						
						int[][] poss = new int[2][];
						int[][] nss = new int[2][];
						
						{
							int[] sa_a = it_a.elsa().sa();
							int[] ns_a = ((INELSA)it_a.elsa()).ns();
							int[] poss_a = new int[it_a.multiplicity()];
							int[] nss_a = new int[it_a.multiplicity()];
							for(int i = it_a.istart(), j=0; i<it_a.iend(); i++, j++){
								poss_a[j] = sa_a[i];
								nss_a[j] = ns_a[i];
							}
							igtools.common.util.Arrays.keyedQsort(poss_a, nss_a);
							poss[0] = poss_a;
							nss[0] = nss_a;
						}
						{
							int[] sa_b = it_b.elsa().sa();
							int[] ns_b = ((INELSA)it_b.elsa()).ns();
							int[] poss_b = new int[it_b.multiplicity()];
							int[] nss_b = new int[it_b.multiplicity()];
							for(int i = it_b.istart(), j=0; i<it_b.iend(); i++, j++){
								poss_b[j] = sa_b[i];
								nss_b[j] = ns_b[i];
							}
							igtools.common.util.Arrays.keyedQsort(poss_b, nss_b);
							poss[1] = poss_b;
							nss[1] = nss_b;
						}
						
						
						if(poss[0].length>0 && poss[1].length>0){
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
								
								if((poss[a_j][poss_i[a_j]]) - (poss[a_i][poss_i[a_i] - 1]) < nss[a_i][poss_i[a_i] - 1])
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
								if((poss[a_j][poss_i[a_j]]) - (poss[a_i][poss_i[a_i] - 1]) < nss[a_i][poss_i[a_i] - 1]){
									dists[i] = (poss[a_j][poss_i[a_j]]) - (poss[a_i][poss_i[a_i] - 1]);
									//System.out.println(poss[a_j][poss_i[a_j]] +"\t"+ poss[a_i][poss_i[a_i] - 1] +"\t"+dists[i]);
									i++;
								}
								a_i = a_i == 0 ? 1 : 0;
					            a_j = a_j == 0 ? 1 : 0;
							}
						
							
							return dists;
							
							
						}
						else{
							int[] dists = {};
							return dists;
						}
					}
				};
			}
			else{
				//overlap noskip
				return new BiMinimalRecurrenceDistancesExtractor() {
					@Override
					public int[] recurrencesDistances(IELSAIterator it_a, IELSAIterator it_b) {
						int[][] poss = new int[2][];
						poss[0] = it_a.sortedPositions();
						poss[1] = it_b.sortedPositions();
						
						
						if(poss[0].length>0 && poss[1].length>0){
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
								//System.out.println(poss[a_j][poss_i[a_j]] +"\t"+ poss[a_i][poss_i[a_i] - 1] +"\t"+dists[i]);
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
					}
				};
			}
		}
		else{
			if(skipNgaps){
				//nooverlap skip
				return new BiMinimalRecurrenceDistancesExtractor() {
					@Override
					public int[] recurrencesDistances(IELSAIterator it_a, IELSAIterator it_b) {
//						System.out.println("biminimal,nooverlap, skip");
						
						int[][] poss = new int[2][];
						int[][] nss = new int[2][];
						
						{
							int[] sa_a = it_a.elsa().sa();
							int[] ns_a = ((INELSA)it_a.elsa()).ns();
							int[] poss_a = new int[it_a.multiplicity()];
							int[] nss_a = new int[it_a.multiplicity()];
							for(int i = it_a.istart(), j=0; i<it_a.iend(); i++, j++){
								poss_a[j] = sa_a[i];
								nss_a[j] = ns_a[i];
							}
							igtools.common.util.Arrays.keyedQsort(poss_a, nss_a);
							poss[0] = poss_a;
							nss[0] = nss_a;
						}
						{
							int[] sa_b = it_b.elsa().sa();
							int[] ns_b = ((INELSA)it_b.elsa()).ns();
							int[] poss_b = new int[it_b.multiplicity()];
							int[] nss_b = new int[it_b.multiplicity()];
							for(int i = it_b.istart(), j=0; i<it_b.iend(); i++, j++){
								poss_b[j] = sa_b[i];
								nss_b[j] = ns_b[i];
							}
							igtools.common.util.Arrays.keyedQsort(poss_b, nss_b);
							poss[1] = poss_b;
							nss[1] = nss_b;
						}
						
						
						if(poss[0].length>0 && poss[1].length>0){
							int[] ks = {it_a.k()-1, it_b.k()-1};
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
								
								//System.out.println(poss[a_j][poss_i[a_j]] +"\t"+ poss[a_i][poss_i[a_i] - 1] +"\t"+
								//		((poss[a_j][poss_i[a_j]]) - (poss[a_i][poss_i[a_i] - 1]))+"\t"+(nss[a_i][poss_i[a_i] - 1]));
								
								if(((poss[a_j][poss_i[a_j]]) - (poss[a_i][poss_i[a_i] - 1]) > ks[a_i]) &&
										((poss[a_j][poss_i[a_j]]) - (poss[a_i][poss_i[a_i] - 1]) < nss[a_i][poss_i[a_i] - 1]))
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
								
								if(((poss[a_j][poss_i[a_j]]) - (poss[a_i][poss_i[a_i] - 1]) > ks[a_i]) &&
										((poss[a_j][poss_i[a_j]]) - (poss[a_i][poss_i[a_i] - 1]) < nss[a_i][poss_i[a_i] - 1])){
									dists[i] = (poss[a_j][poss_i[a_j]]) - (poss[a_i][poss_i[a_i] - 1]);
									//System.out.println(poss[a_j][poss_i[a_j]] +"\t"+ poss[a_i][poss_i[a_i] - 1] +"\t"+dists[i]);
									i++;
								}
								
								a_i = a_i == 0 ? 1 : 0;
					            a_j = a_j == 0 ? 1 : 0;
							}
						
							
							return dists;
							
							
						}
						else{
							int[] dists = {};
							return dists;
						}
					}
				};
			}
			else{
				//nooverlap noskip
				return new BiMinimalRecurrenceDistancesExtractor() {
					@Override
					public int[] recurrencesDistances(IELSAIterator it_a, IELSAIterator it_b) {
//						System.out.println("biminimal,nooverlap, noskip,");
						/*
						int[] pos_a = it_a.sortedPositions();
						int[] pos_b = it_b.sortedPositions();
						
						
						if(pos_a.length>0 && pos_b.length>0){							
							int intersection =  igtools.common.util.Arrays.intersectionSize(pos_a, pos_b);
							
							int t_length = (pos_a.length - intersection) + (pos_b.length - intersection);
							
							int[] e_pos_a = new int[t_length];//e_pos_a[i] = pos_a[x] if pos_a[x] is in pos_b; -1 otherwise
							int[] e_pos_b = new int[t_length];
							
							int a_i = 0, b_i = 0;
							for(int i=0; i<t_length; i++){
								if(pos_a[a_i] == pos_b[b_i]){
									e_pos_a[i] = e_pos_b[i] = pos_a[a_i];
									a_i++;
									b_i++;
								}
								else if(pos_a[a_i] < pos_b[b_i]){
									e_pos_a[i] = pos_a[a_i];
									e_pos_b[i] = -1;
									a_i++;
								}
								else{//pos_a[a_i] > pos_b[b_i]
									e_pos_a[i] = -1;
									e_pos_b[i] = pos_b[b_i];
									b_i++;
								}
							}
							
							pos_a = null;
							pos_b = null;
							
							
							 //   A-AA,--A-
							 //   -B--,BBBB
							 //   
							 //   A-B B-A, A-B B-A A-B
							 //
							
							a_i = b_i = 0;
							if(a_i == -1){

							}
							for(;a_i < t_length; a_i++){
							}
							
						}
						else{
							int[] dists = {};
							return dists;
						}
						*/
						
						
						int[][] poss = new int[2][];
						poss[0] = it_a.sortedPositions();
						poss[1] = it_b.sortedPositions();
						
						
						if(poss[0].length>0 && poss[1].length>0){
							int[] ks = {it_a.k()-1, it_b.k()-1};
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
								//System.out.println(poss[a_j][poss_i[a_j]] +"\t"+ poss[a_i][poss_i[a_i] - 1] +"\t"+((poss[a_j][poss_i[a_j]]) - (poss[a_i][poss_i[a_i] - 1])));
								
								if((poss[a_j][poss_i[a_j]]) - (poss[a_i][poss_i[a_i] - 1]) > ks[a_i])
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
								
								//System.out.println(poss[a_j][poss_i[a_j]] +"\t"+ poss[a_i][poss_i[a_i] - 1] +"\t"+dists[i]);
								
								if((poss[a_j][poss_i[a_j]]) - (poss[a_i][poss_i[a_i] - 1]) > ks[a_i]){
									dists[i] = (poss[a_j][poss_i[a_j]]) - (poss[a_i][poss_i[a_i] - 1]);
									//System.out.println(poss[a_j][poss_i[a_j]] +"\t"+ poss[a_i][poss_i[a_i] - 1] +"\t"+dists[i]);
									i++;
								}
								
								a_i = a_i == 0 ? 1 : 0;
					            a_j = a_j == 0 ? 1 : 0;
							}
						
							
							return dists;
							
							
						}
						else{
							int[] dists = {};
							return dists;
						}
					}
				};
			}
		}
	}

}
