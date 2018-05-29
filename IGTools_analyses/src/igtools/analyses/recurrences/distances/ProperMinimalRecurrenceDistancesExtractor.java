package igtools.analyses.recurrences.distances;


import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.INELSA;


/**
 * 
 * Retrieve proper-recurrences distances of words
 * 
 * @author vbonnici
 *
 */
public abstract class ProperMinimalRecurrenceDistancesExtractor extends ProperRecurrenceDistancesExtractor{
	
	
	/**
	 * 
	 * @param overlapping
	 * @param skipNgaps
	 * @return
	 */
	public static ProperMinimalRecurrenceDistancesExtractor factory(boolean overlapping, boolean skipNgaps){
//		System.out.println("overlapping "+overlapping+";skipNgaps "+skipNgaps);
		if(overlapping){
			if(skipNgaps){
				//overlapping skipping
				return new ProperMinimalRecurrenceDistancesExtractor(){
					@Override
					public int[] recurrencesDistances(IELSAIterator it){
//						System.out.println("recurrencesDistances,overlapping,skip");
						int[] sa = it.elsa().sa();
						int[] ns = ((INELSA)it.elsa()).ns();
						
						int[] poss = new int[it.multiplicity()];
						int[] nss = new int[it.multiplicity()];
						for(int i = it.istart(), j=0; i<it.iend(); i++, j++){
							poss[j] = sa[i];
							nss[j] = ns[i];
						}
						igtools.common.util.Arrays.keyedQsort(poss, nss);
						
//						System.out.println("nof raw positions  "+ poss.length);
						
						int count = 0;
						
//						for(int i=1; i<poss.length; i++){
//							if(nss[i-1] > poss[i]-poss[i-1]){
//								count++;
//							}
//						}
						
						//int p = 0; 
						for(int i=1; i<poss.length; i++){
							//System.out.println("("+i+","+(i-1)+")"+nss[i-1]+"\t"+poss[i]+"\t"+poss[i-1]+"\t"+(poss[i]-poss[i-1]));
							//if(nss[p] > poss[i]-poss[p]){
							if(nss[i-1] > poss[i]-poss[i-1]){
								count++;
								//p = i;
							}
//							else{
//								System.out.println(nss[p]+"\t"+poss[i]+"\t"+poss[p]+"\t"+(poss[i]-poss[p]));
//							}
						}
						
//						System.out.println("count "+count);
						
						
//						int[] recs = new int[count];
//						if(count > 0){
//							int j = 0;
//							
//							for(int i=1; i<poss.length; i++){
//								if(nss[i-1] > poss[i]-poss[i-1]){
//									recs[j] = poss[i]-poss[i-1];
//									j++;
//								}
//							}
//						}
//						
						int[] recs = new int[count];
						if(count > 0){
							int j = 0;
							//p = 0;
							
							for(int i=1; i<poss.length; i++){
								//
								if(nss[i-1] > poss[i]-poss[i-1]){
									recs[j] = poss[i] - poss[i-1];
									j++;
									
									//p = i;
								}
							}
						}
						
						return recs;
					};
				};
			}
			else{
				//overlapping no-skipping
				return new ProperMinimalRecurrenceDistancesExtractor(){
					@Override
					public int[] recurrencesDistances(IELSAIterator it){
//						System.out.println("recurrencesDistances,overlapping");
						int[] poss = it.sortedPositions();
						
						int count = poss.length - 1;
						
						int[] recs = new int[count];
						
						int j = 0;
						if(count > 1){
							for(int i=1; i<poss.length; i++){
								recs[j] = poss[i]-poss[i-1];
								j++;
							}
						}
												
						return recs;
					};
				};
			}
		}
		else{
			if(skipNgaps){
//				System.out.println("overlapping "+overlapping+";skipNgaps "+skipNgaps);
				//no-overlapping skipping
				return new ProperMinimalRecurrenceDistancesExtractor(){
					@Override
					public int[] recurrencesDistances(IELSAIterator it){
//						System.out.println("recurrencesDistances,skip");
						int k = it.k();
						
						int[] sa = it.elsa().sa();
						int[] ns = ((INELSA)it.elsa()).ns();
						
//						Map<Integer,Integer> mm = new HashMap<Integer,Integer>();
						
						int[] poss = new int[it.multiplicity()];
						int[] nss = new int[it.multiplicity()];
						for(int i = it.istart(), j=0; i<it.iend(); i++, j++){
							poss[j] = sa[i];
							nss[j] = ns[i];
							
//							mm.put(poss[j], nss[j]);
						}
						
						System.out.println("-poss "+poss.length);
						
						igtools.common.util.Arrays.keyedQsort(poss, nss);
						
//						for(int i=0; i<poss.length; i++){
//							if(nss[i] != mm.get(poss[i])){
//								System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
//							}
//						}
						
						
//						for(int i=0; i<poss.length; i++){
//							System.out.print("("+poss[i]+","+nss[i]+")");
//						}
//						System.out.println();
						
//						System.out.println("nof occurrences: "+poss.length);
						
						int count = 0;
						
						//int tcount = 0;
						
						int p = 0;
						
						for(int i=1; i<poss.length; i++){
							if(nss[p] < poss[i]-poss[p]){
								p = i; 
							}
							else if(poss[i] - poss[p] >=k){
								count++;
								//if(poss[i] - poss[p] == 154) tcount++;
								p = i;
							}
						}
						
						
						//System.out.println("tcount "+tcount);
//						System.out.println("nof recurreces: "+count);
						
						
						int[] recs = new int[count];
						if(count > 0){
							int j = 0;
							p = 0;
							
							for(int i=1; i<poss.length; i++){
								if(nss[p] < poss[i]-poss[p]){
									p = i; 
								}
								else if(poss[i] - poss[p] >=k){
									recs[j] = poss[i] - poss[p];
									j++;
									
									p = i;
								}
							}
						}
						
						return recs;
					};
				};
			}
			else{
				//no-overlapping no-skipping
				return new ProperMinimalRecurrenceDistancesExtractor(){
					@Override
					public int[] recurrencesDistances(IELSAIterator it){
//						System.out.println("recurrencesDistances");
						
						/*
						int[] poss = it.sortedPositions();
						int k = it.k();
						int count = 0;
						int p = 0;
						for(int i=1; i<poss.length; i++){
							if(poss[i] - poss[p] >=k){
								count++;
								p = i;
							}
						}
						int[] recs = new int[count];
						if(count > 0){
							int j = 0;
							for(int i=1; i<poss.length; i++){
								if(poss[i] - poss[p] >=k){
									recs[j] = poss[i] - poss[p];
									j++;
								
									p = i;
								}
							}
						}
						return recs;
						*/
						
						
						
						int k = it.k();
						
						int[] sa = it.elsa().sa();
						int[] ns = ((INELSA)it.elsa()).ns();
						
						int[] poss = new int[it.multiplicity()];
						int[] nss = new int[it.multiplicity()];
						for(int i = it.istart(), j=0; i<it.iend(); i++, j++){
							poss[j] = sa[i];
							nss[j] = ns[i];
						}
						
						igtools.common.util.Arrays.keyedQsort(poss, nss);

						
						int count = 0;
						
						int p = 0;
						
						for(int i=1; i<poss.length; i++){
							if(poss[i] - poss[p] >=k){
								count++;
								p = i;
							}
						}
						
						int[] recs = new int[count];
						if(count > 0){
							int j = 0;
							p = 0;
							
							for(int i=1; i<poss.length; i++){
								if(poss[i] - poss[p] >=k){
									recs[j] = poss[i] - poss[p];
									j++;
									
									p = i;
								}
							}
						}
						
						return recs;
					};
				};
			}
		}
	}
}
