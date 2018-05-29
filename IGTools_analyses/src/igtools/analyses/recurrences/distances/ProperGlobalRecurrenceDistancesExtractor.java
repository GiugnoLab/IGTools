package igtools.analyses.recurrences.distances;

import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.INELSA;

public abstract class ProperGlobalRecurrenceDistancesExtractor extends ProperRecurrenceDistancesExtractor{
	
	
	/**
	 * 
	 * @param overlapping
	 * @param skipNgaps
	 * @return
	 */
	public static ProperGlobalRecurrenceDistancesExtractor factory(boolean overlapping, boolean skipNgaps){
		if(overlapping){
			if(skipNgaps){
				//overlapping skip
				return new ProperGlobalRecurrenceDistancesExtractor() {
					
					@Override
					public int[] recurrencesDistances(IELSAIterator it) {
						int[] sa = it.elsa().sa();
						int[] ns = ((INELSA)it.elsa()).ns();
						
						int[] poss = new int[it.multiplicity()];
						int[] nss = new int[it.multiplicity()];
						for(int i = it.istart(), j=0; i<it.iend(); i++, j++){
							poss[j] = sa[i];
							ns[j] = ns[i];
						}
						
//						System.out.println("nof pos "+poss.length);
						
						
						int count = 0;
						
						for(int i=0; i<poss.length-1;i++){
							for(int j=i+1;j<poss.length; j++){
								if(nss[i] > poss[j]-poss[i]){
									count++;
								}
								else
									break;
							}
						}
						
						
						
//						System.out.println("count "+count);
						int[] recs = new int[count];

						int m = 0;
						for(int i=0; i<poss.length-1;i++){
							for(int j=i+1;j<poss.length; j++){
								if(nss[i] > poss[j]-poss[i]){
									recs[m] = poss[j] - poss[i];
									m++;
								}
								else
									break;
							}
						}
						return recs;
						
						
					}
				};
			}
			else{
				//overlapping noskip
				return new ProperGlobalRecurrenceDistancesExtractor() {
					
					@Override
					public int[] recurrencesDistances(IELSAIterator it) {
						int[] poss = it.sortedPositions();
						
//						System.out.println("nof pos "+poss.length);
						
						int count = 0;
						for(int i=0; i<poss.length-1;i++){
							count += poss.length-i;
						}
//						System.out.println("count "+count);
						int[] recs = new int[count];
						int p = 0;
						for(int i=0; i<poss.length-1;i++){
							for(int j=i+1; j<poss.length; j++){
								recs[p] = poss[j] - poss[i];
								p++;
							}
						}
						return recs;
						
						
					}
				};
			}
		}
		else{
			if(skipNgaps){
				//nooverlapping skip
				return new ProperGlobalRecurrenceDistancesExtractor() {
					
					@Override
					public int[] recurrencesDistances(IELSAIterator it) {
//						System.out.println("recurrences,global,overlap,skip");
						int[] sa = it.elsa().sa();
						int[] ns = ((INELSA)it.elsa()).ns();
						
						int k = it.k()-1;
						
						int[] poss = new int[it.multiplicity()];
						int[] nss = new int[it.multiplicity()];
						for(int i = it.istart(), j=0; i<it.iend(); i++, j++){
							poss[j] = sa[i];
							nss[j] = ns[i];
						}
						igtools.common.util.Arrays.keyedQsort(poss, nss);
						
//						System.out.println("nof pos "+poss.length);
						
						
						int count = 0;
						
						for(int i=0; i<poss.length-1;i++){
							int j=i+1;
							for(;j<poss.length; j++)
								if((poss[j] - poss[i] >k))
									break;
							
							for(;j<poss.length; j++){
							//for(int j=i+1;j<poss.length; j++){
								if(nss[i] > poss[j]-poss[i]){
									count++;
								}
							}
						}
						
						
						
//						System.out.println("count "+count);
						int[] recs = new int[count];

						int m = 0;
						for(int i=0; i<poss.length-1;i++){
							int j=i+1;
							for(;j<poss.length; j++)
								if((poss[j] - poss[i] >k))
									break;
							
							for(;j<poss.length; j++){
								if(nss[i] > poss[j]-poss[i]){
									recs[m] = poss[j] - poss[i];
									m++;
								}
							}
						}
						return recs;
						
						
					}
				};
			}
			else{
				//nooverlapping noskip
				return new ProperGlobalRecurrenceDistancesExtractor() {
					
					@Override
					public int[] recurrencesDistances(IELSAIterator it) {
						int k = it.k() -1 ;
						int[] poss = it.sortedPositions();
						
//						System.out.println("nof pos "+poss.length);
						
						int count = 0;
						
//						for(int i=0; i<poss.length-1;i++){
//							count += poss.length-i;
//						}
						
						for(int i=0; i<poss.length-1;i++){
							int j=i+1;
							for(;j<poss.length; j++)
								if((poss[j] - poss[i] >k))
									break;
							
							for(; j<poss.length; j++){
								count++;
							}
						}
						
						
//						System.out.println("count "+count);
						int[] recs = new int[count];
						int p = 0;
						for(int i=0; i<poss.length-1;i++){
							int j=i+1;
							for(;j<poss.length; j++)
								if((poss[j] - poss[i] >k))
									break;
							
							for(; j<poss.length; j++){
								recs[p] = poss[j] - poss[i];
								p++;
							}
						}
						return recs;
						
						
					}
				};
			}
		}
	}
}
