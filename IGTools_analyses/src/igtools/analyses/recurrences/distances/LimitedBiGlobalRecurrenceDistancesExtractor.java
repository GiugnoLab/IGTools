package igtools.analyses.recurrences.distances;

import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.INELSA;

public abstract class LimitedBiGlobalRecurrenceDistancesExtractor extends BiRecurrenceDistancesExtractor{

	protected Double maxDist = Double.POSITIVE_INFINITY;
	public void setMaxDistance(Double maxDist){
		this.maxDist = maxDist;
	}
	
	
	public static LimitedBiGlobalRecurrenceDistancesExtractor factory(boolean overlapping, boolean skipNgaps){
//		System.out.println("LimitedBiGlobalRecurrenceExtractor "+overlapping+" "+skipNgaps);
		//TODO
		if(overlapping){
			if(skipNgaps){
				//overlap, skip
				return new LimitedBiGlobalRecurrenceDistancesExtractor() {
					
					@Override
					public int[] recurrencesDistances(IELSAIterator it_a, IELSAIterator it_b) {
//						System.out.println("recurrences,bi,global,overlap,skip");


						int[] sa_a = it_a.elsa().sa();
						int[] ns_a = ((INELSA)it_a.elsa()).ns();
						int[] poss_a = new int[it_a.multiplicity()];
						int[] nss_a = new int[it_a.multiplicity()];
						for(int i = it_a.istart(), j=0; i<it_a.iend(); i++, j++){
							poss_a[j] = sa_a[i];
							nss_a[j] = ns_a[i];
						}
						igtools.common.util.Arrays.keyedQsort(poss_a, nss_a);
					
						int[] sa_b = it_b.elsa().sa();
						int[] ns_b = ((INELSA)it_b.elsa()).ns();
						int[] poss_b = new int[it_b.multiplicity()];
						int[] nss_b = new int[it_b.multiplicity()];
						for(int i = it_b.istart(), j=0; i<it_b.iend(); i++, j++){
							poss_b[j] = sa_b[i];
							nss_b[j] = ns_b[i];
						}
						igtools.common.util.Arrays.keyedQsort(poss_b, nss_b);
						
						int count = 0;
						for(int i=0; i<poss_a.length; i++){
							for(int j=0; j<poss_b.length; j++){
								if(poss_a[i] >= poss_b[j]){
									if(poss_a[i] - poss_b[j] < nss_b[j])
									if(poss_a[i] - poss_b[j] <= maxDist)
										count++;
								}
								else{
									if(poss_b[j] - poss_a[i] < nss_a[i]){
										if(poss_b[j] - poss_a[i] <= maxDist){
											count++;
										}
										else
											break;
									}
								}
							}
						}
						
						int[] dists = new int[count];
						
						
						int d = 0;
						for(int i=0; i<poss_a.length; i++){
							for(int j=0; j<poss_b.length; j++){
								if(poss_a[i] >= poss_b[j]){
									if(poss_a[i] - poss_b[j] < nss_b[j])
									if(poss_a[i] - poss_b[j] <= maxDist){
										dists[d] = poss_a[i] - poss_b[j];
										d++;
									}
								}
								else{
									if(poss_b[j] - poss_a[i] < nss_a[i]){
										if(poss_b[j] - poss_a[i] <= maxDist){
											dists[d] = poss_b[j] - poss_a[i];
											d++;
										}
										else
											break;
									}
								}
							}
						}
						
						
						return dists;
					}
				};
			}
			else{
				//overlap, noskip
				return new LimitedBiGlobalRecurrenceDistancesExtractor() {
					
					@Override
					public int[] recurrencesDistances(IELSAIterator it_a, IELSAIterator it_b) {
//						System.out.println("recurrences,bi,global,overlap,noskip");
						
						int[] pos_a = it_a.sortedPositions();
						int[] pos_b = it_b.sortedPositions();
						
//						System.out.println(pos_a+" "+pos_b);
						
						int count = 0;
						for(int i=0; i<pos_a.length; i++){
							for(int j=0; j<pos_b.length; j++){
								if(pos_a[i] >= pos_b[j]){
									if(pos_a[i] - pos_b[j] <= maxDist)
										count++;
								}
								else{
									if(pos_b[j] - pos_a[i] <= maxDist)
										count++;
									else 
										break;
								}
							}
						}
						
						int[] dists = new int[count];
						
						
						int d = 0;
						for(int i=0; i<pos_a.length; i++){
							for(int j=0; j<pos_b.length; j++){
								if(pos_a[i] >= pos_b[j]){
									if(pos_a[i] - pos_b[j] <= maxDist){
										dists[d] = pos_a[i] - pos_b[j];
										d++;
									}
								}
								else{
									if(pos_b[j] - pos_a[i] <= maxDist){
										dists[d] = pos_b[j] - pos_a[i];
										d++;
									}
									else
										break;
								}
							}
						}
						
						
						return dists;
					}
				};
			}
		}
		else{
			if(skipNgaps){
				//nooverlap, skip
				return new LimitedBiGlobalRecurrenceDistancesExtractor() {
					
					@Override
					public int[] recurrencesDistances(IELSAIterator it_a, IELSAIterator it_b) {
//						System.out.println("recurrences,bi,global,overlap,skip");

						int k_a = it_a.k()-1;
						int k_b = it_b.k()-1;

						int[] sa_a = it_a.elsa().sa();
						int[] ns_a = ((INELSA)it_a.elsa()).ns();
						int[] poss_a = new int[it_a.multiplicity()];
						int[] nss_a = new int[it_a.multiplicity()];
						for(int i = it_a.istart(), j=0; i<it_a.iend(); i++, j++){
							poss_a[j] = sa_a[i];
							nss_a[j] = ns_a[i];
						}
						igtools.common.util.Arrays.keyedQsort(poss_a, nss_a);
					
						int[] sa_b = it_b.elsa().sa();
						int[] ns_b = ((INELSA)it_b.elsa()).ns();
						int[] poss_b = new int[it_b.multiplicity()];
						int[] nss_b = new int[it_b.multiplicity()];
						for(int i = it_b.istart(), j=0; i<it_b.iend(); i++, j++){
							poss_b[j] = sa_b[i];
							nss_b[j] = ns_b[i];
						}
						igtools.common.util.Arrays.keyedQsort(poss_b, nss_b);
						
						int count = 0;
						for(int i=0; i<poss_a.length; i++){
							for(int j=0; j<poss_b.length; j++){
								if(poss_a[i] >= poss_b[j]){
									if(poss_a[i] - poss_b[j] < nss_b[j])
									if(poss_a[i] - poss_b[j] <= maxDist)
									if(poss_a[i] - poss_b[j] > k_b)
										count++;
								}
								else{
									if(poss_b[j] - poss_a[i] < nss_a[i]){
										if(poss_b[j] - poss_a[i] > k_a){
											if(poss_b[j] - poss_a[i] <= maxDist){
												count++;
											}
											else
												break;
										}
									}
								}
							}
						}
						
						int[] dists = new int[count];
						
						
						int d = 0;
						for(int i=0; i<poss_a.length; i++){
							for(int j=0; j<poss_b.length; j++){
								if(poss_a[i] >= poss_b[j]){
									if(poss_a[i] - poss_b[j] < nss_b[j])
									if(poss_a[i] - poss_b[j] <= maxDist)
									if(poss_a[i] - poss_b[j] > k_b){
										dists[d] = poss_a[i] - poss_b[j];
										d++;
									}
								}
								else{
									if(poss_b[j] - poss_a[i] < nss_a[i]){
										if(poss_b[j] - poss_a[i] > k_a){
											if(poss_b[j] - poss_a[i] <= maxDist){
												dists[d] = poss_b[j] - poss_a[i];
												d++;
											}
											else
												break;
										}
									}
								}
							}
						}
						
						
						return dists;
					}
				};
			}
			else{
				//nooverlap noskip
				return new LimitedBiGlobalRecurrenceDistancesExtractor() {
					
					@Override
					public int[] recurrencesDistances(IELSAIterator it_a, IELSAIterator it_b) {
//						System.out.println("recurrences,bi,global,nooverlap,noskip");
						
						int k_a = it_a.k()-1;
						int k_b = it_b.k()-1;
						
						int[] pos_a = it_a.sortedPositions();
						int[] pos_b = it_b.sortedPositions();
						
//						System.out.println(pos_a+" "+pos_b);
						
						int count = 0;
						for(int i=0; i<pos_a.length; i++){
							for(int j=0; j<pos_b.length; j++){
								if(pos_a[i] >= pos_b[j]){
									if(pos_a[i] - pos_b[j] <= maxDist)
									if(pos_a[i] - pos_b[j] > k_b)
										count++;
								}
								else{
									if(pos_b[j] - pos_a[i] <= maxDist){
										if(pos_b[j] - pos_a[i] > k_a){
											count++;
										}
									}
									else
										break;
								}
							}
						}
						
						int[] dists = new int[count];
						
						
						int d = 0;
						for(int i=0; i<pos_a.length; i++){
							for(int j=0; j<pos_b.length; j++){
								if(pos_a[i] >= pos_b[j]){
									if(pos_a[i] - pos_b[j] <= maxDist)
									if(pos_a[i] - pos_b[j] > k_b){
										dists[d] = pos_a[i] - pos_b[j];
										d++;
									}
								}
								else{
									if(pos_b[j] - pos_a[i] <= maxDist){
										if(pos_b[j] - pos_a[i] > k_a){
											dists[d] = pos_b[j] - pos_a[i];
											d++;
										}
									}
									else
										break;
								}
							}
						}
						
						
						return dists;
					}
				};
			}
		}
	}
}
