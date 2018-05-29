package igtools.analyses.elongation;


import igtools.common.sequence.B3Sequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

/**
 * Abstract model for elongate words and select them by a measurement
 * @author vbonnici
 *
 */
public abstract class ElongationAlgorithm {

	
	
	protected NELSA nelsa = null;
	protected B3Sequence b3seq = null;
	
	
	public ElongationAlgorithm(NELSA nelsa){
		this.nelsa = nelsa;
		this.b3seq = nelsa.b3seq();
	}
	
	
	/**
	 * EXTENDED: the elongation is extended
	 * MAXIMAL: the elongation cannot go forward
	 * 
	 * @author vbonnici
	 *
	 */
	public enum ElongationEvent{EXTENDED, MAXIMAL};
	

	public interface ElongationListener{
		public void event(IELSAIterator it, ElongationEvent ev, Double parentValue, Double itValue);
	}
	
	public interface IteratorMeasure{
		/**
		 * 
		 * @param it
		 * @return a value of measure or Double.NAN;
		 */
		public abstract Double value(IELSAIterator it);
	}
	
	public abstract void elongate(int startK, IteratorMeasure im, ElongationListener el);
}
