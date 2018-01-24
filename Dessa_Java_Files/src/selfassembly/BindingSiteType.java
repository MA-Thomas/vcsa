package selfassembly;

import java.util.Set;
import java.util.HashMap;
import javax.vecmath.*;

/**
 * Stores the information which forms a BindingSiteType. These are
 * Position, tolerances, name, partners, and binding angles
 * 
 * @author Rori Rohlfs
 * @author Blake Sweeney
 * @version 1.4
 */
public class BindingSiteType {
	
	/** The position of this BindingSiteType */
	private Vector3d bSTPosition;  
	
	/** The tolerances */
    private double[] tolerances;    
    
    /** The name */
    private String myName;         
   
    private HashMap<String, PartnerProperties> pMap;//Partner Properties
    
    /**
     * Constructs a BindingSiteType from the given parameters
     * 
     * @param bstPos - The position of the BindingSiteType
     * @param tol - The tolerances of the BindingSiteType
     * @param name - The name of the BindingSiteType
     * @param part - A vector of the partners of the BindingSiteType
     * @param ba - A map of the binding angles for the BindingSiteType
     */
	public BindingSiteType(Vector3d bstPos, double[] tol, String name, HashMap<String, PartnerProperties> partnerPropreties) {
		
		bSTPosition = bstPos;
		
		tolerances = new double[3];
		for (int i = 0; i < 3; i++)
			tolerances[i] = tol[i];
		
		myName = name;
		pMap = partnerPropreties;
	}
	
	public HashMap<String, PartnerProperties> getPartnerMap(){ return pMap;}
		
	public final double getBindTime(String partnerName) {
		return pMap.get(partnerName).getBindTime();
	}

	public final double getBreakTime(String partnerName) {
		return pMap.get(partnerName).getBreakTime();
	}

	public final double getFastBindTime(String partnerName) {
		return pMap.get(partnerName).getFastBindTime();
	}

	/**
	 * Returns the Position of this BindingSiteType
	 * 
	 * @return The BindingSiteType's position
	 */
	public Vector3d getBSTPosition() 
	{ return bSTPosition; }
	
	/**
	 * Returns the name of this BindingSiteType
	 * 
	 * @return the name of the BindingSiteType
	 */
	public String getName() 
	{ return myName; }
	
	/**
	 * Returns a vector of the names of all partners this BindingSiteType has
	 * 
	 * @return partners
	 */
	public Set<String> getPartners() 
	{ return pMap.keySet(); }
	
	/**
	 * Returns a double[3] representing the tolerances of this BindingSiteType
	 * 
	 * @return tolerances
	 */
	public double[] getTolerance() 
	{ return tolerances; }
	
	/**
	 * Returns the binding angle for the given partner
	 * 
	 * @param partnerName - A string of the partner name
	 * @return the binding error
	 */
    public double getBindingAngle(String partnerName)    
    { return pMap.get(partnerName).getAngle(); }
    
    
    /**
     * Checks to see if the given BindingSiteType is compatible
     * 
     * @param bst - the BindingSiteType to check
     * @return true if it is compatible
     */
    public boolean isCompatible(String bstName) 
    { return pMap.containsKey(bstName); }
    
    /**
     * This tests if two BindingSiteTypes are equal by testing if they have
     * the same name.
     * 
     * @return boolean - true if the two BindingSiteType names are equal, false otherwise.
     */
    public boolean equals(BindingSiteType bst)
    { return myName.equals(bst.getName()); }
    
    /**
     * A toString method for debugging
     */
    public String toString() 
    {    	
    	StringBuffer buff = new StringBuffer();
    	buff.append("BSTpos ");
    	buff.append(getBSTPosition());
    	buff.append("\n");
    	return buff.toString();
    }
}