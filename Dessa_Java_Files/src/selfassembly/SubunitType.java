package selfassembly;

import java.util.Vector;
import javax.vecmath.Vector3d;;



/**
 * This is a class defining the SubunitType data structure.
 * A SubunitType is a type of Subunits with common Conformations, Domains and
 * BindingSites
 * 
 * @author Rori Rohlfs
 * @author Blake Sweeney
 * @version 1.4
 */
public class SubunitType {
	
	/** Vector of all Domains */
	private Vector<Domain> myDomains; 
	
	/** The name of this SubunitType */
	private String myName;    
	
	/** The mass */
	private final double mymass;    
	
	/** The radius */
	private final double myradius;  
    
	/** The up vector */
	private Vector3d myUpVec; 
    
    
    
    /**
     * Constructs a SubunitType from the given parameters
     * 
     * @param domains - A vector of all domains
     * @param name - A String of the SubunitType name
     * @param mass - The mass of this SubunitType
     * @param radius - The radius of this SubunitType
     * @param up - A Vector3d of the up vector
     */
    public SubunitType (Vector<Domain> domains, String name, double mass, 
    		double radius, Vector3d up) {
		
		myDomains = domains;
		myName = name;
		mymass = mass;
		myradius = radius;
		myUpVec = up;
    }
    
    
    
    /**
     * Constructs a new SubunitType, based off of an old one.
     * 
     * @param other - The other SubunitType to build a new copy of
     */
    public SubunitType (SubunitType other) {
    	
    	Vector<Domain> temp = other.getDomains();
    	myDomains = new Vector<Domain>();
    	int size =  temp.size();
    	for (int i = 0; i  < size; ++i)
    		myDomains.add(new Domain(temp.get(i)));
    	
    	myName = new String(other.getName());
    	mymass = other.getMass();
    	myradius = other.getRadius();
    	myUpVec = new Vector3d(other.getUpVec());
    }
    
    
    
    /**
     * Returns the radius of this SubunitType
     * 
     * @return radius
     */
    public double getRadius() 
    { return myradius; }

    
  
    /**
     * Returns the name of this SubunitType
     * 
     * @return name - The String of this SubunitType name
     */
    public String getName() 
    { return myName; }
    
    
    
    /**
     * Returns the up vector for this SubunitType
     * 
     * @return upVector
     */
    public Vector3d getUpVec() 
    { return myUpVec; }
    
    
    
    /**
     * Checks to see if the name of the SubunitType is the same, if so
     * the two are considered equal
     * 
     * @param s - SubunitType to compare
     * @return boolean - if the two SubunitTypes are equal
     */
    public boolean equals(SubunitType s) 
    { return myName.equals(s.getName()); }
    
    
    
    /**
     * Returns the mass of this SubunitType
     * 
     * @return mass - the Mass of the SubunitType
     */
    public double getMass() 
    { return mymass; }
    
    
    
    /**
     * Returns a Vector of the Domains of this SubunitType
     * 
     * @return domains
     */
    public Vector<Domain> getDomains() 
    { return myDomains; }
    
    
    
    /**
     * Checks to see if conformation c is a part of Domain di in this SubunitType
     * 
     * @param c Conformation
     * @param di Domain id
     * @return true if c is a part of domain di in this SubunitType
     */
    public boolean changeConf(Conformation c, int di) {		
		
		int s = myDomains.size();
		for(int i =0; i < s; ++i)
		{
			Domain d = myDomains.get(i);
			if(d.getDomainId() == di)
			{
				if(!d.getCurConf().bound())
					return d.changeConf(c.getName());
				else
					return false;
			}
		}
		
		return false;
    }
    
    
    
    /**
     * A method to set the up vector of this SubunitType
     * 
     * @param vec - The Vector3d to give this SubunitType
     */
    public void setUpVector(Vector3d vec)
    { myUpVec = vec; }
    
    
    
    /**
     * A toString method for debugging.
     */
    public String toString() {
    	
    	StringBuffer buff = new StringBuffer();
    	
    	buff.append("Name: ");
    	buff.append(myName);
    	return buff.toString();
    }
}
