package selfassembly;



import javax.vecmath.Vector3d;


/**
 * This class stores the information for a <code>BindingSite</code>. This
 * includes the partner, the <code>BindingSiteType</code> the position in
 * the simulation and relative to the <code>Subunit</code>, the ID and the 
 * <code>Subunit</code> this <code>BindingSite</code> is part of. 
 * 
 * @author Rori Rohlfs
 * @author Tiequan Zhang
 * @author Peter Kim
 * @author Blake Sweeney
 * @version 1.4
 */
public class BindingSite {
	
	/** The partner of this <code>BindingSite</code>, null if none */
	private BindingSite partner; 			
	 	
	private String myBSType;
	/** Position in overall simulation */
    private Vector3d posReal; 			
  
    /** The former position of this BindingSite */
    private Vector3d formerPos;
    
    /** The ID of this <code>BindingSite</code> */
    private int myID; 					
    
    /** The <code>Subunit</code> of this <code>BindingSite</code> */
    private Subunit mySub;
    
    
    
    /**
     * A constructor.
     * 
     * @param type - The <code>BindingSiteType</code> of this <code>BindingSite</code>
     * @param id - The ID for this <code>BindingSite</code>
     * @param sub - The <code>Subunit</code> this <code>BindingSite</code> is a part of
     * @param pos - The position of this <code>BindingSite</code>
     */
    public BindingSite(String type, int id, Subunit sub, Vector3d pos) {
    	
    	myBSType = type;
    	myID = id;
    	mySub = sub;
    	partner = null;
    	posReal = new Vector3d(pos);
    	formerPos = new Vector3d();
    }

   
    /**
     * A constructor. The <code>Subunit</code> this is a part of is set to null
     * 
     * @param type - The <code>BindingSiteType</code> of this <code>BindingSite</code>
     * @param id - The ID for this <code>BindingSite</code>
     * @param pos - The position of this <code>BindingSite</code>
     */
    public BindingSite(String type, int id, Vector3d pos) {
    	
    	myBSType = type;
    	myID = id;
    	mySub = null;
    	formerPos = new Vector3d();
    	partner = null;
    	posReal = new Vector3d(pos);
    }
    
    
    
    /**
     * Constructs a copy of a <code>BindingSite</code> based off of another one.
     * 
     * @param bs - The <code>BindingSite</code> to copy
     */
    public BindingSite(BindingSite bs) 
    {
    	myID = bs.getID();
    	mySub = bs.getSubunit();
    	partner = bs.getPartner();
    	myBSType = new String(bs.getBSTName());
    	posReal = new Vector3d(bs.getPosReal());
    	formerPos = new Vector3d(bs.getFormerPos());
    }
    
    public final Vector3d getFormerPos(){ return formerPos;}
    
	/**
	 * Returns the ID of this BindingSite's Subunit
	 * @return int - the ID of this BindingSite's Subunit
	 */
	public int getSubunitID() 
	{ return mySub.getID(); }
	
	/**
	 * Sets the real position of this BindingSite
	 * 
	 * @param pos - the position to set this BindingSite too.
	 */
	public void setPosReal(Vector3d pos) 
	{ posReal.set(pos); }
		
	/**
	 * Returns the real position of this BindingSite.
	 * @return the position of this BindingSite in the overall simulation
	 */
    public Vector3d getRelativePos() 
    {	
    	Vector3d pos = new Vector3d(posReal);
    	pos.sub(getSubunit().getPositionReal());
    	return pos;
    }
	
	/**
	 * Gets the Assembly of this BindingSite.
	 * @return Assembly - The Assembly this BindingSite is part of
	 */
	public Assembly getAssembly() 
	{ return mySub.getAssembly(); }
		
	/**
	 * Returns the Subunit of this BindingSite
	 * @return Subunit - This BindingSite's Subunit
	 */
	public Subunit getSubunit() 
	{ return mySub; }
	
	/**
	 * This breaks the bond between two subunits. Returns true upon success, and 
	 * false upon failure. Fails if the BindingSite is not bound.
	 * 
	 * @return boolean - true if breaking the bond is successful
	 */
	public boolean breakBond() {
		
		if (partner!=null) {//if is bound!			
			partner.partner = null;
			partner = null;			
			return true;
		}
		return false;
	}
	
	/**
	 * Binds this subunit to another. Returns true if binding successful. False
	 * if it fails. This fails if the binding site is already bound.
	 * 
	 * @param - BindingSite - the BindingSite to bind to.
	 * @return boolean - true if binding successful false if it fails.
	 */
	public boolean bindTo(BindingSite bs) {
		
		if (partner==null && !bs.isBound()) {
			
			partner = bs;
			bs.partner = this;
			
			return true;
		}
		return false;
	}
	
	
	/**
	 * Returns the real position of this subunit.
	 * 
	 * @return Vector3d - A vector of the real position of this subunit
	 */
	public Vector3d getPosReal() 
	{ return posReal; }
	
	/**
	 * Returns the ID of this BindingSite
	 * 
	 * @return int - the ID
	 */
	public int getID() 
	{ return myID; }
	
	
	/**
	 * Changes the ID of this <code>BindingSite</code> to the specified one.
	 * 
	 * @param id - the ID to give this <code>BindingSite</code>
	 */
	public void setID(int id) 
	{ myID = id; }
	
	
	
	/**
	 * Returns the BindingSiteType of this BindingSite
	 * 
	 * @return BindingSiteType - the BindingSiteType
	 */
	public String getBSTName() 
	{ return myBSType; }
	
	
	
	/**
	 * Returns the partner of this subunit. Null if none.
	 * 
	 * @return BindingSite - this BindingSite's partner
	 */
	public BindingSite getPartner() 
	{ return partner; }

	
	
	/**
	 * Returns true if this BindingSite is bound.
	 * 
	 * @return boolean - true if the BindingSite is bound.
	 */
	public boolean isBound() 
	{ return partner != null; }
	
	
	
	/**
	 * Sets the SubunitType for this BindingSite
	 * 
	 * @param sub - The Subunit to make this BindingSite part of.
	 */
	public void setSub(Subunit sub) 
	{ mySub = sub; }
	
	/**
	 * An equals method. Tests based upon the BindingSite ID
	 * 
	 * @param bs - The BindingSite to compare to
	 * @return boolean - true of the BindingSite ID's are equal
	 */
	public boolean equals(BindingSite bs)
	{ return myID == bs.getID(); }
	
	
	/**
	 * Gets the energy of this BindingSite based upon it's Translational,
	 * Bending and Rotational Energy
	 * 
	 * @return double - the energy of this BindingSite
	 */
	public double getEnergy() 
	{ 
		if(partner == null)//if not bound return zero
    		return 0;
		else
			return (getTransEnergy() + getBendingEnergy() + getRotationEnergy()); 
	}
	
	
	/**
	 * Gets the Translational energy of this BindingSite
	 * 
	 * @return double - the Translational energy
	 */
    private double getTransEnergy() {
    	
    	Vector3d myPos = new Vector3d(posReal);
    	myPos.sub(partner.posReal);
    	
    	double val = myPos.length();
    	return (val * val);
    }
    
    
    /**
     * Gets the Bending energy for this BindingSite
     * 
     * @return double - the Bending energy
     */
    private double getBendingEnergy() {
    	
    	Vector3d vec = new Vector3d(posReal);
    	vec.sub(mySub.getPositionReal());
    	
    	Vector3d bindVec = new Vector3d(partner.posReal);
    	bindVec.sub(partner.mySub.getPositionReal());
    	
    	double val = Math.sqrt(0.5 + 0.5*(vec.dot(bindVec)));
    	
    	Vector3d cross = new Vector3d();
    	cross.cross(vec, bindVec);
    	cross.scale(val);
    	val = cross.length();

    	return val;
    }
    
    
    
    /**
     * Gets the rotational energy of this BindingSite
     * 
     * @return double - the Rotational energy
     */
    private double getRotationEnergy() {
    	
    	
    	Vector3d upVec = new Vector3d(mySub.getUpRelative());
    	upVec = Test.makePerpendicular(Simulation.getBSTPostion(myBSType), upVec);
    	    	
    	Vector3d bindupVec = new Vector3d(partner.mySub.getUpRelative());
    	bindupVec = Test.makePerpendicular(Simulation.getBSTPostion(partner.getBSTName()), bindupVec);

    	upVec.cross(upVec, bindupVec);
    	Vector3d bs = new Vector3d(posReal);
    	bs.sub(mySub.getPositionReal());
    	bs.scale(bs.dot(upVec));
    	double val = bs.length();
    	
    	return val * 0.0000001;
    }
    
    
    
    /**
     * Saves the current position of this BindingSite
     */
    public void saveCurrentPos()
    { formerPos.set(posReal); }
    
    
    
    /**
     * Restores the position of this BindingSite to the previously saved one
     */
    public void restorePos()
    { posReal.set(formerPos); }
    
    
    
	/**
	 * A toString method useful for debugging
	 */
	public String toString() {
		
		   StringBuffer buff = new StringBuffer();
		   buff.append("\n");
		   buff.append("AName ");
		   buff.append(getAssembly().getID());
		   buff.append("\n");
		   buff.append("posReal() ");
		   buff.append(getPosReal());
		   buff.append("\n");
		   buff.append("relativePos() ");
		   buff.append(getRelativePos());
		   buff.append("\n");
		   
		   return buff.toString();
	}
}
