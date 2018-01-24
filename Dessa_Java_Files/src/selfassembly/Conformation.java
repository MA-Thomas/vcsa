package selfassembly;

import java.util.Vector;


public class Conformation {

    /** BindingSites present */
    private Vector<BindingSite> mySites; 
    
    /** Conformation name */
    private String myName; 
    
    /** Energy in this conformation (not used presently) */
    private double myEnergy;
    
    /**
     * Constructs a Conformation with specified BindingSite Vector bs, energy e,
     * and name. Sets the ID of SubunitTypeID and DomainID associated to -1.
     * 
     * @param bs Vector of BindingSite
     * @param e  double energy
     * @param n  String name of this Conformation
     */
    public Conformation(Vector<BindingSite> bs, double e, String n) {    	
        myName = n;
        mySites = bs;
        myEnergy = e;
    }    
    
    /**
     * Constructs a new conformation as a copy of an old one.
     * 
     * @param conf - the Conformation to copy
     */
    public Conformation(Conformation conf) {
    	
    	myName = new String(conf.getName());
    	myEnergy = conf.getE();
    	
    	mySites = new Vector<BindingSite>();
    	Vector<BindingSite> sites = conf.getBindSites();
    	int size = sites.size();
    	for (int i = 0; i < size; ++i)
    		mySites.add(new BindingSite(sites.get(i)));
    }
    
    
    
    /**
     * Returns the binding sites present in this conformation
     * 
     * @return Vector
     */
    public Vector<BindingSite> getBindSites() 
    { return mySites; }

    
    
    /**
     * Checks if any BindingSite in this Conformation is bound, if so, returns
     * true
     * 
     * @return boolean
     */
    public boolean bound() {
    	int size = mySites.size();
        for (int i = 0; i < size; ++i) {
            if (mySites.get(i).isBound())
                return true;
        }
        return false;
    }

    
    
     
    /**
     * Returns the name of this Conformation
     * 
     * @return String
     */
    public String getName() 
    { return myName; }

    
    
    /**
     * Returns of energy of this Conformation
     * 
     * @return double
     */
    public double getE() 
    { return myEnergy; }

    
    /**
     * Returns a String version of this Conformation
     * 
     * @return String
     */
    public String toString() {
    	StringBuffer buff = new StringBuffer();    	
    	buff.append("Conformation Name: ");
    	buff.append(myName);
    	buff.append(" Energy: ");
    	buff.append(myEnergy);
    	buff.append("\n\t BindingSites:\n");
    	buff.append(mySites);
    	
    	return buff.toString();
    }
}
