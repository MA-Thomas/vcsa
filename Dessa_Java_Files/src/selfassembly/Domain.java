package selfassembly;


import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;
import javax.vecmath.Vector3d;


/**
 * This class represents a Domain. A Domain is made of multiple possible Conformations
 * with one current Conformation. A domain is in a certain position and is identified by
 * a unique ID. 
 * 
 * @author Rori Rolfhs
 * @author Tiequan Zhang
 * @author Blake Sweeney
 * @version 1.4
 */
public class Domain {

	/** Set of possible Conformations */
	private HashMap<String,Conformation> myConfs; 			

	/** Current Conformation for this Domain */
	private Conformation curConf; 				

	/** The position of this Domain. */
	private Vector3d myPos = new Vector3d();		

	/** Domain id */
	private int dID; 								



	/**
	 * Constructs a Domain with specified parameters (SubunitType ID, Vector of
	 * Conformations, Vector of Positions, Domain id)
	 * 
	 * @param c
	 *            Vector of Conformations
	 * @param pos
	 *            Vector3D postion of this Domain in the Subunit
	 * @param di
	 *            Domain id
	 */
	public Domain(HashMap<String,Conformation> c, Vector3d pos, int di) {

		myConfs = c;
		myPos = pos;
		dID = di;
		curConf = null;
	}



	/**
	 * Constructs a Domain with specified parameters (SubunitType ID, Vector of
	 * Conformations, Vector of Positions, Domain id, current Conformation)
	 * 
	 * @param c
	 *            Vector of Conformations
	 * @param pos
	 *            Vector3D postion of this Domain in the Subunit
	 * @param di
	 *            Domain id
	 * @param cc
	 *            Conformation, current Conformation for this Domain
	 */
	public Domain(HashMap<String,Conformation> c, Vector3d pos, int di, Conformation cc) {

		myConfs = c;
		myPos = pos;
		dID = di;
		curConf = cc;
	}



	/**
	 * Constructs a new Domain based off of another one.
	 * 
	 * @param d - the Domain to copy
	 */
	public Domain(Domain d) {

		myPos = new Vector3d(d.getPos());
		dID = d.getDomainId();
		curConf = new Conformation(d.getCurConf());

		myConfs = new HashMap<String,Conformation>();
		HashMap<String,Conformation> confs = d.getConfs();
		Iterator<Conformation> cnfItr = confs.values().iterator();
		while(cnfItr.hasNext())
		{
			Conformation cnf = cnfItr.next();			
			myConfs.put(cnf.getName(),new Conformation(cnf));
		}
	}



	/**
	 * Returns vector of Conformations for this Domain
	 * 
	 * @return Vector of Conformations
	 */
	public HashMap<String,Conformation> getConfs() 
	{ return myConfs; }



	/**
	 * Returns the position of this Domain in the Subunit
	 * 
	 * @return Vector3d
	 */
	public Vector3d getPos() 
	{ return myPos; }



	/**
	 * Returns the current Conformation for this Domain
	 * 
	 * @return Conformation
	 */
	public Conformation getCurConf() 
	{  return curConf; }



	/**
	 * Returns the DomainID of this Domain.
	 * 
	 * @return int - DomainID
	 */
	public int getDomainId() 
	{ return dID; }


	/**
	 * Changes this Domain's current Conformation to newconf, returns true on
	 * success, false on failure
	 * 
	 * @param newconf
	 *            changing this Domain's current Conformation to newconf
	 * @return boolean
	 */
	public boolean changeConf(String confName) { 

		if (myConfs.containsKey(confName)) {
			curConf = myConfs.get(confName);
			return true;
		}
		return false;
	}


	/**
	 * Returns the BindingSites of the current Conformation
	 * 
	 * @return Vector<BindingSite> - a vector of the current BindingSites
	 */
	public Vector<BindingSite> getBindingSites()
	{ return curConf.getBindSites(); }

	/**
	 * Returns a String version of this Domain
	 * 
	 * @return String
	 */
	public String toString() {
		return new String("Domain #" + dID + "  current Conformation: "
				+ curConf.toString());
	}
}
