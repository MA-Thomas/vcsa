package selfassembly;


import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;
import javax.vecmath.*;


/**
 * This class holds the information of <code>Subunit</code>s. Each 
 * <code>Subunit</code> has a type, a set of domains, an ID, and 
 * binding sites.
 * 
 * @author Rori Rolfhs
 * @author Peter Kim
 * @author Tiequan Zhang
 * @author Blake Sweeney
 * @version 1.4
 *
 *
 */
public class Subunit {

	/** The type of this <code>Subunit</code>*/
	private SubunitType mySubunitType;

	/** The ID of this <code>Subunit</code>*/
	private int myID;

	/** The position*/
	private Vector3d myPositionReal;

	/** The up vector */
	private Vector3d upVectorPosReal;

	/** All subunits bound to this on */
	private Vector<Subunit> myPartners;

	/** The assembly this <code>Subunit</code> belongs to */
	private Assembly myAssembly;

	/** The former position of this Subunit */
	private Vector3d formerPos;

	/** The derivative in x for Spring Force models */
	private double xDeriv;

	/** The derivative in y for Spring Force models */
	private double yDeriv;

	/** The derivative in z for Spring Force models */
	private double zDeriv;

	/** The derivative of rotation in x for Spring Force models */
	private double xRotDeriv;

	/** The derivative of rotation in y for Spring Force models */
	private double yRotDeriv;

	/** The derivative of rotation in z for Spring Force models */
	private double zRotDeriv;


	/**
	 * A constructor.
	 * 
	 * @param type - The SubunitType
	 * @param id - The ID
	 * @param pos - The position
	 * @param assem - The Assembly of this Subunit
	 */
	public Subunit (SubunitType type, int id, Vector3d pos, Assembly assem)
	{
		mySubunitType = type;
		upVectorPosReal = new Vector3d(type.getUpVec());
		myID = id;
		myPositionReal = pos;
		myAssembly = assem;
		myPartners = new Vector<Subunit>();
	}

	/**
	 * A constructor. Sets the Assembly of this subunit to null
	 * 
	 * @param type - The SubunitType
	 * @param id - The ID
	 * @param pos - The position of this Subunit
	 */
	public Subunit (SubunitType type, int id, Vector3d pos) {

		mySubunitType = type;
		upVectorPosReal = new Vector3d(type.getUpVec());
		myID = id;
		myPositionReal = new Vector3d(pos);
		myAssembly = null;
		myPartners = new Vector<Subunit>();
	}



	/**
	 * Returns a vector of all unbound BindingSites.
	 * 
	 * @return Vector<BindingSite> - A vector of free BindingSites
	 */
	public Vector<BindingSite> getFreeBindingSites() {

		Vector<BindingSite> temp = new Vector<BindingSite>();
		Vector<Domain> domains = mySubunitType.getDomains();
		int aSize = domains.size();
		int bSize;		
		for (int i = 0; i < aSize; ++i) {
			Vector<BindingSite> sites = domains.get(i).getBindingSites();
			bSize = sites.size();
			for (int j = 0; j < bSize; ++j) {
				BindingSite bs = sites.get(j);
				if (!bs.isBound())
					temp.add(bs);
			}
		}
		return temp;
	}



	/**
	 * Changes the conformation of this <code>Subunit</code> to a new one 
	 * if the <code>Subunit</code>
	 * contains this domain. Currently will not change conformations if
	 * the <code>Subunit</code> is bound.
	 * 
	 * @param dID - the ID of the domain to change
	 * @param conf - the new conformation to change to
	 * @return boolean - true if this <code>Subunit</code> is changed to the new
	 * 	conformation.
	 */
	public final boolean changeConf(int dID, Conformation conf) 
	{
		return mySubunitType.changeConf(conf, dID);
	}



	/**
	 * Returns the Assembly of this <code>Subunit</code>
	 * 
	 * @return Assembly - The assembly of this <code>Subunit</code>
	 */
	public Assembly getAssembly() 
	{  return myAssembly; }



	/**
	 * Sets the Assembly of this <code>Subunit</code> to the specified Assembly. Also
	 * sets the <code>Subunit</code> of all the BindingSites to this to update their 
	 * getAssembly().
	 * 
	 * @param a - The new assembly of this subunit
	 */
	public void setAssembly(Assembly a) { 

		myAssembly = a;

		Vector<BindingSite> sites = getAllBindingSites();
		int size = sites.size();
		for (int i = 0; i < size; ++i)
			sites.get(i).setSub(this);
	}



	/**
	 * Returns the domains of this <code>Subunit</code>.
	 * 
	 * @return Vector<Domain> - A vector of all Domains in this <code>Subunit</code>.
	 */
	public Vector<Domain> getDomains() 
	{ return mySubunitType.getDomains(); }


	/**
	 * Gets all BindingSites in all current Conformations of all Domains
	 * 
	 * @return Vector<BindingSite> - a vector of all available BindingSites
	 */
	public Vector<BindingSite> getBindingSites() {

		Vector<Domain> domains = mySubunitType.getDomains();
		Vector<BindingSite> bs = new Vector<BindingSite>();
		int s = domains.size();
		for (int i = 0; i < s; ++i) {

			Conformation conf = domains.get(i).getCurConf();
			bs.addAll(conf.getBindSites());
		}

		return bs;
	}



	/**
	 * Gets all BindingSites in all Conformations of all Domains
	 * 
	 * @return Vector<BindingSite> - a vector of all available BindingSites
	 */
	public Vector<BindingSite> getAllBindingSites() {

		Vector<Domain> domains = mySubunitType.getDomains();
		Vector<BindingSite> bs = new Vector<BindingSite>();
		int aSize = domains.size();
		for (int i = 0; i < aSize; ++i) {
			HashMap<String,Conformation> confs = domains.get(i).getConfs();
			Iterator<Conformation> cnfItr = confs.values().iterator();
			while(cnfItr.hasNext())
			{
				bs.addAll(cnfItr.next().getBindSites());
			}
		}

		return bs;
	}


	/**
	 * Returns a Vector of all active conformations.
	 * 
	 * @return Vector<Conformation> - all active conformations.
	 */
	public Vector<Conformation> getConfs() {

		Vector<Conformation> confs = new Vector<Conformation>();
		Vector<Domain> myDomains = mySubunitType.getDomains();
		int s = myDomains.size();
		for (int i = 0; i < s; ++i)
			confs.add(myDomains.get(i).getCurConf());

		return confs;
	}



	/**
	 * Returns the ID of this <code>Subunit</code>
	 * 
	 * @return int - the ID
	 */
	public int getID() 
	{ return myID; }



	/**
	 * Returns the name of the SubunitType
	 * 
	 * @return String - The name of the SubunitType
	 */
	public String getSubunitTypeName() 
	{ return mySubunitType.getName(); }



	/**
	 * Returns the SubunitType of this <code>Subunit</code>.
	 * 
	 * @return SubunitType - the <code>SubunitType</code>
	 */
	public SubunitType getSubunitType() 
	{ return mySubunitType; }



	/**
	 * Returns the real position of this <code>Subunit</code>.
	 * 
	 * @return Vector3d - the real position of this subunit
	 */
	public Vector3d getPositionReal() 
	{  return myPositionReal; }

	/** ADDED by MT Feb 24 2016
	 * Returns the real up vector of this <code>Subunit</code>.
	 * 
	 * @return Vector3d - the real position of this subunit
	 */
	public Vector3d getUpVectorPosReal()
	{  return upVectorPosReal; }

	

	/**
	 * Returns a vector of all <code>Subunit</code>s bound to this one.
	 * 
	 * @return Vector<Subunit> - A vector of all bound <code>Subunit</code>s.
	 */
	public Vector<Subunit> getBoundSubunits() 
	{ return myPartners; }



	/**
	 * Returns true if this <code>Subunit</code> has no partners
	 * 
	 * @return boolean - true if unbound
	 */
	public boolean isUnbound() 
	{ return (myPartners.size() == 0); }



	/**
	 * Adds a <code>Subunit</code> as a partner if it is not already
	 * a partner.
	 * 
	 * @param s - <code>Subunit</code> to bind
	 * @return boolean - true if partner is added
	 */
	public boolean addPartner(Subunit s) {

		int i = 0;
		int size = myPartners.size();
		for (; i < size && !myPartners.get(i).equals(s); ++i);

		if (i == myPartners.size()) {

			if (myPartners.size() + 1 > getBindingSites().size()) {
				System.out.println("Error too many partners to this subunit");
				System.out.println(this);
				
				//m for debugging, 7-09-2015
				System.out.println("\n" + Test.seedUsed + "\n");
				//m end debugging
				
				System.exit(-1);
			}

			myPartners.add(s);
			return true;
		}

		return false;
	}



	/**
	 * Returns true if this <code>Subunit</code> contains the <code>Subunit</code> as a partner
	 * and removes it.
	 * 
	 * @param s - the <code>Subunit</code> to remove
	 * @return boolean - true if the <code>Subunit</code> is removed successfully
	 */
	public boolean removePartner(Subunit s) {

		Vector<BindingSite> bss = getBindingSites();//getCurrent Binding Sites
		int size = bss.size();
		for(int i =0; i < size ; ++i)
		{
			BindingSite bs = bss.get(i);
			if(bs.isBound())
			{
				if(bs.getPartner().getSubunitID()==s.getID()){
					/*Some other binding still exist
					 * b/w this and Subunit s,
					 * THIS SHOULD NEVER HAPPEN! 
					 */					 
					System.out.println("Seed:"+Test.seedUsed);
					System.out.println("BS:"+bs.toString());
					System.out.println("PartnerBS:"+bs.getPartner().toString());
					System.out.println("Sub:"+this.toString());
					System.out.println("PartnerSub:"+s.toString());					
					return false;
				}
			}
		}


		for (int i = 0; i < myPartners.size(); ++i) {
			if (myPartners.get(i).equals(s)) {
				myPartners.remove(i);
				return true;
			}
		}

		return false;
	}



	/**
	 * Returns the up vector of the <code>Subunit</code> relative to itself.
	 * 
	 * @return Vector3d - the up vector relative to this <code>Subunit</code>
	 */
	public Vector3d getUpRelative() {

		Vector3d upRelative = new Vector3d();
		upRelative.sub(upVectorPosReal, myPositionReal);

		return upRelative;
	}



	/**
	 * Moves this <code>Subunit</code>.
	 */
	public void subMove(Vector3d ref, Vector3d translation, Quat4d q) {

		Matrix4d m = new Matrix4d();
		m.set(q);
        
		//MT: rotate myPositionReal, i.e. m.transform()
		myPositionReal.sub(ref);
		m.transform(myPositionReal);
		myPositionReal.add(ref);
		
		//MT: now apply translation
		myPositionReal.add(translation);

		//MT: same for upVector
		upVectorPosReal.sub(ref);
		m.transform(upVectorPosReal);
		upVectorPosReal.add(ref);
		upVectorPosReal.add(translation);

		//MT: same for binding sites
		Vector<BindingSite> myBSs = getAllBindingSites();//this.getBindingSites();
		int size = myBSs.size();
		for (int i = 0; i < size; ++i) {
			BindingSite bs = myBSs.get(i);

			Vector3d bsPosReal = bs.getPosReal();

			bsPosReal.sub(ref);

			m.transform(bsPosReal);

			bsPosReal.add(ref);
			bsPosReal.add(translation);
			bs.setPosReal(bsPosReal);

		}
	}


	/**
	 * Gets the energy of this Subunit based upon the sum of the energies
	 * of its BindingSites
	 * 
	 * @return double - the energy of this Subunit
	 */
	public double getEnergy() {

		double energy = 0.0;
		Vector<BindingSite> myBindingSites = getBindingSites();
		int size = myBindingSites.size();
		for(int i = 0; i < size; ++i) 
			energy += myBindingSites.get(i).getEnergy();

		return energy;
	}



	/**
	 * Moves a Subunit and all of its BindingSites
	 * 
	 * @param trans - The Vector to move the Subunit by
	 */
	public void moveSub(Vector3d trans) {

		myPositionReal.add(trans);
		Vector<BindingSite> myBindingSites = getAllBindingSites();
		int size = myBindingSites.size();
		for(int i = 0; i < size; ++i) 
			myBindingSites.get(i).getPosReal().add(trans);
	}



	public void rotateSub(AxisAngle4d rot) {

		Matrix4d m = new Matrix4d();
		m.set(rot);
		Vector<BindingSite> myBindingSites = getAllBindingSites();
		int size = myBindingSites.size();
		for(int i = 0; i < size; ++i) {

			BindingSite bs = myBindingSites.get(i);
			Vector3d bsPos = bs.getRelativePos();

			m.transform(bsPos, bsPos);
			bsPos.add(myPositionReal);

			bs.setPosReal(bsPos);
		}
	}



	/**
	 * Saves the current position of this Subunit as well as the position of
	 * it's BindingSites
	 */
	public void saveCurrentPos() {

		formerPos = new Vector3d(myPositionReal);
		Vector<BindingSite> myBindingSites = getAllBindingSites();
		int size = myBindingSites.size();
		for(int i = 0; i < size; ++i) 
			myBindingSites.get(i).saveCurrentPos();
	}



	/**
	 * Restore the position of this Subunit to it's old position as well
	 * as the position of all of it's BindingSites.
	 */
	public void restorePos() {

		myPositionReal = formerPos;
		Vector<BindingSite> myBindingSites = getAllBindingSites();
		int size = myBindingSites.size();
		for(int i = 0; i < size; ++i)  
			myBindingSites.get(i).restorePos();
	}



	public void applySpringForce(double kVal) {

		Vector3d translation = new Vector3d(xDeriv * kVal, yDeriv * kVal, zDeriv * kVal);

		AxisAngle4d rotation;

		if(xRotDeriv != 0 || kVal != 0) {
			rotation= new AxisAngle4d(1,0,0, xRotDeriv * kVal);
			rotateSub(rotation);
		}

		if(xRotDeriv != 0 || kVal != 0) {
			rotation = new AxisAngle4d(0,1,0, yRotDeriv * kVal);
			rotateSub(rotation);
		}

		if(xRotDeriv != 0 || kVal != 0) {
			rotation = new AxisAngle4d(0,0,1, zRotDeriv * kVal);
			rotateSub(rotation);
		}
	
		moveSub(translation);	
	}



	public void setXDeriv(double x)
	{ xDeriv = x; }



	public void setYDeriv(double y)
	{ yDeriv = y; }



	public void setZDeriv(double z)
	{ zDeriv = z; }



	public void setXRotDeriv(double x) 
	{ xRotDeriv = x; }



	public void setYRotDeriv(double y) 
	{ yRotDeriv = y; }



	public void setZRotDeriv(double z) 
	{ zRotDeriv = z; }



	/**
	 * A toString method for debugging
	 */
	public String toString() {

		StringBuffer buff = new StringBuffer();

		buff.append("ID: ");
		buff.append(myID);
		buff.append(" Type: ");
		buff.append(mySubunitType);
		buff.append(" Assembly: ");
		buff.append(myAssembly.getID());


		buff.append(" Partners: ");
		int size = myPartners.size();
		for (int i = 0; i < size; ++i)
			buff.append(myPartners.get(i).getID()+" ");
		buff.append("\n");

		return buff.toString();
	}
}
