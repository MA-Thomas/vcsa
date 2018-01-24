package selfassembly;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Vector;

import java.util.HashMap;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;


/**
 * Assembly holds an Assembly of Subunits and the AssemblyGraphic. It also
 * tracks the next Event for this Assembly and its validity. Assemblies can be
 * bonded, split, or relaxed.
 * 
 * @author Rori Rohlfs
 * @author Sue Yi Chew
 * @author Tiequan Zhang
 * @author Woo Tae Kim
 * @author Blake Sweeney
 * @version 1.4
 *  
 */


public class Assembly {
	/** A vector of the subunits */
	private Vector<Subunit> mysubs;

	/** The valid time for any events with this Assembly */
	private double validTime;

	/** The ID of this Assembly */
	private int myID;
	
	/** MT - Temporarily store rotation data for the subunits 
	 * to be moved during a bind event. Necessary for
	 * tR in Simulation.java
	 */
	public static HashMap<Integer, Quat4d> tR_temp;

	/**
	 * A constructor
	 * 
	 * @param s - A Vector of the Subunits for this Assembly
	 * @param name - The name of this Assembly
	 */
	public Assembly(Vector<Subunit> s, int id) {

		mysubs = new Vector<Subunit>();
		int size = s.size();
		for (int i = 0; i < size; ++i) {
			mysubs.add(s.get(i));
			mysubs.get(i).setAssembly(this);
		}
		myID = id;
		validTime = -1.0;


	}

	/**
	 * Method binds this assembly to another one using a specified event.
	 * Checks for stearic hindrance and fast binding
	 * 
	 * @param bindasm - Assembly to bind to
	 * @param e - Event to use to specify the binding
	 * @return Assembly - the new Assembly composed of the two previous ones
	 */
	public Assembly bindAssembly(Assembly bindasm, Event e) {


		BindingSite oldbs = e.getBS();
		Subunit oldSub = oldbs.getSubunit();		

		BindingSite bindBS = e.getPartner();
		Subunit bindSub = bindBS.getSubunit();		


		Vector3d oldBSPosRelative = new Vector3d();


		oldBSPosRelative.sub(oldbs.getPosReal(), oldSub.getPositionReal());		
		Vector3d bindBSPosRelative = new Vector3d();
		bindBSPosRelative.sub(bindBS.getPosReal(), bindSub.getPositionReal());
		double angle = oldBSPosRelative.angle(bindBSPosRelative);
		Quat4d firstRotation = new Quat4d();
		AxisAngle4d firstAxisAngle4d = null;
		if (angle < 0.000001) {
			Vector3d firstAxis = new Vector3d();
			double temAngle = bindBSPosRelative.angle(new Vector3d(1, 0, 0));
			if ((temAngle < 0.1) || ((Math.PI - temAngle) < 0.1)) {
				firstAxis.cross(bindBSPosRelative, new Vector3d(0, 1, 0));
			} else {
				firstAxis.cross(bindBSPosRelative, new Vector3d(1, 0, 0));
			}
			firstAxisAngle4d = new AxisAngle4d(firstAxis, Math.PI);
		} else if ((Math.PI - angle) < 0.00001) {
			firstAxisAngle4d = new AxisAngle4d(new Vector3d(1, 0, 0), 0);
		} else {
			Vector3d firstAxis = new Vector3d();
			firstAxis.cross(bindBSPosRelative, oldBSPosRelative);
			firstAxis.normalize();

			firstAxisAngle4d = new AxisAngle4d(firstAxis, angle - Math.PI);
		}

		firstRotation.set(firstAxisAngle4d);
        
		//MT: oldProj is perpendicular to first argument on RHS.
		//MT: second argument on RHS effectively sets the length of oldProj.
		//MT: see my notes in Test.java
		//MT: Also, make new upVector by rotating around axis defined by 
		//MT: the 2 binding sites involved.
		Vector3d oldProj = Test.makePerpendicular(oldBSPosRelative, oldSub
				.getUpRelative());
		Vector3d newUpV = Test.rotateByAxis(bindSub.getUpRelative(),
				firstAxisAngle4d);
		Vector3d newProj = Test.makePerpendicular(oldBSPosRelative, newUpV);

		double currentAngleBetweenProjectedVectors = newProj.angle(oldProj);

		double rotAngle = Simulation.getBSTBindAngle(oldbs.getBSTName(),bindBS.getBSTName());
		Vector3d testVector = new Vector3d();
		if (currentAngleBetweenProjectedVectors>0.00000001) {
			testVector.cross(oldProj, newProj);
			testVector.normalize();
			if (Math.abs(oldBSPosRelative.angle(testVector) - Math.PI) <= 0.000001)
				currentAngleBetweenProjectedVectors = -currentAngleBetweenProjectedVectors;
		}
		
		//MT: An axis angle is a rotation of angle (radians) about the vector (x,y,z)
		AxisAngle4d rot = new AxisAngle4d(oldBSPosRelative, rotAngle
				- currentAngleBetweenProjectedVectors);
		Quat4d totalRotation = new Quat4d();
		
		//MT: see wiki "Quaternions and spatial rotation"
		//MT: Sets the value of this quaternion to the equivalent rotation of the AxisAngle argument.
		totalRotation.set(rot);
        
		//MT: encode application of firstRotation, and then totalRotation.
		//MT: see (32) of http://mathworld.wolfram.com/Quaternion.html
		totalRotation.mul(totalRotation, firstRotation);

		Vector3d translation = new Vector3d();
		translation.sub(oldSub.getPositionReal(), bindSub.getPositionReal());
		Vector3d tem = new Vector3d(oldBSPosRelative);

		double scale = 1 + bindBSPosRelative.length()
		/ oldBSPosRelative.length();
		tem.scale(scale);
		translation.add(tem);

        //MT: collection of subunits already attached to the subunit
		//MT: about to bind. They must all be rotated/translated as well.
		Vector<Subunit> bindSubs = bindasm.getSubunits();

		//MT: each of the already attached subunits must rotate/translate
		//MT: relative to the main binding subunit, which is the refPoint.
		Vector3d refPoint = new Vector3d(bindSub.getPositionReal());		
		int size = bindSubs.size();		
		tR_temp = new HashMap<Integer, Quat4d>(size);
		for (int i = 0; i < size; ++i) {

			Subunit cursub = bindSubs.get(i);
			
			//MT: Keep track of subunits rotations for this binding event. 
			//MT: For use in Simulation.java
			Quat4d prev = new Quat4d();			
			prev.set(Simulation.tR.get(cursub.getID()));
			//System.out.println("pre mul");
			//System.out.println(prev.toString());
			
			//System.out.println("totalRotation is:");
			//System.out.println(totalRotation.toString());
			
			prev.mul(totalRotation,prev);
			//System.out.println("post mul");
			//System.out.println(prev.toString());
			tR_temp.put(cursub.getID(),prev);
			//System.out.println("tR_temp is:");
			//System.out.println(tR_temp.get(cursub.getID()));
			
			cursub.subMove(refPoint, translation, totalRotation);
		}


		if (stericHindrance(bindasm)) {
			return null;
		}

		else {


			processLoop(this, bindasm, oldbs, bindBS);
			Vector<Subunit> bindsubs = bindasm.getSubunits();

			//Adds all Subunit(s) in bindasm to this Assembly
			mysubs.addAll(bindsubs);
			size = bindSubs.size();
			for (int i = 0; i < size; ++i) {
				Subunit cursub = bindsubs.get(i);
				cursub.setAssembly(this);
			}
													
			int sizen = this.numSubunits();
			int sizeo = this.numSubunits() - size;
				
				
			if (sizeo <= size)
				Simulation.bif_out.append(sizeo + "\t" + size + "\t" + this.numSubunits() + "\r\n");
			else
				Simulation.bif_out.append(size + "\t" + sizeo + "\t" + this.numSubunits() + "\r\n");				
			
			
			/*if(Test.springForce)
				updateSpringForce();*/

			return this;
		}
	}



	private void processLoop(Assembly a0, Assembly a1, BindingSite bs0, BindingSite bs1) {
		/*boolean operOK=false;

		operOK = */bs0.bindTo(bs1);
		/*if(!operOK)
			System.out.println("Bind operation failed in process loop");*/
		
		Subunit sub1 = bs1.getSubunit();
		Subunit sub0 = bs0.getSubunit();

		/*operOK = */
		//System.out.println("\n" + "processLoop_sub0" + "\n");
		sub1.addPartner(sub0);

		/*if(!operOK)
			System.out.println("addParter operation failed in process loop");*/

		/*operOK = */
		//System.out.println("\n" + "processLoop_sub1" + "\n");
		sub0.addPartner(sub1);

		/*if(!operOK)
			System.out.println("addParter operation failed in process loop");*/

		NeighborGroups ngloop = new NeighborGroups(a0, true);
		Vector<Subunit> subs = a1.getSubunits();
		int size = subs.size();
		int kSize = 0;
		int jSize = 0;
		for (int i = 0; i < size; ++i) 
		{
			Subunit sub = (Subunit) subs.get(i);
			Vector<BindingSite> bss = sub.getBindingSites();
			kSize = bss.size();
			for (int k = 0; k < kSize; ++k) 
			{
				BindingSite tembs1 = (BindingSite) bss.get(k);
				if (!tembs1.isBound()) {
					Vector3d bsTipTem1 = tembs1.getPosReal();

					Vector<BSA> bsNeighbors = ngloop.findLoopCandidates(bsTipTem1);

					if (bsNeighbors.size() == 0) 
						continue;
					jSize = bsNeighbors.size();
					for (int j = 0; j < jSize; ++j) {

						BSA temBSA = bsNeighbors.get(j);
						BindingSite tembs0 = temBSA.bs;

						if (Simulation.isCompatible(tembs1.getBSTName(),tembs0.getBSTName())) {

							Vector3d temv = new Vector3d(tembs0.getPosReal());

							Vector3d v1 = new Vector3d(temv);
							v1.sub(tembs0.getSubunit().getPositionReal());

							Vector3d v2 = new Vector3d(bsTipTem1);
							v2.sub(tembs1.getSubunit().getPositionReal());                        

							temv.sub(bsTipTem1);

							if (temv.length() <= Test.distanceTolerance ) {
								/*operOK = */tembs1.bindTo(tembs0);
								/*if(!operOK)
									System.out.println("bind operation failed in (for loop) of process loop");

								operOK = */
								//System.out.println("\n" + "processLoop_SUBSUBSUB" + "\n");
								(temBSA.sub).addPartner(sub);
								/*if(!operOK)
									System.out.println("addPartner operation failed in (for loop) of process loop");
                                
								operOK = */
								//System.out.println("\n" + "processLoop_temBSA" + "\n");
								sub.addPartner(temBSA.sub);

								/*if(!operOK)
									System.out.println("addPartner operation failed in (for loop) of process loop");*/
							}
						}
					}
				}
			}
		}
	}

	/**
	 * Performs a FastBindEvent.
	 * 
	 * @param e - The FastBindEvent to perform
	 */
	public void fastBind(Event e) 
	{    	
		e.getBS().bindTo(e.getPartner());
		Subunit sub = (e.getBS()).getSubunit();
		Subunit subPartner = (e.getPartner()).getSubunit();
		
		//System.out.println("\n" + "processLoop_fastbnd_subPartner" + "\n");
		sub.addPartner(subPartner);
		
		//System.out.println("\n" + "processLoop_fastbnd_sub" + "\n");
		subPartner.addPartner(sub);
	}

	/**
	 * This gets the counts of unbound BindingSites in this Assembly. This is
	 * used for faster event sampling in the Simulation. 
	 * 
	 * @return HashMap<String, Integer> - A HashMap of the names of each BindingSiteType
	 * 	and the counts of how many are unbound.
	 */
	public HashMap<String, Integer> getBSCounts() {

		HashMap<String, Integer> tem = new HashMap<String, Integer>();
		int size = mysubs.size();
		int jSize = 0;
		for (int i = 0; i < size; ++i) {

			Subunit s = mysubs.get(i);
			Vector<BindingSite> bss = s.getBindingSites();
			jSize = bss.size();
			for (int j = 0; j < jSize; ++j)
			{            	
				BindingSite b = bss.get(j);

				if (!b.isBound()) {
					String bst = b.getBSTName();

					if (tem.containsKey(bst)) {
						Integer n = (Integer) tem.get(bst);
						tem.put(bst, new Integer(n.intValue() + 1));
					} else 
						tem.put(bst, new Integer(1));
				}
			}
		}
		return tem;
	}

	/**
	 * Returns the number of subunits in this Assembly
	 * 
	 * @return int - the number of Subunits in this Assembly
	 */
	public int numSubunits() 
	{ return mysubs.size(); }

	/**
	 * Given a <code>BreakBondEvent</code> this will split the <code>Assembly</code> 
	 * into two
	 * 
	 * @param e - A <code>BreakBondEvent</code> splitting this <code>Assembly</code>
	 * @param newasmnum - Number to use to create unique name for one of the
	 * new <code>Assembly</code>
	 * @return <code>Assembly</code> - The new split <code>Assembly</code>
	 */
	public Assembly splitAssembly(Event e, int newasmnum) 
	{
		
		int sizet = this.numSubunits();
		
		/*boolean operOK=false;*/
		BindingSite bs1 = e.getBS();
		BindingSite bs2 = e.getPartner();
		Subunit bsSub1 = bs1.getSubunit();
		Subunit bsSub2 = bs2.getSubunit();
		bs1.breakBond();

		/*operOK = */bsSub1.removePartner(bsSub2);
		/*if(!operOK)
			System.out.println("RemovePartner failed in Split Assembly! -BUG EXISTS!");

		operOK = */bsSub2.removePartner(bsSub1);

		/*if(!operOK)
			System.out.println("RemovePartner failed in Split Assembly! -BUG EXISTS!");
*/

		if (!(isConnected(bsSub1, bsSub2))) {

			//second assembly
			Vector<Subunit> detachsubunits = getConnected(bsSub2);
			//detach each SubunitGraphic and update BindingSite Assembly names
			Assembly detachedAssembly = new Assembly(detachsubunits,newasmnum);
			int size = detachsubunits.size();
			for (int j = 0; j < size; ++j) {
				Subunit curSubunit = (Subunit) detachsubunits.get(j);
				mysubs.remove(curSubunit);
				curSubunit.setAssembly(detachedAssembly);
			}
				
			int sizeo = sizet - size;
			if (sizeo <= size)
				Simulation.brf_out.append(sizeo + "\t" + size + "\t" + sizet + "\r\n");
			else
				Simulation.brf_out.append(size + "\t" + sizeo + "\t" + sizet + "\r\n");
			
						
			return detachedAssembly;
		}

		//if bond break didn't make 2 assemblies, return null
		return null;
	}

	private boolean isConnected(Subunit s1, Subunit s2)
	{
		HashSet<Integer> seen = new HashSet<Integer>();
		seen.add(s1.getID());
		return getConnected(s1,seen,s2.getID(),false);    	
	}

	private boolean getConnected(Subunit s, HashSet<Integer> seen, int sId, boolean found)
	{
		if(found) return found;

		Vector<Subunit> boundSubunits = s.getBoundSubunits();
		int size = boundSubunits.size();
		for(int i =0; i < size && (!found); ++i)
		{
			Subunit sub = boundSubunits.get(i);
			int id = sub.getID();
			if(id == sId)
			{
				found = true;
				break;
			}

			if(seen.add(id))
				found = getConnected(sub,seen,sId,found);    		
		}

		return found;
	}

	/**
	 * Returns Vector containing all the subunits are in the same Assembly with
	 * Subunit sub, including sub itself. Similar to Dijkstra algorithm
	 * 
	 * @param sub Subunit
	 * @return Vector of Subunits
	 */
	public Vector<Subunit> getConnected(Subunit sub) {
		Vector<Subunit> ret = new Vector<Subunit>();
		HashSet<Integer> seen = new HashSet<Integer>();
		seen.add(sub.getID());//mark subunit
		ret.add(sub);
		return getConnected(ret, sub,seen);
	}    

	/**
	 * Finds all Subunit objects that are in the same Assembly as Subunit sub
	 * marks them(including Subunit sub) and return a Vector of all Subunits
	 * 
	 * @param v  Vector
	 * @param sub Subunit
	 * @return Vector of all Subunit objects in the same Assembly
	 *  
	 */
	private Vector<Subunit> getConnected(Vector<Subunit> v, Subunit sub,HashSet<Integer> seen) {

		Vector<Subunit> boundSubunits = sub.getBoundSubunits();
		int size = boundSubunits.size();
		for(int i =0; i < size; ++i) {
			Subunit curSub = boundSubunits.get(i);
			if (seen.add(curSub.getID())) {
				v.add(curSub);
				getConnected(v, curSub,seen);
			}
		}
		return v;
	}

	/**
	 * Checks for stearic hindrance of this <code>Assembly</code> with a partner.
	 * 
	 * @param partner - The <code>Assembly</code> to look at for stearic hindrance
	 * @return boolean - True if hindrance is detected, false otherwise
	 */
	private boolean stericHindrance(Assembly partner) {

		NeighborGroups ng = null;
		Vector<Subunit> v = null;

		if (this.numSubunits() > partner.numSubunits()) {
			ng = new NeighborGroups(this);
			v = partner.getSubunits();
		} else {
			ng = new NeighborGroups(partner);
			v = this.getSubunits();
		}
		int size = v.size();
		for (int i = 0; i < size; ++i) {

			if (ng.findHindrance(v.get(i))) 
				return true;
		}

		return false;
	}

	public boolean splitAssemblyInLoop(BindingSite bs1, BindingSite bs2) 
	{
		/*boolean operOK=false;//flags added by RP*/
		Subunit bsSub1 = bs1.getSubunit();
		Subunit bsSub2 = bs2.getSubunit();
		boolean breakLoop = false;
		bs1.breakBond();

		/*operOK = */bsSub1.removePartner(bsSub2);
		/*if(!operOK)
			System.out.println("RemovePartner failed in SplitAssemblyInLoop! -BUG EXISTS!");

		operOK = */bsSub2.removePartner(bsSub1);

		/*if(!operOK)
			System.out.println("RemovePartner failed in SplitAssemblyInLoop! -BUG EXISTS!");
*/

		if ((isConnected(bsSub1, bsSub2))) 
			breakLoop = true;


		/*operOK=false;
		operOK=*/bs1.bindTo(bs2);

		/*if(!operOK)
			System.out.println("bindTo failed in SplitAssemblyInLoop! -BUG EXISTS!");

		operOK = */
		//System.out.println("\n" + "processLoop_bsSub2" + "\n");
		bsSub1.addPartner(bsSub2);

		/*if(!operOK)
			System.out.println("addPartner failed in SplitAssemblyInLoop! -BUG EXISTS!");

		operOK = */
		//System.out.println("\n" + "processLoop_bsSub1" + "\n");
		bsSub2.addPartner(bsSub1);

		/*if(!operOK)
			System.out.println("addPartner failed in SplitAssemblyInLoop! -BUG EXISTS!");*/

		return breakLoop;
	}

	/**
	 * Adds a Subunit to this Assembly
	 * 
	 * @param newSubunit - The Subunit to add
	 */
	public void addSubunit(Subunit newSubunit) 
	{  mysubs.add(newSubunit); }



	/**
	 * Returns the subunits in this Assembly
	 * 
	 * @return Vector<Subunit> - A vector of all subunits in this Assembly
	 */
	public Vector<Subunit> getSubunits() 
	{ return mysubs; }



	/**
	 * Returns the valid time for this Assembly
	 * 
	 * @return double - the valid time
	 */
	public double validTime() 
	{  return validTime; }


	/**
	 * Sets the valid time for this Assembly
	 * 
	 * @param time - the valid time to set to
	 */
	public void setValidTime(double time) 
	{ validTime = time; }

	/**
	 * Returns the name of this Assembly
	 * 
	 * @return String - the Name of this Assembly
	 */
	public int getID() 
	{ return myID; }

	/**
	 * Returns a Vector of all free BindingSites in this Assembly
	 * 
	 * @return Vector<BindingSite> - All free BindingSites in this Assembly
	 */
	public Vector<BindingSite> getFreeBindingSites() {

		Vector<BindingSite> free = new Vector<BindingSite>(); 
		int size = mysubs.size();
		for (int i = 0; i < size; ++i) 
			free.addAll(mysubs.get(i).getFreeBindingSites());

		return free;
	}



	/**
	 * A toString method
	 */
	public String toString() {
		StringBuffer buff = new StringBuffer();

		buff.append("Name: ");
		buff.append(myID);
		buff.append(" ValidTime: ");
		buff.append(validTime);
		buff.append("\n");

		for (int i = 0; i < mysubs.size(); i++)
			buff.append(mysubs.get(i) + "\n");

		buff.append("\n");

		return buff.toString();
	}

	/**
	 * Updates the values the spring force model of this Assembly
	 *//*
	public void updateSpringForce() {

		boolean moreBind = true;

		while(moreBind) {

			double currentEnergy = getAssemblyEnergy();
			updateDeriv(currentEnergy);
			double optimalEnergy = findOptimalState(currentEnergy);
			boolean inside = false;

			while ((currentEnergy - optimalEnergy) > 0.000001) {

				if(currentEnergy < optimalEnergy) {
					System.out.println("Something is wrong for getting currentEnergy, optimal");
					System.exit(-1);
				}

				currentEnergy = optimalEnergy;
				updateDeriv(currentEnergy);
				optimalEnergy = findOptimalState(currentEnergy);
			}

			if(inside)
				System.out.println("end of optimal energy = " + optimalEnergy);

			moreBind = findOtherBond();
		}
	}



	private void updateDeriv(double currentEnergy) {

		double derivStep = 0.000000001;

		for(int i = 0; i < mysubs.size(); i++) {

			double energy;

			Subunit sub = mysubs.get(i);

			sub.saveCurrentPos();

			sub.moveSub(new Vector3d(derivStep, 0, 0));
			energy = getAssemblyEnergy();

			sub.setXDeriv((currentEnergy - energy) / derivStep);
			sub.restorePos();

			sub.saveCurrentPos();    		
			sub.moveSub(new Vector3d(0, derivStep, 0));
			energy = getAssemblyEnergy();
			sub.setYDeriv((currentEnergy - energy) / derivStep);
			sub.restorePos();

			sub.saveCurrentPos();    		
			sub.moveSub(new Vector3d(0, 0, derivStep));

			energy = getAssemblyEnergy();
			sub.setZDeriv((currentEnergy - energy) / derivStep);
			sub.restorePos();

			sub.saveCurrentPos();
			sub.rotateSub(new AxisAngle4d(1,0,0,derivStep));
			energy = getAssemblyEnergy();
			sub.setXRotDeriv((currentEnergy - energy)/derivStep);
			sub.restorePos();

			sub.saveCurrentPos();
			sub.rotateSub(new AxisAngle4d(0,1,0,derivStep));
			energy = getAssemblyEnergy();
			sub.setXRotDeriv((currentEnergy - energy)/derivStep);
			sub.restorePos();

			sub.saveCurrentPos();
			sub.rotateSub(new AxisAngle4d(0,0,1,derivStep));
			energy = getAssemblyEnergy();
			sub.setXRotDeriv((currentEnergy - energy)/derivStep);
			sub.restorePos();
		}
	}



	*//**
	 * Given the current energy this will find the optimal energy for the
	 * Assembly
	 * 
	 * @param currentEnergy - the current energy
	 * @return double - the optimal energy
	 *//*
	public double findOptimalState(double currentEnergy) {

		double optimalEnergy = currentEnergy;
		double energy;
		double upperBound = 100000.0;
		double lowerBound = 0.0;
		double kVal = (upperBound + lowerBound) / 2.0;

		while(true) {

			for(int i = 0; i < mysubs.size(); i++) {

				Subunit sub = mysubs.get(i);
				sub.saveCurrentPos();
				sub.applySpringForce(kVal);
			}

			energy = getAssemblyEnergy();

			if(energy < optimalEnergy) {
				optimalEnergy = energy;
				lowerBound = kVal;
			} else 
				upperBound = kVal;

			kVal = (upperBound + lowerBound) / 2.0;

			if((upperBound - lowerBound) < 0.00000001) {
				if(currentEnergy - optimalEnergy < 0.00000001) {

					for(int i = 0; i < mysubs.size(); i++) 
						mysubs.get(i).restorePos();
				}
				return optimalEnergy;
			}

			for(int i = 0; i < mysubs.size(); i++) 
				mysubs.get(i).restorePos();
		}
	}



	*//**
	 * Gets the total energy for this assembly, but summing the energy
	 * of each subunit
	 * 
	 * @return double - the energy of this Assembly
	 *//*
	private double getAssemblyEnergy() {

		double energy = 0.0;

		for(int i = 0; i < mysubs.size(); i++) 
			energy+= mysubs.get(i).getEnergy();

		return energy;
	}



	*//**
	 * 
	 * @return
	 *//*
	private boolean findOtherBond() {

		boolean found = false;

		for(int i = 0; i < mysubs.size()-1; i++) {

			Subunit subOne = mysubs.get(i);

			for(int j = i+1; j < mysubs.size(); j++) {

				Subunit subTwo = mysubs.get(j);

				Vector<BindingSite> bss1 = subOne.getBindingSites();
				Vector<BindingSite> bss2 = subTwo.getBindingSites();

				for(int k = 0; k < bss1.size(); k++) {

					BindingSite bs1 = bss1.get(k);

					if(!bs1.isBound()) {

						for(int l = 0; l < bss2.size(); l++) {

							BindingSite bs2 = bss2.get(l);

							if(!bs2.isBound() && Simulation.isCompatible(bs1.getBSTName(), bs2.getBSTName())) {

								Vector3d tipPos = new Vector3d(bs1.getPosReal());
								tipPos.sub(bs2.getPosReal());

								if(tipPos.length() <= Test.distanceTolerance) {

									subOne.addPartner(subTwo);
									subTwo.addPartner(subOne);
									bs1.bindTo(bs2);
									found = true;
									System.out.println("found other");
								}
							}
						}
					}
				}
			}
		}
		return found;
	}*/


}