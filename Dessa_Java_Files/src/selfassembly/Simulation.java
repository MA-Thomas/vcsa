package selfassembly;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Set;
import java.util.Vector;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;

import com.sun.xml.internal.bind.v2.schemagen.xmlschema.List;


//TODO Document methods and class variables
//TODO clean up code
//TODO change naming of get/sample methods to be consistent
//TODO understand the the update count methods
//TODO make it so size distribution is constantly tracked not updated at intervals

/**
 * This class runs the simulation. It picks events using a stochastic method,
 * stores them in a queue and then executes them. 
 * 
 * @author Tiequan Zhang
 * @author Blake Sweeney
 * @author Marcus Thomas
 * @version 2017
 */
public class Simulation {
	
	/** MT - */
	

	/** An array of the current size distribution in the simulation */
	private int[] currentSizeDistribution=new int[Test.maxLength];
	
	/** An array of the prev size distribution in the simulation */
	private int[] prevSizeDistribution=new int[Test.maxLength];
	
	/** Needed in printDistribution to keep track of vec info of previous time point. */
	Vector<Assembly> prevallAssemblies;

	/** 
	 * Counter for creating new Assemblies. 
	 */
	private int assemblyNumber;


	/** Current simulation step */
	private static int currentStep = 0;

	/** Current simulation time */
	private double curtime;
    
	/** Previous simulation time */
	private double prevtime;
	
	/** Solution simulation occurs in */
	//private Solution mysoln;

	/** HashMap of Assembly name to Assembly */
	private HashMap<Integer, Assembly> assembliesHashMap;
	
	/** MT - HashMap of subunits and their rotations */
	public static HashMap<Integer, Quat4d> tR;	

	/** MT - If step() produced a valid event, set to 1. Else, set to 0. 
	 * For use in deciding whether to update traj_matrix and timeList (and implicitly, reactionList)*/
	public static int stepValid = 0;
	
	/** The priority queue where events are stored */
	private PriorityQueue pq;

	/** 
	 * A map of the BindingTimes. Maps from name to a HashMap 
	 * of partner names and then a double[3], of bond, break, and fast bind 
	 * average times.
	 */

	private static HashMap<String, BindingSiteType> bstMap;

	/**
	 * A HashMap that maps from Conformation names to a map of
	 * potential switches of Name then Double of average switch time.
	 */
	private HashMap<String, HashMap<String, Double>> confTimes; 

	private HashMap<String, Integer> freeBSs;

	private HashMap<String, Integer> monomerBSs;

	
	private StringBuffer sim_out;
	private StringBuffer scatter_out;
	private StringBuffer vec_out;
	public static StringBuffer bif_out;
	public static StringBuffer brf_out;
	
	// MT and HX: Used to save subset of sim_out and vec_out, at time points we care about.
	private StringBuffer sim_out_reduced_times;
	private StringBuffer vec_out_reduced_times;
	private StringBuffer traj_matrix;
	private StringBuffer timeList;
	private StringBuffer reactionList;

	/**
	 * A constructor. This is meant to be used from the XML loader.
	 * 
	 * @param initialAssemblies - Vector of all initial assemblies
	 * @param startTime - the starting time of this simulation
	 * @param bondTimes - Map of Bind/Break/FastBind times
	 * @param confTimes - Conformational Switching Map
	 * @param bindingPartner - Map of potential binding Partners
	 * @param mySolution - Solution simulation occurs in
	 */
	public Simulation(Vector<Assembly> initialAssemblies, double startTime, 
			HashMap<String, BindingSiteType> bindingSiteTypeMap,
			HashMap<String, HashMap<String, Double>> cnfTimes/*,
			Solution mySolution*/) {

		assemblyNumber = 0;

		assemblyNumber = initialAssemblies.size();
		++assemblyNumber;
		//mysoln = mySolution;

		bstMap = bindingSiteTypeMap;
		curtime = startTime;
		confTimes = cnfTimes;

		pq = new BinaryHeap();

		/** 
		 * Initializes the count HashMaps 
		 */


		freeBSs = new HashMap<String, Integer>(bstMap.size());
		monomerBSs = new HashMap<String, Integer>(bstMap.size());
		Iterator<String> it = bstMap.keySet().iterator();

		while (it.hasNext()) {
			String name = it.next();
			freeBSs.put(name, new Integer(0));
			monomerBSs.put(name, new Integer(0));
		}

		/**
		 * Stores the correct values into the count HashMaps
		 */
		int size = initialAssemblies.size();
		assembliesHashMap = new HashMap<Integer, Assembly>(size);
		for (int i = 0; i < size; ++i) 
		{
			Assembly tem = initialAssemblies.get(i);

			if(tem==null)
				System.out.println("assem is null");
			assembliesHashMap.put(tem.getID(), tem);
			HashMap<String, Integer> a = tem.getBSCounts();


			if (tem.numSubunits() == Test.sizeOfSubunit) 
				updateCounts(monomerBSs, a, true);

			updateCounts(freeBSs, a, true);

		}
	
		initialize();
	}



	/**
	 * This method populates the queue with the initial events.
	 * It considers FormBond, BreakBond and if desired ConfChange events.
	 */
	private void initialize() {

		Vector<Assembly> allAssemblies = new Vector<Assembly>(assembliesHashMap.values());
		int size = allAssemblies.size();
		if (size == 0)
			return;

		Event[] formBondArray = new Event[allAssemblies.size()];

		/**
		 * Populate array with the minimum FormBondEvent for each
		 * assembly
		 */
		for (int i = 0; i < size; ++i) {

			Assembly assememblyFirst = allAssemblies.get(i);

			for (int j = 0; j < i; ++j) {

				Assembly assememblySecond = allAssemblies.get(j);

				Event tem = sampleTwoAssems(assememblyFirst,
						assememblySecond);

				double eventTime;

				if (tem == null)
					eventTime = Double.MAX_VALUE;
				else
					eventTime = tem.getEndTime();

				if (formBondArray[i] == null) 
					formBondArray[i] = tem;

				if (formBondArray[j] == null) 
					formBondArray[j] = tem;

				if (formBondArray[i] != null
						&& eventTime < formBondArray[i].getEndTime()) 
					formBondArray[i] = tem;

				if (formBondArray[j] != null
						&& eventTime < formBondArray[j].getEndTime()) 

					formBondArray[j] = tem;
			}
		}

		//breaking event and conf change event
		for (int k = 0; k < size; ++k) {

			Assembly assemi = (Assembly) allAssemblies.get(k);
			Event minFormBondEvent = formBondArray[k];

			Event brkEvt = sampleBreakBondEvent(assemi);
			Event minEvent = screenEvents(minFormBondEvent, brkEvt);

			//cnfChngEvt could be null
			if (Test.csAllowed) {
				Event cnfChngEvt = getConfChangeEvent(assemi);
				minEvent = screenEvents(cnfChngEvt, minEvent);
			}

			formBondArray[k] = minEvent;
		}
		//	remove duplicated formBondEvent
		for (int k = 0; k < size; ++k) {

			/* for virus shell assembly, there will be different bs, consider it
			 * later
			 */

			if (formBondArray[k] == null)
				continue;

			boolean repeat = false;
			for (int m = k + 1; m < size; ++m) {

				if ((formBondArray[k]) == (formBondArray[m])) {
					repeat = true;
					break;
				}
			}

			if ((repeat == false)
					&& (formBondArray[k].getEndTime() < Double.MAX_VALUE)) {				
				sendEvent(formBondArray[k]);
			}            
		}

		for(int i =0; i < size ; ++i)
			allAssemblies.get(i).setValidTime(curtime);
	}



	/**
	 * This method runs the simulation. It steps for as long as there are
	 * valid events. It will output the distribution and if desired an xml
	 * at specified events. 
	 */
	public void run() {

		int eventsPerPrintout = Test.eventsPerPrint;
		int i = 1;		
		
		Set<Double> expTimeRanges = Test.ranges;
		Double curExpTime = (Double)(Collections.min(expTimeRanges));		
		//System.out.println("Minimum Exp Time is : " + curExpTime);
		
		createOutStringBufferObj();
		
		/** MT - During simulation, maintain this hash map. If a subunit is
		 *  involved in a successful binding event, store the vecmath.Quat4d 
		 *  rotation (representing the difference between its current
		 *  orientation and its orientation at the start of the simulation).
		 *  Print the current value for each subunit on every time step with vec_out.
		 */
		int numSubunits = 720; //only 90 for CCMV 1capsid. But over90 okay.		 
		tR = new HashMap<Integer, Quat4d>(numSubunits);  
		for (int z = 1; z < numSubunits+1; z++){
			Quat4d initQuat = new Quat4d(0,0,0,1);
			tR.put(new Integer(z), initQuat);
		}
		//System.out.println(tR.get(1).toString());
		updateCount();
		printDistribution(); ///matsOut is a trivial 1by1 matrix in Matlabless version.

	
		while (!pq.isEmpty()) {

			step();
			i++;			

			if(curtime >= Test.maxSimulationTime)//stop simulation after max time alloted
				break;

			/** MT - expTimeRanges is a set of experimental times each with 2 decimal places. 
			 * Additional loop condition:
			 * If curtime become larger than curExpTime, we know the previous curtime
			 * is the one at which we care about the state of the system.
			 */
			if( i % eventsPerPrintout == 0 ) {// && expTimeRanges.contains(Math.floor(curtime*100)/100) || i == 2) {
				
				//-------
				if ( curExpTime.doubleValue() <= curtime ) {
					if ( expTimeRanges.isEmpty() ) {
						//System.out.println("Nothing left in set");
					} else {
					    Test.casep = 1;
					    expTimeRanges.remove(curExpTime);
					    if (expTimeRanges.isEmpty() ){
						    curExpTime = Double.MAX_VALUE;
					    } else {
					    curExpTime = (Double)(Collections.min(expTimeRanges));					
					    }
					    //System.out.println("Minimum Exp Time is : " + curExpTime);
			        }
		        }
				//-------
				
			    updateCount();
			    printDistribution(); 	
			}		
		}
		
		updateCount();
		printDistribution();	//print last data in case 'if' never satisfied.


		
		//printOutToFile();
		
		// PRINT TO CONSOLE 
		int lastNewLine_sim = sim_out_reduced_times.lastIndexOf("\n");
		int lastNewLine_vec = vec_out_reduced_times.lastIndexOf("\n");
		int lastNewLine_traj = traj_matrix.lastIndexOf("\n");
		int lastNewLine_timeList = timeList.lastIndexOf("\n");
		int lastNewLine_reactionList = reactionList.lastIndexOf("\n");
		sim_out_reduced_times.deleteCharAt(lastNewLine_sim);
		vec_out_reduced_times.deleteCharAt(lastNewLine_vec);
		traj_matrix.deleteCharAt(lastNewLine_traj);
		timeList.deleteCharAt(lastNewLine_timeList);
		reactionList.deleteCharAt(lastNewLine_reactionList);
		System.out.println(sim_out_reduced_times); //see printDistribution()
		System.out.println("----------------");
		System.out.println(vec_out_reduced_times); //see printDistribution()
		//System.out.println("^^^^^^^^^^^^^^^^");
		//System.out.println(traj_matrix); //see printDistribution()
		//System.out.println("################");
		//System.out.println(timeList);	//see printDistribution()	
		//System.out.println("@@@@@@@@@@@@@@@@");
		//System.out.println(reactionList); // see processEvent()
			
	}

	public void createOutStringBufferObj()
	{
		scatter_out = new StringBuffer();
		sim_out = new StringBuffer();
		vec_out = new StringBuffer();
		bif_out = new StringBuffer();
		brf_out = new StringBuffer();
		
		sim_out_reduced_times = new StringBuffer();
		vec_out_reduced_times = new StringBuffer();
		traj_matrix = new StringBuffer();
		timeList = new StringBuffer();
		reactionList = new StringBuffer();
	}
	
	private void printOutToFile()
	{
		String xmlFile = Test.mode.toString().replace(".xml", "");
		String sim_out_filename = "sim_Dessa2017_"+xmlFile;
		String scatter_filename = "scatter_Dessa2017_"+xmlFile;
		String vector_filename = "vector_Dessa2017_"+xmlFile;
		String bif_out_filename = "bif_Dessa2017_"+xmlFile;
		String brf_out_filename = "brf_Dessa2017_"+xmlFile;
		
		
		File simf = new File(sim_out_filename);
		Integer i = new Integer(0);
		while(simf.exists()){
			++i;
			simf = new File(sim_out_filename+i.toString());
		}
		
		File scatf = new File(scatter_filename);
		i = 0;
		while(scatf.exists()){
			++i;
			scatf = new File(scatter_filename+i.toString());
		}
		
		File vecf = new File(vector_filename);
		i = 0;
		while(vecf.exists()){
			++i;
			vecf = new File(vector_filename+i.toString());
		}
				
		File biff = new File(bif_out_filename);
		i = 0;
		while(biff.exists()){
			++i;
			biff = new File(bif_out_filename+i.toString());
		}
		
		File brff = new File(brf_out_filename);
		i = 0;
		while(brff.exists()){
			++i;
			brff = new File(brf_out_filename+i.toString());
		}		
		
		BufferedWriter simbw = null;
		BufferedWriter scatbw = null;
		BufferedWriter vecbw = null;
		BufferedWriter bifbw = null;
		BufferedWriter brfbw = null;		
		
		try {
			simbw = new BufferedWriter(new FileWriter(simf, false));
			scatbw = new BufferedWriter(new FileWriter(scatf, false));
			vecbw = new BufferedWriter(new FileWriter(vecf, false));
			bifbw = new BufferedWriter(new FileWriter(biff, false));
			brfbw = new BufferedWriter(new FileWriter(brff, false));
			
			simbw.write(sim_out_reduced_times.toString());
			simbw.close();
			
			scatbw.write(scatter_out.toString());
			scatbw.close();
			
			bifbw.write(bif_out.toString());
			bifbw.close();
			
			brfbw.write(brf_out.toString());
			brfbw.close();
			
			vecbw.write(vec_out_reduced_times.toString());
			vecbw.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}


	/**
	 * This keeps track of the counts of each assembly size.
	 */
	private void updateCount() 
	{
		for (int size = 0; size < currentSizeDistribution.length; ++size) 
			currentSizeDistribution[size]=0;

		Iterator<Assembly> aItr = assembliesHashMap.values().iterator();
		while(aItr.hasNext())
		{
			Assembly temA = aItr.next();
			int k = temA.numSubunits();
			currentSizeDistribution[k-1]=currentSizeDistribution[k-1]+1;
		}		
	}



	/**
	 * This method performs the actual step. It gets the next event and then
	 * sends it to be processed if it is still valid.
	 */
	public void step() {

		stepValid = 0;
		
		++currentStep;
		Event ev = (Event) pq.remove();

		curtime = ev.getEndTime();

		double[] validTimes = { -1, -1 };

		Assembly[] assemsInvolved = ev.getAssembliesInvolved();
		Assembly asm1 = assemsInvolved[0];

		//only one Assembly involved and it no longer exists
		if (assemsInvolved[1] == null && !assembliesHashMap.containsKey(asm1.getID())) 
			return;

		else if (assemsInvolved[1] != null) { //two Assemblies involved

			Assembly asm2 = (Assembly) assemsInvolved[1];

			//both Assemblies no longer exist
			if (!assembliesHashMap.containsKey(asm1.getID()) && !assembliesHashMap.containsKey(asm2.getID())) 
				return;

			//only asm1 no longer exists
			if (!assembliesHashMap.containsKey(asm1.getID()) && assembliesHashMap.containsKey(asm2.getID())) {

				//sample new events for asm2 before returning
				Assembly[] tmp = { assemsInvolved[1] };
				validTimes[1] = asm2.validTime();

				if (ev.isValid(validTimes)) 
					getNewEvents(tmp);

				return;
			}

			//only asm2 no longer exists
			if (assembliesHashMap.containsKey(asm1.getID()) && !assembliesHashMap.containsKey(asm2.getID())) {

				Assembly[] tmp = { assemsInvolved[0] };
				validTimes[0] = asm1.validTime();

				if (ev.isValid(validTimes)) 
					getNewEvents(tmp);

				return;
			}
		}

		//track validTimes
		validTimes[0] = asm1.validTime();
		if (assemsInvolved[1] != null)
			validTimes[1] = ((Assembly) assemsInvolved[1]).validTime();

		/*
		 * In the case that Assembly(s) involved in one Event exist, If event is
		 * valid, processes it, generates new events for next iteration. else
		 * picks another event for an Assembly that is valid(try again)
		 *  
		 */
		if (ev.isValid(validTimes)) {
			processEvent(ev);
			stepValid = 1;

			//It is tested when only at most two Assembly(s) are involved
			//One Assembly is invalid, consider the other one in case of
			// FormBondEvent
		} else { 			
			
			//Consider only maximum two Assembly(s) in any event
			if (ev.getEventType() == EventType.formBndEvt) {

				if (ev.getPostTime() >= validTimes[0]) {
					Assembly[] tmp = { assemsInvolved[0] };
					getNewEvents(tmp);

				} else if (ev.getPostTime() >= validTimes[1]) {
					Assembly[] tmp = { assemsInvolved[1] };
					getNewEvents(tmp);
					//MT: //System.out.println("NotValidEvent");
				}
			}
		}
	}



	/**
	 * This method process an event. It will perform the event
	 * and pick a new one for the assemblies involved.
	 * 
	 * @param ev - Event to perform
	 */
	private void processEvent(Event ev) {

		if (ev.getEventType() == EventType.cnfChngEvt) {
			System.out.println("confCHANGE");
			Assembly assem = ev.getAssembliesInvolved()[0];

			updateCounts(freeBSs, assem.getBSCounts(), false);
			if (assem.numSubunits() == Test.sizeOfSubunit) {
				updateCounts(monomerBSs, assem.getBSCounts(), false);

			}

			assem.setValidTime(curtime);
			ev.getSubunit().changeConf(ev.getDomainID(), ev.getNewConf());
			ev.setAssembliesInvolved(new Assembly[] { ev.getAssembliesInvolved()[0] });

			updateCounts(freeBSs, assem.getBSCounts(), true);
			if (assem.numSubunits() == Test.sizeOfSubunit) {
				updateCounts(monomerBSs, assem.getBSCounts(), true);
			}

		} else if (ev.getEventType() == EventType.brkBndEvt) {			
			Assembly oldasm = ev.getBS().getAssembly();
			BindingSite bs = ev.getBS();
			BindingSite partner = ev.getPartner();

			simpleUpdateCountsAll(bs, true);
			simpleUpdateCountsAll(partner, true);

            //----------------------------------------------------------------------------------------
			//MT: keep track of {bound at A or B or C (0,1,2, respectively); bondbreak (0)}      
            //----------------------------------------------------------------------------------------
            String bstName_brkBond_partner = partner.getBSTName();
            String bstName_brkBond = bs.getBSTName();
    		// CCMV binding site types ONLY. Order of binding site types must be same as order of elements in boundAt_ABCD.
            String[] freeBindingSitesCCMVtype0 = new String[4];
            freeBindingSitesCCMVtype0[0]="bst0a";
            freeBindingSitesCCMVtype0[1]="bst0b";
            freeBindingSitesCCMVtype0[2]="bst0c";
            freeBindingSitesCCMVtype0[3]="bst0d";
            String[] freeBindingSitesCCMVtype1 = new String[4];
            freeBindingSitesCCMVtype1[0]="bst1b";
            freeBindingSitesCCMVtype1[1]="bst1d";
            freeBindingSitesCCMVtype1[2]="bst1a";
            freeBindingSitesCCMVtype1[3]="bst1c";
            String boundAt = new String();
            
            if(bstName_brkBond.equals(freeBindingSitesCCMVtype0[2])){
            	if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype0[3])){
            		boundAt = "2"; //bound at C  ( bst0c    bst0d )
            	} 
            }
            else if (bstName_brkBond.equals(freeBindingSitesCCMVtype0[3])){
            	if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype0[2])){
            	    boundAt = "2"; //bound at C  ( bst0d    bst0c )
            		}
            	}            
            
            else if (bstName_brkBond.equals(freeBindingSitesCCMVtype0[0])){
            	if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype1[0])){
            		boundAt = "0"; //bound at A  ( bst0a    bst1b )   		
            	}
            	else if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype1[1])){
            		boundAt = "0"; //bound at A  ( bst0a    bst1d )
            	}
            }
            
            else if (bstName_brkBond.equals(freeBindingSitesCCMVtype1[0])){
            	if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype0[0])){
            		boundAt = "0"; //bound at A  ( bst1b    bst0a )
                }
            }
            
            else if (bstName_brkBond.equals(freeBindingSitesCCMVtype1[1])){
            	if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype0[0])){
            		boundAt = "0"; //bound at A  ( bst1d    bst0a )
                }
            }
            
            else if (bstName_brkBond.equals(freeBindingSitesCCMVtype0[1])){
            	if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype1[2])){
            		boundAt = "1"; //bound at B ( bst0b    bst1a )
            	}
            	else if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype1[3])){
            		boundAt = "1"; //bound at B ( bst0b    bst1c )
            	}
            } 
            
            else if (bstName_brkBond.equals(freeBindingSitesCCMVtype1[2])){
            	if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype0[1])){
            		boundAt = "1"; //bound at B ( bst1a    bst0b )      
            	}
            }
            
            else if (bstName_brkBond.equals(freeBindingSitesCCMVtype1[3])){
            	if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype0[1])){
            		boundAt = "1"; //bound at B ( bst1c    bst0b )
            	}            	
            }   
            
            //System.out.println(boundAt);
            //System.out.println(bstName_brkBond + ' ' + bstName_brkBond_partner);
            //System.out.println('\n');
			reactionList.append(boundAt+' '+0);
            reactionList.append('\n');
            //MT: --------------------------------------------------------------------------
            //    --------------------------------------------------------------------------
            		
			//MT: In general, splitAssembly will produce a new Assembly            
			Assembly newasm = oldasm.splitAssembly(ev, assemblyNumber);

			//MT: if the splitAssembly makes a new assembly
			if (newasm != null) 
			{
				assembliesHashMap.put(newasm.getID(), newasm);

				ev.setAssembliesInvolved(new Assembly[] { oldasm, newasm });
				++assemblyNumber;


				if (newasm.numSubunits() == Test.sizeOfSubunit) 
					updateCounts(monomerBSs, newasm.getBSCounts(), true);

				if (oldasm.numSubunits() == Test.sizeOfSubunit) 
					updateCounts(monomerBSs, oldasm.getBSCounts(), true);

				newasm.setValidTime(curtime);
				oldasm.setValidTime(curtime);

			} else {

				double fastBindingInterval = 0.0;
				oldasm.setValidTime(curtime);
				Assembly[] sameAssembly = { ev.getAssembliesInvolved()[0],
						ev.getAssembliesInvolved()[0] };
				Event fastBindingEvent = new Event(curtime,
						curtime + fastBindingInterval, sameAssembly, ev.getBS(),
						ev.getPartner(),EventType.formBndEvt);
				sendEvent(fastBindingEvent);

				//MT: return to prevent picking of new events for the assembly
				return;
			}
		} else if (ev.getEventType() == EventType.formBndEvt) 
		{
			BindingSite bs = ev.getBS();
			Assembly asm1 = bs.getAssembly();

			BindingSite partner = ev.getPartner();
			Assembly asm2 = partner.getAssembly();

			//    ----------------------------------------------------------------------------------------
			//MT: keep track of {bound at A or B or C (0,1,2, respectively); bondform (1)} 
			//    ----------------------------------------------------------------------------------------
            String bstName_brkBond_partner = partner.getBSTName();
            String bstName_brkBond = bs.getBSTName();
    		// CCMV binding site types ONLY. Order of binding site types must be same as order of elements in boundAt_ABCD.
            String[] freeBindingSitesCCMVtype0 = new String[4];
            freeBindingSitesCCMVtype0[0]="bst0a";
            freeBindingSitesCCMVtype0[1]="bst0b";
            freeBindingSitesCCMVtype0[2]="bst0c";
            freeBindingSitesCCMVtype0[3]="bst0d";
            String[] freeBindingSitesCCMVtype1 = new String[4];
            freeBindingSitesCCMVtype1[0]="bst1b";
            freeBindingSitesCCMVtype1[1]="bst1d";
            freeBindingSitesCCMVtype1[2]="bst1a";
            freeBindingSitesCCMVtype1[3]="bst1c";
            String boundAt = new String();
            
            if(bstName_brkBond.equals(freeBindingSitesCCMVtype0[2])){
            	if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype0[3])){
            		boundAt = "2"; //bound at C  ( bst0c    bst0d )
            	} 
            }
            else if (bstName_brkBond.equals(freeBindingSitesCCMVtype0[3])){
            	if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype0[2])){
            	    boundAt = "2"; //bound at C  ( bst0d    bst0c )
            		}
            	}            
            
            else if (bstName_brkBond.equals(freeBindingSitesCCMVtype0[0])){
            	if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype1[0])){
            		boundAt = "0"; //bound at A  ( bst0a    bst1b )   		
            	}
            	else if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype1[1])){
            		boundAt = "0"; //bound at A  ( bst0a    bst1d )
            	}
            }
            
            else if (bstName_brkBond.equals(freeBindingSitesCCMVtype1[0])){
            	if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype0[0])){
            		boundAt = "0"; //bound at A  ( bst1b    bst0a )
                }
            }
            
            else if (bstName_brkBond.equals(freeBindingSitesCCMVtype1[1])){
            	if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype0[0])){
            		boundAt = "0"; //bound at A  ( bst1d    bst0a )
                }
            }
            
            else if (bstName_brkBond.equals(freeBindingSitesCCMVtype0[1])){
            	if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype1[2])){
            		boundAt = "1"; //bound at B ( bst0b    bst1a )
            	}
            	else if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype1[3])){
            		boundAt = "1"; //bound at B ( bst0b    bst1c )
            	}
            } 
            
            else if (bstName_brkBond.equals(freeBindingSitesCCMVtype1[2])){
            	if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype0[1])){
            		boundAt = "1"; //bound at B ( bst1a    bst0b )      
            	}
            }
            
            else if (bstName_brkBond.equals(freeBindingSitesCCMVtype1[3])){
            	if (bstName_brkBond_partner.equals(freeBindingSitesCCMVtype0[1])){
            		boundAt = "1"; //bound at B ( bst1c    bst0b )
            	}            	
            }       
            //System.out.println(boundAt);
            //System.out.println(bstName_brkBond + ' ' + bstName_brkBond_partner);
            //System.out.println('\n');
			reactionList.append(boundAt+' '+1);
            reactionList.append('\n');
            //MT: --------------------------------------------------------------------------
            //    --------------------------------------------------------------------------
            
			if (asm1 == asm2) {

				asm1.setValidTime(curtime);
				asm1.fastBind(ev);
				ev.setAssembliesInvolved(new Assembly[] { asm1 });

				simpleUpdateCountsAll(bs, false);
				simpleUpdateCountsAll(partner, false);


			} else if ((asm1.numSubunits() + asm2.numSubunits()) <= Test.maxLength) {                
				updateCounts(freeBSs, asm1.getBSCounts(), false);
				updateCounts(freeBSs, asm2.getBSCounts(), false);
				if (asm1.numSubunits() == Test.sizeOfSubunit)
					updateCounts(monomerBSs, asm1.getBSCounts(), false);

				if (asm2.numSubunits() == Test.sizeOfSubunit)
					updateCounts(monomerBSs, asm2.getBSCounts(), false);

				//MT: subunits in asm2 are rotated/translated.
				Assembly tmp = asm1.bindAssembly(asm2, ev);

				if (tmp != null) { //MT: binding event successful
					
					//MT: update HashMap tR by iterating through tR_temp created 
					//    when bingAssembly was called. The same 'keys'(subunit IDs) are  
					//    applicable to both tR_temp and tR. 
					for (Integer key : Assembly.tR_temp.keySet()){
						//System.out.println("tR_temp in Simulation.java is:");
						//System.out.println(Assembly.tR_temp.get(key));
						tR.put(key,Assembly.tR_temp.get(key));
					}
					
					assembliesHashMap.remove(asm2.getID());
					ev.setAssembliesInvolved(new Assembly[] { asm1 });
					asm1.setValidTime(curtime);

					updateCounts(freeBSs, asm1.getBSCounts(), true);

				} else { //MT: binding event unsuccessful
					
					asm1.setValidTime(curtime);
					asm2.setValidTime(curtime);
					updateCounts(freeBSs, asm1.getBSCounts(), true);
					updateCounts(freeBSs, asm2.getBSCounts(), true);
					if (asm1.numSubunits() == Test.sizeOfSubunit) {
						updateCounts(monomerBSs, asm1.getBSCounts(), true);
					}
					if (asm2.numSubunits() == Test.sizeOfSubunit) {
						updateCounts(monomerBSs, asm2.getBSCounts(), true);
					}

				}
			} else { //MT: if the assembly is too big
				asm1.setValidTime(curtime);
				asm2.setValidTime(curtime);
			}
		}

		/*
		 * generates new events for queue based on new current state, only need
		 * to update nextEvent for assemblies that were involved in "ev".if only
		 * one assembly is involved, second element is null.
		 */
		getNewEvents(ev.getAssembliesInvolved());
	}



	/**
	 * Returns the Coefficient 
	 * 
	 * @param primary
	 * @param secondary
	 * @return double - the Coefficient
	 */
	double getCoef(int primary, int secondary) {

		return Math.sqrt((primary + secondary)
				/ (2.0 * primary * secondary));
	}

	public final static Vector3d getBSTPostion(String type)
	{
		return new Vector3d(bstMap.get(type).getBSTPosition());
	}

	public final static double getBSTBindAngle(String type,String partnerType)
	{
		return bstMap.get(type).getBindingAngle(partnerType);
	}

	/**
	 * Picks new events for the two Assemblies. Puts the new Events into
	 * the queue. Events are picked that are either a FormBondEvent, 
	 * BreakBondEvent or if allowed a ConfChangeEvent. 
	 * 
	 * @param assems - Array of the two Assemblies to pick new events for
	 */
	private void getNewEvents(Assembly[] assems) {

		HashMap otherFreeBSs = null;
		HashMap<String, Integer> otherMonomerBSs = null;


		otherFreeBSs = new HashMap<String, Integer>(freeBSs);
		otherMonomerBSs = new HashMap<String, Integer>(monomerBSs);

		HashMap selfBSs0 = null;
		HashMap selfBSs1 = null;
		Event myMinEvent0 = null;
		Event myMinEvent1 = null;

		//Some error checking

		if (assems.length == 1) {

			if (assems[0] == null
					|| (this.assembliesHashMap.get(assems[0].getID()) == null)) {
				System.out
				.println("argument_error_1 in Simulation getNewEvents");
				System.exit(1);
			}

		} else if (assems.length == 2) {
			if ((assems[0] == null)
					|| (this.assembliesHashMap.get(assems[0].getID()) == null)
					|| (assems[1] == null)
					|| (this.assembliesHashMap.get(assems[1].getID()) == null)) {
				System.out
				.println("argument_error_2 in Simulation getNewEvents");
				System.exit(1);
			}
			if (assems[0].equals(assems[1])) {
				System.out
				.println("argument_error_3 in Simulation getNewEvents");
				System.exit(1);
			}
		} else {
			System.out.println("argument_error_4 in Simulation getNewEvents");
			System.exit(1);
		}

		if (assems.length == 1) {

			Assembly a0 = (Assembly) assems[0];
			selfBSs0 = a0.getBSCounts();

			updateCounts(otherFreeBSs, selfBSs0, false);
			a0.setValidTime(curtime);
			myMinEvent0 = getFormBondEvent(a0, otherFreeBSs,
					otherMonomerBSs, null);

			myMinEvent0 = screenEvents(myMinEvent0, getConfChangeEvent(a0));
			Event brkEvt = sampleBreakBondEvent(a0);
			myMinEvent0 = screenEvents(myMinEvent0, brkEvt);


			if (a0.numSubunits() == Test.sizeOfSubunit) {
				if (brkEvt != null)
					System.out.println("error1insimulation");
			}


		} else {

			Assembly a0 = (Assembly) assems[0];
			Assembly a1 = (Assembly) assems[1];
			selfBSs0 = a0.getBSCounts();
			selfBSs1 = a1.getBSCounts();

			updateCounts(otherFreeBSs, selfBSs0, false);
			updateCounts(otherFreeBSs, selfBSs1, false);

			if (a0.numSubunits() == Test.sizeOfSubunit) {
				updateCounts(otherMonomerBSs, selfBSs0, false);
			}
			if (a1.numSubunits() == Test.sizeOfSubunit) {
				updateCounts(otherMonomerBSs, selfBSs1, false);
			}

			a0.setValidTime(curtime);
			a1.setValidTime(curtime);
			Event f0 = null;
			Event f1 = null;

			f0 = getFormBondEvent(a0, otherFreeBSs, otherMonomerBSs, a1);
			f1 = getFormBondEvent(a1, otherFreeBSs, otherMonomerBSs, a0);

			Event c = sampleTwoAssems(a0, a1);
			myMinEvent0 = screenEvents(f0, c);
			myMinEvent0 = screenEvents(myMinEvent0, getConfChangeEvent(a0));
			Event tbe0 = sampleBreakBondEvent(a0);
			myMinEvent0 = screenEvents(myMinEvent0, tbe0);

			if (a0.numSubunits() == Test.sizeOfSubunit) {
				if (tbe0 != null)
					System.out.println(a0.numSubunits() + " " + currentStep
							+ " errorin2simulation");
			}

			myMinEvent1 = screenEvents(f1, c);
			myMinEvent1 = screenEvents(myMinEvent1, getConfChangeEvent(a1));
			Event tbe1 = sampleBreakBondEvent(a1);
			myMinEvent1 = screenEvents(myMinEvent1, tbe1);

			if (a1.numSubunits() == Test.sizeOfSubunit) {
				if (tbe1 != null)
					System.out.println("errorin3simulation");
			}

		}

		if ((myMinEvent0 != null) && (myMinEvent1 == null)) {
			sendEvent(myMinEvent0);
		} else if ((myMinEvent0 == null) && (myMinEvent1 != null)) {
			sendEvent(myMinEvent1);
		} else if ((myMinEvent0 != null) && (myMinEvent1 != null)) {
			if (myMinEvent0 != myMinEvent1) {
				sendEvent(myMinEvent0);
				sendEvent(myMinEvent1);
			} else {
				sendEvent(myMinEvent0);
			}
		} else {
		}
	}



	/**
	 * Using the two assemblies given to picks a new FormBondEvent. This
	 * will be the minimum of the possible events for these Assemblies.
	 * 
	 * @param assem - The fist Assembly to sample
	 * @param otherFreeBSs - A count of the free BindingSites 
	 * 		for the first Assembly
	 * @param otherMonomerBSs - A count of free BindingSites 
	 * 		for the second Assembly
	 * @param assem2 - The second Assembly to sample
	 * @return FormBondEvent - the new minimum bonding event.
	 */
	private Event getFormBondEvent(Assembly assem,
			HashMap otherFreeBSs, HashMap otherMonomerBSs, Assembly assem2) {

		double duration = Double.MAX_VALUE;
		double tem;
		String bstSelfChoice = null;
		String bstOtherChoice = null;
		int countSelf = -1;
		int countOther = -1;

		HashMap<String, Integer> abs = assem.getBSCounts();
		HashMap<String, Integer> other = null;

		boolean monomerOnly = false;

		if ((assem.numSubunits() == Test.sizeOfSubunit) || (!Test.bindMonomerOnly)) 
			other = otherFreeBSs;
		else {
			//current Assembly has more than one monomers and
			// Test.bindMonomerOnly=true
			monomerOnly = true;
			other = otherMonomerBSs;
		}

		Iterator<String> it = abs.keySet().iterator();
		while (it.hasNext()) {

			String bst = it.next();
			BindingSiteType bsType = bstMap.get(bst);
			Iterator<String> partner = bsType.getPartners().iterator();
			while(partner.hasNext())
			{

				String bstP = partner.next();

				int a = abs.get(bst).intValue();
				int b = other.get(bstP).intValue();
				int product = a * b;

				if (product != 0) {

					double bindTime = bsType.getBindTime(bstP);
					tem = getExp(bindTime/product);

					if (tem < duration) {
						duration = tem;
						countSelf = a;
						countOther = b;
						bstSelfChoice = bst;
						bstOtherChoice = bstP;
					}
				}
			}
		}

		Event minFormBondEvent = null;

		if (duration < Double.MAX_VALUE) {
			minFormBondEvent = sampleFormBondEvent(duration, assem, assem2,
					bstSelfChoice, bstOtherChoice,
					(int) (countSelf * getRand()),
					(int) (countOther * getRand()), monomerOnly);
		}
		return minFormBondEvent;
	}

	/**
	 * A Simple method to screen to <code>Event</code>s. Will return the
	 * second if the first is null, first if second is null, or the minimum
	 * of the two if both are not null
	 * 
	 * @param e0 - The first <code>Event</code>
	 * @param e1 - The second <code>Event</code>
	 * @return <code>Event</code> - Whichever <code>Event</code> 
	 * 	is not null and occurs first. 
	 */
	private Event screenEvents(Event e0, Event e1) {

		if (e0 == null) 
			return e1;
		else if (e1 == null) 
			return e0;
		else if (e0.compareTo(e1) == 1) 
			return e0;
		else
			return e1;
	}



	/**
	 * Samples all possible <code>BreakBondEvents</code> for a given
	 * <code>Assembly</code>
	 * 
	 * @param assem - The <code>Assembly</code> to sample
	 * @return <code>BreakBondEvents</code> if an event is found, <code>null</code>
	 * 	otherwise
	 */
	private Event sampleBreakBondEvent(Assembly assem) {

		Vector<Subunit> subs = assem.getSubunits();

		if (subs.size() == Test.sizeOfSubunit) 
			return null;

		double duration = 0;
		double minTime = Double.MAX_VALUE;

		BindingSite temBs, temPartner;
		BindingSite minBs = null;
		BindingSite minPartner = null;
		CheckPairSeen cps = new CheckPairSeen();

		//go thru boundsubs, looking at all possible bond breaking events
		int size = subs.size();
		for (int j = 0; j < size; ++j) {

			Subunit sub = subs.get(j);
			Vector<BindingSite> bindsites = sub.getBindingSites();

			if ((sub.getBoundSubunits().size() > 1) && (Test.breakOnlyEnds)) {
				continue;
			}
			int bsSize = bindsites.size();
			for (int k = 0; k < bsSize; ++k) {

				temBs = bindsites.get(k);

				if (temBs.isBound()) {

					temPartner = temBs.getPartner();

					if (cps.marked(temBs.getID(), temPartner.getID()))
						continue;

					//else mark it seen now
					cps.addPair(temBs.getID(), temPartner.getID());

					if (Test.noLoopOnly) {
						boolean breakLoop = assem.splitAssemblyInLoop(temBs,
								temPartner);
						if (breakLoop)
							continue;
					}

					double brkTime =bstMap.get(temBs.getBSTName()).getBreakTime(temPartner.getBSTName());

					duration = getExp(brkTime);
					//if this is the fastest Event so far, store it
					if (minTime >= duration) {
						minTime = duration;
						minBs = temBs;
						minPartner = temPartner;
					}
				}
			}
		}


		if (minTime == Double.MAX_VALUE)
			return null;

		Assembly[] assems = { minBs.getAssembly(), null };

		return new Event(curtime, minTime + curtime, assems, minBs,
				minPartner,EventType.brkBndEvt);
	}



	private Event sampleFormBondEvent(double duration, Assembly assem,
			Assembly assem2, String bstSelf, String bstOther, int indexSelf,
			int indexOther, boolean monomerOnly) {

		/**
		 * 
		 * Walks through all Assembly(s) in this Simulation, finds and saves all
		 * the other free BindingSite(s) by checking if each one is bound and is
		 * in Assembly(s) specified by assems
		 */
		BindingSite minBs = null;
		BindingSite minPartner = null;
		Vector<Subunit> subs = assem.getSubunits();
		boolean indexSelfFound = false;
		int temIndex = 0;
		int aSize = subs.size();
		for (int i = 0; ((i < aSize) && (!indexSelfFound)); ++i) {

			Subunit sub = subs.get(i);
			Vector<BindingSite> bss = sub.getBindingSites();
			int bSize = bss.size();
			for (int j = 0; j < bSize && !indexSelfFound; ++j) {

				BindingSite bs = bss.get(j);

				if ((!bs.isBound()) && (bs.getBSTName().equals(bstSelf))) {

					if (temIndex == indexSelf) {
						indexSelfFound = true;
						minBs = bs;
					} else {
						temIndex++;
					}
				}
			}
		}

		//find partner BindingSite of type bstOther from all the other
		// Assembly(s) or all the other monomers if argument monomerOnly is true
		boolean indexOtherFound = false;
		temIndex = 0;

		Iterator<Assembly> it = assembliesHashMap.values().iterator();
		Assembly temAssembly;

		while (it.hasNext() && !indexOtherFound) {

			temAssembly = it.next();

			if (temAssembly == assem || temAssembly == assem2)
				continue;
			if (temAssembly.numSubunits() > Test.sizeOfSubunit && monomerOnly) 
				continue;

			Vector<Subunit> temSubs = temAssembly.getSubunits();
			int cSize = temSubs.size();
			for (int j = 0; ((j < cSize) && (!indexOtherFound)); ++j) {

				Subunit temSub = (Subunit) temSubs.get(j);
				Vector<BindingSite> temBss = temSub.getBindingSites();
				int dSize = temBss.size();
				for (int k = 0; ((k < dSize) && (!indexOtherFound)); ++k) 
				{
					BindingSite temBs = temBss.get(k);

					if ((!temBs.isBound()) && (temBs.getBSTName().equals(bstOther))) {

						if (temIndex == indexOther) {
							indexOtherFound = true;
							minPartner = temBs;							
						} else 
							temIndex++;
					}
				}
			}
		}

		Assembly[] assems = { minBs.getAssembly(), 
				minPartner.getAssembly() };

		return new Event(curtime, curtime + duration, assems, minBs,
				minPartner,EventType.formBndEvt);
	}

	public final static boolean isCompatible(String bst1, String bst2)
	{
		return bstMap.get(bst1).isCompatible(bst2);
	}


	private Event sampleTwoAssems(Assembly A1, Assembly A2) {
		//Assume A1 and A2 should be different Assembly(s)
		if (A1.equals(A2)) {
			return null;
		}

		if ((A1.numSubunits() > Test.sizeOfSubunit)
				&& (A2.numSubunits() > Test.sizeOfSubunit)
				&& (Test.bindMonomerOnly))
			return null;

		Vector<BindingSite> bss1 = A1.getFreeBindingSites();
		Vector<BindingSite> bss2 = A2.getFreeBindingSites();

		int primary = -1;
		int partner = -1;

		double minTime = Double.MAX_VALUE;
		double duration = 0;
		int aSize = bss1.size();
		for (int i = 0; i < aSize; ++i) {
			BindingSite temBs = bss1.get(i);
			int bSize = bss2.size();
			for (int j = 0; j < bSize; ++j) {
				BindingSite temPartner = bss2.get(j);
				//check for compatible bindingsiteTypes

				if (!isCompatible(temBs.getBSTName(), temPartner.getBSTName())) {
					continue;
				}


				double bTime = bstMap.get(temBs.getBSTName()).getBindTime(temPartner.getBSTName());


				duration = getExp(bTime);

				//if this is the fastest Event so far, store it
				if (minTime >= duration) {
					minTime = duration;
					primary = i;
					partner = j;
				}
			}
		}
		//To be removed for new initialize()
		Assembly[] assemsInvolved = { A1, A2 };

		if (minTime != Double.MAX_VALUE) {
			return new Event(curtime, curtime + minTime,
					assemsInvolved, (BindingSite) bss1.get(primary),
					(BindingSite) bss2.get(partner),EventType.formBndEvt);
		} else {
			return null;
		}

	}



	/**
	 * This will pick a <code>ConfChangeEvent</code> for the given
	 * <code>Assembly</code>. Assembly should be a monomer, as defined 
	 * by Test.sizeofSubunit.
	 * 
	 * @param assem_i - The <code>Assembly</code> to pick a <code>Event</code>
	 * 	for.
	 * @return <code>ConfChangeEvent</code> if possible, null if none
	 * 	found.
	 */
	private Event getConfChangeEvent(Assembly assem_i) {

		if (!Test.csAllowed || assem_i.numSubunits() != Test.sizeOfSubunit)
			return null;

		Vector<Subunit> subs = assem_i.getSubunits();
		Subunit sub = subs.get(0);

		Vector<Domain> domains = sub.getDomains();
		double minTime = Double.MAX_VALUE;
		int dID = -1;
		Conformation conf = null;

		/**
		 * Go through all domains in the Subunit and look at all possible
		 * conformations
		 */
		int aSize = domains.size();
		for (int i = 0; i < aSize; ++i) {

			Domain curDomain = domains.get(i);
			int id = curDomain.getDomainId();
			Conformation curConf = curDomain.getCurConf();
			HashMap<String,Conformation> allConfs = curDomain.getConfs();
			HashMap<String, Double> switchMap = confTimes.get(curConf.getName());

			/**
			 * For each possible conformation of each possible domain use the
			 * switch map to find the time for the event, storing the smallest
			 */
			Iterator<Conformation> cnfItr = allConfs.values().iterator();
			while(cnfItr.hasNext())
			{
				Conformation testConf = cnfItr.next();

				if (switchMap.containsKey(testConf.getName())) {

					double changeTime = switchMap.get(testConf.getName()).doubleValue();	
					double duration = getExp(changeTime);

					if (duration < minTime) {
						minTime = duration;
						dID = id;
						conf = testConf;
					}
				}
			}
		}

		if (minTime == Double.MAX_VALUE)
			return null;

		return new Event(curtime, curtime + minTime,
				new Assembly[] { assem_i, null}, sub,
				dID, conf);
	}




	private void updateCounts(HashMap<String, Integer> old, HashMap<String, Integer> tem, boolean add) {

		Iterator<String> it = tem.keySet().iterator();

		while (it.hasNext()) {

			String name = it.next();

			int temN = tem.get(name).intValue();
			int oldN = old.get(name).intValue();

			if (add) 
				old.put(name, new Integer(temN + oldN));
			else 
				old.put(name, new Integer(oldN - temN));
		}
	}



	private void simpleUpdateCountsAll(BindingSite b, boolean add) {

		String name = b.getBSTName();
		int n = freeBSs.get(name).intValue();

		if (add)
			freeBSs.put(name, new Integer(n + 1));
		else
			freeBSs.put(name, new Integer(n - 1));
	}



	/**
	 * Returns an random exponentially distributed variable with the given
	 * average value.
	 * 
	 * @param avg - The average to base this around
	 * @return double - a random exponential
	 */
	private double getExp(double avg) {

		double t = -avg * Math.log(getRand());
		return t;
	}



	/**
	 * Returns a random double within [0 , 1)
	 * 
	 * @return double
	 */
	private double getRand() 
	{ return Test.rand.nextDouble(); }





	/**
	 * Adds an <code>Event</code> to the queue
	 * 
	 * @param e - The <code>Event</code> to add
	 */
	private final void sendEvent(Event e) 
	{ pq.add(e); }


	/**
	 * This method writes the size distribution. It will write to 
	 * file if <code>Test.printToScreen</code> = false, otherwise it will print
	 * to screen. The output is formatted as follows:
	 * 		Simulation_Time Count_of_Assemblies_of_Size_1 ... Count_of_Assemblies_of_Size_N
	 * where N is defined by <code>Test.maxOutputSize</code>
	 */
    public void printDistribution() {   
    	

		Vector<Assembly> allAssemblies = new Vector<Assembly>(assembliesHashMap.values());	
		
        // BEGIN trajectory likelihood information SECTION ---------------------------------------------------------------------------------------		
        // If step() produced a valid event, then we know a bond break or bond formation happened (assuming no conf changes).
		// This means the system state has changed and we should append the likelihood info (to timeList, traj_matrix)
		// ---------------------------------------------------------------------------------------------------------------------------------------
	    if (stepValid == 1){
		prevallAssemblies = allAssemblies; 
	    timeList.append(curtime);
	    timeList.append("\n");

	    // Get event type and binding sites involved for the most recent reaction/timestep.
        // SEE processEvent() for 'reactionList' variable.
        
		//MT: Oct23 2016. Append a one-hot (/many-hot) vector indicating which of the following 
		//    is(/are) the case for the current (ccmv) subunit:
		//    monomerType0 monomerType1 boundAtA boundAtB boundAtC boundAtD
		//    Note, if monomerType0 or monomerType1 are the case, each of boundAtA, boundAtB, boundAtC, boundAtD, MUST be 0.
		//    Eventually (in Matlab) this info will be used to calculate the  
		//    likelihood of this trajectory.
		int monomerType0 = 0;
		int monomerType1 = 0;
		int[] boundAt_ABCD = new int[4];
		boundAt_ABCD[0]=0;
		boundAt_ABCD[1]=0;
		boundAt_ABCD[2]=0;
		boundAt_ABCD[3]=0;

		
		// CCMV binding site types ONLY. Order of binding site types must be same as order of elements in boundAt_ABCD.
        String[] freeBindingSitesCCMVtype0 = new String[4];
        freeBindingSitesCCMVtype0[0]="bst0a";
        freeBindingSitesCCMVtype0[1]="bst0b";
        freeBindingSitesCCMVtype0[2]="bst0c";
        freeBindingSitesCCMVtype0[3]="bst0d";
        String[] freeBindingSitesCCMVtype1 = new String[4];
        freeBindingSitesCCMVtype1[0]="bst1b";
        freeBindingSitesCCMVtype1[1]="bst1d";
        freeBindingSitesCCMVtype1[2]="bst1a";
        freeBindingSitesCCMVtype1[3]="bst1c";	    
        
        
        // Iterate through all assemblies at current timepoint.
		int sizep = allAssemblies.size();			
		for (int j = 0; j < sizep; ++j) {
			Assembly assemblyFirst = allAssemblies.get(j);
			Vector<Subunit> allvec = assemblyFirst.getSubunits();
			//System.out.println(allvec.toString());
			
			int k = allvec.size();			
						
	
			// Iterate through all subunits in current assembly.
			for (int l = 0; l < k; ++l) {
				Subunit nxt = allvec.elementAt(l);
				
				//    Append subunit type (for ccmv, there is "only0" and "only1")
				//    Treat only0 as 0, only1 as 1.
				String st = nxt.getSubunitTypeName();
				st = st.replace("only", "");
				
				int numFreeBS = nxt.getFreeBindingSites().size();
				
				// If subunit is monomer.
				if (nxt.isUnbound()){

					if (st.equals("0")){

						monomerType0++;
					} else {

						monomerType1++;
					}					
				
				} 
				
				// Else, if subunit not completely bound.
				else if (numFreeBS > 0) {	
					

                    String[] bstNamesString = new String[3];                    
					for (int nf = 0; nf < numFreeBS; ++nf){
						bstNamesString[nf] = nxt.getFreeBindingSites().get(nf).getBSTName();						
					}
					java.util.List<String> bstNames = Arrays.asList(bstNamesString);
					
					// Check to see which binding site type in bstNames is missing from freeBindingSitesCCMVtype0
					// ( or from freeBindingSitesCCMVtype1 ). ccmv ONLY
					if (st.equals("0")){
						
						for (int i = 0; i < 4; ++i){
						    if (!bstNames.contains(freeBindingSitesCCMVtype0[i])){	
						    	if (i==3){
						    		boundAt_ABCD[i-1]++; // CCMV only has A,B,C binding site types. No boundAtD.
						    	}else {
						    		boundAt_ABCD[i]++; // Increment boundAtA, (boundAtB,boundAtC, respectively)
						    	}
					    	}
						}
						
					}
					else if (st.equals("1")){
						
						for (int i = 0; i < 4; ++i){
						    if (!bstNames.contains(freeBindingSitesCCMVtype1[i])){	
						    	if (i==0 || i==1){
						    		boundAt_ABCD[0]++; // Increment boundAtA
						    	}else if (i==2 || i==3) {
						    		boundAt_ABCD[1]++; // Increment boundAtB
						    	}
					    	}
						}
					}
				} 
				
				// Else, the subunit is completely bound.	
				else {
					boundAt_ABCD[0]++;
					boundAt_ABCD[1]++;
					boundAt_ABCD[2]++;
					//boundAt_ABCD[3]++; //SINCE CCMV HAS NO 'D' BINDING SITES, IT SHOULD ALWAYS BE 0.										
				}				
	
			}

		}
		
		traj_matrix.append(monomerType0+" "+monomerType1+" "+boundAt_ABCD[0]+" "+boundAt_ABCD[1]+" "+boundAt_ABCD[2]+" "+boundAt_ABCD[3]);
		traj_matrix.append("\n");
        }
		// END trajectory likelihood information section -----------------------------------------------------------------------------------
		
		
		/* June 29, 2016. Marcus T and Huiming Xia modification.
		   If prevtime was closest simtime -lt the relevant experimental time 
		   (i.e., if Test.casep=1),
		    append previous sim line and vec line (they'll be printed in run() ).
		*/
		if (Test.casep == 1){
			
			//This block should not exist here. We should be strictly using the prev time and prev state of the system for sim and vec.
			//However, that leads to anomalous results for vec (e.g. at single time point, subunit appearing in more than 1 assembly).
			//So technically this is not the system state we care about, but if the interval (eventsperprint) is 1, the true system state
			//is different by only a single reaction event. With an interval -gt 1, prevtime and the prev state of the
			//system are necessarily not the exact time points we care about anyhow. 
		    prevSizeDistribution = currentSizeDistribution; 
		    prevallAssemblies = allAssemblies; 
		    prevtime = curtime;
		    
		 // SIM SECTION.
			sim_out_reduced_times.append(prevtime + " ");
			for (int i = 0; i < prevSizeDistribution.length && 
					i < Test.maxOutputSize; ++i)
					{						
						sim_out_reduced_times.append(prevSizeDistribution[i] + " ");
					}
			sim_out_reduced_times.append("\n");
			
			// VEC SECTION.
			int sizep = prevallAssemblies.size();			
			for (int j = 0; j < sizep; ++j) {
				Assembly assemblyFirst = prevallAssemblies.get(j);
				Vector<Subunit> allvec = assemblyFirst.getSubunits();
				//System.out.println(allvec.toString());
				
				int k = allvec.size();
				 //Commented out vec_out for loop 10/17/2015 for SLSversion jar file.
				for (int l = 0; l < k; ++l) {
					Subunit nxt = allvec.elementAt(l);
					//vec_out.append("ID ");
					vec_out_reduced_times.append(nxt.getID()+" ");
					//vec_out.append(" Assembly ");
					vec_out_reduced_times.append(assemblyFirst.getID()+" ");
					//vec_out.append(" Location ");
					vec_out_reduced_times.append(nxt.getPositionReal().getX() + " ");
					vec_out_reduced_times.append(nxt.getPositionReal().getY() + " ");
					vec_out_reduced_times.append(nxt.getPositionReal().getZ() + " ");
					//vec_out.append(" Partners ");
					int sze = nxt.getBoundSubunits().size();
					for (int t = 0;t < sze; t++){					
						
						vec_out_reduced_times.append(nxt.getBoundSubunits().get(t).getID()+" ");
					}
					if (sze < 6){
						int remainder = 6 - sze;
						for (int f = 0; f < remainder; f++){
							vec_out_reduced_times.append("0 ");							
						}					
					}
					vec_out_reduced_times.append(prevtime+" ");
					//System.out.println(nxt.toString());
					
					if (Test.caseq == 1){
					//MT: feb26 2016. Append quaternion that will transform initial (t=0) subunit
					//    orientation to current subunit orientation. 
					String s = tR.get(nxt.getID()).toString();
					s = s.replace("(", "");
					s = s.replace(",", "");
					s = s.replace(")", "");
					vec_out_reduced_times.append(s + " "); 
					}
					
					//MT: feb29 2016. Append subunit type (for ccmv, there is "only0" and "only1")
					//    Treat store only0 as 0, only1 as 1.
					String st = nxt.getSubunitTypeName();
					st = st.replace("only", "");
					vec_out_reduced_times.append(st + " ");

					vec_out_reduced_times.append("\n");					
				}											
			}
			
			// Reset Test.casep here since we're finished the Test.casep=1 case.
			Test.casep = 0;
		}
		
		if (curtime == 0.0 && Test.casep == 0){
			
			prevSizeDistribution = currentSizeDistribution.clone();
			prevallAssemblies = (Vector<Assembly>) allAssemblies.clone();
			prevtime = curtime;
		    
			 // SIM SECTION.
				sim_out_reduced_times.append(prevtime + " ");
				for (int i = 0; i < prevSizeDistribution.length && 
						i < Test.maxOutputSize; ++i)
						{						
							sim_out_reduced_times.append(prevSizeDistribution[i] + " ");
						}
				sim_out_reduced_times.append("\n");
				
				// VEC SECTION.
				int sizep = prevallAssemblies.size();			
				for (int j = 0; j < sizep; ++j) {
					Assembly assemblyFirst = prevallAssemblies.get(j);
					Vector<Subunit> allvec = assemblyFirst.getSubunits();
					int k = allvec.size();
					 
					for (int l = 0; l < k; ++l) {
						Subunit nxt = allvec.elementAt(l);
						//vec_out.append("ID ");
						vec_out_reduced_times.append(nxt.getID()+" ");
						//vec_out.append(" Assembly ");
						vec_out_reduced_times.append(assemblyFirst.getID()+" ");
						//vec_out.append(" Location ");
						vec_out_reduced_times.append(nxt.getPositionReal().getX() + " ");
						vec_out_reduced_times.append(nxt.getPositionReal().getY() + " ");
						vec_out_reduced_times.append(nxt.getPositionReal().getZ() + " ");
						//vec_out.append(" Partners ");
						int sze = nxt.getBoundSubunits().size();
						for (int t = 0;t < sze; t++){					
							
							vec_out_reduced_times.append(nxt.getBoundSubunits().get(t).getID()+" ");
						}
						if (sze < 6){
							int remainder = 6 - sze;
							for (int f = 0; f < remainder; f++){
								vec_out_reduced_times.append("0 ");							
							}					
						}
						vec_out_reduced_times.append(prevtime+" ");
						//System.out.println(nxt.toString());
						
						if (Test.caseq == 1){
						//MT: feb26 2016. Append quaternion that will transform initial (t=0) subunit
						//    orientation to current subunit orientation. 
						String s = tR.get(nxt.getID()).toString();
						s = s.replace("(", "");
						s = s.replace(",", "");
						s = s.replace(")", "");
						vec_out_reduced_times.append(s + " "); 
						}
						
						//MT: feb29 2016. Append subunit type (for ccmv, there is "only0" and "only1")
						//    Treat store only0 as 0, only1 as 1.
						String st = nxt.getSubunitTypeName();
						st = st.replace("only", "");
						vec_out_reduced_times.append(st + " ");						
			
						vec_out_reduced_times.append("\n");					
					}											
				}
				
				// Reset Test.casep here since we're finished the Test.casep=1 case.
				Test.casep = 0;
		    
		}
		
		
		// June 29, 2016. Marcus T and Huiming Xia modification.
		// Keep track of state of system for use during next simulation time step.
		// If we cross an experimental time point (i.e. if Test.casep <- 1), we can print
		// the 'prev' state of the system.
		prevSizeDistribution = currentSizeDistribution.clone();
		prevallAssemblies = (Vector<Assembly>) allAssemblies.clone();
		prevtime = curtime;
		
		
		
		// Commented out for reduced times version of Dessa. (Only output given experimental times)
		/*
		float rTheta=(float) (Test.kXcXw *( num / denom));
		scatter_out.append(rTheta);
		
		scatter_out.append("\n");
		sim_out.append("\n");
		output.append("\n");		
		*/
		
		// This print statement prints the sim_out info.
		//System.out.print(output);
 
	}
	
	

	
	/**
	 * A class used when selecting <code>BreakBondEvents</code> to store what
	 * <code>BindingSites</code> have already been seen. This is done by storing
	 * their ID.
	 * 
	 * @author Tiequan Zhang
	 * @author Blake Sweeney
	 * @version 1.4
	 */
	private class CheckPairSeen {

		/** vector of double[2], which are the seen pairs */
		private Vector<double[]> seen; 

		/**
		 * A Constructor
		 */
		public CheckPairSeen() {
			seen = new Vector<double[]>();
		}



		/**
		 * Adds a seen pair
		 * 
		 * @param i - first ID
		 * @param j - second ID
		 */
		public void addPair(int i, int j) { 
			seen.add(new double[] { i, j });
		}

		/**
		 * Checks if a pair has been seen. Both must have been seen together
		 * for it to be considered seen.
		 * 
		 * @param i - First ID to look for
		 * @param j - Second ID to look for
		 * @return <code>boolean</code> - true if both ID's have been seen together
		 */
		public boolean marked(int i, int j) {

			double[] tmp;
			int s = seen.size();
			for (int k = 0; k < s; ++k) {
				tmp = (double[]) (seen.get(k));
				if ((tmp[0] == i && tmp[1] == j)
						|| (tmp[0] == j && tmp[1] == i)) {
					return true;

				}
			}
			return false;
		}


	}
}
