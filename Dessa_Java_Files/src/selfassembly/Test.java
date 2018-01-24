package selfassembly;

import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Random;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

import java.util.Arrays;
import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Array;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * The <code>Test</code> class loads a <code>Simulation</code> from XML, starts
 * it and lists important control variables.
 * 
 * @author Rori Rohlfs
 * @author Tiequan Zhang
 * @author Blake Sweeney
 * @author Marcus Thomas
 * @version 2017
 */
public class Test {

	/* ******************** Output Modifiers ******************* */

	/** Controls the number of <code>Event</code>s between outputs */
	public static int eventsPerPrint = 1000;

	/** Weather or not to print the <code>Conformation</code> distribution */
	public static boolean printConfDist = false;

	/** Maximum size for which distribution statistics are printed */
	public static int maxOutputSize = 8;

	/** MT - Controls which events are printed according to their time (curtime) */
	// public static double[] ranges = new double[2]; // manually input this
	// included are the exp times 0.00086775s-12.952s AS WELL AS each of these
	// times minus 0.001s AS WELL AS each minus 0.002s and minus 0.003s.
	// Why include the second/third/fourth parts? So enough times are
	// available to the PI code, considering that
	// Simulation.java selects from expTimeRanges based on first 3 digits.
	static final Set<Double> ranges = new HashSet<Double>(
			Arrays.asList(new Double[] { 					

                    
                    0.01,0.21,0.41,0.61,0.81,1.01,1.21,
                    1.41,1.61,1.81,2.01,2.21,2.41,2.61,2.81,
                    3.01,3.21,3.41,3.61,3.81,4.01,4.21,4.41,4.61,4.81,
                    5.01,5.21,5.41,5.61,5.81,6.01,6.21,6.41,6.61,6.81,
                    7.01,7.21,7.41,7.61,7.81,8.01,8.21,8.41,8.61,8.81,
                    9.01,9.21,9.41,9.61,9.81,10.01,10.21,10.41,10.61,10.81,
                    11.01,11.21,11.41,11.61,11.81,12.01,12.21,12.41,12.61,12.81,
                    13.01,13.21,13.41,13.61,13.81,14.01,14.21,14.41,14.61,14.81,
                    15.01,15.21,15.41,15.61,15.81,16.4,17.0,17.6,18.2,18.8,19.4,20.0
                    
                    /*                   
                    21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,39.0,40.0,
                    41.0,42.0,43.0,44.0,45.0,46.0,47.0,48.0,49.0,50.0,51.0,52.0,53.0,54.0,55.0,56.0,57.0,58.0,59.0,60.0,
                    61.0,62.0,63.0,64.0,65.0,66.0,67.0,68.0,69.0,70.0,71.0,72.0,73.0,74.0,75.0,76.0,77.0,78.0,79.0,80.0,
                    81.0,82.0,83.0,84.0,85.0,86.0,87.0,88.0,89.0,90.0,91.0,92.0,93.0,94.0,95.0,96.0,97.0,98.0,99.0,
                    100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0, 118.0, 119.0, 120.0,
                    121.0, 122.0, 123.0, 124.0, 125.0, 126.0, 127.0, 128.0, 129.0, 130.0, 131.0, 132.0, 133.0, 134.0, 135.0, 136.0, 137.0, 138.0, 139.0, 140.0, 141.0,
                    142.0, 143.0, 144.0, 145.0, 146.0, 147.0, 148.0, 149.0, 150.0, 151.0, 152.0, 153.0, 154.0, 155.0, 156.0, 157.0, 158.0, 159.0, 160.0, 161.0, 162.0, 
                    163.0, 164.0, 165.0, 166.0, 167.0, 168.0, 169.0, 170.0, 171.0, 172.0, 173.0, 174.0, 175.0, 176.0, 177.0, 178.0, 179.0, 180.0, 181.0, 182.0, 183.0,
                    184.0, 185.0, 186.0, 187.0, 188.0, 189.0, 190.0, 191.0, 192.0, 193.0, 194.0, 195.0, 196.0, 197.0, 198.0, 199.0, 200.0, 201.0, 202.0, 203.0, 204.0, 
                    205.0, 206.0, 207.0, 208.0, 209.0, 210.0, 211.0, 212.0, 213.0, 214.0, 215.0, 216.0, 217.0, 218.0, 219.0, 220.0, 221.0, 222.0, 223.0, 224.0, 225.0, 
                    226.0, 227.0, 228.0, 229.0, 230.0, 231.0, 232.0, 233.0, 234.0, 235.0, 236.0, 237.0, 238.0, 239.0, 240.0, 241.0, 242.0, 243.0, 244.0, 245.0, 246.0, 
                    247.0, 248.0, 249.0, 250.0, 251.0, 252.0, 253.0, 254.0, 255.0, 256.0, 257.0, 258.0, 259.0, 260.0, 261.0, 262.0, 263.0, 264.0, 265.0, 266.0, 267.0, 
                    268.0, 269.0, 270.0, 271.0, 272.0, 273.0, 274.0, 275.0, 276.0, 277.0, 278.0, 279.0, 280.0, 281.0, 282.0, 283.0, 284.0, 285.0, 286.0, 287.0, 288.0, 
                    289.0, 290.0, 291.0, 292.0, 293.0, 294.0, 295.0, 296.0, 297.0, 298.0, 299.0, 300.0, 301.0, 302.0, 303.0, 304.0, 305.0, 306.0, 307.0, 308.0, 309.0, 310.0
                    */

			}));
	
	/* ******************** General Values ******************* */

	/** The name for this Simulation type */
	public static String simType = "xml";

	/** Which xml file to use */
	public static String mode = "";

	/** Mass for a <code>Subunit</code> */
	public static double subMass = 3.4;
	
	/** Case (either 0 don't print or 1 print) to indicate if we print data 
	 * at prev simulation time point. Initialize here for use in Simulation.java */
	public static double casep = 0;

	/** Case (either 0 don't print or 1 print) to indicate if we print quaternions
	 * for subunit rotations in Simulation.java */
	public static double caseq = 0;	

	/**
	 * Ratio of <code>Subunit</code> to <code>BindingSite</code>, should be
	 * <=0.5
	 */
	public static double subunitToBS = 0.5;

	/** Height of a <code>BindingSite</code> */
	public static float bindingSiteHeight = 0.10f;

	/** Radius for a <code>Subunit</code> */
	public static double subRadius = bindingSiteHeight * subunitToBS;

	/* ******************** Simulation Variables ******************* */

	/** Use the diffusion model */
	// public static boolean diffusionEnabled = false;

	/** Pick events that only allow monomer binding */
	public static boolean bindMonomerOnly = false;

	/** only sample breaking event for bonds not involving a loop */
	public static boolean noLoopOnly = true;

	/** For filament assembly, allow only monomers break off */
	public static boolean breakOnlyEnds = false;

	/** Enable conformational switching event sampling */
	public static boolean csAllowed = false;//true;

	/**
	 * The int returned from Assemblies.numbSubunits() that will be treated as a
	 * monomer
	 */
	public static int sizeOfSubunit = 1;

	/** Maximum allowed assembly size */
	public static int maxLength = 360;

	/** The amount of volume to search for nearest neighbors */
	public static double binSize = 0.1;

	/** The distance tolerance for loop detection */
	public static double distanceTolerance = 0.1;

	/** Controls weather or not the Spring Force model is used */
	public static boolean springForce = false;

	/* ******************** Random Number Generator ******************* */

	/** The Random number generator for the simulation */
	public static Random rand;

	/** The random number generator for the seed */
	public static Random seedGenerator;

	/** The seed used for the random number generator */
	public static int seedUsed;

	/**
	 * MAx Simulation time allowed. After which simulation turns off
	 * automatically
	 */
	public static double maxSimulationTime;

	public static double kXcXw;

	/**
	 * main method to start the simulator
	 */
	public static void main(String args[]) {

		/** MT */
		returnRun simulatorRun = null;

		if (args.length == 0) {
			System.out
					.println("usage: xmlFile.xml [Events Per Printout] "
							+ "[Max Simulation Time] [k-Const] [c-Const] [mol. wt Const] [Random Seed] [Max Output Size] ");
			System.exit(-1);
		}

		mode = args[0];
		double k;
		double c;
		double m;
		try {
			eventsPerPrint = Integer.parseInt(args[1]);
		} catch (Exception e) {
			eventsPerPrint = 10000;
		}

		try {
			maxSimulationTime = Double.parseDouble(args[2]);
		} catch (Exception e) {
			maxSimulationTime = Double.MAX_VALUE;
		}
		try {// k const
			k = Double.parseDouble(args[3]);
		} catch (Exception e) {
			k = 2.56e-7;
		}

		try {// c const
			c = Double.parseDouble(args[4]);
		} catch (Exception e) {
			c = 145e-6;
		}

		try {// mol weight const
			m = Double.parseDouble(args[5]);
		} catch (Exception e) {
			m = 250e3;
		}

		kXcXw = k * c * m;// constants for rTheta calculation

		try {
			seedUsed = Integer.parseInt(args[6]);
		} catch (Exception e) {
			seedGenerator = new Random();
			//seedUsed = 821820232; //mt 7/9/2015
			//seedUsed = 822820230;
			seedUsed = seedGenerator.nextInt();
		}

		rand = new Random(seedUsed);

		try {
			maxOutputSize = Integer.parseInt(args[7]);
		} catch (Exception e) {
			maxOutputSize = 200;
		}

		// System.out.println(String.format("Input : %s Seed: %d (Dessa v1.5.7, Bio.Phys.)",mode,
		// seedUsed));
		System.out.println("314 628 999 999 999");
		XMLReader reader = null;
		try {
			reader = new XMLReader(mode);
		} catch (Exception e) {
			e.printStackTrace();
		}

		Simulation sim = reader.getSim();

		sim.run();

		// System.out.println(Arrays.deepToString(m1));
	}

	/**
	 * Rotates a Vector around an axis by the specified angle.
	 * 
	 * @param v
	 *            - The vector to rotate
	 * @param a
	 *            - the AxisAngle to rotate around
	 * @return Vector3d - the rotated Vector3d
	 */
	public static Vector3d rotateByAxis(Vector3d v, AxisAngle4d a) {

		Matrix4d m = new Matrix4d();
		m.set(a);
		Vector3d rv = new Vector3d();
		m.transform(v, rv);
		return rv;
	}

	/**
	 * 
	 * @param axis
	 * @param upV
	 * @return
	 */
	public static Vector3d makePerpendicular(Vector3d axis, Vector3d upV) {

		double angle = axis.angle(upV);

		Vector3d projectedV = new Vector3d();

		if (angle < 0.00001 || (Math.PI - angle) < 0.00001) {
			System.out.println("Bad_choice_in_upVector");
			System.exit(0);
		} else {

			Vector3d tem = new Vector3d(axis);
			tem.scale(axis.dot(upV) / axis.lengthSquared());

			projectedV.sub(upV, tem);
		}
		projectedV.normalize();
        //MT: see wiki on Vector projection. 
		//MT: projectedV ~ a2
		//MT: axis ~ b
		//MT: upV ~ a
		//MT: tem ~ a1
		//MT: NOTE - based on wiki, tem should use axis.length(), not
		//MT: axis.lengthsquared(). But this is only a scalar difference
		return projectedV;
	}
}
