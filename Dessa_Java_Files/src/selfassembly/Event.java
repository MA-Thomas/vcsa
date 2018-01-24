package selfassembly;


/**
 * Event implements Comparable so that two events can be compared based on the
 * end time. The posting time of an Event can be checked against the valid time
 * of Assemblies involved to see if the event has been made invalid by other
 * events on the Assemblies involved.
 * 
 * @author Sue Yi Chew
 * @author Rori Rohlfs
 * @author Tiequan Zhang
 * @author Blake Sweeney
 * @version 1.4
 *  
 */


public class Event implements Comparable<Event> {

	/**
	 * The time this event was posted (when it was created) of this event. It
	 * doesn't matter whether the event will actually occur or not
	 */
	private double postTime;
	/**
	 * the time that this event will be complete (assuming it is not invalidated
	 * before then. Determines sorting order
	 */
	private double endTime; 
	/**
	 * Array of names of Assemblies involved. If this is a ConfChangeEvent or a
	 * BreakBondEvent, only first element is used, second is null. If this a
	 * FormBondEvent, the second element may or may not be used
	 */
	private Assembly[] assembliesInvolved;

	private BindingSite bs;

	private BindingSite partner;

	private Subunit sub;

	private int domainID;

	private Conformation conf;
	
	private EventType evtType;
	
	@SuppressWarnings("unused")
	private Event(){}
	/*
	 * Constructor: FormBondEvent And BreakBondEvent
	 */
	public Event(double curTime, double eTime, Assembly[] assems, 
			BindingSite bindSite, BindingSite partnerBindSite, EventType type)
	{
		postTime = curTime;
		endTime = eTime;
		assembliesInvolved = assems;
		bs = bindSite;
		partner = partnerBindSite;
		sub = null;
		domainID = -1;
		conf = null;
		evtType = type;
	}
	
	public Event(double curTime, double eTime, Assembly[] assems,
			Subunit s, int domID, Conformation cnf)
	{
		postTime = curTime;
		endTime = eTime;
		assembliesInvolved = assems;
		bs = null;
		partner = null;
		sub = s;
		domainID = domID;
		conf = cnf;
		evtType = EventType.cnfChngEvt;
	}

	public final EventType getEventType()
	{
		return evtType;
	}
	public final double getPostTime() {
		return postTime;
	}

	public final void setPostTime(double pTime) {
		postTime = pTime;
	}

	public final double getEndTime() {
		return endTime;
	}

	public final void setEndTime(double eTime) {
		endTime = eTime;
	}

	public final Assembly[] getAssembliesInvolved() {
		return assembliesInvolved;
	}

	public final void setAssembliesInvolved(Assembly[] assemInvolved) {
		assembliesInvolved = assemInvolved;
	}

	public final BindingSite getBS() {
		return bs;
	}

	public final BindingSite getPartner() {
		return partner;
	}

	public final Subunit getSubunit() {
		return sub;
	}

	public final int getDomainID() {
		return domainID;
	}

	public final Conformation getNewConf() {
		return conf;
	}
	
	   /**
     * Compares Object e and this Event based on their end time.
     * 
     * @param e
     *            Object
     * 
     * @return 1 if this Event's end time is earlier than Object e, 0 if they
     *         are the same, -1 if this Event's end time is greater than Object
     *         e
     */
    public int compareTo(Event e) {
    	
    	if ((this == null) && (e == null)) 
    		return 0;
    	else if (this == null) 
    		return -1;
    	else if (e == null) 
    		return 1;
    	else 
    		return -(Double.compare(endTime, e.getEndTime()));
    }



    /**
     * Checks the validity of this Event. Event is valid if it was post after
     * the most recent events for the Assemblies involved. That is, if the
     * posting time of this event is greater than or equal to validTimes of each
     * Assembly involved. validTimes[1] may be null if only one assembly was
     * involved in this event.
     * 
     * @param validTimes
     *            Array of validTimes of all assemblies this event involves
     * @return boolean
     */
    public boolean isValid(double[] validTimes) {
        for (int i = 0; i < validTimes.length; ++i) {
            if (postTime < validTimes[i]) {
                return false;
            }
        }
        return true;
    }

}
