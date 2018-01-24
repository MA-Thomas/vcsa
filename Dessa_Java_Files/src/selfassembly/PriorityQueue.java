package selfassembly;


/**
 * Interface for a priority queue.
 * <p>
 * 
 * The <tt>remove</tt> method returns an item of highest priority from the
 * queue.
 * <p>
 * 
 * If distinct items are of the same priority, an implementation is not obliged
 * to return them in any particular order, unless it declares otherwise.
 * <p>
 * 
 * @author Peter Williams
 * @see <a
 *      href="http://java.sun.com/products/j2se/1.4/docs/api/java/lang/Comparable.html">Comparable
 *      </a>
 */

public interface PriorityQueue {
    /**
     * Indicates the status of the queue.
     * 
     * @return <code>true</code> if the queue is empty.
     */
    public boolean isEmpty();

    /**
     * Returns the current size of the queue
     * 
     * @return the current size of the queue
     */
    public int size();

    /**
     * Adds an item to the queue.
     * 
     * @param item
     *            the item to be added.
     */
    public void add(Comparable item);

    /**
     * Removes an item of highest priority from the queue.
     * 
     * @return an item of highest priority.
     * @exception NoSuchElementException
     *                if the queue is empty.
     */
    public Comparable remove();
}