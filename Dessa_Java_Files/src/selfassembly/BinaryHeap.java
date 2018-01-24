package selfassembly;


import java.util.NoSuchElementException;

/**
 * Binary heap implementation of a priority queue.
 * <p>
 * 
 * @author Peter Williams
 */

public class BinaryHeap implements PriorityQueue {

    private Comparable[] heap; // the heap

    private int size; // number of items in the heap

    /**
     * Constructs the binary heap.
     */
    public BinaryHeap() {
        heap = new Comparable[1];
        size = 0;
    }

    /**
     * Tests if the heap is empty.
     * 
     * @return true if empty, false otherwise.
     */
    public boolean isEmpty() {
        return size == 0;
    }

    /**
     * Returns the current size of the queue
     * 
     * @return the current size of the queue
     */
    public int size() {
        return size;
    }

    /**
     * Adds and item to the heap.
     * 
     * @param item
     *            the item to add.
     */
    public void add(Comparable item) {
        // grow the heap if necessary
        if (size == heap.length) {
            Comparable[] newHeap = new Comparable[2 * heap.length];
            System.arraycopy(heap, 0, newHeap, 0, heap.length);
            heap = newHeap;
        }
        // find where to insert while rearranging the heap if necessary
        int parent, child = size++; // the next available slot in the heap
        while (child > 0 && heap[parent = (child - 1) / 2].compareTo(item) < 0) {
            heap[child] = heap[parent];
            child = parent;
        }
        heap[child] = item;
    }

    /**
     * Removes an item of highest priority from the heap.
     * 
     * @return an item of highest priority.
     * @exception NoSuchElementException
     *                if the heap is empty.
     */
    public Comparable remove() {
        if (size == 0) {
            throw new NoSuchElementException();
        }
        Comparable result = heap[0]; // to be returned
        Comparable item = heap[--size]; // to be reinserted
        int child, parent = 0;
        while ((child = (2 * parent) + 1) < size) {
            // if there are two children, compare them
            if (child + 1 < size && heap[child].compareTo(heap[child + 1]) < 0) {
                ++child;
            }
            // compare item with the larger
            if (item.compareTo(heap[child]) < 0) {
                heap[parent] = heap[child];
                parent = child;
            } else {
                break;
            }
        }
        heap[parent] = item;
        return result;
    }

}