package selfassembly;

import java.util.Vector;
import java.util.HashMap;
import javax.vecmath.*;

/**
 * NeighborGroups stores assembly subunits in neighborgroups for fast collision
 * detection.That is, it divides the subunits in an Assembly into groups of
 * physically close neighbors for that methods which need to compare physically
 * close subunits.
 * 
 * @author Rori Rohlfs
 * @author Tiequan Zhang
 * @version 1.3
 */
public class NeighborGroups {
    // neighborgroups
    public HashMap groups; //keys = x-coords,

    //                        values = hashmap
    //                                 |
    //                                 V
    //                               keys = y-coords,
    //                               values = hashmaps
    //                                        |
    //                                        V
    //                                      keys = z-coords,
    //                                      values = subunit ids in that region

    private double subunitRadius = 0.15;

    private double subdivision = 2 * subunitRadius;

    private double xMin = Double.MAX_VALUE;

    private double yMin = Double.MAX_VALUE;

    private double zMin = Double.MAX_VALUE;

    private double binSize = Test.binSize;

    /**
     * Creates default NeighborGroups
     */
    public NeighborGroups() {
        groups = new HashMap();
    }

    /**
     * Creates NeighborGroups with specified Assembly asm. actually sorts the
     * Subunits in asm into neighborgroups according to a subdivision-spaced
     * cubish grid
     * 
     * @param asm
     */
    public NeighborGroups(Assembly asm) {
        groups = new HashMap();
        Vector subs = asm.getSubunits();
        int size = subs.size();
        for (int i = 0; i < size; ++i) {
            Subunit cursub = (Subunit) subs.get(i);
            Vector3d cursubpos = cursub.getPositionReal();

            Integer x = new Integer((int) Math.round(cursubpos.x / subdivision));
            Integer y = new Integer((int) Math.round(cursubpos.y / subdivision));
            Integer z = new Integer((int) Math.round(cursubpos.z / subdivision));

            if (groups.containsKey(x)) {
                HashMap yhash = (HashMap) groups.get(x);
                if (yhash.containsKey(y)) {
                    HashMap zhash = (HashMap) yhash.get(y);
                    if (zhash.containsKey(z)) {
                        Vector entry = (Vector) zhash.get(z);
                        entry.add(cursub);
                    } else {//zhash doesn't have this z-coord yet
                        Vector entry = new Vector();
                        entry.add(cursub);
                        zhash.put(z, entry);
                    }
                } else { //yhash doesn't have have this y-value yet
                    HashMap newzhash = new HashMap();
                    Vector entry = new Vector();
                    entry.add(cursub);
                    newzhash.put(z, entry);

                    yhash.put(y, newzhash);
                }
            } else { //groups (xhash) doesn't have this x-value yet
                HashMap newzhash = new HashMap();
                Vector entry = new Vector();
                entry.add(cursub);
                newzhash.put(z, (Vector) entry);

                HashMap newyhash = new HashMap();
                newyhash.put(y, newzhash);

                groups.put(x, newyhash);

            }
        }
    }

    /**
     * Constructs NeighborGroups to find near binding sites to bind them
     * 
     * @param asm
     * @param findLoop
     */
    public NeighborGroups(Assembly asm, boolean findLoop) {
        groups = new HashMap();

        Vector subs = asm.getSubunits();
        int size = subs.size();
        int jSize = 0;
        for (int i = 0; i < size; ++i) {//find the position of the left most free BindingSite tip
            Subunit cursub = (Subunit) subs.get(i);
            Vector bss = cursub.getBindingSites();
            jSize = bss.size();
            for (int j = 0; j < jSize; ++j) {
                BindingSite bs = (BindingSite) bss.get(j);
                if (bs.isBound()) continue;
                
                Vector3d TipPos = bs.getPosReal();

                double xPos = TipPos.x;
                if (xPos < xMin)
                    xMin = xPos;
                double yPos = TipPos.y;

                if (yPos < yMin)
                    yMin = yPos;

                double zPos = TipPos.z;
                if (zPos < zMin)
                    zMin = zPos;
            }
        }
        //categorize all free BindingSite(s) in assem based on the positions of
        // their tips
        size = subs.size();
        jSize = 0;
        for (int i = 0; i < size; ++i) {
            Subunit cursub = (Subunit) subs.get(i);
            Vector bss = cursub.getBindingSites();
            jSize = bss.size();
            for (int j = 0; j < jSize; ++j) {
                BindingSite bs = (BindingSite) bss.get(j);
                if (bs.isBound()) {
                    continue;
                }

                Vector3d TipPos =bs.getPosReal();
                double xPos = TipPos.x;
                double yPos = TipPos.y;
                double zPos = TipPos.z;

                Integer x = new Integer((int) convertPos(xPos, xMin));
                Integer y = new Integer((int) convertPos(yPos, yMin));
                Integer z = new Integer((int) convertPos(zPos, zMin));

                if (groups.containsKey(x)) {
                    HashMap yhash = (HashMap) groups.get(x);
                    if (yhash.containsKey(y)) {
                        HashMap zhash = (HashMap) yhash.get(y);
                        if (zhash.containsKey(z)) {
                            Vector entry = (Vector) zhash.get(z);
                            entry.add(new BSA(bs, cursub, asm));
                        } else {//zhash doesn't have this z-coord yet
                            Vector entry = new Vector();
                            entry.add(new BSA(bs, cursub, asm));
                            zhash.put(z, entry);
                        }
                    } else { //yhash doesn't have have this y-value yet
                        HashMap newzhash = new HashMap();
                        Vector entry = new Vector();
                        entry.add(new BSA(bs, cursub, asm));
                        newzhash.put(z, entry);

                        yhash.put(y, newzhash);
                    }
                } else { //groups (xhash) doesn't have this x-value yet
                    HashMap newzhash = new HashMap();
                    Vector entry = new Vector();
                    entry.add(new BSA(bs, cursub, asm));
                    newzhash.put(z, (Vector) entry);

                    HashMap newyhash = new HashMap();
                    newyhash.put(y, newzhash);

                    groups.put(x, newyhash);

                }
            }
        }
    }

    /**
     * Copys every Subunit in newsubs to subunit Vector subs
     * 
     * @param subs
     * @param newsubs
     */
    private void addEntries(Vector subs, Vector newsubs) {
    	int size = newsubs.size();
        for (int i = 0; i < size; ++i) {
            subs.add(newsubs.get(i));
        }
    }

    private int convertPos(double posBefore, double min) {

        return (int) ((posBefore - min) / binSize);

    }

    public Vector<BSA> findLoopCandidates(Vector3d bsTipPos) {

        Vector<BSA> nearBSs = new Vector<BSA>();

        //get the tip of bs absolute pos on the level of assembly or simulation
        double xPos = bsTipPos.x;
        double yPos = bsTipPos.y;
        double zPos = bsTipPos.z;
        int x = convertPos(xPos, xMin);
        int y = convertPos(yPos, yMin);
        int z = convertPos(zPos, zMin);
        for (int i = -1; i < 2; ++i) {

            Integer tx = new Integer(x + i);
            if (groups.containsKey((Object) tx)) {
                HashMap yhash = (HashMap) groups.get(tx);
                for (int j = -1; j < 2; ++j) {
                    Integer ty = new Integer(y + j);
                    if (yhash.containsKey(ty)) {
                        HashMap zhash = (HashMap) yhash.get(ty);

                        for (int k = -1; k < 2; ++k) {

                            Integer tz = new Integer(z + k);
                            if (zhash.containsKey(tz)) {
                                Vector<BSA> entry = (Vector<BSA>) zhash.get(tz);
                                addEntries(nearBSs, entry);

                            }
                        }
                    }
                }
            }
        }

        return nearBSs;

    }

    public boolean findHindrance(Subunit sub) {

        Vector3d subpos = sub.getPositionReal();
       /* double xPos = subpos.x;
        double yPos = subpos.y;
        double zPos = subpos.z;*/

        int x = (int) Math.round(subpos.x / subdivision);
        int y = (int) Math.round(subpos.y / subdivision);
        int z = (int) Math.round(subpos.z / subdivision);
        for (int i = -1; i < 2; ++i) {

            Integer tx = new Integer(x + i);
            if (groups.containsKey((Object) tx)) {
                HashMap yhash = (HashMap) groups.get(tx);
                for (int j = -1; j < 2; ++j) {
                    Integer ty = new Integer(y + j);
                    if (yhash.containsKey(ty)) {
                        HashMap zhash = (HashMap) yhash.get(ty);

                        for (int k = -1; k < 2; ++k) {

                            Integer tz = new Integer(z + k);
                            if (zhash.containsKey(tz)) {
                                Vector entry = (Vector) zhash.get(tz);
                                int tem = entry.size();
                                for (int m = 0; m < tem; ++m) {
                                    Subunit temSub = (Subunit) entry.get(m);
                                    if (temSub.equals(sub)) {
                                        continue;
                                    }
                                    Vector3d temV3d = new Vector3d(temSub
                                            .getPositionReal());
                                    temV3d.sub(subpos);
                                    if (temV3d.length() <= 2 * subunitRadius) {
                                        return true;
                                    }

                                }
                            }
                        }
                    }
                }
            }
        }

        return false;
    }

}
