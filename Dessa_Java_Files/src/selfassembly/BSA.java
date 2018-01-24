package selfassembly;


/**
 * A class to store a <code>BindingSite</code>, <code>Subunit</code> 
 * and an <code>Assembly</code>.
 * 
 * @author Tiequan Zhang
 * @version 1.3
 */
public class BSA {
	
    public BindingSite bs;

    public Subunit sub;

    public Assembly assem;

    public BSA(BindingSite bs0, Subunit sub0, Assembly assem0) {
        bs = bs0;
        sub = sub0;
        assem = assem0;
    }
}