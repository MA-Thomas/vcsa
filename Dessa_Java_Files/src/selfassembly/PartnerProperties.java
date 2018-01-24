package selfassembly;

public class PartnerProperties {
	private double a;//angle
	private double b;//bindTime
	private double bb;//breakTime
	private double fb;//fastBindTime
	
	@SuppressWarnings("unused")
	private PartnerProperties(){}
	
	public PartnerProperties(double angle, double bindTime, double breakTime, double fastBindTime)
	{
		a = angle;
		b = bindTime;
		bb = breakTime;
		fb = fastBindTime;
	}
	
	public final double getAngle(){ return a;}
	public final double getBindTime(){ return b;}
	public final double getBreakTime(){ return bb;}
	public final double getFastBindTime(){ return fb;}

}
