import java.io.*;
import java.util.*;

/* Isotoper
   Initiated August 14, 2002
   David L. Tabb, The Scripps Research Institute

   * Isotoper accepts a molecule's composition and returns an "isotope
   * envelope" which shows the distribution of the molecule into
   * isotopic peaks.  In order to cut down time of calculation, the
   * program accepts a maximum number of peaks to keep in the
   * envelope.
   */

public class Isotoper {
    ElementEnvelopes     Elements = new ElementEnvelopes();
    int                  MaximumDistanceFromMono = 2;
    int                  MaximumPowerForElements = 10;
    IsotopeEnvelope      LastIE = new IsotopeEnvelope();
    
    public Isotoper() {
	this.ReadInitFile();
	this.BuildPowerSeries();
	LastIE.MakeUnity();
    }

    public void ClearComposition() {
	ElementEnvelopes      ERunner = Elements.Next;
	while (ERunner != null) {
	    ERunner.CompositionCount = 0;
	    ERunner = ERunner.Next;
	}
    }

    public void SetComposition(String Element, int Count) {
	// Set the composition of the compound for which isotope
	// distribution is to be generated.
	ElementEnvelopes      ERunner = Elements.Next;
	while ( (ERunner != null) &&
		(!ERunner.Symbol.equals(Element)) )
	    ERunner = ERunner.Next;
	if (ERunner == null) {
	    System.out.println("Element " + Element + " not defined.");
	}
	else {
	    ERunner.CompositionCount = Count;
	}
    }

    public void CalculateEnvelope() {
	ElementEnvelopes  ERunner = Elements.Next;
	IsotopeEnvelope   Builder = new IsotopeEnvelope();
	IsotopeEnvelope   ElementBuilder;
	Builder.MakeUnity();
	while (ERunner != null) {
	    // System.out.println(ERunner.Symbol + "\t" +
	    //    ERunner.CompositionCount);
	    ElementBuilder = ERunner.GetElementalEnvelope();
	    // System.out.println(ERunner.Symbol + "\t" +
	    //    ElementBuilder.ReportHeights());
	    Builder = Builder.SplitBy(ElementBuilder);
	    ERunner = ERunner.Next;
	}
	LastIE = Builder;
	// System.out.println(Builder.ReportHeights());
    }

    public float ReportHeight(int Offset) {
	float      InitialIntensity;
	Intensity  IRunner;
	if (LastIE != null) {
	    IRunner = LastIE.Peaks.Next;
	    InitialIntensity = IRunner.Height;
	    while ( (IRunner != null) &&
		    (Offset > 0) ) {
		IRunner = IRunner.Next;
		Offset--;
	    }
	    if (IRunner != null) {
		return (IRunner.Height / InitialIntensity);
	    }
	    else {
		System.out.println("A");
		return 0.0f;
	    }
	}
	else {
	    System.out.println("B");
	    return 0.0f;
	}
    }

    public void ReadInitFile() {
	ElementEnvelopes    ERunner = Elements;
	Intensity           IRunner;
	try {
	    File            CurrentDirectory = new File(System.getProperty("user.dir"));
	    File            IsotoperIni = new File(CurrentDirectory, "Isotoper.ini");
	    if (IsotoperIni.exists()) {
		FileReader      InputFileReader = new FileReader(IsotoperIni);
		BufferedReader  Incoming = new BufferedReader(InputFileReader);
		String          LineBuffer;
		String          WholeLine;
		StringTokenizer Parser;
		WholeLine = Incoming.readLine();
		while (WholeLine != null) {
		    Parser = new StringTokenizer(WholeLine);
		    if (Parser.hasMoreTokens()) {
			LineBuffer = Parser.nextToken();
			if (LineBuffer.equals("Ele")) {
			    // Add another element to our list
			    ERunner.Next = new ElementEnvelopes();
			    ERunner = ERunner.Next;
			    ERunner.Symbol = Parser.nextToken();
			    IRunner = ERunner.PowerSeries.Peaks;
			    while (Parser.hasMoreTokens()) {
				IRunner.Next = new Intensity();
				IRunner = IRunner.Next;
				IRunner.Height = new Float(Parser.nextToken()).floatValue();
			    }
			    // Normalize sum of peaks to 1.0
			    float   IntensitySum = 0f;
			    IRunner = ERunner.PowerSeries.Peaks.Next;
			    while (IRunner != null) {
				IntensitySum += IRunner.Height;
				IRunner = IRunner.Next;
			    }
			    IRunner = ERunner.PowerSeries.Peaks.Next;
			    while (IRunner != null) {
				IRunner.Height /= IntensitySum;
				IRunner = IRunner.Next;
			    }
			}
		    }
		    WholeLine = Incoming.readLine();
		}
		//		    this.Acids.DebugPrint();
		// Close file
		Incoming.close();
	    }
	    else {
		System.out.println("Could not find Isotoper.ini.");
		System.exit(0);
	    }
	}
	catch (IOException oopsie) {
	    System.out.println("Error while reading Isotoper.ini");
	    System.out.println(oopsie.toString());
	}
    }

    public void BuildPowerSeries() {
	ElementEnvelopes   ERunner = Elements.Next;
	IsotopeEnvelope    IERunner;
	int                Counter;
	while (ERunner != null) {
	    IERunner = ERunner.PowerSeries;
	    for (Counter = 0; Counter < MaximumPowerForElements; Counter++) {
		IERunner.Next = IERunner.SplitBy(IERunner);
		IERunner = IERunner.Next;
	    }
	    ERunner = ERunner.Next;
	}
	ERunner = Elements.Next;
	while (ERunner != null) {
	    IERunner = ERunner.PowerSeries;
	    Counter = 0;
	    while (IERunner != null) {
		/*
		System.out.print(ERunner.Symbol + "\t" + Math.pow(2,Counter));
		System.out.println(IERunner.ReportHeights());
		*/
		IERunner = IERunner.Next;
		Counter++;
	    }
	    ERunner = ERunner.Next;
	}
    }

    public class ElementEnvelopes {
	String           Symbol;
	int              CompositionCount = 0;
	IsotopeEnvelope  PowerSeries = new IsotopeEnvelope();
	ElementEnvelopes Next;

	public IsotopeEnvelope GetElementalEnvelope() {
	    IsotopeEnvelope   Builder = new IsotopeEnvelope();
	    IsotopeEnvelope   IERunner = this.PowerSeries;
	    int               Power = 1;
	    Builder.MakeUnity();
	    if (CompositionCount > Math.pow(2,MaximumPowerForElements)) {
		System.out.println("Too complex a molecule for isotope calculation.");
		System.out.println(Symbol + "\t" + CompositionCount);
	    }
	    while (Power < CompositionCount + 1) {
		if ( (CompositionCount & Power) > 0) {
		    Builder = Builder.SplitBy(IERunner);
		    // System.out.println(Power + "\t" + IERunner.ReportHeights());
		}
		Power *= 2;
		IERunner = IERunner.Next;
	    }
	    return Builder;
	}
    }

    public class IsotopeEnvelope {
	Intensity        Peaks = new Intensity();
	IsotopeEnvelope  Next;

	public void MakeUnity() {
	    Peaks.Next = new Intensity();
	    Peaks.Next.Height = 1.0f;
	}

	public String ReportHeights() {
	    StringBuffer    Returned = new StringBuffer();
	    Intensity       PRunner = this.Peaks.Next;
	    while (PRunner != null) {
		Returned.append("\t");
		Returned.append(PRunner.Height);
		PRunner = PRunner.Next;
	    }
	    return Returned.toString();
	}

	// Split this isotope envelope with the passed one, returning
	// the result
	public IsotopeEnvelope SplitBy (IsotopeEnvelope Other) {
	    IsotopeEnvelope Builder = new IsotopeEnvelope();
	    Intensity       ThisPeak = this.Peaks.Next;
	    Intensity       OtherPeak;
	    Intensity       BuilderPeak = Builder.Peaks;
	    int             ThisOffset = 0;
	    int             OtherOffset;
	    int             Counter;
	    // Generate an empty set of peaks in which to store
	    // isotopic intensities.
	    for (Counter = -1; Counter < MaximumDistanceFromMono; Counter++) {
		BuilderPeak.Next = new Intensity();
		BuilderPeak = BuilderPeak.Next;
		BuilderPeak.Height = 0f;
	    }
	    // Loop once through each peak in this intensity envelope
	    while (ThisPeak != null) {
		OtherPeak = Other.Peaks.Next;
		OtherOffset = 0;
		// Advance to the correct place to be storing intensity
		BuilderPeak = Builder.Peaks.Next;
		for (Counter = 0; Counter < ThisOffset; Counter++) {
		    BuilderPeak = BuilderPeak.Next;
		}
		while ( !(MaximumDistanceFromMono < OtherOffset + ThisOffset) &&
			(OtherPeak != null) ) {
		    BuilderPeak.Height += OtherPeak.Height * ThisPeak.Height;
		    BuilderPeak = BuilderPeak.Next;
		    OtherOffset++;
		    OtherPeak = OtherPeak.Next;
		}
		ThisOffset++;
		ThisPeak = ThisPeak.Next;
	    }
	    /*
	    // Normalize sum of intensities to 1.0;
	    float           IntensitySum = 0f;
	    BuilderPeak = Builder.Peaks.Next;
	    while (BuilderPeak != null) {
	    IntensitySum += BuilderPeak.Height;
	    BuilderPeak = BuilderPeak.Next;
	    }
	    BuilderPeak = Builder.Peaks.Next;
	    while (BuilderPeak != null) {
	    BuilderPeak.Height = BuilderPeak.Height / IntensitySum;
	    BuilderPeak = BuilderPeak.Next;
	    }
	    */
	    return Builder;
	}
    }

    public class Intensity {
	float            Height;
	Intensity        Next;
    }

}
