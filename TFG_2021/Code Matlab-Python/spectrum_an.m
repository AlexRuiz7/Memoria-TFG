    %%Spectrum Analyzer Setting up:
    
    spectrum.center=1679630;    %Center Frequency
    spectrum.units=signal.units;%Frequency units
    spectrum.start=810.000;     %Start Frequency
    spectrum.stop=2549.260;     %Stop Frequency
    spectrum.reflevel=10;       %Reference Level
    spectrum.powerUnits='DBM';  %Power units for Reference Level
    spectrum.points=101;       %Points number for the displayed trace. OPTIONS: 11,21,41,51,101,201,251,401,501,1001,2001,5001,10001  
    
    %Now the function "Anritsu_MS2830A_V01" is used, to send the spectryum analyzer configuration to the case "set_SpectrumAnalyzer"    
    Anritsu_MS2830A_V01(instrument,'set_SpectrumAnalyzer',spectrum);  
    