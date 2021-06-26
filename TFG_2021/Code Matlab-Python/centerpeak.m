       %%Case 'center_peak'. Sets the spectrum analyzer so the main peak
        %%is displayed. It is calculated through a for loop with 3
        %%iterations.
        
        %(PAGES REFERENCE TO SPECTRUM ANALIZER REMOTE MANUAL:)
        
        %INST SPECT - Sets the instrument which is goint to receive the string
        %CALC:MARK:AOFF - P160
        %CALC:MARK:ACT - P115
        %CALC:MARK:RES PEAK - P147
        %CALC:MARK:MAX - P171
        %CALC:PMAR:Y? - P158
        %CALC:MARK:X - P121
        %CALC:MARK:WIDT? - P140
        %\r\n - Compound string is send after them.
        
        case 'center_peak'   
            
            clc
            
            for i=1:3
                
            pause(1);           
            fwrite(serialObject,uint8(sprintf('INST SPECT; CALC:MARK:AOFF\r\n')));
            pause(1);         
            fwrite(serialObject,uint8(sprintf('INST SPECT; CALC:MARK:ACT ON; CALC:MARK:RES PEAK; CALC:MARK:MAX\r\n'))); %; CALC:PMAR:Y?
            pause(1);
            %fwrite(serialObject,uint8(sprintf('INST SIGANA; SYST:COMM:GPIB:SELF:DEL LF\r\n')));
            fwrite(serialObject,uint8(sprintf('INST SPECT; CALC:MARK:X?\r\n'))); 
            fcenter=fscanf(serialObject);     
            pause(1);
            fcenter=str2double(fcenter)            
            fwrite(serialObject,uint8(sprintf('INST SPECT; CALC:MARK:WIDT?\r\n')));
            pause(1);
            marker_span=fscanf(serialObject);
            pause(1);
            marker_span=str2double(marker_span)
            f0=fcenter-marker_span
            f1=fcenter
            f2=fcenter+marker_span
            fwrite(serialObject,uint8(sprintf('INST SPECT; FREQ:START %.0f HZ; FREQ:CENT %.0f HZ; FREQ:STOP %.0f HZ\r\n',f0,f1,f2)));             
            pause(1.5); 
            
            end