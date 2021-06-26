switch command
        %%Case 'open'. Sets the information needed for the GPIB
        %%comunication. Initially, no need to change for this instrument.
        
            case 'open'
                serialPort = instrument.serialPort;
                serialObject = serial(serialPort, 'baud', 115200,'StopBits',1 ,'DataBits', 8,'Parity', 'none', 'InputBufferSize',100000,'OutputBufferSize',1000); %
                set(serialObject,'Terminator','LF');
                set(serialObject,'Timeout',900);
                fopen(serialObject);
                pause(0.5);
                fwrite(serialObject,uint8(sprintf('++addr %d\r\n',instrument.addrGPIB)));
                fwrite(serialObject,uint8(sprintf('*CLS\r\n')));
                pause(0.5);
            
        %%Case 'close'. Close GPIB Comunication.             
            case 'close'
                fclose(serialObject);
                delete(instrfind);