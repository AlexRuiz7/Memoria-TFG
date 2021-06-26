import time
import serial
from datetime import datetime

class RIONNA27:

    #Init object and serial port
    def __init__(self, p, debug):
        self.ser = serial.Serial(
            port=p,
            baudrate=9600,
        )
        self.debug = debug
        self.ser.isOpen()

    #Message encoding based on RION documentation
    def encodeMessage(self, command):
        packet = bytearray()
        packet.extend([0x02, 0x01, 0xFE])           #Header
        packet.extend(command)                      #Payload
        packet.extend([0x1A] * (32 - len(command))) #Padding
        packet.extend([sum(packet[3:]) % 256])      #Checksum
        return packet

    def setParam(self, cmd):
        self.ser.write(self.encodeMessage(cmd))
        while self.ser.inWaiting() < 1:
            pass
        out = self.ser.read(1)
        if self.debug:
            print ''.join('{:02x}'.format(x) for x in bytearray(out))

        #Check if ACK have been received
        if ''.join('{:02x}'.format(x) for x in bytearray(out)) == "06":
            return True
        else:
            return False

    def askParam(self, cmd):
        cmd = cmd + " ?"

        self.ser.write(self.encodeMessage(cmd))
        if self.debug:
            print "Sending " + cmd

        while self.ser.inWaiting() < 1:
            pass

        out = self.ser.read(1)

        s = ''.join('{:02x}'.format(x) for x in bytearray(out))

        #Check if ACK have been received
        if s == '06':
            if self.debug:
                print "ACK Received"
                print "Sending NACK (21)"

            #Send NACK
            self.ser.write(bytearray([21]))

            while self.ser.inWaiting() < 1:
                pass

            #Read packet length
            out = self.ser.read(1)
            tam = ''.join('{:02x}'.format(x) for x in bytearray(out))
            if tam == "01":
                while self.ser.inWaiting() < 131:
                    pass
                out += self.ser.read(self.ser.inWaiting())
            else:
                while self.ser.inWaiting() < 35:
                    pass
                out += self.ser.read(self.ser.inWaiting())

            time.sleep(0.1)


            while True:
                if self.debug:
                    print "Sending ACK (6)"
                #Send ACK
                self.ser.write(bytearray([6]))

                while self.ser.inWaiting() < 1:
                    pass

                outend = ''
                while self.ser.inWaiting() > 0:
                    r = self.ser.read(1)
                    out += r
                    outend += r

                if self.debug:
                    print "Readed: " + str(len(outend))
                    print "String: " + str(outend)
                    print "Hex: " + s

                c = ''.join('{:02x}'.format(x) for x in bytearray(outend))
                #If EOT found, finish communication
                if c == "04":
                    if self.debug:
                        print "EOT Found"
                    break

        #Return payload without padding
        return out[3:].replace("\x1A", "")


#Open device
r = RIONNA27("/dev/tty.usbserial", 0)

# Configure date and time
r.setParam("CLK " + datetime.now().strftime("%y %m %d %H %M %S"))

# Configure SPL mode
r.setParam("IMD 0")

# Configure statistical calculation:
# L10, L90, Lmax and Leq
r.setParam("MSR 1")
r.setParam("STT 0 0 1 0 1 0 0 1 0 1")

# Configure calculation time to 5 min
r.setParam("SMT 5 1")

# Configure group mode
r.setParam("AUT 1 # #")

while true:
    # Start 5 minutes measurement
    r.setParam("SRT 1")

    # Wait until SRT command says: measurement finished
    print "Measuring..."
    while r.askParam("SRT").split(",")[1][0] == "1":
        time.sleep(1)

    # Retrieve the results of the measurement
    ans = r.askParam("DOD 0")

    # Display the results
    print "Waveform peak: " + ans.split(",")[9] + " dB"
    print "LAmax: " + ans.split(",")[10] + " dB"
    print "LAeq: " + ans.split(",")[12] + " dB"
    print "LA10 " + ans.split(",")[14] + " dB"
    print "LA90: " + ans.split(",")[16] + " dB"
