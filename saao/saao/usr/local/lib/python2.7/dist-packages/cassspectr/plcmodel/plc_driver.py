"""This class is a driver for the Cassegrain Spectrograph PLC."""

import time
import serial
import logging as log
from cassspectr.plcmodel.plc_state import PLCState
from cassspectr.interface.spectrograph_interface.PLC import PLCThriftException

class PLCException(PLCThriftException):
    pass
class PLCTimeoutException(PLCThriftException):
    pass

class PLCInvalidResultException(PLCThriftException):
    pass

# The end codes which may be returned by the PLC
end_codes = {"00": "Normal completion",
             "01": "Not executable in RUN mode",
             "02": "Not executable in Monitor",
             "04": "Address over",
             "0B": "Not executable in Program mode",
             "13": "FCS error",
             "14": "Format error",
             "15": "Entry number data error",
             "16": "Command not supported",
             "18": "Frame length error",
             "19": "Not executable",
             "23": "Memory write protected",
             "A3": "FCS error in transmit data",
             "A4": "Format error in transmit data",
             "A5": "Data error in transmit data",
             "A8": "Frame length error in transmit data"
}

# The positions of the control bits. The position of each bit is
# given as a word position, and then a bit position in that word
cbitpos = {"GMCentred":            (0, 0),
           "GMInbeam":             (0, 1),

           "FWInit":               (0, 4),
           "FWMove":               (0, 5),
           "FWReset":              (0, 6),

           "ArcMirrorIB":          (0, 8),
           "Arc1On":               (0, 9),
           "Arc2On":               (0,10),
           "SlitShutter":          (0,11),
           "SlitIllumination":     (0, 12),
           "RearOfSlitMirror":     (0, 13),
           "HartmanA":             (0, 14),
           "HartmanB":             (0, 15),

           "SlitInit":             (1, 0),
           "SlitMove":             (1, 1),
           "SlitReset":            (1, 2),

           "GratingAngleInit":     (1, 4),
           "GratingAngleMoveAbs":  (1, 5),
           "GratingAngleMoveRel":  (1, 6),
           "GratingAngleReset":    (1, 7),

           "CameraFocusInit":      (1, 9),
           "CameraFocusMoveAbs":   (1, 10),
           "CameraFocusMoveRel":   (1, 11),
           "CameraFocusReset":     (1, 12),

           "SlitIlluminationBit0": (2, 0),
           "SlitIlluminationBit1": (2, 1),
           "SlitIlluminationBit2": (2, 2),

           "FilterwheelPosBit0":   (2, 4),
           "FilterwheelPosBit1":   (2, 5),
           "FilterwheelPosBit2":   (2, 6),
           "FilterwheelPosBit3":   (2, 7),

           "SlitWidthValue":       (3, 0),

           "GratingAngleLo":       (4, 0),

           "GratingAngleHi":       (5, 0),

           "CameraFocusLo":        (6, 0),

           "CameraFocusHi":        (7, 0)
}

# Constants for various PLC settings
min_filter_pos = 1
max_filter_pos = 8

min_slit_illumination = 0
max_slit_illumination = 7

min_slit_width = 0
max_slit_width = 28

min_grating_angle = -13
max_grating_angle = 25
minutes_per_degree = 60 # To avoid magic numbers in the code. It never changes!
steps_per_minute = 31.1

min_camera_focus = 2
max_camera_focus = 7

class PLCDriver(object):
    pc2plc = []

    def __init__(self,
                 port='/dev/ttyUSB0',
                 baudrate=9600,
                 parity=serial.PARITY_EVEN,
                 stopbits=serial.STOPBITS_ONE,
                 bytesize=serial.SEVENBITS,
                 timeout=5,
                 test=False):
        self.port = port
        self.baudrate = baudrate
        self.parity = parity
        self.stopbits = stopbits
        self.bytesize = bytesize
        self.timeout = timeout
        self.test = test
        self.ser = None

        if self.test == True:
            print("Driver running in test mode - will write commands to file test.txt")
            self.ser = open("test.txt", "w+b")
        else:
            try:
                self.ser = serial.Serial(port, baudrate, parity=parity, stopbits=stopbits, bytesize=bytesize)
                if not self.ser.isOpen():
                    self.ser.open()
            except ValueError as ve:
                log.error("Value error opening serial port {}".format(port))
                raise PLCException(str(ve))
            except serial.SerialException as se:
                log.error("Serial device could not be found or could not be configured")
                raise PLCException(str(se))
        self.init_pc2plc()
        self.state = PLCState()

    def reconnect_serial(self):
        if self.ser is not None:
            if not self.ser.isOpen():
                self.ser.open()
            else:
                self.ser.close()
                self.ser.open()
            

    def init_pc2plc(self):
        """Initialize the structure sent when a command action is performed."""
        # The structure sent to the PLC is 10 words long, and each
        # word is 16 bits. A 16 bit word will be represented as 4 hex
        # digits. Thus we represent the structure as an array of 40
        # hex characters.
        numwords = 10
        self.pc2plc = ["0000" for i in range(numwords)]

    def set_bit_on(self, pos):
        wordnum = pos[0]
        bitnum = pos[1]
        mask = 1 << bitnum
        # log.debug("Setting bit on at {} {}, using {:016b}".format(wordnum, bitnum, mask))
        word = int(self.pc2plc[wordnum], 16)
        self.pc2plc[wordnum] = "{:04X}".format(word | mask)
        # log.debug("After setting bit {} on: pc2plc[{}] = {:016b} (bin)  {:04X} (hex)".format(bitnum, wordnum, int(self.pc2plc[wordnum], 16), int(self.pc2plc[wordnum], 16)))

    def set_bit_off(self, pos):
        wordnum = pos[0]
        bitnum = pos[1]
        mask = 0xffff - (1 << bitnum)
        word = int(self.pc2plc[wordnum], 16)
        self.pc2plc[wordnum] = "{:04X}".format(word & mask)
        # self.pc2plc[wordnum] &= mask

    def set_raw_value(self, pos, num_nibbles, value):
        """Write the given value directly to the word, starting at the specified position.

        A word is four nibbles long.
        The parameter pos is a tuple of (wordpos, startpos), where
        startpos is the position of the lowest bit to change.
        num_nibbles is a value between 1 and 4
        value is the value to write.

        """
        if len(str(value)) > num_nibbles:
            raise PLCException("Can't fit {} into {} characters".format(value, num_nibbles))
        wordpos = pos[0]
        startpos = pos[1]
        # Get the first nibble to change (lower nibbles may have
        # existing data we don't want to change)
        # These refer to the position from the right
        startnibble = int(startpos/4)
        endnibble = startnibble + num_nibbles - 1
        # Get the current word
        word = "{}".format(self.pc2plc[wordpos])
        # log.debug("Current word = {}".format(word))
        val = "{:0{num_nibbles}d}".format(value, num_nibbles=num_nibbles)

        preword = word[0:4-endnibble-1]
        # log.debug("preword = {}".format(preword))
        postword = word[4-startnibble:]
        # log.debug("postword = {}".format(postword))
        newword = preword + val + postword

        # log.debug("Newword = {}".format(newword.rjust(4, '0')))
        
        self.pc2plc[wordpos] = newword.rjust(4, '0')

    def set_value(self, pos, numbits, value):
        """Set the bits for the value in the specified position in the word.
        The low bit of value is placed at startpos, and the high bit
        no further than startpos + numbits - 1.
        """
        wordpos = pos[0]
        bitpos = pos[1]
        # First check that value will fit in numbits bits
        n, v = numbits, value
        while int(v) > 0:
            v /= 2
            n -= 1
            if n < 0:
                raise PLCException("Can't express {} in {} bits".format(value, numbits))
        # Now ensure that whatever value was previously in this location is zeroed out
        # log.debug("Word (orig) = {:016b}".format(self.pc2plc[wordpos]))
        mask = 0xffff - (2 ** numbits - 1) << bitpos
        self.pc2plc[wordpos] &= mask
        # log.debug("Word (masked) = {:016b}".format(self.pc2plc[wordpos]))
        # Finally, set the bits for the value, starting at bitpos
        self.pc2plc[wordpos] |= (value << bitpos)
        # log.debug("Word (final) = {:016b}".format(self.pc2plc[wordpos]))
        # log.debug("Word (final hex) = {:04X}".format(self.pc2plc[wordpos]))

    def reset_transients(self):
        # log.debug("Resetting transients")        
        self.set_bit_off(cbitpos["GMCentred"])
        self.set_bit_off(cbitpos["GMInbeam"])
        self.set_bit_off(cbitpos["FWInit"])
        self.set_bit_off(cbitpos["FWMove"])
        self.set_bit_off(cbitpos["FWReset"])
        self.set_bit_off(cbitpos["SlitInit"])
        self.set_bit_off(cbitpos["SlitMove"])
        self.set_bit_off(cbitpos["SlitReset"])
        self.set_bit_off(cbitpos["GratingAngleInit"])
        self.set_bit_off(cbitpos["GratingAngleMoveAbs"])
        self.set_bit_off(cbitpos["GratingAngleMoveRel"])
        self.set_bit_off(cbitpos["GratingAngleReset"])
        self.set_bit_off(cbitpos["CameraFocusInit"])
        self.set_bit_off(cbitpos["CameraFocusMoveAbs"])
        self.set_bit_off(cbitpos["CameraFocusMoveRel"])
        self.set_bit_off(cbitpos["CameraFocusReset"])
        self.send_pc2plc()

    def send_pc2plc(self):
        """Generate the command captured in pc2plc, and send it."""
        CMD  = "@00WD0100"
        for word in self.pc2plc:
            CMD += word.zfill(4)

        CMD += self.calculate_fcs(CMD)
        CMD += "*\r"

        # log.debug("Sending command {}".format(CMD))
        # log.debug("Sending command {}".format(self.pretty_print_command(CMD)))
        self.ser.write(bytearray(CMD, 'ascii'))
        time.sleep(0.1)
        if self.test == False:
            result = self.read_back()
            # log.debug("Received result {}".format(result))
            self.validate_result(result)

    def read_back(self):
        """Read the response from the PLC."""
        result = ""
        count = 0
        n = -1
        while True:
            c = self.ser.read(1).decode("utf-8")
            if len(c) < 1:
                # If read returned with no characters, we've hit a time-out.
                raise PLCTimeoutException("Timeout reading from serial port")
            if c >= '0':
                result += c
            if c =='\r':
                break
        return result

    def validate_result(self, result):
        """Check that the result is valid, and that no exception code is returned."""
        # Check start character
        if len(result) < 1:
            log.debug("Nothing to validate (returned string has zero length)")
            return
        
        # Check that the FCS returned is correct
        # # Strip the trailing *CR
        # result = result[:-2:]
        # The FCS is the last two characters
        fcs = result[-2:]
        # The payload is the result without the FCS
        payload = result[:-2]
        fcs_calc = self.calculate_fcs(payload)
        if fcs.lower() != fcs_calc.lower():
            log.error("Returned FCS [{}] does not match calculated FCS [{}]".
                      format(fcs, fcs_calc))
            raise PLCInvalidResultException(
                "Returned FCS [{}] does not match calculated FCS [{}]".
                format(fcs, fcs_calc))                

        # The payload should begin with '@'
        if payload[0] != '@':
            log.error("Result returned does not begin with '@'")
            raise PLCInvalidResultException("Result returned does not begin with '@'")

        # Check if the return string specifies an exception code
        end_code = payload[5:7]
        # Check if an exception code is being returned
        if end_code == "00":
            # Normal completion
            return
        else:
            end_code = int(end_code, 16)
            log.error("End code error {}: {}".format(end_code, end_codes[end_code]))
            raise PLCInvalidResultException(
                "End code error {}: {}".format(end_code, end_codes[end_code]))
            

    def get_status(self):
        log.debug("Getting status. Update state...")
        if not self.test:
            self.update_status()
        state = self.state.get_state()
        log.debug("Got state {}".format(state))
        return state

    def update_status(self):
        """Fetch the status, and update the state object storing this.

        A status request is sent to the PLC. The string sent indicates
        the range of memory locations to be read. The PLC replies with
        the contents of these locations.
        """ 
        # First , construct the string to retrieve the status.
        GET_STAT  = "@"     # Begin with a "@" character
        GET_STAT += "00"    # Node number
        GET_STAT += "RD"    # Request SSC command
        GET_STAT += "0110"  # The starting address of the locations to read
        GET_STAT += "0013"  # The number of words to read
        GET_STAT += self.calculate_fcs(GET_STAT)
        GET_STAT += "*\r"

        # Write the request string to the PLC and fetch the result.
        # Note, we do this twice, as sometimes the status string is truncated.
        log.debug("Writing status request string {}".format(GET_STAT))

        self.ser.write(bytearray(GET_STAT, 'ascii'))
        log.debug("Reading back...")
        result = self.read_back()
        # log.debug("Status request returned {}".format(result))
        log.debug("Status request returned {}".format(self.pretty_print_status(result)))
        
        # Check that the response is valid.
        self.validate_result(result)

        # Update the stored state
        self.process_status(result)

    def process_status(self, result):
        """Update the state of the PLC status."""
        data = result[7:-4]
        self.state.set_state(GMCentred        = int(data[0:4], 16) & 1,
                             GMInbeam         = int(data[0:4], 16) >> 1 & 1,
                             GMMoving         = int(data[0:4], 16) >> 2 & 1,
                             GMFailure        = int(data[0:4], 16) >> 3 & 1,
                             FilterInit       = int(data[0:4], 16) >> 4 & 1,
                             FilterCentred    = int(data[0:4], 16) >> 5 & 1,
                             FilterMoving     = int(data[0:4], 16) >> 6 & 1,
                             FilterFailure    = int(data[0:4], 16) >> 7 & 1,
                             ARCMirror        = int(data[0:4], 16) >> 8 & 1,
                             ARC1             = int(data[0:4], 16) >> 9 & 1,
                             ARC2             = int(data[0:4], 16) >> 10 & 1,
                             SlitShutter      = int(data[0:4], 16) >> 11 & 1,
                             SlitIllumination = int(data[0:4], 16) >> 12 & 1,
                             RearOfSlitMirror = int(data[0:4], 16) >> 13 & 1,
                             HartmanA         = int(data[0:4], 16) >> 14 & 1,
                             HartmanB         = int(data[0:4], 16) >> 15 & 1,

                             SlitWidthInitPos    = int(data[4:8], 16) >> 0 & 1,
                             SlitWidthAtLimits   = int(data[4:8], 16) >> 1 & 1,
                             SlitWidthMoving     = int(data[4:8], 16) >> 2 & 1,
                             SlitWidthFailure    = int(data[4:8], 16) >> 3 & 1,
                             GratingAngleInit    = int(data[4:8], 16) >> 4 & 1,
                             GratingAngleLimit1  = int(data[4:8], 16) >> 5 & 1,
                             GratingAngleLimit2  = int(data[4:8], 16) >> 6 & 1,
                             GratingAngleMoving  = int(data[4:8], 16) >> 7 & 1,
                             GratingAngleFailure = int(data[4:8], 16) >> 8 & 1,
                             CameraFocusInit     = int(data[4:8], 16) >> 9 & 1,
                             CameraFocusLimit1   = int(data[4:8], 16) >> 10 & 1,
                             CameraFocusLimit2   = int(data[4:8], 16) >> 11 & 1,
                             CameraFocusMoving   = int(data[4:8], 16) >> 12 & 1,
                             CameraFocusFailure  = int(data[4:8], 16) >> 13 & 1,
                             SlitWidthInitReq    = int(data[4:8], 16) >> 14 & 1,
                             AngleInitReq        = int(data[4:8], 16) >> 15 & 1,

                             FilterwheelPosition = int(data[11],16),
                             GratingID           = int(data[9:11], 16) & 0b11111, 

                             GratingInserted     = int(data[8:12], 16) >> 9 & 1,
                             GratingHatchClosed  = int(data[8:12], 16) >> 10 & 1,
                             SlitShutterFailure  = int(data[8:12], 16) >> 12 & 1,
                             ARCMirrorFailure    = int(data[8:12], 16) >> 13 & 1,
                             RoSMirrorFailure    = int(data[8:12], 16) >> 14 & 1,
                             HartmanFailure      = int(data[8:12], 16) >> 15 & 1,

                             # The following results from:
                             # Interpret the two words (8 hex nibbles) as two numbers directly.
                             # The number is sign * (10000 * high word + low word )
                             # The MSB of the two words gives the sign.
                             # So if data[19] = 80 (=1000) => (int(data[19]) - 4)/(-1/4) = -1
                             # So if data[19] = 00 (=0000) => (int(data[19]) - 4)/(-1/4) = +1
                             GratingAngleSteps   = ((int(data[12:16]) + 10000 * int(data[19:20]))
                                                  * (int(data[16]) - 4)*(-1/4.0)),
                             
                             FocusPositionPot    = int(data[20:24]),

                             FocusPosition       = int(data[24:28])/10000.0 + int(data[31]),

                             SlitWidthPosition   = int(data[32:36]),

                             TopCrateInterlock     = int(data[36:40], 16) >> 0 & 1,
                             FilterInterlock       = int(data[36:40], 16) >> 1 & 1,
                             PneumaticsInterlock   = int(data[36:40], 16) >> 2 & 1,
                             ARCInterlock          = int(data[36:40], 16) >> 3 & 1,
                             SlitWidthInterlock    = int(data[36:40], 16) >> 4 & 1,
                             BottomSignalInterlock = int(data[36:40], 16) >> 5 & 1,
                             BottomDriveInterlock  = int(data[36:40], 16) >> 6 & 1,

                             SlitIlluminationValue = int(data[36]),

                             SlitWidthErrorCodes   = int(data[40:44]),

                             GratingAngleErrorCodes = int(data[44:48]),

                             FocusError            = int(data[48:52], 16) >> 0 & 1,
                             
                             FocusMotorOn          = int(data[48:52], 16) >> 4 & 1,
                             FocusMoving           = int(data[48:52], 16) >> 5 & 1,
                             FocusReferencing      = int(data[48:52], 16) >> 6 & 1,
                             FocusAtPosition       = int(data[48:52], 16) >> 7 & 1,
                             FocusNegativeLimit    = int(data[48:52], 16) >> 8 & 1,
                             FocusReference        = int(data[48:52], 16) >> 9 & 1,
                             FocusPositiveLimit    = int(data[48:52], 16) >> 10 & 1,

                             FocusIn1              = int(data[48:52], 16) >> 12 & 1,
                             FocusIn2              = int(data[48:52], 16) >> 13 & 1,
                             FocusIn3              = int(data[48:52], 16) >> 14 & 1,
                             FocusIn4              = int(data[48:52], 16) >> 15 & 1
                         )

    def get_errors(self):
        if not self.test:
            self.update_status()
        state = self.state.get_state()
        failure_states = ["GMFailure", 
                          "FilterFailure", 
                          "SlitWidthFailure", 
                          "GratingAngleFailure",
                          "CameraFocusFailure",
                          "SlitShutterFailure",
                          "ARCMirrorFailure",
                          "RoSMirrorFailure",
                          "HartmanFailure"]
        failed_states = []
        for failure_code in failure_states:
            if state[failure_code] == 1:
                failed_states.append(failure_code)
        return failed_states

    def calculate_fcs(self, payload):
        """Calculate the LRC of the string.

        The Frame Check Sequence is calculated as follows: For each
        byte of the payload, interpret this as an ascii character. For
        each ascii character, get the two-character hex
        representation. the FCS is the XOR of the leftmost of these
        characters, followed by the XOR of the rightmost of these
        characters.

        """
        left = 0
        right = 0
        for byte in payload:
            hexcode = hex(ord(byte))
            left ^= int(hexcode[2],16)
            right ^= int(hexcode[3],16)
        fcs = "{:X}{:X}".format(left, right)
        return fcs

    def set_gm_centred(self):
        # Set GM Centred bit on, GM inbeam bit off
        log.info("Setting guide mirror to centred")
        self.reset_transients()
        self.set_bit_on(cbitpos["GMCentred"])
        self.send_pc2plc()
        if self.test:
            self.state.state["GMCentred"] = 1
            self.state.state["GMInbeam"] = 0
            self.state.state["GMMoving"] = 0

    def set_gm_inbeam(self):
        # Set GM Centred bit off, GM inbeam bit on
        log.info("Setting guide mirror to in-beam")
        self.reset_transients()
        self.set_bit_on(cbitpos["GMInbeam"])
        self.send_pc2plc()
        if self.test:
            self.state.state["GMCentred"] = 0
            self.state.state["GMInbeam"] = 1
            self.state.state["GMMoving"] = 0

    def set_fw_init(self):
        #Set FWInit bit on, FWMove and FWReset bits off
        log.info("Initializing filterwheel")
        self.reset_transients()
        self.set_bit_on(cbitpos["FWInit"])
        self.set_bit_off(cbitpos["FWMove"])
        self.set_bit_off(cbitpos["FWReset"])
        self.send_pc2plc()
        if self.test:
            self.state.state["FilterInit"] = 1
            self.state.state["FilterCentred"] = 1
            self.state.state["FilterMoving"] = 0
            self.state.state["FilterwheelPosition"] = -1

    def set_fw_reset(self):
        #Set FWReset bit on, FWMove and FWInit bits off
        log.info("Resetting filterwheel")
        self.reset_transients()
        self.set_bit_on(cbitpos["FWReset"])
        self.send_pc2plc()
        if self.test:
            self.state.state["FilterInit"] = 0
            self.state.state["FilterCentred"] = 0
            self.state.state["FilterMoving"] = 0
            self.state.state["FilterwheelPosition"] = -1

    def set_fw_move(self, pos):
        log.info("Moving filterwheel to {}".format(pos))
        if pos < min_filter_pos or pos > max_filter_pos:
            raise PLCException("The filter position entered ({}) "
                               "is outside the allowed range ({} - {})"
                               .format(pos, min_filter_pos, max_filter_pos))
        # Set FWMove bit on, FWReset and FWInit bits off
        self.reset_transients()
        self.set_bit_on(cbitpos["FWMove"])
        # Set bits for move position
        self.set_raw_value(cbitpos["FilterwheelPosBit0"], 1, pos)
        self.send_pc2plc()
        if self.test:
            self.state.state["FilterInit"] = 1
            self.state.state["FilterCentred"] = 1
            self.state.state["FilterMoving"] = 0
            self.state.state["FilterwheelPosition"] = pos

    def set_arc_mirror(self, in_beam = True):
        log.info("Setting arc mirror to {}".format(in_beam))
        if in_beam:
            self.set_bit_on(cbitpos["ArcMirrorIB"])
            if self.test:
                self.state.state["ARCMirror"] = 1
        else:
            self.set_bit_off(cbitpos["ArcMirrorIB"])
            if self.test:
                self.state.state["ARCMirror"] = 0
        self.send_pc2plc()

    def set_arc_lamp_1(self, lamp_on=True):
        log.info("Setting arc lamp 1 to {}".format(lamp_on))
        lampstr = "Arc1On"
        if lamp_on:
            self.set_arc_mirror(in_beam = True)
            self.set_bit_on(cbitpos[lampstr])
            if self.test:
                self.state.state["ARC1"] = 1
        else:
            self.set_bit_off(cbitpos[lampstr])
            if self.test:
                self.state.state["ARC1"] = 0
        self.send_pc2plc()

    def set_arc_lamp_2(self, lamp_on=True):
        log.info("Setting arc lamp 2 to {}".format(lamp_on))
        lampstr = "Arc2On"
        if lamp_on:
            self.set_arc_mirror(in_beam = True)
            self.set_bit_on(cbitpos[lampstr])
            if self.test:
                self.state.state["ARC2"] = 1
        else:
            self.set_bit_off(cbitpos[lampstr])
            if self.test:
                self.state.state["ARC2"] = 0
        self.send_pc2plc()

    def set_slit_init(self):
        log.info("Initializing slit")
        self.reset_transients()
        self.set_bit_on(cbitpos["SlitInit"])
        self.send_pc2plc()
        if self.test:
            self.state.state["SlitWidthInitPos"] = 1
            self.state.state["SlitWidthMoving"] = 0
            self.state.state["SlitWidthPosition"] = 1

    def set_slit_reset(self):
        log.info("Resetting slit")
        self.reset_transients()
        self.set_bit_on(cbitpos["SlitReset"])
        self.send_pc2plc()
        if self.test:
            self.state.state["SlitWidthInitPos"] = 0
            self.state.state["SlitWidthMoving"] = 0
            self.state.state["SlitWidthPosition"] = -1


    def set_slit_width(self, width):
        log.info("Setting slit width to {}".format(width))
        if width < min_slit_width or width > max_slit_width:
            raise PLCException("Invalid slit [{}] width chosen."
                               "Value must be in range [{} - {}]"
                               .format(width, min_slit_width, max_slit_width))
        self.reset_transients()
        self.set_bit_on(cbitpos["SlitMove"])
        self.set_raw_value(cbitpos["SlitWidthValue"], 4, width)
        self.send_pc2plc()
        if self.test:
            self.state.state["SlitWidthPosition"] = width
            self.state.state["SlitWidthInitPos"] = 0
            self.state.state["SlitWidthMoving"] = 0
                         
    def set_slit_shutter(self, open = True):
        log.info("Setting slit shutter to {}".format(open))
        if open:
            self.set_bit_on(cbitpos["SlitShutter"])
            if self.test:
                self.state.state["SlitShutter"] = 1
        else:
            self.set_bit_off(cbitpos["SlitShutter"])
            if self.test:
                self.state.state["SlitShutter"] = 0
        self.send_pc2plc()

    def set_slit_illumination_status(self, illumination_on = True):
        log.info("Setting slit illumination to {}".format(illumination_on))
        if illumination_on == True:
            self.set_bit_on(cbitpos["SlitIllumination"])
            if self.test:
                self.state.state["SlitShutter"] = 0
                self.state.state["SlitIllumination"] = 1
        else:
            self.set_bit_off(cbitpos["SlitIllumination"])
            if self.test:
                self.state.state["SlitIllumination"] = 0
        self.send_pc2plc()

    def set_slit_illumination_value(self, value):
        if value < min_slit_illumination or value > max_slit_illumination:
            raise PLCException("Invalid slit illumination value [{}]. "
                               "Value must be be in range [{} - {}]"
                               .format(value, min_slit_illumination, 
                                       max_slit_illumination))
#        self.set_slit_illumination_status(True)
        self.set_raw_value(cbitpos["SlitIlluminationBit0"], 1, value)
        self.send_pc2plc()

    def set_rear_of_slit_mirror(self, in_beam = True):
        log.info("Setting rear of slit mirror to {}".format(in_beam))
        if in_beam:
            self.set_bit_on(cbitpos["RearOfSlitMirror"])
            if self.test:
                self.state.state["RearOfSlitMirror"] = 1
        else:
            self.set_bit_off(cbitpos["RearOfSlitMirror"])
            if self.test:
                self.state.state["RearOfSlitMirror"] = 0
        self.send_pc2plc()

    def set_grating_angle_init(self):
        log.info("Initializing grating")
        self.reset_transients()
        self.set_bit_on(cbitpos["GratingAngleInit"])
        self.send_pc2plc()
        if self.test:
            self.state.state["GratingAngleInit"] = 1
            self.state.state["GratingAngle"] = 0

    def set_grating_angle_reset(self):
        log.info("Resetting grating")
        self.reset_transients()
        self.set_bit_on(cbitpos["GratingAngleReset"])
        self.send_pc2plc()
        if self.test:
            self.state.state["GratingAngleInit"] = 0
            self.state.state["GratingAngle"] = -1

    def grating_angle_move_abs(self, angle):
        log.info("Absolute move of grating to {}".format(angle))
        if angle < min_grating_angle or angle > max_grating_angle:
            raise PLCException("Invalid grating angle chosen [{}]. "
                               "Must be in range[{} - {}]"
                               .format(angle, min_grating_angle, max_grating_angle))
        self.set_grating_angle_value(angle)
        self.set_grating_angle_abs()

    def grating_angle_move_rel(self, angle):
        log.info("Relative move of grating by {}".format(angle))
        self.update_status()
        current_angle = self.state.state["GratingAngle"]
        log.info("Current angle = {}. Relative move = {}. New angle will be {}".
                  format(current_angle, angle, current_angle + angle))
        if (current_angle + angle < min_grating_angle or 
            current_angle + angle > max_grating_angle):
            raise PLCException("Invalid grating angle chosen [{}]. "
                               "Must be in range[{} - {}]"
                               .format(current_angle + angle, min_grating_angle,max_grating_angle))
        self.set_grating_angle_value(angle)
        self.set_grating_angle_rel()

    def set_grating_angle_abs(self):
        self.reset_transients()
        self.set_bit_on(cbitpos["GratingAngleMoveAbs"])
        self.send_pc2plc()
        
    def set_grating_angle_rel(self):
        self.reset_transients()
        self.set_bit_on(cbitpos["GratingAngleMoveRel"])
        self.send_pc2plc()

    def set_grating_angle_value(self, angle):
        # TODO: Check if we lose precision when setting steps
        # Take an angle, convert to steps
        if angle < min_grating_angle or angle > max_grating_angle:
            raise PLCException("Invalid grating angle chosen [{}]. "
                               "Must be in range[{} - {}]"
                               .format(angle, min_grating_angle, max_grating_angle))
        # Convert the angle to steps
        steps = angle * minutes_per_degree * steps_per_minute
        self.set_grating_angle_steps(steps)

    def set_grating_angle_steps(self, steps):
        # The angle is expressed in steps. It must correspond to an
        # angle in the range [-30 - 30]. 
        min_steps = int(min_grating_angle * minutes_per_degree * steps_per_minute)
        max_steps = int(max_grating_angle * minutes_per_degree * steps_per_minute)
        if steps < min_steps or steps > max_steps:
            raise PLCException("Invalid grating angle steps chosen [{}]. "
                               "Must be in range[{} - {}]"
                               .format(angle, min_steps, max_steps))
        # TODO: do this properly!        
        steps = int(steps)
        if len(str(abs(steps))) <= 4:
            loval = abs(steps)
            hival = 0
        else:
            loval = int(str(steps)[-4:])
            hival = str(int(abs(steps)))
            if len(hival) > 4:
                hival = hival[:-4]
            hival = int(hival)
        self.reset_transients()
        self.set_raw_value(cbitpos["GratingAngleLo"], 4, loval)
        self.set_raw_value(cbitpos["GratingAngleHi"], 1, hival)
        if steps < 0:
            self.set_bit_on((cbitpos["GratingAngleHi"][0], 15))
        else:
            self.set_bit_off((cbitpos["GratingAngleHi"][0], 15))
        # log.debug("Sending  to PLC")
        self.send_pc2plc()
        if self.test:
            self.state.state["GratingAngleInit"] = 0
            self.state.state["GratingAngle"] = steps / (steps_per_minute * minutes_per_degree)

    def set_camera_focus_init(self):
        log.info("Initializing focus")
        self.reset_transients()
        self.set_bit_on(cbitpos["CameraFocusInit"])
        self.send_pc2plc()
        if self.test:
            self.state.state["CameraFocusInit"] = 1

    def set_camera_focus_reset(self):
        log.info("Resetting focus")
        self.reset_transients()
        self.set_bit_on(cbitpos["CameraFocusReset"])
        self.send_pc2plc()
        if self.test:
            self.state.state["CameraFocusInit"] = 0

    def camera_focus_move_abs(self, pos):
        log.info("Moving camera focus to absolute position {}".format(pos))
        if pos < min_camera_focus or pos > max_camera_focus:
            raise PLCException("Invalid camera focus position chosen [{}]. "
                               "Must be in range[{} - {}]"
                               .format(pos, min_camera_focus, max_camera_focus))
        self.set_camera_focus_value(pos)
        self.set_camera_focus_abs()

    def camera_focus_move_rel(self, val):
        log.info("Moving camera focus to relative value {}".format(val))
        self.update_status()
        current_pos = self.state.state["FocusPosition"]
        log.info("Current pos = {}. Relative move = {}. New pos will be {}".
                  format(current_pos, val, current_pos + val))
        if (current_pos + val < min_camera_focus or 
            current_pos + val > max_camera_focus):
            raise PLCException("Invalid camera focus position chosen [{}]. "
                               "Must be in range[{} - {}]"
                               .format(current_pos + val, min_camera_focus, max_camera_focus))
        self.set_camera_focus_value(val)
        self.set_camera_focus_rel()

    def set_camera_focus_abs(self):
        self.reset_transients()
        self.set_bit_on(cbitpos["CameraFocusMoveAbs"])
        self.send_pc2plc()

    def set_camera_focus_rel(self):
        self.reset_transients()
        self.set_bit_on(cbitpos["CameraFocusMoveRel"])
        self.send_pc2plc()

    def set_camera_focus_value(self, val):
        log.info("Setting focus move value to {}".format(val))
        self.reset_transients()
        sign = 1
        if val < 0:
            sign = -1
        intval, fracval = str(float(val)).split('.')
        intval = int(intval)
        fracval = int(fracval.ljust(4, '0'))
        self.set_raw_value(cbitpos["CameraFocusLo"], 4, fracval)
        self.set_raw_value(cbitpos["CameraFocusHi"], 1, abs(intval))
        if sign < 0:
            self.set_bit_on((cbitpos["CameraFocusHi"][0], 15))
        else:
            self.set_bit_off((cbitpos["CameraFocusHi"][0], 15))
        self.send_pc2plc()
        if self.test:
            steps = "{}.{}".format(intval, fracval)
            steps = float(steps)
            if sign == -1:
                # Deal with int = 0 but val negative
                steps = -1 * abs(steps)
            self.state.state["FocusPosition"] = steps

    def set_hartman_a(self, in_beam = True):
        log.info("Set Hartmann A to {}".format(in_beam))
        self.reset_transients()
        if in_beam:
            self.set_bit_on(cbitpos["HartmanA"])
            self.set_bit_off(cbitpos["HartmanB"])
            if self.test:
                self.state.state["HartmanA"] = 1
                self.state.state["HartmanB"] = 0
        else:
            self.set_bit_off(cbitpos["HartmanA"])
            if self.test:
                self.state.state["HartmanA"] = 0
        self.send_pc2plc()
    
    def set_hartman_b(self, in_beam = True):
        log.info("Set Hartmann B to {}".format(in_beam))
        self.reset_transients()
        if in_beam:
            self.set_bit_on(cbitpos["HartmanB"])
            self.set_bit_off(cbitpos["HartmanA"])
            if self.test:
                self.state.state["HartmanA"] = 0
                self.state.state["HartmanB"] = 1
        else:
            self.set_bit_off(cbitpos["HartmanB"])
            if self.test:
                self.state.state["HartmanB"] = 0
        self.send_pc2plc()
    
    def pretty_print_status(self, status):
        pp = "[{} {} {} {}] {} {} {} {} {} {} {} {} {} {} {} {} {} [{}]".format(
            status[0], status[1:3], status[3:5], status[5:7], 
            status[7:11], status[11:15], status[15:19], status[19:23],
            status[23:27], status[27:31], status[31:35], status[35:39],
            status[39:43], status[43:47], status[47:51], status[51:55],
            status[55:59], status[59:61])
        return pp

    def pretty_print_command(self, cmd):
        pp = "[{} {} {} {}] {} {} {} {} {} {} {} {} {} {} [{}]".format(
            cmd[0], cmd[1:3], cmd[3:5], cmd[5:9], 
            cmd[9:13], cmd[13:17], cmd[17:21], cmd[21:25],
            cmd[25:29], cmd[29:33], cmd[33:37], cmd[37:41],
            cmd[41:45], cmd[45:49], cmd[49:51])
        return pp

if __name__ == '__main__':
    logfile = __file__.split('.py')[0] + ".log"
    loglevel = "debug"
    log.basicConfig(filename = logfile, 
                    format='%(asctime)s: %(message)s', 
                    level=getattr(log, loglevel.upper()))
    port = "/dev/ttyUSB0"
    plcd = PLCDriver(port, test=True)
