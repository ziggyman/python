import logging as log

minutes_per_degree = 60 # To avoid magic numbers in the code. It never changes!
steps_per_minute = 31.1

class PLCState():
    state = {}
    def __init__(self):
        self.set_state()

    def set_state(self,
                  GMCentred=1,
                  GMInbeam=0,
                  GMMoving=0,
                  GMFailure=0,
                  FilterInit=0,
                  FilterCentred=0,
                  FilterMoving=0,
                  FilterFailure=0,
                  ARCMirror=0,
                  ARC1=0,
                  ARC2=0,
                  SlitShutter=0,
                  SlitIllumination=0,
                  RearOfSlitMirror=0,
                  HartmanA=0,
                  HartmanB=0,

                  SlitWidthInitPos=0,
                  SlitWidthAtLimits=0,
                  SlitWidthMoving=0,
                  SlitWidthFailure=0,
                  GratingAngleInit=0,
                  GratingAngleLimit1=0,
                  GratingAngleLimit2=0,
                  GratingAngleMoving=0,
                  GratingAngleFailure=0,
                  CameraFocusInit=0,
                  CameraFocusLimit1=0,                  
                  CameraFocusLimit2=0,
                  CameraFocusMoving=0,
                  CameraFocusFailure=0,
                  SlitWidthInitReq=0,
                  AngleInitReq=0,

                  FilterwheelPosition=1,
                  GratingID=12345,
                  GratingInserted=0,
                  GratingHatchClosed=0,
                  SlitShutterFailure=0,
                  ARCMirrorFailure=0,
                  RoSMirrorFailure=0,
                  HartmanFailure=0,

                  GratingAngleSteps=0,

                  FocusPositionPot=0,
                  FocusPosition=0,
                  SlitWidthPosition=0,

                  TopCrateInterlock=0,
                  FilterInterlock=0,
                  PneumaticsInterlock=0,
                  ARCInterlock=0,
                  SlitWidthInterlock=0,
                  BottomSignalInterlock=0,
                  BottomDriveInterlock=0,

                  SlitIlluminationValue=0,

                  SlitWidthErrorCodes=0,
                  GratingAngleErrorCodes=0,

                  FocusError=0,
                  FocusMotorOn=0,
                  FocusMoving=0,
                  FocusReferencing=0,
                  FocusAtPosition=0,
                  FocusNegativeLimit=0,
                  FocusReference=0,
                  FocusPositiveLimit=0,

                  FocusIn1=0,
                  FocusIn2=0,
                  FocusIn3=0,
                  FocusIn4=0
):
        self.state["GMCentred"] = GMCentred
        self.state["GMInbeam"] = GMInbeam
        self.state["GMMoving"] = GMMoving
        self.state["GMFailure"] = GMFailure
        self.state["FilterInit"] = FilterInit
        self.state["FilterCentred"] = FilterCentred
        self.state["FilterMoving"] = FilterMoving
        self.state["FilterFailure"] = FilterFailure
        self.state["ARCMirror"] = ARCMirror
        self.state["ARC1"] = ARC1
        self.state["ARC2"] = ARC2
        self.state["SlitShutter"] = SlitShutter
        self.state["SlitIllumination"] = SlitIllumination
        self.state["RearOfSlitMirror"] = RearOfSlitMirror
        self.state["HartmanA"] = HartmanA
        self.state["HartmanB"] = HartmanB

        self.state["SlitWidthInitPos"] = SlitWidthInitPos
        self.state["SlitWidthAtLimits"] = SlitWidthAtLimits
        self.state["SlitWidthMoving"] = SlitWidthMoving
        self.state["SlitWidthFailure"] = SlitWidthFailure
        self.state["GratingAngleInit"] = GratingAngleInit
        self.state["GratingAngleLimit1"] = GratingAngleLimit1
        self.state["GratingAngleLimit2"] = GratingAngleLimit2
        self.state["GratingAngleMoving"] = GratingAngleMoving
        self.state["GratingAngleFailure"] = GratingAngleFailure
        self.state["CameraFocusInit"] = CameraFocusInit
        self.state["CameraFocusLimit1"] = CameraFocusLimit1                  
        self.state["CameraFocusLimit2"] = CameraFocusLimit2
        self.state["CameraFocusMoving"] = CameraFocusMoving
        self.state["CameraFocusFailure"] = CameraFocusFailure
        self.state["SlitWidthInitReq"] = SlitWidthInitReq
        self.state["AngleInitReq"] = AngleInitReq

        self.state["FilterwheelPosition"] = FilterwheelPosition
        self.state["GratingID"] = GratingID
        self.state["GratingInserted"] = GratingInserted
        self.state["GratingHatchClosed"] = GratingHatchClosed
        self.state["SlitShutterFailure"] = SlitShutterFailure
        self.state["ARCMirrorFailure"] = ARCMirrorFailure
        self.state["RoSMirrorFailure"] = RoSMirrorFailure
        self.state["HartmanFailure"] = HartmanFailure

        self.state["GratingAngleSteps"] = GratingAngleSteps
        self.state["GratingAngle"] = GratingAngleSteps / (minutes_per_degree * steps_per_minute)

        self.state["FocusPositionPot"] = FocusPositionPot

        self.state["FocusPosition"] = FocusPosition
        self.state["SlitWidthPosition"] = SlitWidthPosition

        self.state["TopCrateInterlock"] = TopCrateInterlock
        self.state["FilterInterlock"] = FilterInterlock
        self.state["PneumaticsInterlock"] = PneumaticsInterlock
        self.state["ARCInterlock"] = ARCInterlock
        self.state["SlitWidthInterlock"] = SlitWidthInterlock
        self.state["BottomSignalInterlock"] = BottomSignalInterlock
        self.state["BottomDriveInterlock"] =  BottomDriveInterlock

        self.state["SlitIlluminationValue"] = SlitIlluminationValue

        self.state["SlitWidthErrorCodes"] = SlitWidthErrorCodes
        self.state["GratingAngleErrorCodes"] = GratingAngleErrorCodes

        self.state["FocusError"] = FocusError
        self.state["FocusMotorOn"] = FocusMotorOn
        self.state["FocusMoving"] = FocusMoving
        self.state["FocusReferencing"] = FocusReferencing
        self.state["FocusAtPosition"] = FocusAtPosition
        self.state["FocusNegativeLimit"] = FocusNegativeLimit
        self.state["FocusReference"] = FocusReference
        self.state["FocusPositiveLimit"] = FocusPositiveLimit

        self.state["FocusIn1"] = FocusIn1
        self.state["FocusIn2"] = FocusIn2
        self.state["FocusIn3"] = FocusIn3
        self.state["FocusIn4"] = FocusIn4

    def get_state(self):
        return self.state
