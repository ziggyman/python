"""This is the top-level controller. It uses the Composite pattern to
contain controllers for each of the instrument sub-systems.
"""

import logging as log
from cassspectr.controller.detector_controller import DetectorController
from cassspectr.controller.plc_controller import PLCController
from cassspectr.utils.tcs_pull import read_remote_file, get_last_line, get_tcs_info
from astropy.io import fits
from datetime import datetime
from collections import OrderedDict

try:
    from tcsview.client.tcsview_client import TCSViewClient
except ImportError as e:
    log.warn("No such module: {}".format(e))
try:
    from weather.client.weather_client import WeatherClient
except ImportError as e: 
    log.warn("No such module: {}".format(e))

class Controller:
    def __init__(self, 
                 plc_host="localhost",
                 plc_port=9090, 
                 detector_host="localhost",
                 detector_port=9091,
                 tcsview_host="localhost",
                 tcsview_port=9093,
                 weather_host="localhost",
                 weather_port=9094):
        # self.plcc = PLCController(plc_host, plc_port)
        # self.dc = DetectorController(tcp_host=detector_host, tcp_port=detector_port)
        self.fits_info = OrderedDict({})
        self.plc_host = plc_host
        self.plc_port = plc_port
        self.detector_host = detector_host
        self.detector_port = detector_port
        self.tcsview_host=tcsview_host
        self.tcsview_port=tcsview_port
        self.weather_host = weather_host
        self.weather_port = weather_port

    @property
    def plcc(self):
        if not getattr(self, '_plcc', False):
            self._plcc = PLCController(
                tcp_host=self.plc_host, tcp_port=self.plc_port)
        return self._plcc

    @property
    def dc(self):
        if not getattr(self, '_dc', False):
            self._dc = DetectorController(
                tcp_host=self.detector_host, tcp_port=self.detector_port)
        return self._dc

    @property
    def tcsvc(self):
        if not getattr(self, '_tcsvc', False):
            self._tcsvc = TCSViewClient(
                tcp_host=self.tcsview_host, tcp_port=self.tcsview_port)
        return self._tcsvc

    @property
    def weathc(self):
        if not getattr(self, '_weathc', False):
            self._weathc = WeatherClient(
                tcp_host=self.weather_host, tcp_port=self.weather_port)
        return self._weathc

    def get_tcs_info(self, fits_info):
        # Get the TCS information 
        fits_info["TELRA"] = ("NA", "The telescope right ascension")
        fits_info["TELDEC"] = ("NA", "The telescope declination")
        fits_info["AIRMASS"] = ("NA", "The airmass (sec(z))")
        fits_info["ZD"] = ("NA", "The telescope zenith distance")
        fits_info["HA"] = ("NA", "The telescope hour angle")
        fits_info["TELFOCUS"] = ("NA", "The telescope focus")
        fits_info["INSTANGL"] = ("NA", "The instrument angle")
        fits_info["DOMEPOS"] = ("NA", "The dome position")
        try:
            tcs = self.tcsvc.get_tcs_info()
            if len(tcs) > 0:
                fits_info["TELRA"] = (tcs["TELRA"], "The telescope right ascension")
                fits_info["TELDEC"] = (tcs["TELDEC"], "The telescope declination")
                fits_info["AIRMASS"] = (tcs["AIRMASS"], "The airmass (sec(z))")
                fits_info["ZD"] = (tcs["ZD"], "The telescope zenith distance")
                fits_info["HA"] = (tcs["HA"], "The telescope hour angle")
                fits_info["TELFOCUS"] = (tcs["TELFOCUS"], "The telescope focus")
                fits_info["INSTANGL"] = (tcs["INSTANGL"], "The instrument angle")
                fits_info["DOMEPOS"] = (tcs["DOMEPOS"], "The dome position")
        except Exception as e:
            log.error("Unable to get TCS info: {}".format(str(e)))
        
    def get_weather_info(self, fits_info):
        # Get the weather and seeing information
        fits_info["TMTDEW"] = ("NA", "Difference between current and dew-point temps")
        fits_info["HUMIDIT"] = ("NA", "The average humidity")
        fits_info["RELSKYT"] = ("NA", "The relative sky temperature (avg cloud cover)")
        fits_info["WIND"] = ("NA", "The average wind speed")
        fits_info["ENVTEM"] = ("NA", "The average ambient temperature")
        fits_info["SEEING"] = ("NA", "The curent seeing")
        try:
            weather_info = self.weathc.get_weather_info()
            if len(weather_info) > 0:
                fits_info["TMTDEW"] = (weather_info["avg_t_min_tdew"], "Difference between current and dew-point temps")
                fits_info["HUMIDIT"] = (weather_info["avg_hum"], "The average humidity")
                fits_info["RELSKYT"] = (weather_info["avg_cloud"], "The relative sky temperature (avg cloud cover)")
                fits_info["WIND"] = (weather_info["avg_wind"], "The average wind speed")
                fits_info["ENVTEM"] = (weather_info["avg_temp"], "The average ambient temperature")
                fits_info["SEEING"] = (weather_info["seeing"], "The curent seeing")
        except Exception as e:
            log.warn("Unable to get weather and seeing information: {}".format(e))

    def add_supp_info(self, supp_fits_info):
        if not supp_fits_info:
            return
        for it in supp_fits_info.items():
            k = it[0]
            v = it[1]
            # print("key = {}  val = {}".format(k, v))
            self.fits_info[k] = v

    def start_exposure(self, exp_time, 
                       rows, cols, 
                       rowbin, colbin, 
                       rowcen, colcen, 
                       supp_fits_info = {}):
        """Start an exposure.
        
        Calls the detector driver start_exposure method, and then
        ensures that the FITS information is up to date.
        """
        self.dc.start_exposure(exp_time, rows, cols, rowbin, colbin, rowcen, colcen)
        self.add_supp_info(supp_fits_info)
        self.dc.add_detector_info(self.fits_info)
        # self.plcc.add_plc_info(fits_info)
        self.get_tcs_info(self.fits_info)
        self.get_weather_info(self.fits_info)
 
    def gen_fits(self, arr, filename):

        """Generate a simple fits file.

        Only the basic data is included (no extra keywords).
        """
        hdu = fits.PrimaryHDU()
        hdu.data = arr
        hdu.writeto(filename, clobber=True)

    def gen_fits_with_keys(self, arr, filename):
        """Generate a complete fits file.

        The keys parameter is a dictionary with keywords and
        values. These are included in the header.
        """
        hdu = fits.PrimaryHDU()
        hdu.data = arr
#        print("Fits info: {}".format(self.fits_info))
        for key, val in self.fits_info.items():
            print("{} = {}".format(key, val))
            hdu.header[key] = val
        hdu.writeto(filename, clobber=True)

