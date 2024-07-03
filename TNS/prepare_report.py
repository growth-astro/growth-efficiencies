# /usr/bin/env python
#
# Copyright
#
# Distributed under terms of the intermediate Palomar Transient Factory license
#
# Author: Yi Cao <ycao@astro.caltech.edu>
# History
#     - July 10: Create this file
#     - Sep 5: Iair Arcavi - added the internal name fields (prefix,year_format,postfix)
#     - Apr 25 2018: C. Fremling - adapted to ZTF
"""
Report a transient to TNS. Then the TNS returns a name and a flag that
indicates whether it is a new or existing object.
"""

import datetime
import requests
import json
import time
#import api_bulk_report as api_bulk_report


# Sandbox
BASE_URL = "http://sandbox-tns.weizmann.ac.il/api"
API_KEY = "54916f1700966b3bd325fc1189763d86512bda1d"

# Real website
#BASE_URL = "https://wis-tns.weizmann.ac.il/api"
#API_KEY = "54916f1700966b3bd325fc1189763d86512bda1d"
AT_REPORT_FORM = "bulk-report"
AT_REPORT_REPLY = "bulk-report-reply"


class OptionsError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class Options(object):
    def __init__(self, **kwargs):
        _ZTF_start_date = datetime.datetime(2018, 3, 4)
        _now = datetime.datetime.utcnow()
        self.ra = self._getvalue(kwargs, "ra")
        self.ra_error = self._getvalue(kwargs, "ra_error", default_value=0.2)
        self.dec = self._getvalue(kwargs, "dec")
        self.dec_error = self._getvalue(kwargs, "dec_error", default_value=0.2)
        self.groupid = self._getvalue(kwargs, "groupid", default_value=48)
        self.discovery_datetime =\
            self._getvalue(kwargs, "discovery_datetime",
                           default_value=_now, default_type=datetime.datetime)
        self.reporter = self._getvalue(kwargs, "reporter",
                                       default_value="C. Fremling (Caltech) on behalf of the Zwicky Transient Facility (ZTF) collaboration")
        self.at_type = self._getvalue(kwargs, "at_type", default_value=1)
        #self.internal_name_format_prefix =\
        #    self._getvalue(kwargs,"internal_name_format_prefix",
        #                   default_value="ZTF")
        #self.internal_name_format_year_format =\
        #    self._getvalue(kwargs,"internal_name_format_year_format",
        #                   default_value="YY")
        #self.internal_name_format_postfix =\
        #    self._getvalue(kwargs,"internal_name_format_postfix",
        #                   default_value="")
        self.internal_name =\
            self._getvalue(kwargs, "internal_name",
                           default_value="")                   
        self.proprietary_period_groups =\
            self._getvalue(kwargs, "proprietary_period_groups",
                           default_value=[48])
        self.proprietary_period_value =\
            self._getvalue(kwargs, "proprietary_period_value",
                           default_value=0)
        self.proprietary_period_units =\
            self._getvalue(kwargs, "proprietary_period_units",
                           default_value="days")

        self.non_detection_obsdate =\
            self._getvalue(kwargs,
                           "non_detection_obsdate",
                           default_value=_ZTF_start_date)
        self.non_detection_limiting_flux =\
            self._getvalue(kwargs, "non_detection_limiting_flux",
                           default_value=21.5)
        self.non_detection_filter_name =\
            self._getvalue(kwargs, "non_detection_filter_name",
                           default_value="r")
        self.non_detection_exptime =\
            self._getvalue(kwargs, "non_detection_exptime",
                           default_value=30)
        self.non_detection_observer =\
            self._getvalue(kwargs, "non_detection_observer",
                           default_value="ZTF")
        self.non_detection_instrument =\
            self._getvalue(kwargs, "non_detection_instrument",
                           default_value=196)

        self.photometry_obsdate =\
            self._getvalue(kwargs, "photometry_obsdate",
                           default_type=datetime.datetime)
        self.photometry_flux =\
            self._getvalue(kwargs, "photometry_flux")
        self.photometry_flux_error =\
            self._getvalue(kwargs, "photometry_flux_error")
        self.photometry_filter_name =\
            self._getvalue(kwargs, "photometry_filter_name")
        self.photometry_instrument =\
            self._getvalue(kwargs, "photometry_instrument",
                           default_value=196)
        self.photometry_exptime =\
            self._getvalue(kwargs, "photometry_exptime",
                           default_value=30)
        self.photometry_observer =\
            self._getvalue(kwargs, "photometry_observer",
                           default_value="ZTF")
        
        self.photometry2_obsdate =\
            self._getvalue(kwargs, "photometry2_obsdate",
                           default_type=datetime.datetime)
        self.photometry2_flux =\
            self._getvalue(kwargs, "photometry2_flux")
        self.photometry2_flux_error =\
            self._getvalue(kwargs, "photometry2_flux_error")
        self.photometry2_filter_name =\
            self._getvalue(kwargs, "photometry2_filter_name")
        self.photometry2_instrument =\
            self._getvalue(kwargs, "photometry2_instrument",
                           default_value=196)

        self._prepare_report()

        self.base_url =\
            self._getvalue(kwargs, "base_url", BASE_URL)
        self.at_report_form = self._getvalue(kwargs, "bulk_report_form",
                                             AT_REPORT_FORM)
        self.at_report_reply = self._getvalue(kwargs, "bulk_report_reply",
                                              AT_REPORT_REPLY)
        self.api_key = self._getvalue(kwargs, "api_key", API_KEY)

    def _getvalue(self, kwargs, key, default_value=None, default_type=None):
        if default_value is None and key not in kwargs:
            raise OptionsError("%s is required" % key)
        value = kwargs[key] if key in kwargs else default_value
        if default_type is not None and not(value, default_type):
            raise OptionsError("%s needs be an object of class %s" %
                               (key, default_type))
        return value

    def _filter_id(self, filter_name):
        if filter_name == 'r':
            return 111
        elif filter_name == 'g':
            return 110
        elif filter_name == "i":
            return 112
        else:
            error_message =\
                "Filter has to be either r or g: current value " + filter_name
            raise OptionsError(error_message)

    def _prepare_report(self):
        self.proprietary_period =\
            {"proprietary_period_value": "%i" % self.proprietary_period_value,
             "proprietary_period_units": self.proprietary_period_units}
        #self.internal_name_format =\
        #    {"prefix": "ZTF",
        #     "year_format": self.internal_name_format_year_format,
        #     "postfix": self.internal_name_format_postfix}
        if self.non_detection_limiting_flux == 99:
           self.non_detection =\
               {"archiveid": "0",
                "archival_remarks": "Non existent in SDSS/PS1"}
        else:            
           self.non_detection =\
               {"obsdate": self._convert_time(self.non_detection_obsdate),
                "limiting_flux": self.non_detection_limiting_flux,
                "flux_units": "1",
                "filter_value": self._filter_id(self.non_detection_filter_name),
                "instrument_value": self.non_detection_instrument,
                "exptime": self.non_detection_exptime,
                "observer": self.non_detection_observer}
        self.photometry =\
            {"obsdate": self._convert_time(self.photometry_obsdate),
             "flux": self.photometry_flux,
             "flux_err": self.photometry_flux_error,
             "flux_units": "1",
             "filter_value": self._filter_id(self.photometry_filter_name),
             "instrument_value": self.photometry_instrument,
             "exptime": self.photometry_exptime,
             "observer": self.photometry_observer}
        self.photometry2 =\
            {"obsdate": self._convert_time(self.photometry2_obsdate),
             "flux": self.photometry2_flux,
             "flux_err": self.photometry2_flux_error,
             "flux_units": "1",
             "filter_value": self._filter_id(self.photometry2_filter_name),
             "instrument_value": self.photometry_instrument,
             "exptime": self.photometry_exptime,
             "observer": self.photometry_observer}
        self.dictionary =\
            {"at_report": {"0": {"ra": {"value": self.ra,
                                        "error": self.ra_error,
                                        "units": "arcsec"},
                                 "dec": {"value": self.dec,
                                         "error": self.dec_error,
                                         "units": "arcsec"},
                                 "groupid": self.groupid,
                                 "internal_name_format": {"prefix": "ZTF",
                                                          "year_format": "YY",
                                                          "postfix": ""},
                                 "internal_name": self.internal_name,
                                 "reporter": self.reporter,
                                 "discovery_datetime":
                                     self._convert_time(self.discovery_datetime),
                                 "at_type": self.at_type,
                                 "proprietary_period_groups":
                                     self.proprietary_period_groups,
                                 "proprietary_period": self.proprietary_period,
                                 "non_detection": self.non_detection,
                                 "photometry": {"photometry_group":
                                                {"0": self.photometry, "1": self.photometry2}}
                                 }
                           }
             }

    def _convert_time(self, obsdate):
        frac_of_day = (float(obsdate.hour) +
                       (float(obsdate.minute) +
                        float(obsdate.second) / 60) / 60) / 24
        return "%04i-%02i-%06.3f" % (obsdate.year,
                                     obsdate.month,
                                     obsdate.day + frac_of_day)


class TNSClient(object):
    def __init__(self, opts):
        self._options = opts
        self._api_key = API_KEY

    def _get_url(self, components):
        return "/".join(components)

    def send_bulk_report(self):
        opts = self._options
        print(json.dumps(opts.dictionary))
        with open('data.txt', 'w') as outfile:
           json.dump(opts.dictionary, outfile, indent=2, sort_keys=True)
        report_filename="data.txt"
        report_type="json"
        import api_bulk_report
        #feed_url = self._get_url([opts.base_url, opts.at_report_form])
        #fp = cStringIO.StringIO()
        #fp.write(json.dumps(opts.dictionary))
        #fp.seek(0)
        #params = {"api_key": self._api_key,
        #          "data": fp}
        #r = requests.post(feed_url, data=params)
        #fp.close()
#
#        d = r.json()
#        if d["id_code"] != 200:
#            raise TNSClientError(("Unable to report to TNS,"
#                                  "the id_code is %i" % d["id_code"]))
#        self._report_id = d["data"]["report_id"]
#
    def get_reply_report(self):
        opts = self._options
        feed_url = self._get_url([opts.base_url, opts.at_report_reply])
        params = {"api_key": self._api_key, "report_id": self._report_id}
        r = requests.post(feed_url, data=params)

        new = True
        d = r.json()
        feedback = d["data"]["feedback"]["at_report"][0]
        if "100" in feedback:
            # new object submited
            name = feedback["100"]["objname"]
            prefix = "AT"
        elif "101" in feedback:
            # existing object
            name = feedback["101"]["objname"]
            prefix = feedback["101"]["prefix"]
            new = False
        else:
            name = None
            prefix = None
        return new, name, prefix


class TNSClientError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message


def report(internal_name, ra, dec, mag, mag_err,lim_mag,time_lim, filter_name,filt_lim, obsdate, time_last, mag_last, filt_last, magerr_last, reporter):
    report = Options(internal_name=internal_name, ra=ra, dec=dec,
                     photometry_flux=round(mag,2),
                     photometry_flux_error=round(mag_err,2),
                     photometry_filter_name=filter_name,
                     photometry_obsdate=obsdate,
	                 photometry2_flux=round(mag_last,2),
                     photometry2_flux_error=round(magerr_last,2),
                     photometry2_filter_name=filt_last,
                     photometry2_obsdate=time_last,
                     discovery_datetime=obsdate,
                     non_detection_limiting_flux=round(lim_mag,2),
                     non_detection_obsdate=time_lim,
                     non_detection_filter_name=filt_lim,
                     reporter=reporter
                     )
    connection = TNSClient(report)
    connection.send_bulk_report()
#    time.sleep(3)
#    object_flag, object_name, prefix = connection.get_reply_report()
#    return object_flag, object_name, prefix
#
#
# def test():
#obsdate = datetime.datetime(2015, 7, 3, 3, 45, 15)
#ra = "01:49:12"
#dec = "30:10:21"
#mag = 17.
#mag_err = 0.1
#filter_name = "r"
#internal_name= "ZTF18aaacwya"

#object_flag, object_name = report(internal_name, ra, dec, mag, mag_err, filter_name, obsdate)
#print object_flag, object_name

#if object_name:
#    if object_flag:
#        print "new transient: %s" % object_name
#    else:
#        print "old transient: %s" % object_name[2:]
#else:
#    print "Bad Request"
#
# if __name__ == "__main__":
#     test()
