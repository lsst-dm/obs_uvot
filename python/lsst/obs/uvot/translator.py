
# This file is part of astro_metadata_translator.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the LICENSE file at the top-level directory of this distribution
# for details of code ownership.
#
# Use of this source code is governed by a 3-clause BSD-style
# license that can be found in the LICENSE file.

"""Metadata translation code for uvot FITS headers."""

from __future__ import annotations

from lsst.utils import getPackageDir

__all__ = ("UvotTranslator",)

import os
import logging
import posixpath
import re
from collections.abc import Iterator, Mapping, MutableMapping
from typing import TYPE_CHECKING, Any

import astropy.time
import astropy.units as u
from astropy.coordinates import Angle, EarthLocation
from astropy.io import fits

from astro_metadata_translator import cache_translation, FitsTranslator
from astro_metadata_translator.translators.helpers import altaz_from_degree_headers, is_non_science, tracking_from_degree_headers

if TYPE_CHECKING:
    import astropy.coordinates

log = logging.getLogger(__name__)


class UvotTranslator(FitsTranslator):
    """Metadata translator for uvot standard headers."""

    name = "UVOT"
    """Name of this translation class"""

    supported_instrument = "UVOT"
    """Supports the uvot instrument."""

    default_resource_root = os.path.join(getPackageDir("obs_uvot"),"corrections")
    """Default resource path root to use to locate header correction files."""

    DETECTOR_MAX = 1


    _const_map = {"detector_num": 0,
        "science_program": "None", #is this the appropriate name?
        "instrument": "UVOT",
        "pressure": None,
        "observation_type": "science",
        "relative_humidity": None,
        "pressure": None,
        "temperature": None,
        "focus_z": None,
        "boresight_airmass": None,
        "location": None,
        "boresight_rotation_coord": "sky",
        "detector_serial": None,
        "detector_num": 0,
        "altaz_begin": None,
        }

    _trivial_map: dict[str, str | list[str] | tuple[Any, ...]] = {
        "object": "OBJECT",
        "observation_id": "OBS_ID",
        "telescope": ("TELESCOP", dict(default="SWIFT")),
        "physical_filter": "FILTER",
        "detector_name": "DETNAM",
        }


    @classmethod
    def can_translate(cls, header: Mapping[str, Any], filename: str | None = None) -> bool:
        """Indicate whether this translation class can translate the
            supplied header.

            There is no ``INSTRUME`` header in early HSC files, so this method
            looks for HSC mentions in other headers.  In more recent files the
            instrument is called "Hyper Suprime-Cam".

            Parameters
            ----------
            header : `dict`-like
                Header to convert to standardized form.
            filename : `str`, optional
                Name of file being translated.

            Returns
            -------
            can : `bool`
                `True` if the header is recognized by this class. `False`
                otherwise.
            """
        if "INSTRUME" in header:
            return header["INSTRUME"] == "UVOTA"

        return False

    @cache_translation
    def to_exposure_id(self) -> int:
        """Calculate exposure ID.

        Returns
        -------
        id : `int`
            ID of exposure.
        """
        value = self._header["OBS_ID"]
        self._used_these_cards("OBS_ID")
        return int(value)

    @cache_translation
    def to_observation_counter(self) -> int:
        """Return the lifetime exposure number.

            Returns
            -------
            sequence : `int`
                The observation counter.
            """
        return self.to_exposure_id()

    @cache_translation
    def to_visit_id(self) -> int:
            # Docstring will be inherited. Property defined in properties.py
        return self.to_exposure_id()

    #taken from translators/sdss.py
    @cache_translation
    def to_tracking_radec(self) -> astropy.coordinates.SkyCoord:
        # Docstring will be inherited. Property defined in properties.py
        radecsys = None
        radecpairs = (("RA_PNT", "DEC_PNT"),)
        import pdb
        pdb.set_trace()
        tfdh =  tracking_from_degree_headers(self, radecsys, radecpairs, unit=u.deg)

    @cache_translation
    def to_boresight_rotation_angle(self):
        #double check whether this angle should be 90 or 270 or ?
        angle = Angle(90.*u.deg) - Angle(self.quantity_from_card("PA_PNT", u.deg))
        angle = angle.wrap_at("360d")
        return angle

    @cache_translation
    def to_detector_exposure_id(self) -> int | None:
        # Docstring will be inherited. Property defined in properties.py
        exposure_id = self.to_exposure_id()
        if exposure_id is None:
            return None
        return int(exposure_id)

    @cache_translation
    def to_datetime_begin(self):
        self._used_these_cards("DATE-OBS")
        return Time(self._header["DATE-OBS"], scale="TT", format="isot")

    @cache_translation
    def to_datetime_end(self):
        self._used_these_cards("DATE-END")
        return Time(self._header["DATE-END"], scale="TT", format="isot")

    @cache_translation
    def to_exposure_time(self):
        self._used_these_cards("TSTOP")
        self._used_these_cards("TSTART")
        return self.quantity_from_card("TSTOP", u.s) - self.quantity_from_card("TSTART", u.s)

    @cache_translation
    def to_dark_time(self):
        return self.to_exposure_time()

'''
    @cache_translation
    def to_datetime_end(self) -> astropy.time.Time:
        # Docstring will be inherited. Property defined in properties.py
        # Instcals have no DATE-END or DTUTC
        datetime_end = self._from_fits_date("DTUTC", scale="utc")
        if datetime_end is None:
            datetime_end = self.to_datetime_begin() + self.to_exposure_time()
        return datetime_end

    def _translate_from_calib_id(self, field: str) -> str:
        """Fetch the ID from the CALIB_ID header.

        Calibration products made with constructCalibs have some metadata
        saved in its FITS header CALIB_ID.

        Parameters
        ----------
        field : `str`
            Field to extract from the ``CALIB_ID`` header.

        Returns
        -------
        value : `str`
            The value extracted from the calibration header for that field.
        """
        data = self._header["CALIB_ID"]
        match = re.search(r".*%s=(\S+)" % field, data)
        if not match:
            raise RuntimeError(f"Header CALIB_ID with value '{data}' has not field '{field}'")
        self._used_these_cards("CALIB_ID")
        return match.groups()[0]

    @cache_translation
    def to_physical_filter(self) -> str | None:
        """Calculate physical filter.

        Return `None` if the keyword FILTER does not exist in the header,
        which can happen for some valid Community Pipeline products.

        Returns
        -------
        filter : `str`
            The full filter name.
        """
        if self.is_key_ok("FILTER"):
            value = self._header["FILTER"].strip()
            self._used_these_cards("FILTER")
            return value
        elif self.is_key_ok("CALIB_ID"):
            return self._translate_from_calib_id("filter")
        else:
            return None

    @cache_translation
    def to_location(self) -> astropy.coordinates.EarthLocation:
        """Calculate the observatory location.

        Returns
        -------
        location : `astropy.coordinates.EarthLocation`
            An object representing the location of the telescope.
        """
        if self.is_key_ok("OBS-LONG"):
            # OBS-LONG has west-positive sign so must be flipped
            lon = self._header["OBS-LONG"] * -1.0
            value = EarthLocation.from_geodetic(lon, self._header["OBS-LAT"], self._header["OBS-ELEV"])
            self._used_these_cards("OBS-LONG", "OBS-LAT", "OBS-ELEV")
        else:
            # Look up the value since some files do not have location
            value = EarthLocation.of_site("ctio")

        return value

    @cache_translation
    def to_observation_type(self) -> str:
        """Calculate the observation type.

        Returns
        -------
        typ : `str`
            Observation type. Normalized to standard set.
        """
        if not self.is_key_ok("OBSTYPE"):
            return "none"
        obstype = self._header["OBSTYPE"].strip().lower()
        self._used_these_cards("OBSTYPE")
        if obstype == "object":
            return "science"
        return obstype

    @cache_translation
    def to_altaz_begin(self) -> astropy.coordinates.AltAz:
        # Docstring will be inherited. Property defined in properties.py
        return altaz_from_degree_headers(self, (("ZD", "AZ"),), self.to_datetime_begin(), is_zd={"ZD"})

    @cache_translation
    def to_detector_group(self) -> str:
        # Docstring will be inherited. Property defined in properties.py
        name = self.to_detector_unique_name()
        return name[0]

    @cache_translation
    def to_detector_name(self) -> str:
        # Docstring will be inherited. Property defined in properties.py
        name = self.to_detector_unique_name()
        return name[1:]

    @cache_translation
    def to_focus_z(self) -> u.Quantity:
        # Docstring will be inherited. Property defined in properties.py
        # ``TELFOCUS`` is a comma-separated string with six focus offsets
        # (fx, fy, fz, tx, ty, tz) recorded in units of microns.
        tel_focus_list = self._header["TELFOCUS"].split(",")
        return float(tel_focus_list[2]) * u.um

    @classmethod
    def fix_header(
        cls, header: MutableMapping[str, Any], instrument: str, obsid: str, filename: str | None = None
    ) -> bool:
        """Fix DECam headers.

        Parameters
        ----------
        header : `dict`
            The header to update.  Updates are in place.
        instrument : `str`
            The name of the instrument.
        obsid : `str`
            Unique observation identifier associated with this header.
            Will always be provided.
        filename : `str`, optional
            Filename associated with this header. May not be set since headers
            can be fixed independently of any filename being known.

        Returns
        -------
        modified = `bool`
            Returns `True` if the header was updated.

        Notes
        -----
        Fixes the following issues:

        * If OBSTYPE contains "zero" or "bias",
          update the FILTER keyword to "solid plate 0.0 0.0".

        Corrections are reported as debug level log messages.
        """
        modified = False

        # Calculate the standard label to use for log messages
        log_label = cls._construct_log_prefix(obsid, filename)

        obstype = header.get("OBSTYPE", "unknown")

        if "bias" in obstype.lower() or "zero" in obstype.lower():
            header["FILTER"] = "solid plate 0.0 0.0"
            modified = True
            log.debug("%s: Set FILTER to %s because OBSTYPE is %s", log_label, header["FILTER"], obstype)

        return modified

    @classmethod
    def determine_translatable_headers(
        cls, filename: str, primary: MutableMapping[str, Any] | None = None
    ) -> Iterator[MutableMapping[str, Any]]:
        """Given a file return all the headers usable for metadata translation.

        DECam files are multi-extension FITS with a primary header and
        each detector stored in a subsequent extension.  DECam uses
        ``INHERIT=T`` and each detector header will be merged with the
        primary header.

        Guide headers are not returned.

        Parameters
        ----------
        filename : `str`
            Path to a file in a format understood by this translator.
        primary : `dict`-like, optional
            The primary header obtained by the caller. This is sometimes
            already known, for example if a system is trying to bootstrap
            without already knowing what data is in the file. Will be
            merged with detector headers if supplied, else will be read
            from the file.

        Yields
        ------
        headers : iterator of `dict`-like
            Each detector header in turn. The supplied header will be merged
            with the contents of each detector header.

        Notes
        -----
        This translator class is specifically tailored to raw DECam data and
        is not designed to work with general FITS files. The normal paradigm
        is for the caller to have read the first header and then called
        `determine_translator()` on the result to work out which translator
        class to then call to obtain the real headers to be used for
        translation.
        """
        # Circular dependency so must defer import.
        from ..headers import merge_headers

        # This is convoluted because we need to turn an Optional variable
        # to a Dict so that mypy is happy.
        primary_hdr = primary if primary else {}

        # Since we want to scan many HDUs we use astropy directly to keep
        # the file open rather than continually opening and closing it
        # as we go to each HDU.
        with fits.open(filename) as fits_file:
            # Astropy does not automatically handle the INHERIT=T in
            # DECam headers so the primary header must be merged.
            first_pass = True

            for hdu in fits_file:
                if first_pass:
                    if not primary_hdr:
                        primary_hdr = hdu.header
                    first_pass = False
                    continue

                header = hdu.header
                if "CCDNUM" not in header:  # Primary does not have CCDNUM
                    continue
                if header["CCDNUM"] > 62:  # ignore guide CCDs
                    continue
                yield merge_headers([primary_hdr, header], mode="overwrite")

    @classmethod
    def observing_date_to_offset(cls, observing_date: astropy.time.Time) -> astropy.time.TimeDelta | None:
        """Return the offset to use when calculating the observing day.

        Parameters
        ----------
        observing_date : `astropy.time.Time`
            The date of the observation. Unused.

        Returns
        -------
        offset : `astropy.time.TimeDelta`
            The offset to apply. The offset is always 12 hours. DECam has
            no defined observing day concept in its headers. To ensure that
            observations from a single night all have the same observing_day,
            adopt the same offset used by the Vera Rubin Observatory of
            12 hours.
        """
        return cls._observing_day_offset
'''
