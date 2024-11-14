"""This module provides convenient tools for loading netcdf files"""
import glob
import logging
import os
import re
from dataclasses import dataclass
from datetime import datetime
from typing import List, Tuple

import numpy as np
import xarray as xr
import dask.array

SWOT_EXPERT_PATTERN = re.compile(
    r"SWOT_L2_LR_SSH_Expert_(?P<ncycle>\d+)_(?P<ntrack>\d{3})_"
    r"(?P<date_start>\d{8}T\d{6})_(?P<date_end>\d{8}T\d{6})_DG10_01.nc")

logger = logging.getLogger(__name__)


def open_netcdf(files):
    raise NotImplementedError("Functional interface to do")

    # if wanted_variables is not None:
    #     files = glob.glob(regex)
    #     ds = xr.open_mfdataset(files[:1], **reading_options)
    #     reading_options["drop_variables"] = [
    #         v for v in ds if v not in wanted_variables]



@dataclass
class NetcdfFilesUtilities:
    """ Utilities centered around the information of the netcdf file names.

    Given the regex pattern of the netcdf files names, one can extract
    information about the cycle, track, and date range one file covers.

    This class has been built with the assumption that a file covers multiple
    dates, and that it is reflected in its name.

    Attributes
    ----------
    regex
        The regular expression pattern of the file names.
    date_start_name
        Group name for date range start in the regex
    date_end_name
        Group name for date range end in the regex
    date_fmt
        Date format for converting the date string of the name into a numpy
        datetime64
    cycle_name
        Group name for the cycle number in the regex
    track_name
        Group name for the track number name in the regex
    """
    regex: re.Pattern = SWOT_EXPERT_PATTERN
    date_start_name: str = "date_start"
    date_end_name: str = "date_end"
    cycle_name: str = "ncycle"
    track_name: str = "ntrack"
    date_fmt: str = "%Y%m%dT%H%M%S"

    def __post_init__(self):
        # Check groups are in the pattern
        self._check_groups_in_pattern(self.regex, [
            self.date_start_name, self.date_end_name, self.cycle_name,
            self.track_name
        ])

    def filter_date_range(self, path, date_start: np.datetime64,
                          date_end: np.datetime64) -> List[str]:
        """ Filter input files names in a given time interval.

        Parameters
        ----------
        path
            glob expression for listing the netcdf files
        date_start
            Starting date for the filtering
        date_end
            Ending date for the filtering

        Returns
        -------
        :
            The filtered files list
        """
        files = glob.glob(path)
        return [
            f for f in files if self._date_range_filter(
                os.path.basename(f), date_start, date_end)
        ]

    def _date_range_filter(self, file_name: str, date_start: np.datetime64,
                           date_end: np.datetime64) -> bool:
        # File selection
        match_object = self.regex.search(file_name)
        if match_object is None:
            return False

        # In case there is no match
        try:
            date_start_file = np.datetime64(
                datetime.strptime(match_object.group(self.date_start_name),
                                  self.date_fmt))
            date_end_file = np.datetime64(
                datetime.strptime(match_object.group(self.date_end_name),
                                  self.date_fmt))

            if date_start_file <= date_end and date_end_file >= date_start:
                return True
        except ValueError:
            # In case the integer conversion failed. This should not happen if
            # the input regex is properly configured (with groups defined with
            # \d)
            raise Exception(
                f"File name '{file_name}' matched but the extracted dates"
                f" '{match_object.group(self.date_start_name)}' and"
                f" '{match_object.group(self.date_end_name)}' could not be"
                f" converted to a numpy datetime.")

        return False

    def _check_groups_in_pattern(self, regex: re.Pattern,
                                 group_names: List[str]):
        # Check that the pattern has the necessary groups. It will raise an
        # exception if there is a missing group
        regex_groups = regex.groupindex.keys()

        missing_groups = [g for g in group_names if g not in regex_groups]
        if len(missing_groups) > 0:
            raise Exception(
                f"Input file name pattern '{regex.pattern}' misses one of the"
                f" regex groups: '{missing_groups}'")

    def _track_cycle_from_filename(self, file_name: str) -> Tuple[int, int]:
        # Match the file name
        match_object = self.regex.search(file_name)
        try:
            ncycle = int(match_object.group(self.cycle_name))
            ntrack = int(match_object.group(self.track_name))
        except AttributeError:
            # In case the returned match object is None
            raise Exception(f"File name '{file_name}' did not match pattern"
                            f"'{self.regex.pattern}'")
        except ValueError:
            # In case the integer conversion failed. This should not happen if
            # the input regex is properly configured (with groups defined with
            # \d)
            raise Exception(
                f"File name '{file_name}' matched but the extracted cycle"
                f" '{match_object.group(self.cycle_name)}' or track"
                f" '{match_object.group(self.track_name)}' could not be"
                f" converted to a valid number.")

        return ncycle, ntrack

    def preprocess_track_cycle(self, ds: xr.Dataset) -> xr.Dataset:
        """ Extract track and cycle from the file name put it in the dataset.

        Parameters
        ----------
        ds:
            Input dataset for which we wish to add the ntrack and ncycle fields

        Returns
        -------
        : xr.Dataset
            Dataset with tracks and cycles numbers added in dataset
        """
        # Retrieve the file name. This means the source of the dataset should
        # be a netcdf
        try:
            file_name = ds.encoding["source"]
        except KeyError:
            raise Exception(
                "Dataset did not contained the filename, check that the input"
                " dataset has been opened from a file set")

        ncycle, ntrack = self._track_cycle_from_filename(file_name)

        # Put the cycle and track numbers in the dataset
        ds[self.cycle_name] = (ds.num_lines.dims, ncycle*dask.array.ones(ds.num_lines.size))
        ds[self.track_name] = (ds.num_lines.dims, ntrack*dask.array.ones(ds.num_lines.size))
        return ds

    def preprocess_groupby_cycle(self, ds):
        """ Add the cycle dimension from the file name.

        Parameters
        ----------
        ds:
            Input dataset for which we wish to group the cycle

        Returns
        -------
        : xr.Dataset
            Dataset with (cycle, num_lines, num_pixels) dimensions
        """
        # Retrieve the file name. This means the source of the dataset should
        # be a netcdf
        try:
            file_name = ds.encoding["source"]
        except KeyError:
            raise Exception(
                "Dataset did not contained the filename, check that the input"
                " dataset has been opened from a file set")

        ncycle, _ = self._track_cycle_from_filename(file_name)

        return ds.expand_dims(cycle=[ncycle])

    def preprocess_groupby_track(self, ds: xr.Dataset) -> xr.Dataset:
        """ Add the track and cycle dimensions from the file name.

        Parameters
        ----------
        ds:
            Input dataset for which we wish to group the track and cycle

        Returns
        -------
        : xr.Dataset
            Dataset with (cycle, track, num_lines, num_pixels) dimensions
        """
        # Retrieve the file name. This means the source of the dataset should
        # be a netcdf
        try:
            file_name = ds.encoding["source"]
        except KeyError:
            raise Exception(
                "Dataset did not contained the filename, check that the input"
                " dataset has been opened from a file set")

        ncycle, ntrack = self._track_cycle_from_filename(file_name)

        return ds.expand_dims(track=[ntrack], cycle=[ncycle])

    def nested_files(self, path: str, crop_if_incomplete: bool = True):
        """Scan a folder and organize the netcdf files between cycle and tracks.

        This function is useful to preprocess the file names and prepare a list
        compatible with xr.open_mfdataset(combine="nested") with multiple
        dimensions.

        Parameters
        ----------
        path
            Folder containing the netcdf files
        crop_if_incomplete
            In case the last cycle does not have the same number of tracks,
            whether to keep it or not

        Returns
        -------
        :
            A list containing lists of netcdf files. The first level of nesting
            is for cycle and the second one for tracks.
        """

        netcdf_files = sorted(glob.glob(path))

        output = []
        last_cycle = 0
        for file_name in netcdf_files:
            ncycle, _ = self._track_cycle_from_filename(file_name)

            if ncycle != last_cycle:
                output.append([])
            output[-1].append(file_name)

            last_cycle = ncycle

        nested_files = output
        if len(output[-1]) != len(output[0]):
            if crop_if_incomplete:
                logger.warn("Last cycle is incomplete and has been left out")
                nested_files = output[:-1]
            else:
                logger.warn(
                    "Last cycle is incomplete but still in the list, the"
                    "resulting nested list will be malformed"
                )

        return nested_files
