import datetime
import logging
from datetime import timedelta

import numpy as np
import pandas 
import xarray as xr


def retrieve_list_of_files_from_url(path_catalog, path_data, prefix3='dt_',read_online=True):
    """
    Retrieve a list of files from a URL catalog and construct download links.

    Parameters
    ----------
    path_catalog : str
        The URL to the catalog page containing the list of files to retrieve.
    path_data : str
        The base URL path to the directory where the files are hosted.
    prefix3 : str, optional
        A prefix filter to select files with a specific prefix (default is 'dt*').
    read_online : bool, optional
        If True, constructs online-readable URLs (default is True).

    Returns
    -------
    list
        A sorted list of file URLs, ready for download or online access.

    Notes
    -----
    This function reads the contents of a catalog page hosted at the provided URL (`path_catalog`). It then extracts
    the filenames and constructs file URLs based on the `path_data`. The constructed URLs can be either for online
    reading (default) or direct download, depending on the value of the `read_online` parameter.

    The optional `prefix3` parameter allows filtering files by a specific prefix. Only files with names starting with
    the specified prefix will be included in the list.

    Examples
    --------
    >>> catalog_url = "https://example.com/catalog/"
    >>> data_base_url = "https://example.com/data/"
    >>> file_list = retrieve_list_of_files_from_url(catalog_url, data_base_url)
    >>> for file_url in file_list:
    ...     print(file_url)  # Print the constructed file URLs.
    """
    
    from urllib.request import Request, urlopen, urlretrieve
    from bs4 import BeautifulSoup
    def read_url(url):
        allfiles = list()
        url = url.replace(" ","%20") 
        req = Request(url)
        a = urlopen(req).read()
        soup = BeautifulSoup(a, 'html.parser')
        x = (soup.find_all('a'))
        for i in x:
            file_name = i.extract().get_text() 
            if prefix3 != None:
                if file_name[:3] == prefix3:
                    allfiles.append(file_name)  
            else:
                allfiles.append(file_name)

        return(allfiles)

    list_of_files=read_url(path_catalog)

    for i in range(np.shape(list_of_files)[0]):
        if read_online:
            list_of_files[i]=path_data+list_of_files[i]+"#mode=bytes"
        else:
            list_of_files[i]=path_data+list_of_files[i]
        
    return sorted(list_of_files)
