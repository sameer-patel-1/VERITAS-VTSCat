
import git
import yaml
import glob
import shutil
import ast
import os
from astropy.table import QTable, Table
from astropy.io import fits
import pandas as pd
import numpy as np
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.nddata import NDData
from collections import OrderedDict
import warnings
warnings.simplefilter("ignore", UserWarning)
import re

def flatten_dicts_lists(filename): # has to be read in first
    file = Table.read(filename, format='ascii.ecsv')
    meta_copy = file.meta.copy()
    for meta in file.meta:
        if type(file.meta[meta]) is list:
           meta_copy[meta] = ",".join(file.meta[meta])

        if type(file.meta[meta]) is dict:
            for key in file.meta[meta]:
                meta_copy[str(meta)+"_"+str(key)] = file.meta[meta][key]
            del meta_copy[meta]

    file.meta = meta_copy
    return file

def fix_line_break(file):
    for meta in file.meta:
       if type(file.meta[meta]) is str:
           if "\n" in file.meta[meta]:
               file.meta[meta] = file.meta[meta].replace('\n', '')
           file.meta[meta] = file.meta[meta].replace('–', '-') # this is annoying!!
    return file

def assert_str(file):
    for meta in file.meta:
        file.meta[meta] = str(file.meta[meta])
    return file

def merge_fits(filetype, src_id, paper_ref):
    print("Currently in %s/%s" % (src_id, paper_ref))

    try:
        fixed_files = []

        filelst = [i for i in sorted(glob.glob("*.ecsv")) if filetype in i]

        filename_lst = [i.replace('.ecsv', '') for i in filelst]

        if filelst == []:
            print("No %s files found in %s/%s" % (filetype, src_id, paper_ref))
            return None

        if len(filelst) > 1:
            filelst = sorted(filelst, key=lambda x:float(re.findall("(\d+)",x)[1]))
            filename_lst = sorted(filename_lst, key=lambda x:float(re.findall("(\d+)",x)[1]))

        for file in filelst:
            file = flatten_dicts_lists(file)
            file = fix_line_break(file)
            file = assert_str(file)
            fixed_files.append(file)


        hdr = fits.Header()
        hdr['DATATYPE'] = str(filetype).strip()
        hdr['SRC_ID'] = str(src_id).strip()
        hdr['REFERENC'] = str(paper_ref).strip()

        primary_hdu = fits.PrimaryHDU(header=hdr)
        hdu_list = [primary_hdu]

        for i in range(len(fixed_files)):
            table_name = filename_lst[i]
            file_data = fixed_files[i]
            hdu_list.append(fits.BinTableHDU(data=file_data,header=file_data.meta, name=table_name))

        hdul = fits.HDUList(hdu_list)
        hdul.writeto('%s-%s.fits' % (src_id, filetype), overwrite=True)

    except IndexError:
        print("No %s files found in %s/%s" % (filetype, src_id, paper_ref))
    return None

def merge_main(data_dir, filetype):
    os.chdir(data_dir)
    dir_list = sorted([i for i in sorted(os.listdir()) if os.path.isdir(i)])

    for ver_x in dir_list:
        os.chdir(data_dir + ver_x)

        paper_lst = sorted([i for i in sorted(os.listdir()) if os.path.isdir(i)])

        for paper_dir in paper_lst:
            try:
                os.chdir(data_dir + ver_x + "/" + paper_dir)
                merge_fits(filetype, ver_x, paper_dir)
                # merge_fits('lc')

                # write_tables(obs_data, '%s-obs_data' % ver_x)
            except FileNotFoundError or IndexError:
                print("No %s files found in %s/%s" % (filetype, ver_x, paper_dir))
                pass
    return None

# =============================================================================
# Main Program
# =============================================================================
if __name__ == "__main__":
    repo = git.Repo(".", search_parent_directories=True)
    repo_dir = repo.working_tree_dir + "/" # establish pwd as the git repo
    heasarc_dir = repo_dir+"heasarc/" # base dir for heasarc files/folders
    source_dir = heasarc_dir+"sources/"

    data_dir = heasarc_dir+'detected_data/'