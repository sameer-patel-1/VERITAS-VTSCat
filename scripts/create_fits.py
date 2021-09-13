
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
# from astropy.nddata import NDData
from collections import OrderedDict
import warnings
warnings.simplefilter("ignore", UserWarning)
import re
from astropy.time import Time

def sort_source(name):
    basename = os.path.basename(name).split('.')[0]
    num = basename.split('-')[1]
    return int(num)

def sort_data(name):
    num = name.split('/')[4] # sorted by year
    return int(num)

def list_dir(directory): # default os.listdir returns only basepath
    return [os.path.join(directory, file) for file in os.listdir(directory)]

def get_src_data(src_id, field, source_dir):
    stream = open(source_dir+'/VER-%d' % src_id, 'r')
    src_file = yaml.load(stream, Loader=yaml.FullLoader)
    df = pd.json_normalize(src_file, sep='_')
    val = df[field]
    return None

def get_main_header(type_data, src_id, field, paper_ref, source_dir):
    headers={'TELESCOP': ('VERITAS','Telescope Name'),
    'OBS_ID': ('','Placeholder'),
    'OBJECT': (get_src_data(src_id, 'veritas_name_name', source_dir),'Object Name'),
    'RA_OBJ': (get_src_data(src_id, 'pos_ra', source_dir), 'Object Right Ascension'),
    'DEC_OBJ': (get_src_data(src_id, 'pos_dec', source_dir), 'Object Declination'),
    'DATE-OBS': ('Start Date'),
    'DATE-END': ('Stop Date'),
    'SRC_ID': (src_id, 'VERITAS Catalog Source ID'),
    'REFERENC': (os.path.basename(paper_ref), 'ADS Bibcode for data in paper'),
    'DATATYPE': (type_data.upper(), 'Type of Data File'),
    'DATE': (Time(Time.now()).fits, 'File Creation Date')}
    return headers[field]

def get_data_headers(type_data, src_id, field, paper_ref, source_dir):
    headers={'DATA_TYPE': (type_data.upper(), 'type of data'),
    'SRC_ID': (src_id, 'VERITAS Catalog Source ID'),
    'REFERENC': (os.path.basename(paper_ref), 'Reference for the data'),
    'REFNOTES': (os.path.basename(paper_ref), 'Reference for the data'),
    'TELESCOP': ('Telescope Name'),
    'UL_CONF': ('upper level confidence'),
    'COMMENTS': ('reference notes'),
    'SED_TYPE': ('type of SED'),
    'DATA_SOURCE': ('where/how the data was obtained'),
    'MJD_MAX': ('MJD start time'),
    'MJD_MIN': ('MJD end time'),
    'DESCRIPTION': ('misc. comments pertaining to data file')}
    return None

def get_data_comments():
    data={'e_ref': 'reference value for energy bin of spectral measurements',
    'e_min': 'min value of energy bin',
    'e_max': 'max value of energy bin',
    'e_diff': 'width of energy bin',
    'dnde': 'differential photon flux',
    'dnde_err': 'average error in differential photon flux',
    'dnde_errn': 'lower error in differential photon flux',
    'dnde_errp': 'upper error in differential photon flux',
    'significance': 'statistical significance of measurement of gamma-ray excess in this energy bin',
    'dnde_ul': 'upper limit on differential photon flux',
    'e2dnde': 'differential energy flux',
    'e2dnde_err': 'error in differential energy flux',
    'dnde_syserr': 'average systematic error in differential photon flux',
    'dnde_syserrn': 'lower systematic error in differential photon flux',
    'dnde_syserrp': 'upper systematic error in differential photon flux',
    'e2dnde_ul': 'upper limit on differential energy flux',
    'non': 'number of counts in on region',
    'alpha': 'normalisation factor between on and off regions',
    'noff': 'number of counts in off region',
    'e_min': 'minimum energy for flux integration',
    'e_max': 'maximum energy for flux integration',
    'e_cut': 'maximum energy for flux integration',
    'time_min': 'lower edge of time bin',
    'time_max': 'upper edge of time bin',
    'time_err': 'width of time bin',
    'flux': 'integral energy flux above e_min',
    'flux_err': 'average error in integral photon flux',
    'flux_errn': 'lower error in integral photon flux',
    'flux_errp': 'upper error in integral photon flux',
    'flux_ul': 'upper limit on photon energy flux',
    'significance': 'statistical significance of measurement of gamma-ray excess in this bin',
    'live_time': 'observation time',
    'time': 'time bin center',
    'phase_min': 'lower value of orbital phase',
    'phase_max': 'upper value of orbital phase',
    'index': 'index of the power law fit',
    'index_err': 'error in the index of the power law fit'
    'time_bin': 'width of time bin',
    'fit_e_min': 'lower energy bound of power-law function fit range',
    'fit_e_max': 'upper energy bound of power-law function fit range',
    'f_var': 'fractional variability',
    'f_var_err': 'error in fractional variability',
    'non': 'number of counts in on region'',
    'noff': 'number of counts in off region',
    'alpha': 'normalisation factor between on and off regions',
    'nex': '(number) excess'}
    return None

def get_col_names(data_dir):
    lc_cols = []
    sed_cols = []

    data_directories = sorted([d for d in list_dir(data_dir) if os.path.isdir(d)], key=sort_source)
    for src_dir in data_directories:
        paper_directories = sorted([d for d in list_dir(src_dir) if os.path.isdir(d)])
        for paper_dir in paper_directories:
            # print(paper_dir)
            sed_lst = glob.glob(paper_dir+"/*-sed.ecsv")
            lc_lst = glob.glob(paper_dir+"/*-lc.ecsv")

            for sed_file in sed_lst:
                sed_table = Table.read(sed_file, format='ascii.ecsv')
                # sed_col_lst = list(sed_table.columns)
                sed_col_lst = list(sed_table.meta.keys())
                for col in sed_col_lst:
                    if col not in sed_cols:
                        sed_cols.append(col)
                        # print(col, sed_file)

            for lc_file in lc_lst:
                lc_table = Table.read(lc_file, format='ascii.ecsv')
                # lc_col_lst = list(lc_table.columns)
                lc_col_lst = list(lc_table.meta.keys())
                if 'SED_TYPE' in lc_col_lst:print(lc_file)
                for col in lc_col_lst:
                    if col not in lc_cols:
                        lc_cols.append(col)
                        # print(col, lc_file)

    return lc_cols, sed_cols

def flatten_dicts_lists(filename): # has to be read in first
    file = Table.read(filename, format='ascii.ecsv')
    meta_copy = file.meta.copy()
    for meta in file.meta:
        if type(file.meta[meta]) is list:meta_copy[meta] = ",".join(file.meta[meta])

        if type(file.meta[meta]) is dict:
            for key in file.meta[meta]:
                meta_copy[str(meta)+"_"+str(key)] = file.meta[meta][key]
            del meta_copy[meta]

    file.meta = meta_copy
    return file

def fix_line_break(file):
    for meta in file.meta:
       if type(file.meta[meta]) is str:
           if "\n" in file.meta[meta]:file.meta[meta] = file.meta[meta].replace('\n', '')
           file.meta[meta] = file.meta[meta].replace('–', '-') # this is annoying!!
    return file

def assert_str(file):
    for meta in file.meta:
        file.meta[meta] = str(file.meta[meta])
    return file

def get_val_from_file(file, field):

    return None

def get_main_headers(filetype, src_id, paper_ref, file):
    # hdr = fits.Header()
    # c1 = fits.Card('TELESCOP', 'VERITAS', get_main_header('TELESCOP'))
    # c1 = fits.Card('OBS_ID', '', get_header_comments('OBS_ID'))

    
    hdr['DATATYPE'] = str(filetype).strip()
    hdr['SRC_ID'] = str(src_id).strip()
    hdr['REFERENC'] = str(paper_ref).strip()
    return hdr

def merge_fits(filetype, src_id, paper_ref):
    print("Currently in %s/%s" % (src_id, paper_ref))

    try:
        fixed_files = []

        filelst = [i for i in sorted(glob.glob("*.ecsv")) if filetype in i]

        filename_lst = [i.replace('.ecsv', '') for i in filelst]

        if not filelst:
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

        main_hdr = get_main_headers(filetype, src_id, paper_ref, file)
        hdu_list = [fits.PrimaryHDU(header=main_hdr)]

        for i in range(len(fixed_files)):
            table_name = filename_lst[i]
            file_data = fixed_files[i]
            hdu_list.append(fits.BinTableHDU(data=file_data,header=file_data.meta, name=table_name))

        hdul = fits.HDUList(hdu_list)
        hdul.writeto('%s-%s.fits' % (src_id, filetype), overwrite=True, checksum=True)

    except IndexError:
        print("No %s files found in %s/%s" % (filetype, src_id, os.path.basename(paper_ref)))
    return None

def merge_main(data_dir, filetype):
    os.chdir(data_dir)
    dir_list = sorted([d for d in list_dir(data_dir) if os.path.isdir(d)], key=sort_source)


    for ver_x in dir_list:
        os.chdir(ver_x)
        src_id = int(os.path.basename(ver_x).split('-')[-1])
        paper_lst = sorted([i for i in sorted(os.listdir()) if os.path.isdir(i)], sort=sort_data)
        paper_directories = sorted([d for d in list_dir(ver_x) if os.path.isdir(d)])

        for paper_dir in paper_directories:
            try:
                os.chdir(paper_dir)
                merge_fits(filetype, src_id, paper_dir)
                # merge_fits('lc')

                # write_tables(obs_data, '%s-obs_data' % ver_x)
            except FileNotFoundError or IndexError:
                print("No %s files found in %s" % (filetype, ver_x))
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

    os.chdir(repo_dir+'scripts')
    # lc_cols, sed_cols = get_col_names(data_dir)