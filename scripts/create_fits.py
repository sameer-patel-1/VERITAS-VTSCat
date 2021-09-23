import git
import yaml
import glob
import os
import tarfile
import pandas as pd
import numpy as np

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import QTable, Table
from astropy.io import fits
import warnings
warnings.simplefilter("ignore", UserWarning)

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
    stream = open(source_dir+'VER-%d.yaml' % src_id, 'r')
    src_file = yaml.load(stream, Loader=yaml.FullLoader)
    df = pd.json_normalize(src_file, sep='_')
    try:
        df[field][0]
    except KeyError:
        if src_id == 175: # special exception
            if field == 'pos_ra':return Angle('17h46m41.01487s').degree
            if field == 'pos_dec':return Angle('-28d38m30.5345s').degree
            print('Setting RA/DEC for src_id = 175 manually...')

        if src_id>999:return None # undetected sources have no VERITAS name
        else:
            print('Check source file VER-%d for %s entry' % (src_id, field))
    val = df[field][0]
    return val

def get_obs_time(data_type, file=None):
    if file is None:data_files = sorted(glob.glob('*-%s.ecsv' % data_type)) # will iterate over all files in the paper
    else:data_files = [file] # only get values from the single SED/LC file

    date_min = np.inf # set a global min
    date_max = -np.inf # set a global max

    for data_file in data_files:
        data = Table.read(data_file, format='ascii.ecsv')
        df_meta = pd.json_normalize(data.meta, sep='_')
        meta_cols = list(df_meta.columns)

        if 'mjd_min' in meta_cols and 'mjd_max' in meta_cols: # meta has mjd_min and mjd_max; first priority
            if df_meta['mjd_min'][0] < date_min:date_min = df_meta['mjd_min'][0] # global min
            if df_meta['mjd_max'][0] > date_max:date_max = df_meta['mjd_max'][0] # global max
        else: # not in meta; try in data cols
            data_cols = list(data.columns)
            if 'time' in data_cols:
                if min(data['time'].data) < date_min:date_min = min(data['time'].data) # global min
                if max(data['time'].data) > date_max:date_max = max(data['time'].data) # global max
            if 'time_min' in data_cols and 'time_max' in data_cols:
                if min(data['time_min'].data) < date_min:date_min = min(data['time_min'].data) # global min
                if max(data['time_max'].data) > date_max:date_max = max(data['time_max'].data) # global max

    if date_min == np.inf and date_max == -np.inf: # if date_min/max not found anywhere, set to None
        return None,None,None,None
    else:
        return date_min, date_max, Time(date_min,format='mjd').fits, Time(date_max,format='mjd').fits

def get_main_header(data_type, src_id, ads_paper, source_dir):
    min_mjd,max_mjd,min_fits,max_fits = get_obs_time(data_type)

    hdr_lst = [
    fits.Card('TELESCOP','VERITAS','Telescope name'),
    fits.Card('OBS_ID',None,None),
    fits.Card('SRC_ID', src_id, 'VERITAS catalog source ID'),
    fits.Card('OBJECT', get_src_data(src_id, 'veritas_name_name', source_dir), 'Object name'),
    fits.Card('RA_OBJ', get_src_data(src_id, 'pos_ra', source_dir), 'Object right ascension, in deg'),
    fits.Card('DEC_OBJ', get_src_data(src_id, 'pos_dec', source_dir), 'Object declination, in deg'),
    fits.Card('DATE-OBS', min_fits, 'Start date'),
    fits.Card('DATE-END', max_fits, 'Stop date'),
    fits.Card('REFERENC', os.path.basename(ads_paper), 'ADS bibcode where data appears'),
    fits.Card('DATATYPE', data_type.upper(), 'Type of data file'),
    fits.Card('DATE',Time(Time.now()).fits, 'File creation date')
    ]
    return hdr_lst

def fix_line_break(file): # to get rid of escape chars and non-standard strings
    for meta in file.meta:
       if type(file.meta[meta]) is str:
           if "\n" in file.meta[meta]:file.meta[meta] = file.meta[meta].replace('\n', '')
           file.meta[meta] = file.meta[meta].replace('–', '-') # this is annoying!!
    return file

def fix_mjd_unit(file): # because Unit 'MJD' is not supported by the FITS standard
    for column in file.columns:
        if str(file[column].unit).lower() == 'mjd':
            file[column].unit = u.day
    return file

def get_meta_data(data_file, field):
    data = Table.read(data_file, format='ascii.ecsv')
    file = fix_line_break(data) # first fix all the annoying/non-standard things
    file = fix_mjd_unit(file) # fix mjd units

    meta = pd.json_normalize(file.meta, sep='_')
    try:
        if type(meta[field][0]) is list: # because we only want strings
            meta[field][0] = '; '.join(meta['field'][0])
        return meta[field][0]
    except KeyError:
        return None

def get_file_data(data_file, field, arg=None):
    data = Table.read(data_file, format='ascii.ecsv')
    try:
        val = data[field].data
        if arg=='min':val = min(val)
        if arg=='max':val = max(val)
        if arg=='unit':val = data[field].unit
    except KeyError:
        val = None
    return val

def get_data_header(data_type, data_file, src_id, ads_paper, source_dir):
    min_mjd,max_mjd,min_fits,max_fits = get_obs_time(data_type,file=data_file)

    main_lst = [
    fits.Card('TELESCOP','VERITAS','Telescope name'),
    fits.Card('OBS_ID',None,None),
    fits.Card('SRC_ID', src_id, 'VERITAS catalog source ID'),
    fits.Card('OBJECT', get_src_data(src_id, 'veritas_name_name', source_dir), 'Object name'),
    fits.Card('RA_OBJ', get_src_data(src_id, 'pos_ra', source_dir), 'Object right ascension, in deg'),
    fits.Card('DEC_OBJ', get_src_data(src_id, 'pos_dec', source_dir), 'Object declination, in deg'),
    fits.Card('DATATYPE', data_type.upper(), 'Type of data file'),
    fits.Card('DATE-OBS', min_fits, 'Start date'),
    fits.Card('DATE-END', max_fits, 'Stop date'),
    fits.Card('REFERENC', os.path.basename(ads_paper), 'ADS bibcode where data appears'),
    fits.Card('DATE',Time(Time.now()).fits, 'File creation date')
    ]

    meta_lst = [
    fits.Card('REFNOTES', get_meta_data(data_file, 'comments'), 'Reference notes'),
    fits.Card('DATASRCE', get_meta_data(data_file, 'data_source'), 'Where/how the data was obtained'),
    fits.Card('DESCRIPT', get_meta_data(data_file, 'description'), 'Misc. comments pertaining to data file'),
    fits.Card('UL_CONF', get_meta_data(data_file, 'UL_CONF'), 'Upper limit confidence')
    ]

    if data_type=='sed':
        meta_lst.append(fits.Card('MJD_MIN', min_mjd, 'MJD start time'))
        meta_lst.append(fits.Card('MJD_MAX', max_mjd, 'MJD end time'))

    if data_type=='lc':
        meta_lst.append(fits.Card('TIMESYS', 'MJD', 'Reference time'))
        meta_lst.append(fits.Card('TIMEUNIT', str(u.Unit('day')), 'Unit for TSTART, TSTOP, TIMEZER'))
        meta_lst.append(fits.Card('TSTART', min_mjd, 'Start time counted from TIMESYS'))
        meta_lst.append(fits.Card('TSTOP', max_mjd, 'Stop time counted from TIMESYS'))
        meta_lst.append(fits.Card('E_MIN', get_file_data(data_file, 'e_min', arg='min'), 'Min. of the energy range'))
        meta_lst.append(fits.Card('E_MAX', get_file_data(data_file, 'e_max', arg='max'), 'Max. of the energy range'))
        meta_lst.append(fits.Card('E_UNIT', str(get_file_data(data_file, 'e_min', arg='unit')), 'Physical unit for E_MIN, E_MAX'))

    return main_lst+meta_lst

def chk_header(hdu):
    hdr = hdu.header
    for hdr_key in list(hdr):
        if hdr_key != 'OBS_ID': # this wil be empty since Antara will fill out later
            try:
                assert hdr.comments[hdr_key] != '' # should not be blank
            except AssertionError:
                print('Warning: %s has no comment card' % hdr_key)
    return None

def add_hdr_comments(data_hdu, data_type):
    
    ttype_lc = {
    'e_min': 'min. energy for flux integration',
    'e_max': 'max. energy for flux integration',
    'e_cut': 'cut-off energy for flux integration',
    'time_min': 'lower edge of time bin',
    'time_max': 'upper edge of time bin',
    'time_err': 'width of time bin',
    'flux': 'integral energy flux above e_min',
    'flux_err': 'average error in integral photon flux',
    'flux_errn': 'lower error in integral photon flux',
    'flux_errp': 'upper error in integral photon flux',
    'flux_ul': 'upper limit on photon energy flux',
    'significance': 'stat. significance of gamma-ray excess in bin',
    'live_time': 'observation time',
    'time': 'time bin center',
    'phase_min': 'lower value of orbital phase',
    'phase_max': 'upper value of orbital phase',
    'index': 'index of the power law fit',
    'index_err': 'error in the index of the power law fit',
    'time_bin': 'width of time bin',
    'fit_e_min': 'lower energy bound of power-law function range',
    'fit_e_max': 'upper energy bound of power-law function range',
    'f_var': 'fractional variability',
    'f_var_err': 'error in fractional variability',
    'non': 'number of counts in on region',
    'noff': 'number of counts in off region',
    'alpha': 'normalisation factor between on and off regions',
    'nex': 'number of excess counts'}

    ttype_sed = {
    'e_ref': 'ref. value for energy bin',
    'e_min': 'min. value of energy bin',
    'e_max': 'max. value of energy bin',
    'e_diff': 'width of energy bin',
    'dnde': 'diff. photon flux',
    'dnde_err': 'average error in diff. photon flux',
    'dnde_errn': 'lower error in diff. photon flux',
    'dnde_errp': 'upper error in diff. photon flux',
    'significance': 'stat. significance of gamma-ray excess in bin',
    'dnde_ul': 'upper limit on diff. photon flux',
    'e2dnde': 'diff. energy flux',
    'e2dnde_err': 'error in diff. energy flux',
    'dnde_syserr': 'average systematic error in diff. photon flux',
    'dnde_syserrn': 'lower systematic error in diff. photon flux',
    'dnde_syserrp': 'upper systematic error in diff. photon flux',
    'e2dnde_ul': 'upper limit on diff. energy flux',
    'non': 'number of counts in on region',
    'alpha': 'normalisation factor between on and off regions',
    'noff': 'number of counts in off region'
    }

    # taken from 2010A&A...524A..42P, Table 18
    tform={
    'L': 'field data format: Logical',
    'X': 'field data format: Bit',
    'B': 'field data format: Unsigned byte',
    'I': 'field data format: 16-bit integer',
    'J': 'field data format: 32-bit integer',
    'K': 'field data format: 64-bit integer',
    'A': 'field data format: Character',
    'E': 'field data format: Single precision float',
    'D': 'field data format: Double precision float',
    'C': 'field data format: Single precision complex',
    'M': 'field data format: Double precision complex',
    'P': 'field data format: Array Descriptor (32-bit)',
    'Q': 'field data format: Array Descriptor (64-bit)',
    }

    hdr = data_hdu.header
    for hdr_key in list(hdr):
        val = hdr[hdr_key]
        if 'TTYPE' in hdr_key and data_type=='lc':hdr.comments[hdr_key] = str.capitalize(ttype_lc[val])
        if 'TTYPE' in hdr_key and data_type=='sed':hdr.comments[hdr_key] = str.capitalize(ttype_sed[val])
        if 'TFORM' in hdr_key:hdr.comments[hdr_key] = str.capitalize(tform[val])
        if 'TUNIT' in hdr_key:hdr.comments[hdr_key] = 'Physical unit of field'
    return data_hdu

def get_data_hdu(data_type, data_file, src_id, ads_paper, source_dir):
    data_hdr = get_data_header(data_type, data_file, src_id, ads_paper, source_dir)
    data = Table.read(data_file, format='ascii.ecsv')
    data.meta = {} # remove the meta since otherwise there will be duplicate entries in the header
    data = fix_mjd_unit(data) # fix mjd_unit
    num_file = int(data_file.split('.')[0].split('-')[-2])
    ext_name = 'VER-%d-%d-%s' % (src_id, num_file, data_type.upper())
    data_hdu = fits.BinTableHDU(data=data, header=fits.Header(cards=data_hdr), name=ext_name)
    data_hdu = add_hdr_comments(data_hdu, data_type)
    data_hdu.add_checksum() # adds checksum and datasum
    chk_header(data_hdu)
    return data_hdu

def get_main_hdu(data_type, src_id, ads_paper, source_dir):
    main_hdr = get_main_header(data_type, src_id, ads_paper, source_dir)
    main_hdu = fits.PrimaryHDU(header=fits.Header(cards=main_hdr)) # no data, just header
    main_hdu.add_checksum() # adds checksum and datasum
    chk_header(main_hdu)
    return main_hdu

def get_col_names(data_dir): # just to get a list of all SED and LC datafile column names
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
                        # if col=='SED_TYPE':print(col, sed_file)
                        # print(col, sed_file)

            for lc_file in lc_lst:
                lc_table = Table.read(lc_file, format='ascii.ecsv')
                # lc_col_lst = list(lc_table.columns)
                lc_col_lst = list(lc_table.meta.keys())
                # if 'SED_TYPE' in lc_col_lst:print(lc_file)
                if 'e_max' in lc_col_lst:print(lc_file)
                for col in lc_col_lst:
                    if col not in lc_cols:
                        lc_cols.append(col)
                        # print(col, lc_file)

    return lc_cols, sed_cols

def get_data_files(data_type, file_type): # should be inside the folder
    data_files = [i for i in sorted(glob.glob("*.%s" % file_type)) if data_type in i]

    if not data_files:return []
    else: return data_files

def merge_fits(data_type, src_id, ads_paper, source_dir):
    # print("Currently in %s/%s" % (src_id, os.path.basename(ads_paper)))

    data_files = get_data_files(data_type, 'ecsv')

    if not data_files:print("No %s files found in VER-%s/%s" % (data_type, src_id, os.path.basename(ads_paper)))

    else:
        main_hdu = get_main_hdu(data_type, src_id, ads_paper, source_dir)
        hdu_list = [main_hdu]

        for data_file in data_files:
            data_hdu = get_data_hdu(data_type, data_file, src_id, ads_paper, source_dir)
            hdu_list.append(data_hdu)

        hdul = fits.HDUList(hdu_list)
        hdul.writeto('VER-%s-%s.fits' % (src_id, data_type), overwrite=True, checksum=True)

    return None

def cleanup_datafiles(data_type, file_type):
    data_files = get_data_files(data_type, file_type)
    if data_files:[os.remove(i) for i in data_files]
    return None

def cleanup_ql_no_data(data_type):
    data_files = sorted(glob.glob("*.ecsv"))
    data_files = [i for i in data_files if 'sed' in i or 'lc' in i]
    if not data_files:
        ql_files = sorted(glob.glob("*.tar.gz"))
        ql_files = [i for i in ql_files if data_type in i]
        c = [os.remove(i) for i in ql_files]
    return None

def compress(data_type, output_file="archive.tar.gz", output_dir='', root_dir='.', items=[]):
    """compress dirs.

    KWArgs
    ------
    output_file : str, default ="archive.tar.gz"
    output_dir : str, default = ''
        absolute path to output
    root_dir='.',
        absolute path to input root dir
    items : list
        list of dirs/items relative to root dir

    """
    # print(items)
    # os.chdir(root_dir)
    if items == []:
        # print("No %s images in %s" % (data_type, root_dir))
        return None
    else:
        with tarfile.open(os.path.join(output_dir, output_file), "w:gz") as tar:
            for item in items:
                tar.add(item, arcname=item)
        return None

def merge_main(data_type, data_dir, source_dir):
    dir_list = sorted([d for d in list_dir(data_dir) if os.path.isdir(d)], key=sort_source)

    for ver_x in dir_list:
        os.chdir(ver_x)
        src_id = int(os.path.basename(ver_x).split('-')[-1])
        paper_dirs = sorted([d for d in list_dir(ver_x) if os.path.isdir(d)])

        for ads_paper in paper_dirs:
            os.chdir(ads_paper)
            merge_fits(data_type, src_id, ads_paper, source_dir)
            cleanup_datafiles(data_type, 'ecsv')

            compress(data_type, output_file="%s-%s_ql.tar.gz" % (os.path.basename(ver_x),data_type), output_dir=ads_paper, root_dir=ads_paper, items=get_data_files(data_type, 'png'))
            cleanup_datafiles(data_type, 'png')
            cleanup_ql_no_data(data_type)


    print('All %s FITS files successfully created!' % data_type)
    # for ads_paper in paper_dirs:
    return None

def make_fits(source_dir, data_dir):
    # lc_cols, sed_cols = get_col_names(data_dir)
    merge_main('lc', data_dir, source_dir)
    merge_main('sed', data_dir, source_dir)
    
    if not 'undetected_data' in data_dir:
        os.remove(data_dir+'VER-146/2014ApJ...783...16A/VER-146-1-extension-table.ecsv')
        os.remove(data_dir+'VER-146/2018ApJ...867L..19A/VER-146-1-spectralFits-table.ecsv')
        os.remove(data_dir+'VER-53/2020ApJ...891..170V/VER-53-1-spectralFits-table.ecsv')
    return None