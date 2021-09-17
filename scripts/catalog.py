#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 09:31:05 2021

@author: sam
"""
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

# =============================================================================
# General functions for loading & processing yaml files
# =============================================================================
def sort_source(name):
    basename = os.path.basename(name).split('.')[0]
    num = basename.split('-')[1]
    return int(num)

def sort_data(name):
    num = name.split('/')[4] # sorted by year
    return int(num)

def get_lists(repo_dir, source_dir, undetected=False):
    source_glob = sorted(glob.glob(source_dir+"VER-*.yaml"), key=sort_source)
    data_glob = sorted(glob.glob(repo_dir+"/20*/*/VER-*.yaml"), key=sort_data)

    if undetected: # although we don't require VER names for undetected sources
        source_lst = [source_glob[i] for i in range(len(source_glob)) if 300000 > float(os.path.basename(source_glob[i]).split('.')[0].split('-')[1]) > 999] # only non-detections which have 300000 > source id > 999

        data_lst = [data_glob[i] for i in range(len(data_glob)) if 300000 > float(os.path.basename(data_glob[i]).split('.')[0].split('-')[1]) > 999] # only 300000 > src_ids > 999

    else:
        source_lst = [source_glob[i] for i in range(len(source_glob)) if float(os.path.basename(source_glob[i]).split('.')[0].split('-')[1]) < 999] # ignore non-detections which have source id > 999

        data_lst = [data_glob[i] for i in range(len(data_glob)) if float(os.path.basename(data_glob[i]).split('.')[0].split('-')[1]) < 999] # remove src_ids > 999

    return source_lst, data_lst

def list_dir(directory): # default os.listdir returns only basepath
    return [os.path.join(directory, file) for file in os.listdir(directory)]

def load_yaml(file):
    stream = open(file, 'r')
    file = yaml.load(stream, Loader=yaml.FullLoader)
    return file

def none_field_dtype(field_name):
    str_fields = ['reference_id', 'description', 'morph_type', 'pa_frame', 'spec_model_type']
    int_fields = ['source_id', 'n_dof']
    float_fields = ['excess',
    'significance',
    'significance_post_trial',
    'livetime',
    'n_on',
    'n_off',
    'alpha',
    'ra',
    'ra_err',
    'ra_sys_err',
    'dec',
    'dec_err',
    'dec_sys_err',
    'glon',
    'glon_err',
    'glon_sys_err',
    'glat',
    'glat_err',
    'glat_sys_err',
    'sigma',
    'sigma_err',
    'sigma_sys_err',
    'sigma2',
    'sigma2_err',
    'sigma2_sys_err',
    'pa',
    'pa_err',
    'pa_sys_err',
    'mjd0',
    'period',
    'period_err',
    'period_err_n',
    'period_err_p',
    'flux',
    'flux_err',
    'flux_sys_err',
    'eflux',
    'flux_ul',
    'conf',
    'e_min',
    'erange_min',
    'erange_max',
    'mjd_min',
    'mjd_max',
    'phase_min',
    'phase_max',
    'theta',
    'spec_norm',
    'spec_norm_err',
    'spec_norm_sys_err',
    'spec_index',
    'spec_index_err',
    'spec_index_sys_err',
    'spec_e_ref',
    'chi2',
    'spec_flux',
    'spec_flux_err',
    'spec_flux_sys_err',
    'spec_e_min',
    'spec_e_max',
    'spec_e_cut',
    'spec_e_cut_err',
    'spec_e_cut_sys_err',
    'spec_alpha',
    'spec_alpha_err',
    'spec_alpha_sys_err',
    'spec_beta',
    'spec_beta_err',
    'spec_beta_sys_err']

    if field_name in str_fields:
        return np.str_(np.ma.masked)
    elif field_name in int_fields:
        return np.int64(np.ma.masked)
    elif field_name in float_fields:
        return np.float64(np.ma.masked)
    return np.ma.masked

def get_field(df, field_name, field_id):
    try:
        # print(field_name)
        if type(df.iloc[0][field_name]) is list and len(df.iloc[0][field_name])==0:
            return none_field_dtype(field_name) # Needed for empty lists in source yaml files. This is really annoying!
        return df.iloc[0][field_name]
    except KeyError:
        # return np.ma.masked
        return none_field_dtype(field_id)
    return None

def get_unit(df, unit_field):
    try:
        return u.Unit(df.iloc[0][unit_field].replace(" ","."), format="cds") # returns CDS units
    except KeyError:
        # return 'Problem - No Units!'
        return None # gets converted to dimensionless units
    return None

def process(df, fields_tup):
    data_dict = {}
    for i in range(len(fields_tup)):
        data_dict[fields_tup[i][0]] = get_field(df,fields_tup[i][1],fields_tup[i][0])

        if type(data_dict[fields_tup[i][0]]) is list:
            data_dict[fields_tup[i][0]] = ",".join(data_dict[fields_tup[i][0]]) # convert list to 1 comma separated string
    return data_dict

def get_src_list(data_dir):
    # os.chdir(data_dir)
    directories = sorted([d for d in list_dir(data_dir) if os.path.isdir(d)])
    src_list = sorted([int(os.path.basename(i).split("-")[1]) for i in directories])
    filenames = ["VER-%d.yaml" % i for i in src_list]
    dirnames = ["VER-%d" % i for i in src_list]
    return filenames, dirnames, src_list

def get_paper_dir(dir_name): # returns a list of paper directories in each directory
    paper_directories = sorted([d for d in list_dir(dir_name) if os.path.isdir(d)])
    return paper_directories

def get_yaml_list(dir_name): # returns a list of paper directories in each directory
    yaml_lst = sorted(glob.glob(dir_name+"/*.yaml"))
    return yaml_lst
# =============================================================================
# Functions for Creating & Writing Tables
# =============================================================================

def create_table(merged_dict):
    t = QTable(rows=list(merged_dict), masked=True)
    return t

def write_tables(t, name):
    ascii.write(t, name + ".ecsv", format='ecsv', fast_writer=True, overwrite=True)
    t.write(name + ".html", format='jsviewer')
    # t.write(name + ".fits", format='fits', overwrite=True)
    return None

# =============================================================================
# Functions for Source Catalog
# =============================================================================

def get_obs_count(data_dir, dirname):
    print(dirname)
    return len(next(os.walk(data_dir + dirname))[1])

def create_src_dict(repo_dir, source_dir, df, src_cnt, undetected=False): # create a source catalog dictionary for a single entry
    src_tup = [('source_id', 'source_id'),
               ('veritas_name', 'veritas_name_name'),
               ('component', 'veritas_name_components'),
               ('n_obs', 'n_obs'),
               ('common_name', 'common_name'),
               ('simbad_id', 'pos_simbad_id'),
               ('other_names', 'other_names'),
               ('location_flag', 'where'),
               ('simbad_ra', 'pos_ra'),
               ('simbad_dec', 'pos_dec'),
               ('ref_sources', 'reference_id')] # for now, assume they are only VERITAS papers

    src_dict = process(df, src_tup)

    if src_dict['veritas_name'] is None: src_dict['veritas_name'] = np.ma.masked
    if src_dict['component'] is None: src_dict['component'] = np.ma.masked

    src_dict['n_obs'] = src_cnt

    if src_dict['other_names']: # this is to remove duplicates in other_names
        other_names = src_dict['other_names'].split(',')
        for other_name in other_names:
            if other_name == src_dict['common_name']:
                other_names.remove(other_name)
        src_dict['other_names'] =  ','.join(other_names)

    src_dict['simbad_ra'] = src_dict['simbad_ra']*u.deg
    src_dict['simbad_dec'] = src_dict['simbad_dec']*u.deg
    if src_dict['simbad_ra'] == 0.0: # so that entries don't have Simbad RA & Dec get converted gracefully
        src_dict['simbad_ra'] = src_dict['simbad_dec'] = np.nan*u.deg
    return src_dict

def create_src_cat(repo_dir, source_dir, data_dir):
    if 'undetected_data' in data_dir:undetected=True
    else:undetected=False

    src_files, data_files = get_lists(repo_dir, source_dir, undetected)
    lst_dicts = []

    for file in src_files:
        data = load_yaml(file)
        df = pd.json_normalize(data, sep='_')
        src_cnt = get_obs_count(data_dir,os.path.basename(file).split(".")[0])
        src_dict = create_src_dict(repo_dir, source_dir, df, src_cnt, undetected)
        lst_dicts.append(src_dict)
    t = create_table(lst_dicts)
    write_tables(t, data_dir+'src_cat')
    print("Created source catalog table sucessfully!")
    return None

# =============================================================================
# Units and Fixes
# Note: The unit conversion will need to be done at this level since masked
# entries get converted to 0 when multiplied by astropy units
# =============================================================================

def fix_livetime(val): # default units in hours
    if pd.isnull(val): # masked value workaround since val*u.hour = 0 and [val]*u.hour = [NaN]
        nan = [val]*u.hour
        fixed_val = nan.value[0]*nan.unit
        return fixed_val

    if 'hour' in val:
        return float(val.replace('hour',''))*u.hour
    elif 'h' in val:
        return float(val.replace('h',''))*u.hour
    elif 'min' in val:
        return float(float(val.replace('min',''))/60.0)*u.hour
    elif 'm' in val:
        return float(float(val.replace('m',''))/60.0)*u.hour
    elif 'sec' in val:
        return float(float(val.replace('sec',''))/3600.0)*u.hour
    elif 's' in val:
        return float(float(val.replace('s',''))/3600.0)*u.hour
    # elif val is np.ma.masked: # masked value workaround since val*u.hour = 0 and [val]*u.hour = [NaN]

    return 'Incorrect value for livetime'

def fix_coords(val, *args): # default units in degrees; works for ra, dec, glon, glat, sigma, sigma2, pa

    if pd.isnull(val):
        nan = [val]*u.degree
        fixed_val = nan.value[0]*nan.unit
        return fixed_val
    elif args[0] == 'ra' and 'h' not in val and 'm' in val:
        return Angle('0h' + str(val)).degree
    elif args[0] == 'dec' and 'd' not in val and 'm' in val:
        return Angle('0d' + str(val)).degree
    return Angle(val).degree*u.degree

def fix_orbit(val, *args): # default units in days
    if args:
        if args[0] == 'mjd0' and not pd.isnull(val) and 'd' in str(val):
            return float(val.replace('d',''))*u.day
    if pd.isnull(val):
        nan = [val]*u.day
        fixed_val = nan.value[0]*nan.unit
        return fixed_val
    return val*u.day

def fix_flux(param, val, unit_str): # default units specified
    try:

        if param == 'flux':
            def_unit = u.Unit("cm-2.s-1", format="cds")
        elif param == 'eflux':
            def_unit = u.Unit("erg.cm-2.s-1", format="cds")
        elif param == 'energy':
            def_unit = u.Unit("TeV", format="cds")
            
        if pd.isnull(val) or unit_str is None:
            nan = [val]*def_unit
            fixed_val = nan.value[0]*nan.unit
            return fixed_val
        else:
            value = val * unit_str
            return value.to_value(def_unit)*def_unit
    except AttributeError:
        return "Problem in fix_flux"
        # return np.ma.masked
    return None

def fix_spec(param, val, unit_str, *scale): # default units specified
    try:
        if param == 'theta':
            if not pd.isnull(val):val = float(val.replace("d", ""))
            def_unit = u.Unit("degree")
        elif param == 'norm': # norm needs a scale
            if not pd.isnull(val):val = val*float(scale[0])
            def_unit = u.Unit("cm-2.s-1.TeV-1", format="cds")
        elif param == 'energy':
            def_unit = u.Unit("TeV", format="cds")
        elif param == 'flux': # flux needs a scale
            if not pd.isnull(val):val = val*float(scale[0])
            def_unit = u.Unit("cm-2.s-1", format="cds")

        if pd.isnull(val) or unit_str is None:
            nan = [val]*def_unit
            fixed_val = nan.value[0]*nan.unit
            return fixed_val
        else:
            value = val * unit_str
            return value.to_value(def_unit)*def_unit
    except AttributeError:
        return "Problem in fix_spec"
    return None

# =============================================================================
# Functions for converting data yaml files to ecsv inside each VER-x/* folder
# =============================================================================
def create_id_dict(df):
    id_tup = [('source_id', 'source_id'),
              ('reference_id', 'reference_id'),
              ('description', 'description')]
    id_dict = process(df, id_tup)
    return OrderedDict(id_dict)

def create_obs_dict(df):
    obs_tup = [('excess', 'data_excess'),
              ('significance', 'data_significance'),
              ('significance_post_trial', 'data_significance_post_trial'),
              ('livetime', 'data_livetime'),
              ('n_on', 'data_non'),
              ('n_off', 'data_noff'),
              ('alpha', 'data_alpha')]
    obs_dict = process(df, obs_tup)

    obs_dict['livetime'] = fix_livetime(obs_dict['livetime'])
    return OrderedDict(obs_dict)

def create_coords_dict(df):
    coords_tup = [('ra', 'pos_ra_val'),
                  ('ra_err', 'pos_ra_err'),
                  ('ra_sys_err', 'pos_ra_err_sys'),
                  ('dec', 'pos_dec_val'),
                  ('dec_err', 'pos_dec_err'),
                  ('dec_sys_err', 'pos_dec_err_sys'),
                  ('glon', 'pos_glon_val'),
                  ('glon_err', 'pos_glon_err'),
                  ('glon_sys_err', 'pos_glon_err_sys'),
                  ('glat', 'pos_glat_val'),
                  ('glat_err', 'pos_glat_err'),
                  ('glat_sys_err', 'pos_glat_err_sys')]
    coords_dict = process(df, coords_tup)

    coords_dict['ra'] = fix_coords(coords_dict['ra'],'ra')
    coords_dict['ra_err'] = fix_coords(coords_dict['ra_err'],'ra')
    coords_dict['ra_sys_err'] = fix_coords(coords_dict['ra_sys_err'],'ra')

    coords_dict['dec'] = fix_coords(coords_dict['dec'],'dec')
    coords_dict['dec_err'] = fix_coords(coords_dict['dec_err'],'dec')
    coords_dict['dec_sys_err'] = fix_coords(coords_dict['dec_sys_err'],'dec')

    coords_dict['glon'] = fix_coords(coords_dict['glon'],'glon')
    coords_dict['glon_err'] = fix_coords(coords_dict['glon_err'],'glon')
    coords_dict['glon_sys_err'] = fix_coords(coords_dict['glon_sys_err'],'glon')

    coords_dict['glat'] = fix_coords(coords_dict['glat'],'glat')
    coords_dict['glat_err'] = fix_coords(coords_dict['glat_err'],'glat')
    coords_dict['glat_sys_err'] = fix_coords(coords_dict['glat_sys_err'],'glat')
    return OrderedDict(coords_dict)

def create_morph_dict(df):
    morph_tup = [('sigma', 'morph_sigma_val'),
                  ('sigma_err', 'morph_sigma_err'),
                  ('sigma_sys_err', 'morph_sigma_err_sys'),
                  ('sigma2', 'morph_sigma2_val'),
                  ('sigma2_err', 'morph_sigma2_err'),
                  ('sigma2_sys_err', 'morph_sigma2_err_sys'),
                  ('pa', 'morph_pa_val'),
                  ('pa_err', 'morph_pa_err'),
                  ('pa_sys_err', 'morph_pa_err_sys'),
                  ('morph_type', 'morph_type'),
                  ('pa_frame', 'morph_pa_frame')]
    morph_dict = process(df, morph_tup)

    morph_dict['sigma'] = fix_coords(morph_dict['sigma'], 'sigma')
    morph_dict['sigma_err'] = fix_coords(morph_dict['sigma_err'], 'sigma')
    morph_dict['sigma_sys_err'] = fix_coords(morph_dict['sigma_sys_err'], 'sigma')
    morph_dict['sigma2'] = fix_coords(morph_dict['sigma2'], 'sigma')
    morph_dict['sigma2_err'] = fix_coords(morph_dict['sigma2_err'], 'sigma')
    morph_dict['sigma2_sys_err'] = fix_coords(morph_dict['sigma2_sys_err'], 'sigma')
    morph_dict['pa'] = fix_coords(morph_dict['pa'], 'pa')
    morph_dict['pa_err'] = fix_coords(morph_dict['pa_err'], 'pa')
    morph_dict['pa_sys_err'] = fix_coords(morph_dict['pa_sys_err'], 'pa')
    return OrderedDict(morph_dict)

def create_orbit_dict(df):
    orbit_tup = [('mjd0', 'orbit_mjd0'),
                  ('period', 'orbit_period_val'),
                  ('period_err', 'orbit_period_err'),
                  ('period_err_n', 'orbit_period_errn'),
                  ('period_err_p', 'orbit_period_errp')]
    orbit_dict = process(df, orbit_tup)

    orbit_dict['mjd0'] = fix_orbit(orbit_dict['mjd0'],'mjd0')
    orbit_dict['period'] = fix_orbit(orbit_dict['period'])
    orbit_dict['period_err'] = fix_orbit(orbit_dict['period_err'])
    orbit_dict['period_err_n'] = fix_orbit(orbit_dict['period_err_n'])
    orbit_dict['period_err_p'] = fix_orbit(orbit_dict['period_err_p'])

    return OrderedDict(orbit_dict)

def create_flux_dict(df):
    flux_tup = [('flux', 'flux_flux_val'),
                  ('flux_err', 'flux_flux_err'),
                  ('flux_sys_err', 'flux_flux_err_sys'),
                  ('eflux', 'flux_eflux_val'),
                  ('flux_ul', 'flux_flux_ul_val'),
                  ('conf', 'flux_conf_val'),
                  ('e_min', 'flux_e_min_val')]
    flux_dict = process(df, flux_tup)

    flux_dict['flux'] = fix_flux('flux', flux_dict['flux'], get_unit(df, "flux_flux_unit"))
    flux_dict['flux_err'] = fix_flux('flux', flux_dict['flux_err'], get_unit(df, "flux_flux_unit"))
    flux_dict['flux_sys_err'] = fix_flux('flux', flux_dict['flux_sys_err'], get_unit(df, "flux_flux_unit"))
    flux_dict['flux_ul'] = fix_flux('flux', flux_dict['flux_ul'], get_unit(df, "flux_flux_ul_unit"))

    flux_dict['eflux'] = fix_flux('eflux', flux_dict['eflux'], get_unit(df, "flux_eflux_unit"))

    flux_dict['e_min'] = fix_flux('energy', flux_dict['e_min'], get_unit(df, "flux_e_min_unit"))
    return OrderedDict(flux_dict)

def create_spec_dict(df):
    spec_tup = [('erange_min', 'spec_erange_min'),
                ('erange_max', 'spec_erange_max'),
                  ('mjd_min', 'spec_mjd_min'),
                  ('mjd_max', 'spec_mjd_max'),
                  ('phase_min', 'spec_phase_min'),
                  ('phase_max', 'spec_phase_max'),
                  ('theta', 'spec_theta'),
                  ('spec_model_type', 'spec_model_type'),
                  ('spec_norm', 'spec_model_parameters_norm_val'),
                  ('spec_norm_err', 'spec_model_parameters_norm_err'),
                  ('spec_norm_sys_err', 'spec_model_parameters_norm_err_sys'),
                  ('spec_index', 'spec_model_parameters_index_val'),
                  ('spec_index_err', 'spec_model_parameters_index_err'),
                  ('spec_index_sys_err', 'spec_model_parameters_index_err_sys'),
                  ('spec_e_ref', 'spec_model_parameters_e_ref_val'),
                  ('chi2', 'spec_model_parameters_chi2_val'), # changed 23/6/21
                  ('n_dof', 'spec_model_parameters_n_dof_val'), # changed 23/6/21
                  ('spec_flux', 'spec_model_parameters_flux_val'),
                  ('spec_flux_err', 'spec_model_parameters_flux_err'),
                  ('spec_flux_sys_err', 'spec_model_parameters_flux_err_sys'),
                  ('spec_e_min', 'spec_model_parameters_e_min_val'),
                  ('spec_e_max', 'spec_model_parameters_e_max_val'),
                  ('spec_e_cut', 'spec_model_parameters_e_cut_val'),
                  ('spec_e_cut_err', 'spec_model_parameters_e_cut_err'),
                  ('spec_e_cut_sys_err', 'spec_model_parameters_e_cut_err_sys'),
                  ('spec_alpha', 'spec_model_parameters_alpha_val'),
                  ('spec_alpha_err', 'spec_model_parameters_alpha_err'),
                  ('spec_alpha_sys_err', 'spec_model_parameters_alpha_err_sys'),
                  ('spec_beta', 'spec_model_parameters_beta_val'),
                  ('spec_beta_err', 'spec_model_parameters_beta_err'),
                  ('spec_beta_sys_err', 'spec_model_parameters_beta_err_sys')]
    spec_dict = process(df, spec_tup)

    spec_dict['erange_min'] = fix_spec('energy', spec_dict['erange_min'], get_unit(df, "spec_erange_unit"))
    spec_dict['erange_max'] = fix_spec('energy', spec_dict['erange_max'], get_unit(df, "spec_erange_unit"))

    spec_dict['theta'] = fix_spec('theta', spec_dict['theta'], u.degree)

    spec_dict['spec_norm'] = fix_spec('norm', spec_dict['spec_norm'], get_unit(df, "spec_model_pl_norm_unit"), get_field(df, "spec_model_pl_norm_scale", "spec_norm"))
    spec_dict['spec_norm_err'] = fix_spec('norm', spec_dict['spec_norm_err'], get_unit(df, "spec_model_pl_norm_unit"), get_field(df, "spec_model_pl_norm_scale", "spec_norm_err"))
    spec_dict['spec_norm_sys_err'] = fix_spec('norm', spec_dict['spec_norm_sys_err'], get_unit(df, "spec_model_pl_norm_unit"), get_field(df, "spec_model_pl_norm_scale", "spec_norm_sys_err"))

    spec_dict['spec_e_ref'] = fix_spec('energy', spec_dict['spec_e_ref'], get_unit(df, "spec_model_parameters_e_ref"))

    spec_dict['spec_flux'] = fix_spec('flux', spec_dict['spec_flux'], get_unit(df, "spec_model_parameters_flux_unit"), get_field(df, "spec_model_parameters_flux_scale", "spec_flux"))
    spec_dict['spec_flux_err'] = fix_spec('flux', spec_dict['spec_flux_err'], get_unit(df, "spec_model_parameters_flux_unit"), get_field(df, "spec_model_parameters_flux_scale", "spec_flux_err"))
    spec_dict['spec_flux_sys_err'] = fix_spec('flux', spec_dict['spec_flux_sys_err'], get_unit(df, "spec_model_parameters_flux_unit"), get_field(df, "spec_model_parameters_flux_scale", "spec_flux_sys_err"))

    spec_dict['spec_e_min'] = fix_spec('energy', spec_dict['spec_e_min'], get_unit(df, "spec_model_parameters_e_min_unit"))
    spec_dict['spec_e_max'] = fix_spec('energy', spec_dict['spec_e_max'], get_unit(df, "spec_model_parameters_e_max_unit"))

    spec_dict['spec_e_cut'] = fix_spec('energy', spec_dict['spec_e_cut'], get_unit(df, "spec_model_parameters_e_cut_unit"))
    spec_dict['spec_e_cut_err'] = fix_spec('energy', spec_dict['spec_e_cut_err'], get_unit(df, "spec_model_parameters_e_cut_unit"))
    spec_dict['spec_e_cut_sys_err'] = fix_spec('energy', spec_dict['spec_e_cut_sys_err'], get_unit(df, "spec_model_parameters_e_cut_unit"))

    return OrderedDict(spec_dict)

def paper_dir_cat(yaml_lst):
    '''
    Creates (for a single paper_dir), an obs_data dictionary

    '''
    lst_dicts = []

    for yaml_file in yaml_lst:
        data = load_yaml(yaml_file)
        df = pd.json_normalize(data, sep='_')

        id_dict = create_id_dict(df)
        merged_dict = id_dict

        obs_dict = create_obs_dict(df)
        merged_dict.update(obs_dict)

        coords_dict = create_coords_dict(df)
        merged_dict.update(coords_dict)

        morph_dict = create_morph_dict(df)
        merged_dict.update(morph_dict)

        orbit_dict = create_orbit_dict(df)
        merged_dict.update(orbit_dict)

        flux_dict = create_flux_dict(df)
        merged_dict.update(flux_dict)

        spec_dict = create_spec_dict(df)
        merged_dict.update(spec_dict)

        lst_dicts.append(merged_dict)

    return lst_dicts

def create_data_cat(data_dir):
    '''
    Creates a data_catalog in each paper_dir and removes the data yaml files
    '''
    # filenames, dirnames, src_list = get_src_list()
    data_directories = sorted([d for d in list_dir(data_dir) if os.path.isdir(d)], key=sort_source)

    for src_dir in data_directories:
        paper_dir_lst = get_paper_dir(src_dir) # contains paper directories in src_dir
        for paper_dir in paper_dir_lst:
            print(paper_dir.split('/')[-2]+"/"+paper_dir.split('/')[-1])
            yaml_lst = get_yaml_list(paper_dir) # contains yaml files in paper_dir
            if not yaml_lst: # no data tables to create
                continue
            lst_dict = paper_dir_cat(yaml_lst)
            t = create_table(lst_dict)
            name_no_ext = os.path.basename(yaml_lst[0]).split(".")[0]
            src_id = int(name_no_ext.split("-")[1])
            table_name = "VER-" + str(src_id) + "-obs_data"
            write_tables(t, paper_dir+"/"+table_name)
            for yaml_file in yaml_lst:
                os.remove(yaml_file)
    return None

def build_catalog(repo_dir, source_dir, data_dir):

    create_src_cat(repo_dir, source_dir, data_dir)
    print("Succesfully built source catalog!")
    create_data_cat(data_dir)
    print("Succesfully built data catalogs!")
    
    return None
