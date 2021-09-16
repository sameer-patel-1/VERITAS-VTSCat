# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 18:19:41 2019

@author: Sam
"""
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import yaml
from pathlib import Path
import glob
import git
import os
from astroquery.simbad import Simbad
from collections import Counter

def get_lists(repo_dir, source_dir, undetected=False):
        source_glob = glob.glob(source_dir+"VER-*.yaml")
        data_glob = sorted(glob.glob(repo_dir+"/20*/*/VER-*.yaml"), key=len)

        if undetected: # although we don't require VER names for undetected sources
            source_lst = sorted([source_glob[i] for i in range(len(source_glob)) if 300000 > float(os.path.basename(source_glob[i]).split('.')[0].split('-')[1]) > 999], key=len) # only non-detections which have 300000 > source id > 999

            data_lst = sorted([data_glob[i] for i in range(len(data_glob)) if 300000 > float(os.path.basename(data_glob[i]).split('.')[0].split('-')[1]) > 999], key=len) # only 300000 > src_ids > 999

        else:
            source_lst = sorted([source_glob[i] for i in range(len(source_glob)) if float(os.path.basename(source_glob[i]).split('.')[0].split('-')[1]) < 999], key=len) # ignore non-detections which have source id > 999

            data_lst = sorted([data_glob[i] for i in range(len(data_glob)) if float(os.path.basename(data_glob[i]).split('.')[0].split('-')[1]) < 999], key=len) # remove src_ids > 999

        return source_lst, data_lst


def sort_dict(dict, reverse=False):
    """Sorts a dictionary in order of its keys

    Args:
        dict ([type]): [description]

    Returns:
        [type]: [description]
    """
    d = {}
    for key in sorted(dict.keys()):
        if reverse:
            d.update({key:sorted(dict[key],key=len,reverse=True)})
        else:
            d.update({key:dict[key]})
    return d

def create_data_dict(source_lst, data_lst):
    """Creates a dict with src_id:list(all data yaml files) for src_id<999

    Args:
        source_lst ([type]): [description]
        data_lst ([type]): [description]

    Returns:
        [type]: [description]
    """
    data_dict = {}

    for entry in data_lst:
        src_id = int(os.path.basename(entry).split('.')[0].split('-')[1])
        src_id_lst = [entry]
        for src_entry in data_lst:
            src_entry_id = int(os.path.basename(src_entry).split('.')[0].split('-')[1])
            if src_id == src_entry_id:
                if src_entry not in src_id_lst:
                    src_id_lst.append(src_entry)

        data_dict.update({src_id:src_id_lst}) #dict with file path:src_id (src_id contains duplicates)


    for source in source_lst: # use source yaml files for missing data yaml files
        src_id = int(os.path.basename(source).split('.')[0].split('-')[1])
        if src_id in list(data_dict.keys()):
            data_dict.update({src_id:[source]+sorted(data_dict[src_id])})

        elif src_id not in list(data_dict.keys()):
            data_dict.update({src_id:[source]})

    data_dict = sort_dict(data_dict)

    return data_dict

def load_yaml(file):
    stream = open(file, 'r')
    file = yaml.load(stream, Loader=yaml.FullLoader)
    return file

def sign(dec_dms):
    if dec_dms<0:
        return '-'
    else:
        return '+'

def get_source_name(data):
    try:
        gamma_names = data['gamma_names']
        for name in gamma_names:
            if 'VER' in name:
                return name
    except:
        return None

def get_simbad_cat_name(data):
    try:
        # common_name = data['common_name']
        simbad_id = data['pos']['simbad_id']
        obj_ids = Simbad.query_objectids(simbad_id)
        for name in range(len(obj_ids)):
            if 'VER' in obj_ids[name][0]:
                name = obj_ids[name][0]
                return name
    except:
        return None

def get_simbad_gen_name(data):
    try:
        c = SkyCoord(data['pos']['ra'], data['pos']['dec'], unit = "deg", frame='icrs')
        simbad_gen_name = 'VER J'+str("%.2d" % c.ra.hms[0])+str("%.2d" % c.ra.hms[1])+sign(c.dec.dms[0])+str("%.2d" % abs(c.dec.dms[0]))+str(int(np.floor(abs(c.dec.dms[1]/60*10))))
        return simbad_gen_name
    except:
        return None

def get_ver_coord_gen_name(data):
    try:
        c = SkyCoord(data['pos']['ra']['val'], data['pos']['dec']['val'], unit = "deg", frame='icrs')
        ver_gen_name = 'VER J'+str("%.2d" % c.ra.hms[0])+str("%.2d" % c.ra.hms[1])+sign(c.dec.dms[0])+str("%.2d" % abs(c.dec.dms[0]))+str(int(np.floor(abs(c.dec.dms[1]/60*10))))
        return ver_gen_name
    except:
        return None

def merge_simbad(data_dict):
    ver_name_dict = {}
    data_dict_upd = {}

    for source in data_dict:
        for data_entry in data_dict[source]:
            data = load_yaml(data_entry)
            source_name = get_source_name(data)
            if bool(source_name) is True: # found ver_name from source file
                ver_name_dict.update({source:source_name})
                data_dict_upd.update({source:data_entry})
                break
            elif bool(source_name) is False: # could not find ver_name from source file -> check Simbad catalog
                simbad_cat_name = get_simbad_cat_name(data)
                if bool(simbad_cat_name) is True:
                    ver_name_dict.update({source:simbad_cat_name})
                    data_dict_upd.update({source:data_entry})
                    break

    all(map(data_dict.pop, data_dict_upd)) # now data_dict only contains src_ids for which ver_names don't exist
    return data_dict, ver_name_dict

def generate_names(data_dict, ver_name_dict):

    data_dict = sort_dict(data_dict, reverse=True)

    for source in data_dict:
        src_file = data_dict[source][-1]
        if len(data_dict[source]) == 1:
            data = load_yaml(data_dict[source][0])
            simbad_gen_name = get_simbad_gen_name(data)
            ver_name_dict.update({source:simbad_gen_name})
        else:
            data_dict[source].pop()
            for data_entry in data_dict[source]:
                data = load_yaml(data_entry)
                ver_gen_name = get_ver_coord_gen_name(data)
                if bool(ver_gen_name) is True:
                    ver_name_dict.update({source:ver_gen_name})
                    break
                else:
                    data = load_yaml(src_file)
                    simbad_gen_name = get_simbad_gen_name(data)
                    ver_name_dict.update({source:simbad_gen_name})

    ver_names = sort_dict(ver_name_dict)
    return ver_names

def chk_names(source_lst, ver_name_dict):
    for source in source_lst:
        src_id = int(os.path.basename(source).split('.')[0].split('-')[1])
        src_name = load_yaml(source)['veritas_name']['name']
        gen_name = ver_name_dict[src_id]
        if src_name == gen_name: # everything is in order
            pass
        else:
            textstr = '\n'.join((
                    "Naming inconsistency found in src_id = %d" % src_id,
                    "Generated name = %s" % gen_name,
                    "Source name = %s" % src_name,
                    "Setting the generated name to source name"))
            print(textstr+"\n")

            ver_name_dict.update({src_id:src_name})
    return ver_name_dict

def get_ver_names(repo_dir, source_dir):

    src_lst, data_lst = get_lists(repo_dir, source_dir)
    data_dict = create_data_dict(src_lst, data_lst)
    merged_data_dict, ver_name_dict = merge_simbad(data_dict)
    ver_names = generate_names(merged_data_dict, ver_name_dict)
    ver_names = chk_names(src_lst, ver_names)
    print("Generated VERITAS source names successfully!")
    return ver_names

def gen_ver_names():
    repo = git.Repo(".", search_parent_directories=True)
    repo_dir = repo.working_tree_dir + "/" # establish pwd as the git repo
    heasarc_dir = repo_dir+"heasarc/" # base dir for heasarc files/folders
    source_dir = heasarc_dir+"sources/"

    ver_names = get_ver_names(repo_dir, source_dir)
    return ver_names