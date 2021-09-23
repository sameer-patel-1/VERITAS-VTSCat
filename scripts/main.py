import restructure
import ver_name
import catalog
import create_fits

import os
import git
import glob
import yaml
import shutil

undetected=False # set true to include undetected sources.

#NB: undetected sources don't play well with the scripts because of different naming convention of data files. So, don't use it (ie. set it to False always; we're not providing HEASARC with undetected VERITAS data anyway)

repo = git.Repo(".", search_parent_directories=True)
repo_dir = repo.working_tree_dir + "/" # establish pwd as the git repo
heasarc_dir = repo_dir+"heasarc/" # base dir for heasarc files/folders
if not os.path.isdir(heasarc_dir):os.mkdir(heasarc_dir) # if folder doesn't exist, create it!
os.chdir(repo_dir+'scripts')

restructure.restructure_to_heasarc(repo_dir, heasarc_dir)
source_dir = heasarc_dir+"sources/"
# ver_names = ver_name.get_ver_names(repo_dir, source_dir) # no longer needed

data_dirs = [heasarc_dir+'detected_data/']
if undetected:data_dirs.append(heasarc_dir+'undetected_data/')

for data_dir in data_dirs:
    # lc_cols, sed_cols = create_fits.get_col_names(data_dir)

    catalog.build_catalog(repo_dir, source_dir, data_dir)

    create_fits.make_fits(source_dir, data_dir)

    print("Everything is done. Don\'t forget to compress the detected data and send it to HEASARC!")