import os
import git
import glob
import yaml
import shutil
import pandas as pd
from astropy.io import ascii
from astropy.table import QTable, Table

def load_yaml(file):
    stream = open(file, 'r')
    file = yaml.load(stream, Loader=yaml.FullLoader)
    return file

def list_dir(directory): # default os.listdir returns only basepath
    return [os.path.join(directory, file) for file in os.listdir(directory)]

def make_copy(repo_dir, heasarc_dir): # make a copy of data directories & source directory to heasarc

    unwanted_dirs = ['.git', '.ipynb_checkpoints', 'heasarc', 'scripts'] # unwanted directories in repo_dir
    directories=sorted([d for d in list_dir(repo_dir) if os.path.isdir(d) and (os.path.basename(d) not in unwanted_dirs)])

    if list_dir(heasarc_dir): # /heasarc directory is not empty
        existing_directories=sorted([d for d in list_dir(heasarc_dir) if os.path.isdir(d)])
        a = [shutil.rmtree(d) for d in existing_directories] # remove all existing directories first
        b = [shutil.copytree(d, heasarc_dir+os.path.basename(d)) for d in directories]
    else:
        c = [shutil.copytree(d, heasarc_dir+os.path.basename(d)) for d in directories]
    return None

def rename_sources(heasarc_dir):
    source_dir = heasarc_dir+'sources/'
    yaml_files = sorted(glob.glob(source_dir + "*.yaml"))

    for file in yaml_files:
        fnm = os.path.basename(file).split('.')[0]
        split_lst = fnm.split("-")
        split_lst[1] = split_lst[1].split(".")[0]
        fnm = 'VER-' + str(int(split_lst[1])) + ".yaml"
        shutil.move(file,source_dir+fnm)
    return None

def rem_paper(heasarc_dir): # remove specific paper directory or subdirectory
    paths = []
    paths.append('2021/2021arXiv210601386A') # contains non-standard data; remove for now
    paths.append('2017/2017PhRvD..95h2001A') # contains non data for designated sources in info.yaml; remove for now
    # paths
    # paths.append('2011/2011ApJ...738..25A/data') # no longer required
    # paths.append('2013/2013ApJ...779..150A/data') # no longer required
    # paths.append('2016/2016ApJ...819..156B/data') # no longer required
    # paths.append('2017/2017ICRC...35..712K/data') # no longer required
    a = [shutil.rmtree(heasarc_dir+dir, ignore_errors=True) for dir in paths]
    return None

def isint(fnm):
    '''
    Return True if fnm can be coverted to an integer else return False
    '''
    try:
        int(fnm)
        return True
    except ValueError:
        return False
    return None

def rem_non_ver_files():
    '''
    Remove files which don't contain 'VER' in them; but keep fits files since they may contain veritas skymaps
    '''
    files = sorted(glob.glob("./**/*.*", recursive=True))
    non_ver_files = [i for i in files if 'VER-' not in i]
    non_ver_files = [i for i in non_ver_files if 'info.yaml' not in i] # don't remove info.yaml files
    non_ver_files = [i for i in non_ver_files if '.fits' not in i or '.FITS' not in i] # don't remove fits files

    a = [os.remove(file) for file in non_ver_files if non_ver_files]

    return None

def rem_non_skymap_fits():
    '''
    Remove fits files which are not skymaps
    '''
    files = sorted(glob.glob("./**/*.fits", recursive=True))
    non_skymaps = [i for i in files if 'skymap' not in i]

    a = [os.remove(file) for file in non_skymaps if non_skymaps]
    return None

def rem_nosrc_ecsv():
    '''
    Remove ecsv files which don't have src_ids'. Should only be run after removing non_ver files!
    '''
    files = sorted(glob.glob("./**/*.ecsv", recursive=True))

    if files:
        for file in files:
            split_lst = file.split("-")
            try:
                int(split_lst[1])
            except ValueError:
                os.remove(file)
    return None

def rem_non_ver_data():
    '''
    Remove datafiles (LC/SED) which contain data from multiple telescopes in them since HEASARC needs only VERITAS data
    '''
    files = sorted(glob.glob("*.ecsv"))
    files = [i for i in files if 'sed' in i or 'lc' in i]

    if files:
        for file in files:
            data = Table.read(file, format='ascii.ecsv')
            meta = pd.json_normalize(data.meta, sep='_')
            try:
                telescope = meta['telescope'][0]
                if len(telescope.split(',')) > 1: # meaning data from more than one telescope
                    os.remove(file)
            except KeyError:
                pass
    return None

def rename_figures(paper_dir):
    '''
    Removes the paper name prefix from figures
    Converts %06d to %0d int and renames figures
    Note: Gernot has already added suffixes in figures which are the observation number
    Note: It removes the pesky arxiv figures in 2017ICRC...3S..729M since the published paper files are already present
    '''
    paper_name = os.path.basename(paper_dir)

    png_files = sorted(glob.glob("./figures/*.png")) # round 1 - take care of paper name prefixes
    if png_files:
        for file in png_files:
            if 'arXiv' in file:
                os.remove(file)
                continue
            if paper_name in os.path.basename(file):
                new_fnm = os.path.basename(file).replace(paper_name+"-","") # remove the paper name prefix
                shutil.move(file, './figures/'+new_fnm)

    png_files = sorted(glob.glob("./figures/*.png")) # round 2 - rename
    if png_files:
        for file in png_files:
            split_lst = os.path.basename(file).split("-")
            shutil.move(file, './figures/'+os.path.basename(file).replace(split_lst[1],str(int(split_lst[1]))))
    return None

def rename_yaml():
    '''
    Converts %06d to %0d int and renames data yaml files
    Note: It also adds a "-1" at the end of the file name if it is the only one in the folder
    '''
    yaml_files = sorted(glob.glob("*.yaml"))
    yaml_files = [i for i in yaml_files if 'VER' in i]

    for file in yaml_files:
        split_lst = file.split("-")
        if '.yaml' in split_lst[1]:
            split_lst[1] = split_lst[1].replace(".yaml", "")
            split_lst.append(".yaml")
            shutil.move(file, file.replace(split_lst[1],str(int(split_lst[1]))+"-1"))
        else:
            shutil.move(file, file.replace(split_lst[1],str(int(split_lst[1]))))
    return None

def rename_skymap():
    fits_files = sorted(glob.glob("*.fits"))
    fits_files = [i for i in fits_files if 'VER' in i]

    for file in fits_files:
        fnm = file.split('.')[0]
        split_lst = fnm.split("-")
        split_lst[1] = split_lst[1].split(".")[0]
        shutil.move(file, file.replace(split_lst[1],str(int(split_lst[1]))))
    return None

def change_ecsv_format():
    '''
    This changes from VER-X-Z-Y to VER-X-Y-Z; where X = src_id, Y = obs_id & Z is obs_type (SED or LC)
    '''
    ecsv_files = sorted(glob.glob("*.ecsv"))
    ecsv_files = [i for i in ecsv_files if 'VER' in i]

    for file in ecsv_files:
            fnm = file.split('.')[0]
            split_lst = fnm.split("-")
            try:
                int(split_lst[-1])
                new_name = file.replace(split_lst[1],str(int(split_lst[1])) + "-" + str(int(split_lst[-1])))
                new_name = new_name.replace("-%d.ecsv" % int(split_lst[-1]), ".ecsv")
                shutil.move(file, new_name)
            except ValueError:
                try: # somehow reqd. for /2017/2017ICRC...35..712K
                    int(split_lst[2])
                    shutil.move(file, file.replace(split_lst[1], "%d" % int(split_lst[1])))
                    # pass # do nothing because in this case, it is already in the VER-XX-Y-Z format
                except ValueError:
                    shutil.move(file, file.replace(split_lst[1],str(int(split_lst[1]))+"-1"))
    return None

def chk_update_src(src_dict, src_id):
    try:
        src_cnt = src_dict[src_id]
        src_dict[src_id] = src_cnt+1
    except KeyError:
        src_dict[src_id] = 1
    return src_dict

def load_sources():
    '''
    Returns a (non-repeated) list of sources in the paper directory
    '''
    try:
        src_lst = load_yaml('info.yaml')['source_id']
    except KeyError:
        src_lst = []
        print('%s will not be translated into HEASARC as its info.yaml contains no source_ids' % os.path.basename(os.getcwd()))
    
    return sorted(src_lst)

def get_src_files(src_id):
    file_list = glob.glob("./**/*.*", recursive=True)
    file_list = [i for i in file_list if "VER-%d" % src_id in i]
    return file_list

def move_dirs(heasarc_dir, paper_dir, src_lst):
    if src_lst: # only if len(src_lst)>0
        for src_id in src_lst:
            if src_id < 999: dest_path = heasarc_dir+"detected_data/VER-%d/%s" % (src_id,os.path.basename(paper_dir))
            else: dest_path = heasarc_dir+"undetected_data/VER-%d/%s" % (src_id,os.path.basename(paper_dir))

            os.makedirs(dest_path, exist_ok=True) # create the src_id/paper_dir directory
            file_lst = glob.glob("./**/*.*", recursive=True) # so the destination dir is flattened
            file_lst = [i for i in file_lst if "VER-%d" % src_id in i] # for that source_id
            cp_files = [shutil.copy(file, dest_path) for file in file_lst]
    return None

def process_files(heasarc_dir):
    directories=sorted([d for d in list_dir(heasarc_dir) if os.path.isdir(d)])
    data_dir = [i for i in directories if isint(os.path.basename(i))] # get data diectories sorted by increasing years

    for year_dir in data_dir:
        papers_dir =  sorted([d for d in list_dir(year_dir) if os.path.isdir(d)])
        for paper_dir in papers_dir:
            # print(paper_dir) # for debugging
            os.chdir(paper_dir) # should be common for all following functions
            rem_non_ver_files()
            rem_non_skymap_fits()
            rem_nosrc_ecsv()
            print(paper_dir)
            rem_non_ver_data()
            rename_figures(paper_dir)
            rename_skymap()
            rename_yaml()
            change_ecsv_format()

            src_lst = load_sources()
            move_dirs(heasarc_dir, paper_dir, src_lst)
    
    a = [shutil.rmtree(d) for d in data_dir] # remove data directories (except /sources) after restructring is complete
    return None

def restructure_to_heasarc(repo_dir, heasarc_dir):

    make_copy(repo_dir, heasarc_dir)
    rename_sources(heasarc_dir)
    rem_paper(heasarc_dir)
    process_files(heasarc_dir)
    print("Restructuring completed successfully!")
    return None