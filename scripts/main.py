from restructure import run_restructure
from ver_name import gen_ver_names
from catalog import build_catalog
from create_fits import make_fits

run_restructure()
ver_names = gen_ver_names()
build_catalog()
make_fits()