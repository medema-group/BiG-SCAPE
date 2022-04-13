#!/usr/bin/env python

import os
from distutils.core import setup

def generate_data_files():
    data_files = []
    for directory in ['html_template','Annotated_MIBiG_reference']:
        for path, dirs, files in os.walk(directory):
            install_dir = os.path.join('bin',path)
            list_entry = (install_dir, [os.path.join(path, f) for f in files ])
            data_files.append(list_entry)
    return data_files

setup(
     name='BiG-SCAPE',
     version='1.1.3',
     py_modules=['functions','ArrowerSVG'],
     scripts=['bigscape.py'],
     data_files=generate_data_files(),
     description='Biosynthetic Genes Similarity Clustering and Prospecting Engine.',
     keywords=['secondary metabolites', 'natural products', 'genome mining'],
     url='https://github.com/medema-group/BiG-SCAPE',
     author='Jorge Navarro',
     author_email='jorge.c.navarro.munoz@gmail.com',
     license='GNU Affero v3'
     )


