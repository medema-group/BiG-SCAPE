#!/usr/bin/env python

import os
from distutils.core import setup

def generate_package_data():
    data_files = []
    for directory in ['html_template','Annotated_MIBiG_reference']:
        for path, dirs, files in os.walk(directory):
            for f in files:
                data_files.append(os.path.join(path, f))
    data_files.append('anchor_domains.txt')
    data_files.append('domain_includelist.txt')
    data_files.append('domains_color_file.tsv')
    package_data={'BiG-SCAPE': data_files}
    return package_data

setup(
     name='BiG-SCAPE',
     version='1.1.4',
     py_modules=['functions','ArrowerSVG'],
     packages=['BiG-SCAPE'],
     package_dir={'BiG-SCAPE': '.'},
     package_data=generate_package_data(),
     scripts=['bigscape.py'],
     description='Biosynthetic Genes Similarity Clustering and Prospecting Engine.',
     keywords=['secondary metabolites', 'natural products', 'genome mining'],
     url='https://github.com/medema-group/BiG-SCAPE',
     author='Jorge Navarro',
     author_email='jorge.c.navarro.munoz@gmail.com',
     license='GNU Affero v3'
     )
