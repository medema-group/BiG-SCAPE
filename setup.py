#!/usr/bin/env python

from distutils.core import setup

setup(
     name='BiG-SCAPE',
     version='1.1.3',
     py_modules=['functions','ArrowerSVG'],
     scripts=['bigscape.py'],
     package_data=['': ['html_template/*',
                        'html_template/*/*',
                        'html_template/*/*/*/*',
                        '*domain*'
                        ]],
     include_package_data=True,
     description='Biosynthetic Genes Similarity Clustering and Prospecting Engine.',
     keywords=['secondary metabolites', 'natural products', 'genome mining'],
     url='https://github.com/medema-group/BiG-SCAPE',
     author='Jorge Navarro',
     author_email='jorge.c.navarro.munoz@gmail.com',
     license='GNU Affero v3'
     )


