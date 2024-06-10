#!/usr/bin/env python3

"""
BiG-SCAPE

PI: Marnix Medema                   marnix.medema@wur.nl

Maintainers:
Jorge Navarro                       jorge.navarromunoz@wur.nl
Arjan Draisma                       arjan.draisma@wur.nl
Catarina Loureiro                   catarina.salesesantosloureiro@wur.nl
Nico Louwen                         nico.louwen@wur.nl

Developers:
Jorge Navarro                       jorge.navarromunoz@wur.nl
Satria Kautsar                      sakautsar@lbl.gov
Emmanuel (Emzo) de los Santos       E.De-Los-Santos@warwick.ac.uk
Arjan Draisma                       arjan.draisma@wur.nl
Catarina Loureiro                   catarina.salesesantosloureiro@wur.nl
Nico Louwen                         nico.louwen@wur.nl


Usage:   Please see `python bigscape.py -h`
                    `python bigscape.py cluster -h`

Example: python bigscape.py cluster -i ./inputfiles -o ./results -c 8 -p ./Pfam-A.hmm


Official repository:
https://github.com/medema-group/BiG-SCAPE


License: GNU Affero General Public License v3 or later
A copy of GNU AGPL v3 should have been included in this software package in
LICENSE.txt.
"""

from big_scape import main

if __name__ == "__main__":
    main.main()
