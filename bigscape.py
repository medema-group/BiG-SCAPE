#!/usr/bin/env python3

"""
BiG-SCAPE

PI: Marnix Medema                   marnix.medema@wur.nl

Maintainers:
Jorge Navarro                       jorge.navarromunoz@wur.nl
Arjan Draisma                       arjan.draisma@wur.nl
Catarina Sales E Santos Loureiro    catarina.salesesantosloureiro@wur.nl
Nico Louwen                         nico.louwen@wur.nl

Developers:
Jorge Navarro                       jorge.navarromunoz@wur.nl
Satria Kautsar                      sakautsar@lbl.gov
Emmanuel (Emzo) de los Santos       E.De-Los-Santos@warwick.ac.uk
Arjan Draisma                       arjan.draisma@wur.nl
Catarina Sales E Santos Loureiro    catarina.salesesantosloureiro@wur.nl
Nico Louwen                         nico.louwen@wur.nl


Usage:   Please see `python bigscape.py -h`

Example: python bigscape.py -c 8 --pfam_dir ./ -i ./inputfiles -o ./results


Official repository:
https://github.com/medema-group/BiG-SCAPE


# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in
# LICENSE.txt.
"""

from big_scape.main import run_bigscape

if __name__ == "__main__":
    run_bigscape()
