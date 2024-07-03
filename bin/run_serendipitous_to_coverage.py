
import os, sys

graceids = ['S190408an', # BBH
            'S190412m', # BBH
            'S190421ar', # BBH
            'S190425z', # BNS
            'S190426c', # NSBH
            'S190503bg', # BBH
            'S190510g', # BNS
            'S190512at', # BBH
            'S190521g', # IMBH
            'S190627e', # BBH
            'S190701ah', # BBH
            'S190718y', # BNS 
            'S190720a', # BBH 
            'S190727h' # BBH 
            'S190728q', # BBH
            'S190814bv', # NSBH
            'S190816u', # NSBH 
            'S190828j', # BBH
            'S190901ap', # BNS 
            'S190910d', # BBH 
            'S190910h', # BBH 
            'S190915ak', # BBH 
            'S190923y', # NSBH
            'S190924h', # BBH
            'S190928c', # cosmic string
            'S190930s', # massgap
            'S190930t' # NSBH
            ]

graceids = ['S190816u', # NSBH 
            'S190828j', # BBH
            'S190901ap', # BNS 
            'S190910d', # BBH 
            'S190910h', # BBH 
            'S190915ak', # BBH 
            'S190923y', # NSBH
            'S190924h', # BBH
            'S190928c', # cosmic string
            'S190930s', # massgap
            'S190930t' # NSBH
            ]

for graceid in graceids:
    system_command = "python serendipitous_to_coverage.py --doGraceDB -f ../data/serendipitous/ztf_obs_history_2018_03_17.txt -t allobsfile -g %s -o ../data/%s/ZTF_fields_MSIP.dat --doMSIPOnly" % (graceid, graceid)
    os.system(system_command)

    system_command = "python serendipitous_to_coverage.py --doGraceDB -f ../data/serendipitous/ztf_obs_history_2018_03_17.txt -t allobsfile -g %s -o ../data/%s/ZTF_fields.dat" % (graceid, graceid)
    os.system(system_command)
