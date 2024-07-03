from astropy.io import ascii
from astropy.table import Table

def get_candidate_dict(events, candidate_file):
    """file should be in the format, [event, name, ra, dec, ...]"""
    new_cands = Table(ascii.read(candidate_file))
    candidate_dict = {}
    for event in events:
        mask = new_cands['event'] == event
        candidate_dict[event] = [new_cands[mask]['name']]
#file = 'results_kowalski_queries_not_in_paper.csv'
#events = ['S190425z', 'S190426c', 'S190901ap', 'S190910d', 'S190910h', 'S190923y', 'S190930t']
#candidate_dict = get_candidate_dict(events, file)
