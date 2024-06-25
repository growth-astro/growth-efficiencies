import fastavro
import numpy as np


def read_avro(alert_name):
	'''
	read avro data into dictionary
	'''
	with open(alert_name, 'rb') as fo:		
		reader = fastavro.reader(fo)
		candidate = next(reader, None)
	return candidate


def get_arrays(candidate, keys):
	'''
	>>> arr1, arr2 = get_arrays(candidate, ["key1", "key2"])
	get numpy array for keys

	written by Sjoert van Velzen, May 2018
	'''
	return [np.array(arr) for arr in [np.array([can[key] for can in candidate['prv_candidates']] + [candidate['candidate'][key]]) for key in keys]]
