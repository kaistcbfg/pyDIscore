
import numpy as np
import argparse
import pickle
import gzip
import sys

def str2bool(v):
	if isinstance(v, bool): return v
	if v.lower() in ['yes', 'true', 't', 'y', '1']: return True
	elif v.lower() in ['no', 'false', 'f', 'n', '0']: return False
	else: raise argparse.ArgumentTypeError('Boolean value expected.')
#

def gen_chrsize_num_dict(genome_info, resolution):

	chrsize_dict = {}
	chrnum_dict = {}

	idx = 1
	f = open(genome_info)
	for line in f:
		size = int(int(line.split()[1])/resolution) + 1
		targetchrname = line.split()[0]
		chrsize_dict[targetchrname] = size
		chrnum_dict[targetchrname] = str(idx)
		idx+=1
	#
	
	return chrnum_dict,chrsize_dict
#

def make_map_from_covnormfile(covnormfile, bin_length, resolution):
	
	contact_map = np.zeros((bin_length, bin_length))
	
	f = gzip.open(covnormfile)
	f.readline()
	for line in f:
		line = line.rstrip()
		linedata = line.split('\t')
		bin1 = int(int(linedata[0].split('.')[1])/resolution)
		bin2 = int(int(linedata[1].split('.')[1])/resolution)
		freq = float(linedata[8])

		contact_map[bin1][bin2] += freq
		contact_map[bin2][bin1] += freq
	#
	f.close()

	return contact_map
#

def find_starter(contact_map):

	rowsum = np.sum(contact_map, axis=1)
	starter = np.argwhere(rowsum).squeeze()[0]

	return starter
#

def calc_DI(contact_map, starter, window_bin):

	rowsize, colsize = contact_map.shape
	DIdict = {}

	for rowidx in range(starter, rowsize):

		row = contact_map[rowidx]
		A = 0
		B = 0

		for z in range(rowidx-window_bin, rowidx):
			if not ((z < starter) or (z >= rowsize)): A += row[z]
		#
		for z in range(rowidx+1, rowidx+window_bin):
			if not ((z < starter) or (z >= rowsize)): B += row[z]
		#

		E = (A+B)/2
		if E == 0: DI = 0
		elif A == B: DI = 0
		else: DI = ((B-A)/np.abs(B-A)) * ( (((A-E)**2)/E) + (((B-E)**2)/E) )

		DIdict[rowidx] = DI
	#

	return DIdict
# 

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='python DI score calc.')
	parser.add_argument('--input-file', type=str, help='input *.gz file (required)', required=True)
	parser.add_argument('--input-format', type=str, default='cov', help='default cov, format: cov (covNorm) or pkl (numpy pickled array)')
	parser.add_argument('--chrname', type=str, help='target chromosome (required)', required=True)
	parser.add_argument('--fai-file', type=str, help='FAI file for chr size (required)', required=True)
	parser.add_argument('--resolution', type=int, default=40000, help='bin resolution (default 40kb)')
	parser.add_argument('--window-size', type=int, default=2000000, help='window size (default 2Mb)')
	parser.add_argument('--double-count-flag', type=str2bool, default=False, help='default False, if True, apply /2 to matrix')
	parser.add_argument('--chrname-number-flag',type=str2bool, default=False, help='default False, if True chr1 -> 1 (X:23, Y:24, M:25')
	parser.add_argument('--fullbin-output-flag', type=str2bool, default=True, help='default True, if False, print from starter bin')
	parser.add_argument('--output-file', type=str, help='output DI score bedgraph file. if None, print')
	args = parser.parse_args()

	chrnum_dict,chrsize_dict = gen_chrsize_num_dict(args.fai_file, args.resolution)
	bin_length = chrsize_dict[args.chrname]
	window_bin = int(args.window_size/args.resolution)

	if args.input_format == 'cov': contact_map = make_map_from_covnormfile(args.input_file, bin_length, args.resolution)
	elif args.input_format == 'pkl': contact_map = pickle.load(gzip.open(args.input_file,'rb'))
	else: sys.exit('Error: wrong --input-format value') 

	if args.double_count_flag == True: contact_map = contact_map/2
	starter_bin = find_starter(contact_map)

	DIdict = calc_DI(contact_map, starter_bin, window_bin)
	DIbins = list(DIdict.keys())
	DIbins.sort()

	if args.chrname_number_flag: printchrname = chrnum_dict[args.chrname]
	else: printchrname = args.chrname

	if args.output_file!=None: o = open(args.output_file,'w')
	if args.fullbin_output_flag == True:
		for i in range(bin_length):
			if i in DIbins: DIval = DIdict[i]
			else: DIval = 0
			outstr = "{}\t{}\t{}\t{}\n".format(printchrname, i*args.resolution, (i+1)*args.resolution, DIval)
			if args.output_file!=None: o.write(outstr)
			else: print(outstr[:-1])
		# 
	else:
		for i in DIbins:
			outstr = "{}\t{}\t{}\t{}\n".format(printchrname, i*args.resolution, (i+1)*args.resolution, DIdict[i])
			if args.output_file!=None: o.write(outstr)
			else: print(outstr[:-1])
		#
	#
	if args.output_file!=None: o.close()

##
