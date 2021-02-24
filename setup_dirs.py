import os

f = open('prefixes.txt')
line = f.readline()
while line:
	prefix = line.rstrip('\n')
	src = prefix[prefix.index('GLEAM'):].split('.ms')[0]
	os.system('mkdir '+src)
	os.chdir(src)
	os.system('ln -s ../alpha/* .')
	g = open('run_'+src+'.sh','w')
	g.write('#!/bin/bash\n')
	g.write('python 01_get_images.py '+prefix+'\n')
	g.write('python 02_convolve_restore_pbcor.py '+prefix.split('/')[-1]+'\n')
	g.write('python 03_make_alpha_maps.py '+prefix.split('/')[-1]+'\n')
	g.close()
	os.chdir('../')
	print(prefix,src)
	line = f.readline()
f.close()
