#!/usr/bin/python3

import cgi,cgitb
import shutil
import requests
import re
import json
import sys
import random
import datetime
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib as mpl

cgitb.enable()

def get_peptides(_l):
	url = 'http://gpmdb.thegpm.org/protein/model/%s&excel=1' % (_l)
	session = requests.session()
	try:
		r = session.get(url,timeout=20)
	except requests.exceptions.RequestException as e:
		print(e)
		return None

	text = re.sub('\r\n','\n',r.text)
	return text.splitlines()

def get_protein(_l):
	url = 'http://gpmdb.thegpm.org/1/protein/sequence/acc=%s' % (_l)
	session = requests.session()
	try:
		r = session.get(url,timeout=20)
	except requests.exceptions.RequestException as e:
		print(e)
		return None
	try:
		values = json.loads(r.text)
	except:
		return None
	return values[0]

def get_description(_l):
	a = _l
	if re.search(r'A.+\.\d+$',_l):
		a = re.sub(r'\.\d+$','',_l)
	url = 'http://gpmdb.thegpm.org/1/protein/description/acc=%s' % (a)
	session = requests.session()
	try:
		r = session.get(url,timeout=20)
	except requests.exceptions.RequestException as e:
		print(e)
		return None
	try:
		values = json.loads(r.text)
	except:
		return None
	return values[0]

def make_peptide_png(_l,_plength,_lines,_title,_file):
	start = {}
	end = {}
	res = {}
	max = 0
	for line in _lines[1:]:
		vs = line.split('\t')
		vs[0] = int(vs[0])
		vs[1] = int(vs[1])
		vs[2] = int(vs[2])
		if vs[0] in start:
			start[vs[0]] = start[vs[0]] + vs[2]
		else:
			start[vs[0]] = vs[2]
		if vs[1] in end:
			end[vs[1]] = end[vs[1]] + vs[2]
		else:
			end[vs[1]] = vs[2]
		if vs[1] > max:
			max = vs[1]
		for a in range(vs[0],vs[1]+1):
			if a in res:
				res[a] = res[a] + vs[2]
			else:
				res[a] = vs[2]

	max = _plength
	max_res = 0
	for a in res:
		if res[a] > max_res:
			max_res = res[a]
	a = 1
	xs = []
	ys = []
	ok = False
	for a in range(1,max+1):
		xs.append(a)
		if a in res:
			ys.append(res[a])
			ok = True
		else:
			ys.append(0)
	if not ok:
		print('<br />')
		print('<p>No peptides observed for "%s"</p>' % (_l))
		return 0
	mpl.style.use('seaborn-notebook')
	plt.xlim(0,int(1.02*_plength))
	ms = 10
	plt.plot(xs,ys,color=(0.25,0,.75,.8),marker='',linestyle='solid',linewidth=1)
#	plt.yscale('log')
	plt.ylabel('PSM (tabbs)')
	plt.xlabel('residue')
	plt.grid(True, lw = 1, ls = '--', c = '.9')
#	plt.axvline(x=_plength,color=(.2,.2,.2,.5),linestyle='dotted',linewidth=1)
	desc = get_description(_l)
	gene = ''
	if desc.find(':p') != -1:
		gene = re.sub(r'\:p.+','',desc)
		gene += ':p'
	if gene:
		plt.title(_title+' %s'% (gene))
	else:
		plt.title(_title)
	fig = plt.gcf()
	fig.set_size_inches(10, 5)
	cl = re.sub('[\|\:]','_',_file)
	plt.gca().get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
	plt.gca().get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
	fig.savefig('/var/www/intrinsicdisorder/ptm_png/%s_peps.png' % (cl), dpi=200, bbox_inches='tight')
	try:
		shutil.copy2('/var/www/intrinsicdisorder/ptm_png/%s_peps.png' % (cl),'/mnt/Actinium/ptm_png_a')
	except:
		pass
	nl = ''
	up = re.findall(r'UP\:(\w+)',desc)
	nl = '<a class="bluesq" href="/a/ptm_png.py?l=%s" target="_nl" title="Check for PTMs">&#x1F384;</a>&nbsp;' % (_l)
	nl += '<a class="bluesq" href="/a/nl_png.py?l=%s" target="_nl" title="Check for N-linked glycosylation">&#x1F36C;</a>&nbsp;' % (_l)
	nl += '<a class="bluesq" href="/a/seq.py?l=%s" target="_seq" title="Sequence display">&#x270E;</a>' % (_l)
	if len(up) and up[0] != 'NA':
		nl += '&nbsp;<a class="bluesq" href="https://alphafold.ebi.ac.uk/entry/%s" target="_af" title="AF structure prediction">&#129526;</a>' % (up[0])
	print('''<div class="pic"><img src="/ptm_png/%s_peps.png" height="400" width="800"/></div>''' % (cl))
	print('<p class="desc">%s (%s)</p>' % (desc,nl))
	print('''<p class="con">This Westmore-Standing (W-S) diagram shows the number of times a residue has been observed in a peptide-to-spectrum match in GPMDB, as a function of the residue's position in 
	the corresponding protein sequence.</p>''')
	return

print('Content-type: text/html\n\n')
form = cgi.FieldStorage()
try:
	label = form['l'].value
except:
	print('There must be a protein accession value specified')
	exit()

filename = label
title = '%s W-S diagram' % (label)
sys.stdout.flush()

protein = get_protein(label)

ls = get_peptides(label)
make_peptide_png(label,len(protein),ls,title,filename)

