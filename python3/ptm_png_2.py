#!/usr/bin/python3

import cgi,cgitb

import shutil
import sys
import requests
import re
import json
import random
import datetime
import copy
import codecs
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib as mpl

cgitb.enable()

motifs = [
	('R..[ST]P',3,'RxxΦP'),
	('R..[ST]F',3,'RxxΦF'),
	('G..[ST]P',3,'GxxΦP'),
	('P.[ST]P',2,'PxΦP'),
	('[ST]P.K',0,'ΦPxK'),
	('[ST]P',0,'ΦP'),
	('R..S..[ST]',6,'RxxSxxΦ'),
	('L.R..[ST]',5,'LxRxxΦ'),
	('S..[ST]L',3,'SxxΦL'),
	('R..[ST][LIV]',3,'RxxΦL/I/V'),
	('R..[ST]',3,'RxxΦ'),
	('[ST]..[DE]',0,'ΦxxD/E'),
	('RR.[ST]',3,'RRxΦ'),
	('N.[ST]',2,'NxΦ'),
	('[RK]R..[ST]',4,'R/KRxxΦ'),
	('G[ST]',1,'GΦ'),
	('[ST][DE].[DE]',0,'ΦD/ExD/E'),
	('RS.[ST]',3,'RSxΦ'),
	('[ST].E.[LIV]',0,'ΦxExL/I/V'),
	('[ST]F',0,'ΦF'),
	('N.[ST][LIV]',2,'NxΦL/I/V'),
	('[ST]Q',0,'ΦQ'),
	('[ST][LIV]',0,'ΦL/I/V'),
	('Y..P',0,'YxxP: Src')
]

def get_motif(_a,_pm):
	for i,p in enumerate(_pm):
		if _a in p:
			return motifs[i][2]
	return ''

def load_motifs(_seq):
	values = []
	for motif in motifs:
		r = [m.start()+motif[1]+1 for m in re.finditer(motif[0],_seq)]
		values.append(set(r))
	return values

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

def create_plot(plt,xs,ys,ptype):
	ms = 10
	ms8 = 8
	clrs = {		'Azure blue':		'#0007FFAA',
				'British racing green':	'#004225AA',
				'Honey':		'#EC9702AA',
				'India green':		'#138808AA',
				'Lime green':		'#32CD32AA',
				'Purple':		'#6C0BA9AA',
				'Red':			'#FF0000AA',
				'Turquoise':		'#00FFEFAA',
				'Violet':		'#710193AA',
				'White':		'#EEEEEEAA',
				'Khaki':		'#F0E68CAA',
				'Dark grey':		'#060606AA'
	}

	if not len(ptype) or ptype.find('N') > -1:
		plt.plot(xs['N-acetyl'],ys['N-acetyl'],color=clrs['Purple'],markersize=ms,marker='o',linestyle='None',label='n-acetyl')
	if not len(ptype) or ptype.find('A') > -1:
		plt.plot(xs['acetyl'],ys['acetyl'],color=clrs['Azure blue'],markersize=ms,marker='o',linestyle='None',label='K-acetyl')
	if not len(ptype) or ptype.find('S') > -1:
		plt.plot(xs['succinyl'],ys['succinyl'],color=clrs['India green'],markersize=ms,marker='o',linestyle='None',label='K-succyl')
	if not len(ptype) or ptype.find('P') > -1:
		plt.plot(xs['S-phosphoryl'],ys['S-phosphoryl'],markersize=ms,color=clrs['Red'],marker='v',linestyle='None',label='S-phos')
		plt.plot(xs['SP-phosphoryl'],ys['SP-phosphoryl'],markersize=ms,color=clrs['Honey'],marker='v',linestyle='None',label='SP-phos')
		plt.plot(xs['T-phosphoryl'],ys['T-phosphoryl'],markersize=ms,color=clrs['Red'],marker='^',linestyle='None',label='T-phos')
		plt.plot(xs['TP-phosphoryl'],ys['TP-phosphoryl'],markersize=ms,color=clrs['Honey'],marker='^',linestyle='None',label='TP-phos')
		plt.plot(xs['Y-phosphoryl'],ys['Y-phosphoryl'],markersize=ms,color=clrs['Red'],marker='o',linestyle='None',label='Y-phos')
	if not len(ptype) or ptype.find('s') > -1:
		plt.plot(xs['S-phosphoryl'],ys['S-phosphoryl'],markersize=ms,color=clrs['Red'],marker='v',linestyle='None',label='S-phos')
		plt.plot(xs['SP-phosphoryl'],ys['SP-phosphoryl'],markersize=ms,color=clrs['Honey'],marker='v',linestyle='None',label='SP-phos')
	if not len(ptype) or ptype.find('t') > -1:
		plt.plot(xs['T-phosphoryl'],ys['T-phosphoryl'],markersize=ms,color=clrs['Red'],marker='^',linestyle='None',label='T-phos')
		plt.plot(xs['TP-phosphoryl'],ys['TP-phosphoryl'],markersize=ms,color=clrs['Honey'],marker='^',linestyle='None',label='TP-phos')
	if not len(ptype) or ptype.find('y') > -1:
		plt.plot(xs['Y-phosphoryl'],ys['Y-phosphoryl'],markersize=ms,color=clrs['Red'],marker='o',linestyle='None',label='Y-phos')
	if not len(ptype) or ptype.find('U') > -1:
		plt.plot(xs['ubiquitinyl'],ys['ubiquitinyl'],markersize=ms,color=clrs['Lime green'],marker='v',linestyle='None',label='K-GG')
		plt.plot(xs['K-sumoyl'],ys['K-sumoyl'],markersize=ms,color=clrs['Lime green'],marker='X',linestyle='None',label='K-sumo')
	if not len(ptype) or ptype.find('M') > -1:
		plt.plot(xs['R-dimethyl'],ys['R-dimethyl'],markersize=ms,color=clrs['British racing green'],marker='d',linestyle='None',label='R-dimet')
	if not len(ptype) or ptype.find('O') > -1:
		plt.plot(xs['K-oxidation'],ys['K-oxidation'],markersize=ms,color=clrs['Turquoise'],marker='*',linestyle='None',label='K-oxy')
		plt.plot(xs['P-oxidation'],ys['P-oxidation'],markersize=ms,color=clrs['Violet'],marker='*',linestyle='None',label='P-oxy')
	if not len(ptype) or ptype.find('C') > -1:
		plt.plot(xs['citrulline'],ys['citrulline'],markersize=ms,color=clrs['British racing green'],marker='.',linestyle='None',label='R-citr')
	if not len(ptype) or ptype.find('G') > -1:
		plt.plot(xs['ST-glyco'],ys['ST-glyco'],markersize=ms,markerfacecolor=clrs['White'],markeredgewidth=1,markeredgecolor=clrs['Dark grey'],marker='s',linestyle='None',label='O-gly')
	if not len(ptype) or ptype.find('L') > -1:
		plt.plot(xs['gla'],ys['gla'],markersize=ms,markerfacecolor=clrs['Khaki'],markeredgewidth=1,markeredgecolor=clrs['Dark grey'],marker='o',linestyle='None',label='Gla')

def make_ptm_csv(_l,_plength,_title,_protein,_file,_y,_s,_ltype,_ptype):
#	print('<div id="diagram" class="pic"><center><img src="/pics/loading%i.gif"/></center></div>' % (random.randint(2,5)))
#	sys.stdout.flush()
	session = requests.session()
	seq = list(_protein)

	url = 'http://gpmdb.thegpm.org/1/peptide/pf/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	values = {'acetyl':None,'phosphoryl':None,'ubiquitinyl':None}
	try:
		r = session.get(url,timeout=20)
	except:
		print('Could not connect to GPMDB')
		return None
	try:
		values['phosphoryl'] = json.loads(r.text)
	except:
		print('JSON error: problem with phosphorylation data')
		return None

	url = 'http://gpmdb.thegpm.org/1/peptide/af/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except:
		print('Could not connect to GPMDB')
		return None
	try:
		values['acetyl'] = json.loads(r.text)
	except:
		print('JSON error: problem with acetylation data')
		return None
	try:
		values['acetyl'] = json.loads(r.text)
	except:
		print('JSON error: problem with acetylation data')
		return None

	url = 'http://gpmdb.thegpm.org/1/peptide/gl/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except:
		print('Could not connect to GPMDB')
		return None
	try:
		values['gla'] = json.loads(r.text)
	except:
		print('JSON error: problem with GLA data')
		return None

	url = 'http://gpmdb.thegpm.org/1/peptide/uf/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except:
		print('Could not connect to GPMDB')
		return None
	try:
		values['ubiquitinyl'] = json.loads(r.text)
	except:
		print('JSON error: problem with ubiquitinylation data')
		return None
	#formulate a URL to request information about sumoylation for the protein identified by _l
	url = 'http://gpmdb.thegpm.org/1/peptide/su/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except requests.exceptions.RequestException as e:
		print(e)
		return None
	try:
		values['K-sumoyl'] = json.loads(r.text)
	except:
		return None

	#formulate a URL to request information about succinylation for the protein identified by _l
	url = 'http://gpmdb.thegpm.org/1/peptide/sc/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except requests.exceptions.RequestException as e:
		print(e)
		return None
	try:
		values['succinyl'] = json.loads(r.text)
	except:
		return None

	#formulate a URL to request information about R-dimethylation for the protein identified by _l
	url = 'http://gpmdb.thegpm.org/1/peptide/di/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except:
		print('Could not connect to GPMDB') 
		return None
	try:
		values['dimethyl'] = json.loads(r.text)
	except:
		print('JSON error: problem with dimethylation data')
		return None
	#formulate a URL to request information about O-linked glycosylation for the protein identified by _l
	url = 'http://gpmdb.thegpm.org/1/peptide/ol/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except:
		print('Could not connect to GPMDB')
		return None
	try:
		values['ST-glyco'] = json.loads(r.text)
	except:
		print('JSON error: problem with dimethylation data')
		return None
	#formulate a URL to request information about R-citrullination for the protein identified by _l
	url = 'http://gpmdb.thegpm.org/1/peptide/ct/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except requests.exceptions.RequestException as e:
		print(e)
		return None
	try:
		values['citrulline'] = json.loads(r.text)
	except:
		return None
	#formulate a URL to request information about KP-oxidation for the protein identified by _l
	url = 'http://gpmdb.thegpm.org/1/peptide/ox/acc=%s&pos=1-%i&w=n' % (_l,_plength)
	try:
		r = session.get(url,timeout=20)
	except:
		print('Could not connect to GPMDB')
		return None
	try:
		values['oxidation'] = json.loads(r.text)
	except:
		print('JSON error: problem with oxidation data')
		return None
	a = 1;
	xs = {	'citrulline':[],
			'succinyl':[],
			'acetyl':[],
			'N-acetyl':[],
			'S-phosphoryl':[],
			'SP-phosphoryl':[],
			'T-phosphoryl':[],
			'TP-phosphoryl':[],
			'Y-phosphoryl':[],
			'ubiquitinyl':[],
			'K-sumoyl':[],
			'R-dimethyl':[],
			'P-oxidation':[],
			'K-oxidation':[],
			'ST-glyco':[],
			'gla':[]
		}
	ys = copy.deepcopy(xs)
	ts = {'acetyl':0,'phosphoryl':0,'oxidation':0,'ubiquitinyl':0,'K-sumoyl':0,'dimethyl':0,'citrulline':0,'ST-glyco':0,'succinyl':0,'gla':0}
	min_obs = 5
	lines = []
	href = 'http://gpmdb.thegpm.org/thegpm-cgi/dblist_label_model.pl?label=%s&residues=' % (_l)
	ms = re.finditer(r'(?=N[^P][ST])',_protein)
	nl = ''
	if True or len([m for m in ms]) > 0:
		nl = '<a class="bluesq" href="/a/nl_png.py?l=%s" target="_nl" title="Check for N-linked glycosylation">&#127852;</a>&nbsp;' % (_l)
		nl += '<a class="bluesq" href="/a/peptides_png.py?l=%s" target="_nl" title="Check for observable peptides">&#x1F527;</a>&nbsp;' % (_l)
		nl += '<a class="bluesq" href="/a/seq.py?l=%s" target="_seq" title="Sequence display">&#x270E;</a>&nbsp;' % (_l)

	phospho_motifs = load_motifs(_protein)

	phospho = {'S':0,'T':0,'Y':0}
	hydroxy = {'K':0,'P':0}
	oglyco = {'S':0,'T':0}
	kms = [_protein.find('M',1)]
	if kms[0] != -1:
		kms[0] = kms[0] + 1
		kms.append(kms[0]+1)
	Ktypes = {}
	Ktype = ''
	Ktabbs = {}
	Ptypes = {}
	Ptabbs = {}
	for a in range(1,_plength+1):
		b = str(a)
		if b in values['acetyl']:
			if values['acetyl'][b] is None:
#				print('<script>document.getElementById("diagram").style="display: none;"</script>')
				print('<br />')
				print('<div><p>No PTMs detected for "%s".</p></div>' % (_l))
				return
		pre = ''
		if a - 4 >= 0:
			pre = seq[a-4].lower()
		elif a - 4 == -1:
			pre = '['
		if a - 3 >= 0:
			pre += seq[a-3].lower()
		elif a - 3 == -1:
			pre = '['
		if a - 2 >= 0:
			pre += seq[a-2].lower()
		elif a - 2 == -1:
			pre = '['
		post = ''
		if a < len(seq):
			post = seq[a].lower()
			if a+1 < len(seq):
				post += seq[a+1].lower()
			else:
				post += ']'
			if a+2 < len(seq):
				post += seq[a+2].lower()
			elif post.find(']') == -1:
				post += ']'
		elif a == len(seq):
			post = ']'
		if len(lines) % 2 != 0:
			line = '<tr><td><a href="%s%i" target="_s">%i</a></td><td>%s</td>' % (href,a,a,'<i style="font-size: 10pt">%s&middot;</i>%s<i style="font-size: 10pt">&middot;%s</i>' % (pre,seq[a-1],post))
		else:
			line = '<tr class="alt"><td><a href="%s%i" target="_s">%i</a></td><td>%s</td>' % (href,a,a,'<i style="font-size: 10pt">%s&middot;</i>%s<i style="font-size: 10pt">&middot;%s</i>' % (pre,seq[a-1],post))
		ok = False
		mname = 'acetyl'
		Ktype = ''
		tabbs = 0
		if(b in values[mname]):
			tabbs = 0
			if values[mname][b] >= min_obs and (seq[a-1] == 'K' or a < 4):
				if seq[a-1] == 'K':
					xs[mname].append(a)
					ys[mname].append(values[mname][b])
					ts[mname] += 1
					Ktype += '+acetyl'
					tabbs += values[mname][b]
				else:
					xs['N-acetyl'].append(a)
					ys['N-acetyl'].append(values[mname][b])
					ts[mname] += 1
				line += '<td>%i</td>' % (values['acetyl'][b])
				ok = True
			elif(b in values[mname]) and a > 3 and a in kms:
				if values[mname][b] >= min_obs:
					xs['N-acetyl'].append(a)
					ys['N-acetyl'].append(values['acetyl'][b])
					ts[mname] += 1
					line += '<td>%i</td>' % (values[mname][b])
					tabbs += values[mname][b]
					ok = True
				else:
					line += '<td></td>'
			else:
				line += '<td></td>'
		else:
			line += '<td></td>'

		mname = 'succinyl'
		if(b in values[mname]) and seq[a-1] == 'K':
			if values[mname][b] >= 2:
				xs[mname].append(a)
				ys[mname].append(values[mname][b])
				ts[mname] += 1
				line += '<td>%i</td>' % (values[mname][b])
				ok = True
			else:
				line += '<td></td>'
		else:
			line += '<td></td>'

		mname = 'ubiquitinyl'
		if(b in values[mname]) and seq[a-1] == 'K':
			if values[mname][b] >= min_obs:
				xs[mname].append(a)
				ys[mname].append(values[mname][b])
				ts[mname] += 1
				line += '<td>%i</td>' % (values[mname][b])
				tabbs += values[mname][b]
				Ktype += '+GGyl'
				ok = True
			else:
				line += '<td></td>'
		else:
			line += '<td></td>'

		mname = 'K-sumoyl'
		if(b in values[mname]):
			if values[mname][b] >= 2 and seq[a-1] == 'K':
				xs[mname].append(a)
				ys[mname].append(values['K-sumoyl'][b])
				ts[mname] += 1
				Ktype += '+SUMOyl'
				tabbs += values[mname][b]
				line += '<td>%i</td>' % (values[mname][b])
				ok = True
			else:
				line += '<td></td>'
		else:
			line += '<td></td>'
		if Ktype:
			if Ktype in Ktypes:
				Ktypes[Ktype] += 1
			else:
				Ktypes[Ktype] = 1
			if Ktype in Ktabbs:
				Ktabbs[Ktype] += tabbs
			else:
				Ktabbs[Ktype] = tabbs
		mname = 'phosphoryl'
		if(b in values[mname]):
			if values[mname][b] >= min_obs:
				if seq[a-1] == 'S':
					if a < _plength and seq[a] == 'P':
						xs['SP-phosphoryl'].append(a)
						ys['SP-phosphoryl'].append(values[mname][b])
					else:
						xs['S-phosphoryl'].append(a)
						ys['S-phosphoryl'].append(values[mname][b])
					ts[mname] += 1
					phospho['S'] += 1
					mv = get_motif(a,phospho_motifs)
					if mv:
						line += '<td>%i</td><td>%s</td>' % (values[mname][b],mv)
						if mv in Ptypes:
							Ptypes[mv] += 1
						else:
							Ptypes[mv] = 1
						if mv in Ptabbs:
							Ptabbs[mv] += values[mname][b]
						else:
							Ptabbs[mv] = values[mname][b]
					else:
						line += '<td>%i</td><td></td>' % (values[mname][b])
					ok = True
				elif seq[a-1] == 'T':
					if a < _plength and seq[a] == 'P':
						xs['TP-phosphoryl'].append(a)
						ys['TP-phosphoryl'].append(values[mname][b])
					else:
						xs['T-phosphoryl'].append(a)
						ys['T-phosphoryl'].append(values[mname][b])
					ts[mname] += 1
					phospho['T'] += 1
					mv = get_motif(a,phospho_motifs)
					if mv:
						if mv in Ptypes:
							Ptypes[mv] += 1
						else:
							Ptypes[mv] = 1
						if mv in Ptabbs:
							Ptabbs[mv] += values[mname][b]
						else:
							Ptabbs[mv] = values[mname][b]
						line += '<td>%i</td><td>%s</td>' % (values[mname][b],mv)
					else:
						line += '<td>%i</td><td></td>' % (values[mname][b])
					ok = True
				elif seq[a-1] == 'Y':
					xs['Y-phosphoryl'].append(a)
					ys['Y-phosphoryl'].append(values[mname][b])
					ts[mname] += 1
					phospho['Y'] += 1
					mv = get_motif(a,phospho_motifs)
					if mv:
						if mv in Ptypes:
							Ptypes[mv] += 1
						else:
							Ptypes[mv] = 1
						if mv in Ptabbs:
							Ptabbs[mv] += values[mname][b]
						else:
							Ptabbs[mv] = values[mname][b]
						line += '<td>%i</td><td>%s</td>' % (values[mname][b],mv)
					else:
						line += '<td>%i</td><td></td>' % (values[mname][b])
					ok = True
				else:
					line += '<td></td><td></td>'
			else:
				line += '<td></td><td></td>'
		else:
			line += '<td></td><td></td>'
		mname = 'dimethyl'
		if(b in values[mname]):
			if values[mname][b] >= min_obs and seq[a-1] == 'R':
					xs['R-dimethyl'].append(a)
					ys['R-dimethyl'].append(values[mname][b])
					ts[mname] += 1
					line += '<td>%i</td>' % (values[mname][b])
					ok = True
			else:
				line += '<td></td>'
		else:
			line += '<td></td>'

		mname = 'oxidation'
		if(b in values[mname]):
			if values[mname][b] >= min_obs:
				if seq[a-1] == 'K':
					xs['K-oxidation'].append(a)
					ys['K-oxidation'].append(values[mname][b])
					line += '<td>%i</td>' % (values[mname][b])
					ts[mname] += 1
					hydroxy['K'] += 1
					ok = True
				elif seq[a-1] == 'P':
					xs['P-oxidation'].append(a)
					ys['P-oxidation'].append(values[mname][b])
					ts['oxidation'] += 1
					hydroxy['P'] += 1
					line += '<td>%i</td>' % (values[mname][b])
					ok = True
				else:
					line += '<td></td>'
			else:
				line += '<td></td>'
		else:
			line += '<td></td>'

		mname = 'citrulline'
		if(b in values[mname]):
			if values[mname][b] >= min_obs and seq[a-1] == 'R':
					xs[mname].append(a)
					ys[mname].append(values[mname][b])
					ts[mname] += 1
					line += '<td>%i</td>' % (values[mname][b])
					ok = True
			else:
				line += '<td></td>'

		else:
			line += '<td></td>'
		mname = 'ST-glyco'
		if(b in values[mname]):
			if values[mname][b] >= min_obs and seq[a-1] == 'S':
				xs[mname].append(a)
				ys[mname].append(values[mname][b])
				ts[mname] += 1
				oglyco['S'] += 1
				line += '<td>%i</td>' % (values[mname][b])
				ok = True
			elif values[mname][b] >= min_obs and seq[a-1] == 'T':
				xs[mname].append(a)
				ys[mname].append(values[mname][b])
				ts[mname] += 1
				oglyco['T'] += 1
				line += '<td>%i</td>' % (values[mname][b])
				ok = True
			else:
				line += '<td></td>'
		mname = 'gla'
		if(b in values[mname]) and seq[a-1] == 'E':
			if values[mname][b] >= 2:
				xs[mname].append(a)
				ys[mname].append(values[mname][b])
				ts[mname] += 1
				line += '<td>%i</td>' % (values[mname][b])
				ok = True
			else:
				line += '<td></td>'
		else:
			line += '<td></td>'

		if ok:
			lines.append(line + '</tr>')
	if(_s == 'xkcd'):
		plt.xkcd()
	else:
		if len(_s) > 0:
			try:
				plt.style.use(_s)
				mpl.style.use(_s)
			except:
				mpl.style.use('seaborn-notebook')
		else:
			mpl.style.use('seaborn-notebook')
	plt.xlim(0,int(1.02*_plength))
	create_plot(plt,xs,ys,_ptype)
	if _ltype == 'linear':
		plt.yscale('linear')
	else:
		plt.yscale('log')
	plt.ylabel('PSM (tabbs)')
	plt.xlabel('residue')
	plt.legend(loc='best')
	plt.grid(True, lw = 1, ls = '--', c = '.8')
#	plt.axvline(x=_plength,color=(.2,.2,.2,.5),linestyle='dotted',linewidth=1)
	desc = re.sub(r' \[',r'<br />[&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;',get_description(_l))
	desc = re.sub(r'[\[\]]',r'',desc)
	gene = ''
	if desc.find(':p') != -1:
		gene = re.sub(r'\:p.+','',desc)
		gene += ':p'
	if gene:
		plt.title(_title+' %s'% (gene))
	else:
		plt.title(_title)
	ax = plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
	ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
	fig = plt.gcf()
	fig.set_size_inches(10, 5)
	cl = re.sub('[\|\:]','_',_file)
	plt.gca().get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
	if _y is not None:
		print(_y)
		plt.ylim(1,_y)
	else:
		plt.ylim(1,None)
	ps = ax.get_ylim()
	if _ltype != 'linear':
		if ps[1] < 2000:
			ax.set_ylim([1,2000])
	if len(lines) == 0:
		plt.ylim(0,1000)
#		print('<br />')
#		print('<div><p>No PTMs detected for "%s". (<a href="http://gpmdb.thegpm.org/~/dblist_label/label=%s" title="main GPMDB page for this sequence" target="_blank">more</a>)</p></div>' % (_l,_l))
#		print('<p class="desc">%s %s</p>' % (re.sub(r'alt\:',r'<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;alt:',desc),nl))
#		return
	fname = '%s_ptms.png' % (cl)
	if _ptype:
		fname = '%s_%s_ptms.png' % (cl,_ptype)

	fig.savefig('/var/www/intrinsicdisorder/ptm_png/%s' % (fname), dpi=200, bbox_inches='tight')
	try:
		shutil.copy2('/var/www/intrinsicdisorder/ptm_png/%s' % (fname),'/mnt/Actinium/ptm_png_a')
	except:
		pass
	script = "<div id='diagram' class='pic'><center><img src='/ptm_png/%s' height='400' width='800' /></center></div>" % (fname)
	print(script)
	up = re.findall(r'UP\:(\w+)',desc)
	link = ''
	if len(up) and up[0] != 'NA':
		nl += '<a href="https://uniprot.org/uniprot/%s#ptm_processing" target="_blank" title="Un!Pr*t entry">%s</a> ' % (up[0],'&#x1F5D1;')
		nl += '<a href="https://glygen.org/protein/%s-1#Glycosylation" target="_blank" title="GlyGen entry">%s</a>' % (up[0],'&#x2600;')
		desc = re.sub(r'UP\:\w+',link,desc)
		nl += '<a class="bluesq" href="https://alphafold.ebi.ac.uk/entry/%s" target="_af" title="AF structure prediction">&#129526;</a>' % (up[0])
	else:
		up = None
	if _l.find('ENS') == 0:
		nl += '&nbsp;<a class="bluesq" href="/a/exon.py?l=%s" target="_ex" title="exon structure">&#x2612;</a>' % (_l)
	if desc.find('Source:SGD;Acc:') != -1:
		desc = re.sub(r'Acc:(S00\d+)',r'<a href="https://www.yeastgenome.org/locus/\1" target="_blank">\1</a>',desc)
	if desc.find('; MGI:') != -1:
		desc = re.sub(r'MGI:(\d+)',r'<a href="http://www.informatics.jax.org/marker/MGI:\1" target="_blank">MGI:\1</a>',desc)
	print('<p class="desc">%s (%s), %i aa</p>' % (re.sub(r'alt\:',r'<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;alt:',desc),nl,len(seq)))
	ktext = 'title="'
	ktotal = float(sum(Ktabbs.values()))
	for k in sorted(Ktypes.items(), key=lambda x:x[1]):
		if ktotal > 0.0:
			ktext += '%i× K%s (%i Ta, %.1f%%)\n' % (k[1],re.sub(r'(\w)\+(\w)',r'\1/\2',k[0]),Ktabbs[k[0]],100.0*Ktabbs[k[0]]/ktotal)
		else:
			ktext += '%i× K%s (%i Ta)\n' % (k[1],re.sub(r'(\w)\+(\w)',r'\1/\2',k[0]),Ktabbs[k[0]])
	ktext += '"'
	ptext = ''
	ptotal = float(sum(Ptabbs.values()))
	for p in sorted(Ptypes.items(), key=lambda x:x[1]):
		if ptotal > 0.0:
			ptext += '%i× %s (%i Ta, %.1f%%)\n' % (p[1],p[0],Ptabbs[p[0]],100.0*Ptabbs[p[0]]/ptotal)
		else:
			ptext += '%i× %s (%i Ta)\n' % (p[1],p[0],Ptabbs[p[0]])

	create_table(lines,ts,_l,phospho,hydroxy,oglyco,ktext,ptext)
	return

def create_table(_lines,_ts,_l,_sty,_pk,_og,_kt,_pt):
	print('''
	<p class="con">This Modification-Abundance (M-A) diagram shows the number of times a residue has been observed with a particular PTM in a peptide-to-spectrum match in GPMDB, 
	as a function of the residue's position in the corresponding protein sequence. 
	Note that the Y-axis has a log scale. 
	The particular residues and observation numbers are detailed in the following table:</p>''')
	if len(_lines) == 0:
		print('<div><p style="margin: auto;border-collapse: collapse;">No PTMs detected</p></div>')
		return
	print('<div><table style="margin: auto;border-collapse: collapse;">')
	head = '''<tr class="heads"><td><span title="Location of modification in protein coordinates">pos</span></td>
	<td><span title="Amino acid residue with the modification">res</span></td>
	<td>&nbsp;&nbsp;&nbsp;<span title="PSMs with N-terminal or K-acetylation at a site (tabbs)">acetyl</span>&nbsp;&nbsp;&nbsp;</td>
	<td><span title="PSMs with K-succinylation at a site (tabbs)">succyl</span></td>
	<td><span title="PSMs with  K-GG conjugation at a site (tabbs)">GG</span></td>
	<td><span title="PSMs with K-sumoylation at a site (tabbs)">sumo</span></td>
	<td><span title="PSMs with S/T/Y-phosphorylation at a site (tabbs)">phos</span></td>
	<td><span title="phosphorylation sequence motifs">p-motif</span></td>
	<td><span title="PSMs with R-dimethylation at a site (tabbs)">R-dimet</span></td>
	<td><span title="PSMs with P/K-oxidation at a site (tabbs)">P/K-oxy</span></td>
	<td><span title="PSMs with citrulline at a site (tabbs)">R-citr</span></td>
	<td><span title="PSMs with O-glycosylation at a site (tabbs)">O-gly</span></td>
	<td><span title="PSMs with gamma-carboxylation at a site (tabbs)">Gla</span></td>
	</tr>'''
	print(head)
	for l in _lines:
		print(l)
	if len(_lines) > 10:
		print(head)
	pline = 'title="%i× S+phosphoryl\n%i× T+phosphoryl\n%i× Y+phosphoryl\n%s"' % (_sty['S'],_sty['T'],_sty['Y'],_pt)
	kline = 'title="%i× K+oxide\n%i× P+oxide"' % (_pk['K'],_pk['P'])
	gline = 'title="%i× S+glycosyl\n%i× T+glycosyl"' % (_og['S'],_og['T'])
	print('''<tr class="tots"><td style="text-align: right">sites:</td>
	<td>%i</td>
	<td>%i</td>
	<td>%i</td>
	<td>%i</td>
	<td %s>%i</td>
	<td %s>%i</td>
	<td> </td>
	<td>%i</td>
	<td %s>%i</td>
	<td>%i</td>
	<td %s>%i</td>
	<td>%i</td>
	</tr>''' % (len(_lines),_ts['acetyl'],_ts['succinyl'],_ts['ubiquitinyl'],_kt,_ts['K-sumoyl'],pline,_ts['phosphoryl'],_ts['dimethyl'],kline,_ts['oxidation'],_ts['citrulline'],gline,_ts['ST-glyco'],_ts['gla']))
	print('</table></div>')
	print('''<div id="content"><ol>
	<li><b>pos:</b> location of the modification in protein coordinates;</li>
	<li><b>res:</b> amino acid residue with the modification and flanking residues;</li>
	<li><b>acetyl:</b> PSMs with N-terminal or K-acetylation at <i>pos</i> (tabbs);</li>
	<li><b>sucyl:</b> PSMs with K-succinylation at <i>pos</i> (tabbs);</li>
	<li><b>phos:</b> PSMs with S/T/Y-phosphorylation at <i>pos</i> (tabbs);</li>
	<li><b>p-motif:</b> assigned <a href=" https://doi.org/10.1371/journal.pbio.3000341" target="_blank" title="Evolution of protein kinase substrate recognition at the active site">phosphorylation motif pattern</a> where Φ is a modified S/T;</li>
	<li><b>GG:</b> PSMs with K-GG conjugation at <i>pos</i> (tabbs);</li>
	<li><b>sumo:</b> PSMs with K-sumoylation at <i>pos</i> (tabbs);</li>
	<li><b>R-dimeth:</b> PSMs with R-dimethylation at <i>pos</i> (tabbs);</li>
	<li><b>P/K-oxy:</b> PSMs with 4-hydroxyproline or 5-hydroxylysine at <i>pos</i> (tabbs);</li>
	<li><b>R-citr:</b> PSMs with citrulline at <i>pos</i> (tabbs);</li>
	<li><b>O-gly</b> PSMs with O-linked glycosylation at <i>pos</i> (tabbs); &amp;</li>
	<li><b>Gla</b> PSMs with E+gamma-carboxyl at <i>pos</i> (tabbs).</li>
	</ol></div>''')
	print('''<div id="content">p-code values (<i>e.g.</i>, <a href="/a/ptm_png.py?l=%s&amp;p=NAS" target="_blank">&amp;p=NAS</a>):
	<ol>
	<li><a href="/a/ptm_png.py?l=%s&amp;p=N" target="_blank">&amp;p=N</a>: N-terminal acetylation;</li>
	<li><a href="/a/ptm_png.py?l=%s&amp;p=A" target="_blank">&amp;p=A</a>: K-acetylation;</li>
	<li><a href="/a/ptm_png.py?l=%s&amp;p=S" target="_blank">&amp;p=S</a>: K-succinylation;</li>
	<li><a href="/a/ptm_png.py?l=%s&amp;p=P" target="_blank">&amp;p=P</a>: S/T/Y-phosphorylation;</li>
	<li><a href="/a/ptm_png.py?l=%s&amp;p=U" target="_blank">&amp;p=U</a>: K-GG or K-SUMO conjugation;</li>
	<li><a href="/a/ptm_png.py?l=%s&amp;p=M" target="_blank">&amp;p=M</a>: R-dimethylation;</li>
	<li><a href="/a/ptm_png.py?l=%s&amp;p=O" target="_blank">&amp;p=O</a>: 4-hydroxyproline or 5-hydroxylysine;</li>
	<li><a href="/a/ptm_png.py?l=%s&amp;p=C" target="_blank">&amp;p=C</a>: citrulline;</li>
	<li><a href="/a/ptm_png.py?l=%s&amp;p=G" target="_blank">&amp;p=G</a>: O-linked glycosylation; &amp;</li>
	<li><a href="/a/ptm_png.py?l=%s&amp;p=L" target="_blank">&amp;p=L</a>: Gla.</li>
	</ol></div>''' % (_l,_l,_l,_l,_l,_l,_l,_l,_l,_l,_l))

form = cgi.FieldStorage()
print('Content-type: text/html\n\n')
style = ''
try:
	style = form['s'].value
except:
	style = ''

ltype = 'log'
try:
	ltype = form['t'].value
except:
	ltype = 'log'
ptype = ''
try:
	ptype = form['p'].value
except:
	ptype = ''

label = ''
try:
	label = form['l'].value
except:
	print('There must be a protein accession value specified')
	exit()
filename = label
title = 'ω-mod %s PTMs (%s)' % (label,ltype)
y_axis = None
protein = get_protein(label)
make_ptm_csv(label,len(protein),title,protein,filename,y_axis,style,ltype,ptype)


