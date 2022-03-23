#!/usr/bin/python3

import cgi,cgitb
import re
import requests
import json
import random
import sys
import subprocess
import time
import datetime
import re
import hashlib

proteome_AAA = {
'A':6.919235,
'C':2.183718,
'D':4.811060,
'E':7.131376,
'F':3.593881,
'G':6.522858,
'H':2.609194,
'I':4.328273,
'K':5.712237,
'L':9.873164,
'M':2.202668,
'N':3.584903,
'O':0.000000,
'P':6.297458,
'Q':4.770981,
'R':5.639044,
'S':8.385303,
'T':5.461958,
'U':0.000264,
'V':6.014760,
'W':1.229189,
'Y':2.634526,
}

def AAA(_seq):
	aaa = {}
	res = 'ACDEFGHIKLMNPOQRSTUVWY'
	for a in list(res):
		aaa[a] = 0
	for a in _seq:
		if a in aaa:
			aaa[a] += 1
		else:
			aaa[a] = 1
	return aaa

def trypsin(_seq):
	# a list of elegible sites
	ps = [-1]
	# residues cut by trypsin
	ss = set(['K','R'])
	# residues that block cleavage when C-terminal to site
	cbad = set(['P'])
	# length of protein
	lseq = len(_seq)
	seqs = list(_seq)
	peps = []
	# iterate through the sequence
	for i,res in enumerate(seqs):
		if i == lseq - 1:
			continue
		# if you just passed a cleavage site, act
		if i == 0 or (seqs[i-1] in ss and seqs[i] not in cbad):
			j = i + 1
			# is the next residue a cleavage site too
			if j < lseq-1 and seqs[j] in ss and seqs[j+1] not in cbad:
				peps.append({'seq':'%s' % (_seq[i:j+1]),'f':i+1,'l':j+1})
				if i == 0:
					j += 1
				else:
					if j == lseq - 1:
						peps.append({'seq':'%s' % (_seq[i:j+1]),'f':i+1,'l':j+1})
					continue
			# find the next cleavage site
			while j < lseq-1 and not (seqs[j] in ss and seqs[j+1] not in cbad):
				j += 1
			if j < lseq -2:
				peps.append({'seq':'%s' % (_seq[i:j+1]),'f':i+1,'l':j+1})
			j += 1
			# deal with the last residue cleavage problem
			if j < lseq-1 and seqs[j] in ss and seqs[i+1] not in cbad:
				peps.append({'seq':'%s' % (_seq[i:j+1]),'f':i+1,'l':j+1})
			elif j >= lseq - 1:
				peps.append({'seq':'%s' % (_seq[i:j+1]),'f':i+1,'l':lseq})
		else:
			pass
	# make sure everything is in order
	peps = [p for p in sorted(peps, key=lambda k: k['f'])]
	return peps

def get_description(_l):
	url = 'http://gpmdb.thegpm.org/1/protein/description/acc=%s' % (_l)
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

def get_interpro(_l):
	ds = list()
	if not _l:
		return ds
	url = 'http://gpmdb.thegpm.org/1/protein/domains/acc=%s' % (_l)
	session = requests.session()
	try:
		r = session.get(url,timeout=20)
	except requests.exceptions.RequestException as e:
		print(e)
		return None
	try:
		ds = json.loads(r.text)
	except:
		return None
	return ds


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

def get_domains(_s,_hex):
	url = 'http://gpmdb.thegpm.org/thegpm-cgi/tm_seq.py'
	session = requests.session()
	try:
		r = session.post(url,params={'s':_s},timeout=20)
	except requests.exceptions.RequestException as e:
		print(e)
		return None
	try:
		values = json.loads(r.text)
	except:
		return None
	return set(values[0]),values[1]

def get_highlights(_h,_s):
	rs = set()
	if not _h:
		return rs
	h = _h.strip()
	vs = h.split(',')
	for v in vs:
		v = v.strip()
		if v.find('re:') == 0:
			p = re.compile(v[3:])
			for m in p.finditer(_s):
				for j in range(m.start()+1,m.end()+1):
					rs.add(j)
		elif v.find('-') != -1:
			gs = v.split('-')
			for j in range(int(gs[0]),int(gs[1])+1):
				rs.add(j)
		else:
			try:
				rs.add(int(v))
			except:
				try:
					p = re.compile('[%s]' % v)
					for m in p.finditer(_s):
						for j in range(m.start()+1,m.end()+1):
							rs.add(j)
				except:
					continue
	return rs

def get_marked(_h,_s):
	rs = set()
	if not _h:
		return rs
	h = _h.strip()
	vs = h.split(',')
	for v in vs:
		v = v.strip()
		if v.find('re:') == 0:
			p = re.compile(v[3:])
			for m in p.finditer(_s):
				for j in range(m.start()+1,m.end()+1):
					rs.add(j)
		elif v.find('-') != -1:
			gs = v.split('-')
			for j in range(int(gs[0]),int(gs[1])+1):
				rs.add(j)
		else:
			p = re.compile('[%s]' % v)
			for m in p.finditer(_s):
				for j in range(m.start()+1,m.end()+1):
					rs.add(j)
	return rs

def print_seq(_s,_m,_b,_l,_h,_fe):
	if len(_s) == 0:
		return
	m = hashlib.sha3_256()
	m.update(_s.encode())
	hexdigest = m.hexdigest()
	ds,text = get_domains(_s,hexdigest)
	hs = get_highlights(_h,_s)
	red = get_marked(_m,_s)
	blue = get_marked(_b,_s)
	interpro = get_interpro(_l)
	display = '<div class="ex1"><hr width="650" style="margin-left: -20px;"/>'
	mem = ''
	highlight = ''
	marks = {}
	for i,r in enumerate(_s):
		if i != 0 and i % 50 == 0:
			display += '&nbsp;&nbsp;<span class="num">%i</span></div>\n<div class="ex1">' % (i)
		elif i % 10 == 0:
			display += ''
		mem = ''
		if i+1 in ds and _fe['grey']:
			mem = ' mem'
			if i+1 not in marks:
				marks[i+1] = {'grey':1,'green':0,'red':0,'blue':0}
			else:
				marks[i+1]['grey'] = 1
		highlight = ''
		if i+1 in hs and _fe['green']:
			highlight = ' highlight'			
			if i+1 not in marks:
				marks[i+1] = {'grey':0,'green':1,'red':0,'blue':0}
			else:
				marks[i+1]['green'] = 1
		if i+1 in red and _fe['red']:
			display += '<span class="red%s%s" title="%s %i">%s</span>' % (highlight,mem,r,i+1,r)
			if i+1 not in marks:
				marks[i+1] = {'grey':0,'green':0,'red':1,'blue':0}
			else:
				marks[i+1]['red'] = 1
		elif i+1 in blue and _fe['blue']:
			display += '<span class="blue%s%s" title="%s %i">%s</span>' % (highlight,mem,r,i+1,r)
			if i+1 not in marks:
				marks[i+1] = {'grey':0,'green':0,'red':0,'blue':1}
			else:
				marks[i+1]['blue'] = 1
		else:
			display += '<span class="unmarked%s%s" title="%s %i">%s</span>' % (highlight,mem,r,i+1,r)
	display = re.sub(r'\<div class=\"ex1\"\>$',r'',display)
	display += '</div>\n'
	print(display)
	print('<div class="ex1"><hr width="650" style="margin-left: -20px;"/></div>\n')
	if text:
		print('<div class="ex1"><u>TM domains:</u></div>\n')
		print(text)
	if interpro:
		print('<div class="ex1"><u>Interpro domains:</u></div>\n')
		for v in interpro:
			print('<div class="ex1">(%i-%i) %s </div>' % (v['b'],v['e'],v['ldesc']))
	print('<div class="ex1"><u>SHA3 256:</u><br>\n%s\n</div>\n' % (hexdigest))
	return marks

def print_form(_s,_m,_b,_l,_h,_fe):
	but = '<input type="submit" class="button" value="&#8635;" title="refresh display" />'
	print('<div class="ex1">')
	print('<form style="display: inline;" name="seq_form" action="/a/seq.py" METHOD="POST" ENCTYPE="multipart/form-data">')
	print('<hr width="650" style="margin-left: -20px;"/>')
	print('<input type="hidden" value="no" name="red" />')
	print('<input type="hidden" value="no" name="blue" />')
	print('<input type="hidden" value="no" name="green" />')
	print('<input type="hidden" value="no" name="grey" />')
	if _fe['grey']:
		print('<span class="mem"><input type="checkbox" id="grey_box" name="grey" value="yes" CHECKED/>&nbsp;TM domains</span>&nbsp;%s<br />' % (but))
	else:
		print('<span class="mem"><input type="checkbox" id="grey_box" name="grey" value="yes" />&nbsp;TM domains </span>&nbsp;%s<br />' % (but))
	if _fe['blue']:
		print('<span class="blue"><input type="checkbox" id="blue_box" name="blue" value="yes" CHECKED/>&nbsp;&nbsp;residues:</span>&nbsp;<input id="blue" name="b" size="20" value="%s" placeholder="ED" />&nbsp;%s<br/>' % (_b,but))
	else:
		print('<span class="blue"><input type="checkbox" id="blue_box" name="blue" value="yes" />&nbsp;&nbsp;residues:</span>:&nbsp;<input id="blue" name="b" size="20" value="%s" placeholder="ED" />&nbsp;%s<br/>' % (_b,but))
	if _fe['red']:
		print('<span class="red"><input type="checkbox" id="red_box" name="red" value="yes" CHECKED/>&nbsp;&nbsp;residues:</span>&nbsp;<input id="red" name="m" size="20" value="%s" placeholder="KR" />&nbsp;%s<br/>' % (_m,but))
	else:
		print('<span class="red"><input type="checkbox" id="red_box" name="red" value="yes"/>&nbsp;&nbsp;residues:</span>&nbsp;<input id="red" name="m" size="20" value="%s" placeholder="KR" />&nbsp;%s<br/>' % (_m,but))
	if _fe['green']:
		print('<span class="highlight"><input type="checkbox" id="green_box" name="green" value="yes" CHECKED/>&nbsp;&nbsp;&nbsp;&nbsp;ranges:</span>&nbsp;<input id="blue" name="h" size="20" value="%s" placeholder="1-20,25 or re:N[^P][ST]" />&nbsp;%s<br/>' % (_h,but))
	else:
		print('<span class="highlight"><input type="checkbox" id="green_box" name="green" value="yes" />&nbsp;&nbsp;&nbsp;&nbsp;ranges:</span>&nbsp;<input id="blue" name="h" size="20" value="%s" placeholder="1-20,25 or re:N[^P][ST]" />&nbsp;%s<br/>' % (_h,but))
	print('<span class="accession">&nbsp;&nbsp;&nbsp;accession:</span>&nbsp;<input id="label" name="l" size="20" value="%s" onChange="clearSeq();" placeholder="ENSP00000242786" />&nbsp;%s<br/>' % (_l,but))
	print('<span class="accession">&nbsp;&nbsp;&nbsp;&nbsp;sequence:</span><br/><textarea rows="10" cols="50" id="seq" name="s" placeholder="protein sequence">%s</textarea>&nbsp;%s' % (_s,but))
	print('</form></div>\n')
	
cgitb.enable()
form = cgi.FieldStorage()
print('Content-type: text/html\n\n')
seq = ''
form_entries = dict()
try:
	seq = form['s'].value.upper()
except:
	seq = ''
seq = re.sub(r'[^A-Z]+',r'',seq)
mark = ''
try:
	mark = form['m'].value
except:
	mark = ''
form_entries['red'] = True
blue = ''
try:
	blue = form['b'].value
except:
	blue = ''
form_entries['blue'] = True
label = ''
try:
	label = form['l'].value.upper()
except:
	label = ''
highlight = ''
try:
	highlight = form['h'].value
except:
	highlight = ''
boxes = ['grey','red','green','blue']
for b in boxes:
	try:
		if form.getvalue(b) == 'yes':
			form_entries[b] = True
		else:
			form_entries[b] = False
	except:
		form_entries[b] = True
marks = {}
if seq:
	marks = print_seq(seq,mark,blue,label,highlight,form_entries)
if not seq and label:
	seq = get_protein(label)
	marks = print_seq(seq,mark,blue,label,highlight,form_entries)
print_form(seq,mark,blue,label,highlight,form_entries)
peps = trypsin(seq)
if len(peps) == 0:
	exit()
print('<hr width="650" style="margin-left: -20px;"/><div id="tryptic_peptides">')
print('<table cellspacing="1" cellpadding="4">')
print('<tr><td width="50">start</td><td width="50">end</td><td>1&deg; tryptic peptide</td></tr>')
for p in peps:
	line = '<tr><td>%i</td><td>%i</td><td>' %  (p['f'],p['l'])
	for i,a in enumerate(p['seq']):
		c = p['f']+i
		if c in marks:
			if marks[c]['grey']:
				line += '<span class="mem">%s</span>' % (a)
			elif marks[c]['red']:
				line += '<span class="red">%s</span>' % (a)
			elif marks[c]['blue']:
				line += '<span class="blue">%s</span>' % (a)
			elif marks[c]['green']:
				line += '<span class="highlight">%s</span>' % (a)
		else:
			line += '%s' % (a)		
	line += '</td></tr>\n'
	print(line)
print('</table></div>')
aaa = AAA(seq)
taaa = sum(aaa.values())
print('<hr width="650" style="margin-left: -20px;"/><div id="aaa"><table cellspacing="1" cellpadding="2">')
print('<tr align="center" valign="top"><td width="100">Residue</td><td width="100">Count</td><td width="100">Fraction</td><td width="100">Proteome</td></tr>')
for a in sorted(aaa):
	print('<tr align="center" valign="top"><td>%s</td><td>%i</td><td>%.2f</td><td>%.2f</td></tr>' % (a,aaa[a],100.0*aaa[a]/taaa,proteome_AAA[a]))
print('</table></div>')

