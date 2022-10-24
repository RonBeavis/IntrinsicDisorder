#!/usr/bin/python3

#!c:/python36/python.exe


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
import math

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

def get_exons(_l,_s):
	url = 'https://rest.ensembl.org/map/translation/%s/1..%i?content-type=application/json' % (_l,len(_s))
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
	return values

def get_sequence(_l):
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

def get_domains(_e,_s):
	s = 1
	e = 1
	d = 0
	vs = []
	bases = 0
	for ex in _e['mappings']:
		js = {}
		bases += ex['end']-ex['start']+1
		e = math.ceil(bases/3)
		js['start'] = s
		js['end'] = e
		js['first'] = _s[s-1]
		js['last'] = _s[e-1]
		js['aa'] = bases/3
		js['bases'] = bases
		js['b_start'] = ex['start']
		js['b_end'] = ex['end']
		vs.append(js)
		if bases % 3 == 0:
			s = e + 1
		else:
			s = e
	return vs

def print_seq(_l,_s):
	if len(_s) == 0:
		return
	js = get_exons(_l,_s)
	if js is None:
		print('<p>Error: exon information not available</p>')
		return None
	exons = get_domains(js,_s)
	display = '<div class="ex1"><hr width="650" style="margin-left: -20px;"/>'
	ex = 0
	toggle = 0
	for i,r in enumerate(_s):
		if i != 0 and i % 50 == 0:
			display += '&nbsp;&nbsp;<span class="num">%i</span></div>\n<div class="ex1">' % (i)
		elif i % 10 == 0:
			display += ''
		if i+1 >= exons[ex]['start'] and i+1 < exons[ex]['end']:
			if toggle == 0:
				display += '<span class="red" title="%s %i, %i">%s</span>' % (r,i+1,ex+1,r)
			else:
				display += '<span class="blue" title="%s %i, %i">%s</span>' % (r,i+1,ex+1,r)
		elif i+1 == exons[ex]['end']:
			if ex < len(exons) - 1:
				if exons[ex]['end'] == exons[ex+1]['start']:
					display += '<span class="unmarked" title="%s %i, %i+%i">%s</span>' % (r,i+1,ex+1,ex+2,r)
				else:
					if toggle == 0:
						display += '<span class="xred" title="%s %i, %i">%s</span>' % (r,i+1,ex+1,r)
					else:
						display += '<span class="xblue" title="%s %i, %i">%s</span>' % (r,i+1,ex+1,r)
				if toggle == 0:
					toggle = 1
				else:
					toggle = 0
				ex += 1
			elif ex == len(exons) - 1:
				if toggle == 0:
					display += '<span class="xred" title="%s %i, %i">%s</span>' % (r,i+1,ex+1,r)
				else:
					display += '<span class="xblue" title="%s %i, %i">%s</span>' % (r,i+1,ex+1,r)
				ex += 1
	display = re.sub(r'\<div class=\"ex1\"\>$',r'',display)
	display += '</div>\n'
	print(display)
	return exons

def print_form(_l,_x):
	but = '<input type="submit" class="button" value="&#8635;" title="refresh display" />'
	print('<div class="ex2">')
	print('<form style="display: inline;" name="exon_form" action="/a/exon.py" METHOD="POST" ENCTYPE="multipart/form-data">')
	print('<hr width="650" style="margin-left: -20px;"/>')
	if _x is not None:
		desc = get_description(_l)
		print('<p class="desc">%s. (<a href="/a/exonj.py?l=%s" target="_blank" title="Exon boundary information in JSON format">JSON</a>)</p>' % (desc,_l))
		print('''<p class="desc">Exons are indicated by alternating <span class="red">RED</span> and <span class="blue">BLUE</span> residues.<br/>
			<span class="unmarked">X</span> indicates a boundary involving bases from both exons, otherwise<br/>
			<span class="xblue">X</span> or <span class="xred">X</span> indicates the last residue in an exon.</p>''')
	print('<span class="accession">&nbsp;&nbsp;&nbsp;accession:</span>&nbsp;<input id="label" name="l" size="20" value="%s" placeholder="ENSP00000242786" />&nbsp;%s<br/>' % (_l,but))
	print('</form></div>\n')

cgitb.enable()
form = cgi.FieldStorage()
print('Content-type: text/html\n\n')
label = ''
try:
	label = form['l'].value.upper()
except:
	label = ''
xs = None
if label:
	seq = get_sequence(label)
	xs = print_seq(label,seq)
print_form(label,xs)
