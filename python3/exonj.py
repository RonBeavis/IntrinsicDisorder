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
	s = 0
	d = 0
	vs = []
	for e in _e['mappings']:
		js = {}
		v = e['end']-e['start']+1
		x = v/3
		d += math.ceil(x) - x
		js['start'] = int(math.floor(s))+1
		js['end'] = int(math.ceil(s+x))
		js['first'] = _s[int(math.floor(s))]
		js['last'] = _s[int(math.ceil(s+x)-1)]
		vs.append(js)
		s += v/3
	return vs

def print_seq(_l):
	seq = get_sequence(_l)
	if not seq:
		return None
	js = get_exons(_l,seq)
	if not js:
		return None
	exons = get_domains(js,seq)
	if not js:
		return None
	desc = get_description(_l)
	if not desc:
		return None
	jv = {"info":{"accession":"%s" % (_l)}}
	jv['info']['sequence'] = seq
	jv['info']['description']= desc
	jv['info']['assembly'] = js['mappings'][0]['assembly_name']
	jv['info']['chromosome'] = js['mappings'][0]['seq_region_name']
	jv['info']['strand'] = js['mappings'][0]['strand']
	jv['map'] = []
	for i,j in enumerate(js['mappings']):
		ex = {}
		ex['order'] = i+1
		ex['gstart'] = j['start']
		ex['gend'] = j['end']
		ex['pstart'] = exons[i]['start']
		ex['pend'] = exons[i]['end']
		ex['pfirst'] = exons[i]['first']
		ex['plast'] = exons[i]['last']
		jv['map'].append(ex)
	print(json.dumps(jv))
#	print(js)
#	print(exons)
#	print(desc)
	return True

cgitb.enable()
form = cgi.FieldStorage()
print('Content-type: application/json\n\n')
label = ''
try:
	label = form['l'].value.upper()
except:
	label = ''
if label:
	if print_seq(label) is None:
		print('{"error":"required information not available"}')
else:
	print('{"error":"no label (parameter: l) supplied"}')

