#!/usr/bin/python3

import cgi,cgitb
import sys
import requests
import re
import json
import datetime

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

def get_peptides(_e):
	url = 'http://gpmdb.thegpm.org/protein/model/%s&excel=1' % (_e)
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

def print_top(_l,_l2):
	desc = "Sequence overlap %s ∩ %s" % (_l,_l2)

	print('''<!DOCTYPE html>
		<html lang="en" class="no-js">
		<head>
		<meta http-equiv="X-UA-Compatible" content="IE=edge">
		<meta charset="utf-8">
		<title>Sequence overlap/intersection display</title>
		<meta name="viewport" content="width=device-width,initial-scale=1" />
		<meta name="robots" content="index,nofollow,noarchive">''')
	print('''
		<meta property="og:locale" content="en_EN" />
		<meta property="og:type" content="website" />
		<meta property="og:title" content="GPMDB Sequence overlap/intersection display" />
		<meta property="og:description" content="%s" />
		<meta property="og:url" content="https://intrinsicdisorder.com" />
		<meta property="og:image:width" content="800" />
		<meta property="og:image:height" content="400" />
		<meta property="og:image" content="https://intrinsicdisorder.com/pics/ov.png" />
		<meta property="og:image:secure_url" content="https://intrinsicdisorder.com/pics/ov.png" />
		''' % (desc))
	v = re.sub(r'[\|\:]',r'_',_l)
	print('''
		<meta name="twitter:url" content="https://intrinsicdisorder.com/a/overlap.py?l=%s&l2=%s">
		<meta name="twitter:domain" content="intrinsicdisorder.com">
		<meta name="twitter:card" content="summary_large_image" />
		<meta name="twitter:site" content="@norsivaeb" />
		<meta name="twitter:description" content="%s" />
		<meta name="twitter:title" content="Intrinsic Disorder observed residue overlap - %s" />
		<meta name="twitter:image" content="https://intrinsicdisorder.com/pics/ov.png" />
		'''  % (_l,_l2,desc,desc))
	print('''
		<style media="screen" type="text/css">
		@font-face	{
			font-family: 'Anonymous Pro';
			font-style: normal;
			font-weight: 400;
			src: local('Anonymous Pro'), local('Anonymous Pro-Regular'), url('/fonts/AnonymousPro-Regular.ttf');
			format('ttf');
		}
		@font-face	{
			font-family: 'Anonymous Pro';
			font-style: normal;
			font-weight: 700;
			src: local('Anonymous Pro-Bold'), url('/fonts/AnonymousPro-Bold.ttf');
			format('ttf');
		}
		@font-face	{
			font-family: 'Anonymous Pro';
			font-style: italic;
			font-weight: 400;
			src: local('Anonymous Pro-Italic'), url('/fonts/AnonymousPro-Italic.ttf');
			format('ttf');
		}
		body {
			color: #000000;
			background-color: #FFFFFF;
			font-weight: normal;
			font-family: "Anonymous Pro",serif;
			font-size: 13pt;
			margin: auto;
		}
		.cdiv	{
			  display: table;
			  margin: 0 auto;
		}
		.num	{
			color: grey;
			font-size: 11pt;
		}
		.red	{
			background-color: #ff6666;
			color: white;
			border: 1px solid white;
			border-radius: 5px;
			cursor: pointer;
		}
		.blue	{
			background-color: #6666ff;
			color: white;
			border: 1px solid white;
			border-radius: 5px;
			cursor: pointer;
		}
		.accession	{
			background-color: #996633;
			color: white;
			border: 1px solid white;
			border-radius: 5px;
			cursor: pointer;
		}
		.unmarked	{
			color: grey;
			border: 1px solid white;
			border-radius: 5px;
			cursor: pointer;
		}
		.mem	{
			background-color: #aaaaaa;
			color: #FFFFFF;
		}
		.highlight	{
			background-color: #00cc99;
			color: white;
		}
		div.ex1	{
			margin: 3px 3px 3px 3px;
		}
		.button {
		  font-weight: normal;
		  position: relative;
		  background-color: #4CAF50;
		  border: none;
		  font-size: 16px;
		  color: #FFFFFF;
		  padding: 3px;
		  width: 30px;
		  text-align: center;
		  -webkit-transition-duration: 0.4s; /* Safari */
		  transition-duration: 0.4s;
		  text-decoration: none;
		  overflow: hidden;
		  cursor: pointer;
		  margin-bottom: 8px;
		  margin-top: 8px;

		}

		.button:after {
			font-weight: normal;
		  content: "";
		  background: #f1f1f1;
		  display: block;
		  position: absolute;
		  padding-top: 300%;
		  padding-left: 10%;
		  margin-left: -20px !important;
		  margin-top: -120%;
		  opacity: 0;
		  transition: all 0.8s;
		}

		.button:active:after {
			font-weight: normal;
		  padding: 0;
		  margin: 0;
		  opacity: 1;
		  transition: 0s
		}
	</style>
	</head>''')
	print('''\n<body>
	<div id="main"><div id="main_body" class="cdiv">\n''')
	return

def print_bottom():
	print('<p id="copyright">%s Intrinsic Disorder</p>' % (datetime.datetime.now()))

	t = '''</div></div></body></html>\n'''
	print(t)
	return

def print_form(_l1,_l2):
	but = '<input type="submit" class="button" value="&#8635;" title="refresh display" />'
	print('<div class="ex1">')
	print('<form style="display: inline;" name="seq_form" action="/a/overlap.py/" METHOD="GET">')
	print('<input id="red" name="l" size="20" value="%s" placeholder="ENSP0...." /> overlaps on <input name="l2" size="20" value="%s" placeholder="ENSP0...." />' % (_l1,_l2))
	print(but)
	print('</form>')

cgitb.enable()

form = cgi.FieldStorage()
print('Content-type: text/html\n\n')
label1 = None
label2 = None
try:
	label1 = form['l'].value
except:
	label1 = ''
try:
	label2 = form['l2'].value
except:
	label2 = ''

print_top(label1,label2)
print_form(label1,label2)

if label1 == '' or label2 == '':
	print('</div>\n</div>\n</body>\n</html>\n')
	exit()

ls = get_peptides(label1)
protein_1 = {}
for l in ls:
	v = l.strip()
	if v.find('Sequence') != -1:
		continue
	vs = v.split('\t')
	if len(vs) < 5:
		continue
	sj = re.sub(r'[LI]',r'J',vs[5])
	if sj not in protein_1:
		protein_1[sj] = [(int(vs[0]),int(vs[1]))]
	else:
		protein_1[sj].append((int(vs[0]),int(vs[1])))

protein_2 = {}
ls = get_peptides(label2)

mset = set()

for l in ls:
	v = l.strip()
	if v.find('Sequence') != -1:
		continue
	vs = v.split('\t')
	if len(vs) < 5:
		continue
	sj = re.sub(r'[LI]',r'J',vs[5])
	if sj not in protein_1:
		continue
	if sj not in protein_2:
		protein_2[sj] = [(int(vs[0]),int(vs[1]))]
	else:
		protein_2[sj].append((int(vs[0]),int(vs[1])))

for p in protein_2:
	ls = protein_2[p]
	for l in ls:
		for i in range(l[0],l[1]+1):
			mset.add(i)

seq = get_protein(label2)
des = re.sub('\[.+',r'',get_description(label1))
print('<p><a href="/protein/model/%s" target="_blank">%s</a>&mdash;%s <br/>' % (label1,label1,des))
print('overlap with observed peptides shown on <br/>')
des = re.sub('\[.+',r'',get_description(label2))
print('<a href="/protein/model/%s" target="_blank">%s</a>&mdash;%s</p>' % (label2,label2,des))
print('<p><span class="red">overlapping residues</span>: %i<br/>total residues: %i<br/><span class="red">overlap</span>: %.1f%%</p>' % (len(mset),len(seq),100*(len(mset)/len(seq))))
#print('<pre>')
display = '<div class="ex1">\n'
for i,s in enumerate(seq):
	if i != 0 and i % 50 == 0:
			display += '&nbsp;&nbsp;<span class="num">%i</span></div>\n<div class="ex1">' % (i)
	if i+1 in mset:
		display += '<span class="red" title="%s %i">%s</span>' % (s,i+1,s)
	else:
		display += '<span class="unmarked" title="%s %i">%s</span>' % (s,i+1,s)
display += '</div>\n'
print(display)
print_form(label1,label2)
peps = trypsin(seq)
if len(peps) == 0:
	exit()
print('<hr width="650" style="margin-left: -20px;"/><div id="tryptic_peptides">')
print('<table cellspacing="1" cellpadding="2">')
print('<tr><td width="50">start</td><td width="50">end</td><td>1&deg; tryptic peptide</td></tr>')
full = 0
full_aa = 0
for p in peps:
	line = '<tr><td>%i</td><td>%i</td><td>' %  (p['f'],p['l'])
	cf = 0
	for i,a in enumerate(p['seq']):
		c = p['f']+i
		if c in mset:
			line += '<span class="red">%s</span>' % (a)
			cf += 1
		else:
			line += '%s' % (a)
	if cf == len(p['seq']):
		full += 1
		full_aa += len(p['seq'])
	line += '</td></tr>\n'
	print(line)
print('<tr><td></td><td></td><td>peptides: %i/%i (%.1f%%), residues: %i/%i (%.1f%%)</td>' % (full,len(peps),100.0*full/len(peps),full_aa,len(seq),100*full_aa/len(seq)))
print('</table></div>')

print_bottom()

