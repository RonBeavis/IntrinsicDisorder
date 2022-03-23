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

def print_top(_l,_url):
	desc = "Sequence viewer"
	if _l:
		desc = get_description(_l)
	max = 100
	if desc:
		desc =  '%s' % (desc)
		if len(desc) > max:
			desc = re.sub(r'\s+\[.+','',desc)
		if len(desc) > max:
			desc = desc[0:max-2] + '…'

	print('''<!DOCTYPE html>
		<html lang="en" class="no-js">
		<head>
		<meta http-equiv="X-UA-Compatible" content="IE=edge">
		<meta charset="utf-8">
		<title>☙ 1&deg; sequence display</title>
		<meta name="viewport" content="width=device-width,initial-scale=1" />
		<meta name="robots" content="index,nofollow,noarchive">''')
	print('''
		<meta property="og:locale" content="en_EN" />
		<meta property="og:type" content="website" />
		<meta property="og:title" content="Intrinsic Disorder sequence display" />
		<meta property="og:description" content="%s" />
		<meta property="og:url" content="https://intrinsicdisorder.com" />
		<meta property="og:image:width" content="800" />
		<meta property="og:image:height" content="400" />
		<meta property="og:image" content="https://intrinsicdisorder.com/pics/sq.png" />
		<meta property="og:image:secure_url" content="https://intrinsicdisorder.com/pics/sq.png" />
		''' % (desc))
	v = re.sub(r'[\|\:]',r'_',_l)
	print('''
		<meta name="twitter:url" content="https://intrinsicdisorder.com/a/seq.py?l=%s">
		<meta name="twitter:domain" content="intrinsicdisorder.com">
		<meta name="twitter:card" content="summary_large_image" />
		<meta name="twitter:site" content="@norsivaeb" />
		<meta name="twitter:description" content="%s" />
		<meta name="twitter:title" content="Intrinsic Disorder sequence display - %s" />
		<meta name="twitter:image" content="https://intrinsicdisorder.com/pics/sq.png" />
		'''  % (re.sub(r'\|',r'~',_l),desc,_l))
	print('''<script type="text/javascript">
		<!--
		function clearSeq()
		{
			var obj = document.getElementById("seq");
			obj.value="";
		}
		-->\n</script>''')
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
<SCRIPT LANGUAGE="JavaScript" type="text/javascript"> 
function fill_main(url)
{
	fetch("/a/seq2.py?"+url)
		.then(function(response) {
		 return response.text();
	})
  .then(function(text) {
    document.getElementById("main_body").innerHTML =  text;
  })
		.catch(function() {
	  console.log(url);
	});
}
</script>
	</head>''')
	print('''\n<body onLoad="fill_main('%s');">
	<div id="main"><div id="main_body" class="cdiv">\n''' % (_url))
	print('<img src="/pics/loading%i.gif"/>' %  (random.randint(2,5)))
	return

def print_bottom():
	print('</div><p id="copyright">%s Intrinsic Disorder</p>' % (datetime.datetime.now()))

	t = '''</div></body></html>\n'''
	print(t)
	return


cgitb.enable()
form = cgi.FieldStorage()
print('Content-type: text/html\n\n')
seq = ''
url = ''
try:
	seq = form['s'].value.upper()
except:
	seq = ''
seq = re.sub(r'[^A-Z]+',r'',seq)
url += 's=' +seq
mark = ''
try:
	mark = form['m'].value
except:
	mark = ''
url += '&m=' + mark
blue = ''
try:
	blue = form['b'].value
except:
	blue = ''
url += '&b=' + blue
label = ''
try:
	label = re.sub('~',r'|',form['l'].value.upper())
except:
	label = ''
url += '&l=' + label
highlight = ''
try:
	highlight = form['h'].value
except:
	highlight = ''
url += '&h=' + highlight
boxes = ['grey','red','green','blue']
for b in boxes:
	try:
		if form[b].value == 'yes':
			url += '&' + b + '=yes'
		else:
			url += '&' + b + '='
	except:
			url += '&' + b + '=yes'
print_top(label,url)
print_bottom()
