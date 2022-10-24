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
	desc = "Exon viewer"
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
		<title>Exon sequence display</title>
		<meta name="viewport" content="width=device-width,initial-scale=1" />
		<meta name="robots" content="index,nofollow,noarchive">''')
	print('''
		<meta property="og:locale" content="en_EN" />
		<meta property="og:type" content="website" />
		<meta property="og:title" content="Intrinsic Disorder exon display" />
		<meta property="og:description" content="%s" />
		<meta property="og:url" content="https://intrinsicdisorder.com" />
		<meta property="og:image:width" content="800" />
		<meta property="og:image:height" content="400" />
		<meta property="og:image" content="https://intrinsicdisorder.com/pics/ex.png" />
		<meta property="og:image:secure_url" content="https://intrinsicdisorder.com/pics/ex.png" />
		''' % (desc))
	v = re.sub(r'[\|\:]',r'_',_l)
	print('''
		<meta name="twitter:url" content="https://intrinsicdisorder.com/a/exon.py?l=%s">
		<meta name="twitter:domain" content="intrinsicdisorder.com">
		<meta name="twitter:card" content="summary_large_image" />
		<meta name="twitter:site" content="@norsivaeb" />
		<meta name="twitter:description" content="%s" />
		<meta name="twitter:title" content="Intrinsic Disorder exon display - %s" />
		<meta name="twitter:image" content="https://intrinsicdisorder.com/pics/ex.png" />
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
			font-family: 'Amino New';
			font-style: normal;
			font-weight: 400;
			src: local('Amino New'), url('/fonts/AMINONEW.TTF');
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
			background-color: #ffffff;
			color: red;
			border: 1px solid white;
			border-radius: 5px;
			cursor: pointer;
		}
		.blue	{
			background-color: #ffffff;
			color: blue;
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
			background-color: #ff00ff;
			color: white;
			border: 1px solid #ff00ff;
			border-radius: 5px;
			cursor: pointer;
		}
		.xblue	{
			background-color: #0000ff;
			color: white;
			border: 1px solid #0000ff;
			border-radius: 5px;
			cursor: pointer;
		}
		.xred	{
			background-color: #ff0000;
			color: white;
			border: 1px solid #ff0000;
			border-radius: 5px;
			cursor: pointer;
		}
		.mem	{
			background-color: #aaaaaa;
			color: #FFFFFF;
		}
		.bar	{
			background-color: #dddddd;
			height: 10px;
			width: 600px;
			display: inline;
			margin: 0 auto;
			font-size: 8px;
		}
		p {
			display: block;
			margin-top: 0.5em;
			margin-bottom: 0.5em;
			margin-left: 5em;
			margin-right: auto;
		}
		.desc {
	  		width: 650px ;
			font-size: 12pt;
	  		margin-left: auto ;
	  		margin-right: auto ;
	  		background-color: #f9eee5;
	  		font-family:  "Merriweather", serif;
		}
		.highlight	{
			background-color: #00cc99;
			color: white;
		}
		div.ex1	{
			font-family: "Amino New";
			margin: 3px 3px 3px 3px;
		}
		div.ex2	{
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
	fetch("/a/exon2.py?"+url)
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
label = ''
try:
	label = re.sub('~',r'|',form['l'].value.upper())
except:
	label = ''
url += '&l=' + label
print_top(label,url)
print_bottom()
