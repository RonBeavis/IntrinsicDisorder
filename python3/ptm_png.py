#!/usr/bin/python3

import cgi,cgitb

import sys
import requests
import re
import json
import random
import datetime
from urllib.parse import quote

cgitb.enable()

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

def start_page(_l = '', _s = '',_t = '',_p = ''):
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
	<title>ID ω-mod diagram</title>
	<meta name="viewport" content="width=device-width,initial-scale=1" />
	<meta name="robots" content="index,nofollow,noarchive">''')
	print('''
		<meta property="og:locale" content="en_EN" />
		<meta property="og:type" content="website" />
		<meta property="og:title" content="Intrinsic Disorder PTMs" />
		<meta property="og:description" content="%s" />
		<meta property="og:url" content="https://intrinsicdisorder.com" />
		<meta property="og:site_name" content="Intrinsic Disorder" />
		<meta property="og:image:width" content="800" />
		<meta property="og:image:height" content="400" />''' % (desc))
	v = re.sub(r'[\|\:]',r'_',_l)
	print('''
		<meta property="og:image" content="https://intrinsicdisorder.com/ptm_png/%s_ptms.png" />
		<meta property="og:image:secure_url" content="https://intrinsicdisorder.com/ptm_png/%s_ptms.png" />
		'''  % (quote(v),quote(v)))
	print('''
		<meta name="twitter:url" content="https://intrinsicdisorder.com/a/ptm_png.py?l=%s&s=%s&p=%s">
		<meta name="twitter:domain" content="intrinsicdisorder.com">
		<meta name="twitter:card" content="summary_large_image" />
		<meta name="twitter:description" content="%s" />
		<meta name="twitter:title" content="Intrinsic Disorder PTMs - %s" />
		<meta name="twitter:image" content="https://intrinsicdisorder.com/ptm_png/%s_ptms.png" />
		'''  % (quote(re.sub(r'\|',r'~',_l)),quote(_s),quote(_p),desc,_l,quote(v)))
	print('''
<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Merriweather+Sans:300,300italic,regular,italic,600,600italic,700,700italic,800,800italic&amp;subset=latin,latin-ext" type="text/css" />
<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Merriweather:300,300italic,regular,italic,600,600italic,700,700italic,800,800italic&amp;subset=latin,latin-ext" type="text/css" />

<style media="screen" type="text/css">
	#content {
	  width: 800px ;
	  margin-left: auto ;
	  margin-right: auto ;
	}
	.pic {
	  width: 800px;
	  height: 400px;
	  margin-left: auto;
	  margin-right: auto;
	}
	.con {
	  width: 800px ;
	  margin-left: auto ;
	  margin-right: auto ;
	}
	.desc {
	  width: 800px ;
	  margin-left: auto ;
	  margin-right: auto ;
	  background-color: #f9eee5;
	  font-family:  "Merriweather", serif;
	}
	#copyright {
		font-size: 10pt;
	    margin-top: auto ;
	    margin-bottom: auto ;
	}
	ol	{ 
		padding-left:0; 
		list-style-position:inside; 
		line-height:140%;
	}
	body {
		color: #000000;
		background-color: #FFFFFF;
		font-weight: normal;
		font-family:  "Merriweather Sans", Helvetica, Arial, sans-serif;
		font-size: 12pt;
		margin: 3px 3px 3px 3px;
		text-align:left;
	}
	a {
		color: #369;
		text-decoration: none;
	}
	.bluesq {
		color: #369;
		text-decoration: none;
	}
	.tots {
		background:#000000;
		color: #FFFFFF;
	}
	.heads {
		background:#f9eee5;
	}
	.alt {
		background:#e5ecf9;
	}
	p {
		display: block;
		margin-top: 0.5em;
		margin-bottom: 0.5em;
		margin-left: 5em;
		margin-right: auto;
	}
	table { 
		border: none;
		border-collapse: collapse;
		font-family:  "Merriweather Sans", Helvetica, Arial, sans-serif;
		font-size: 12pt;
	}
	table td { border-left: 1px solid #000; text-align: center;	vertical-align: top; padding: 3px 10px 3px 10px;}
	table td:first-child { border-left: none; }</style>
<SCRIPT LANGUAGE="JavaScript" type="text/javascript"> 
function fill_main(url,style,ltype,ptype)
{
	fetch("/a/ptm_png_2.py?l="+url+"&s="+style+"&t="+ltype+"&p="+ptype)
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
	print('''

<body onLoad="fill_main('%s','%s','%s','%s')">
<div id="main_body">'''  % (_l,_s,_t,_p))
	print('<div id="diagram" class="pic"><center><img src="/pics/loading%i.gif"/></center></div>' % (random.randint(2,5)))

def end_page():
	print('''</div><br /><p id="copyright">%s Intrinsic Disorder</p>''' % (datetime.datetime.now()))
	print('''</body>\n</html>\n''')

form = cgi.FieldStorage()
print('Content-type: text/html\n\n')
style = ''
try:
	style = form['s'].value
except:
	style = ''
ltype = ''
try:
	ltype = form['t'].value
except:
	ltype = ''
ptype = ''
try:
	ptype = form['p'].value
except:
	ptype = ''
try:
	label = re.sub(r'~',r'|',form['l'].value)
	start_page(label,style,ltype,ptype)
except:
	start_page('',style,ltype,ptype)
	print('There must be a protein accession value specified')
	end_page()
	exit()
end_page()

