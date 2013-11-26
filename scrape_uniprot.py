#!/usr/bin/env python
"""
Taken from: https://bitbucket.org/dalloliogm/scrape-uniprot/src/8a10e6a4e2beab4045af2cf764998bf6a750f0ce/README.TXT?at=default

Given an Uniprot ID, scrape the Uniprot server and return a list of infos

$: python scrape_uniprot.py P12314
Name    HugoName    ...

"""

import mechanize
from html2text import html2text
import re
import time
import cookielib
import sys

# DEFINE YOUR EMAIL ADDRESS HERE
user_mail = 'hannahaf@hotmail.com'
if not user_mail:
    sys.exit("ERROR: Define your email in scrape_uniprot.py")

def initialize_browser():

    br = mechanize.Browser()
    # Cookie Jar
    cj = cookielib.LWPCookieJar()
    br.set_cookiejar(cj)

    # Browser options
    br.set_handle_equiv(True)
    br.set_handle_gzip(True)
    br.set_handle_redirect(True)
    br.set_handle_referer(True)
    br.set_handle_robots(False)

    # Follows refresh 0 but not hangs on refresh > 0
    br.set_handle_refresh(mechanize._http.HTTPRefreshProcessor(), max_time=1)

    # Want debugging messages?
#    br.set_debug_http(True)
#    br.set_debug_redirects(True)
#    br.set_debug_responses(True)

    # User-Agent (this is cheating, ok?)
    br.addheaders = [('User-agent', 'Mechanize Scraper. Contact: %s')]
    br.addheaders.append(('email', '%s'))

    return br

br = initialize_browser()

def scrape_uniprot(uniprot_id):
    baseurl = "http://www.uniprot.org/uniprot/"
    page = br.open(baseurl + uniprot_id)
    html = page.read()
    infos = []

    # get the line containing most information.
    # This will look like:
    # ['UniProt identifier</acronym></p><p><textarea name="query" rows="6" cols="61">&gt;sp|P12314|FCGR1_HUMAN High affinity immunoglobulin gamma Fc receptor I OS=Homo sapiens GN=FCGR1A PE=1 SV=2\n']   
 
    uniprot_info_re = re.compile('UniProt identifier</acronym>.*\n')
    uniprot_info = uniprot_info_re.findall(html)[0]

#    print uniprot_info
    parse_uniprot_info_re = re.compile('UniProt identifier</acronym></p><p><textarea name="query" rows="\d" cols="\d+">&gt;sp\|.*\|(.*_HUMAN) (.*) OS=Homo sapiens GN=(.*) PE=\d+ SV=\d+\n')
    print '!!!!!!!!!!!!!! prase: ' + str(parse_uniprot_info_re.findall(uniprot_info))
    info = parse_uniprot_info_re.findall(uniprot_info)
    if len(info) is not 0:
        (uniprot_name, full_name, hugo_name) = parse_uniprot_info_re.findall(uniprot_info)[0]
    else:
        full_name = " "
        uniprot_name = " "
        hugo_name = " "
    found = re.findall('''DR-Ensembl-Gene\'\);">(.*?)</a>''', html)
    if len(found) is not 0:
        ensembl_id = found[0]
    else: 
        ensembl_id = " "

#    function = re.findall('Function<\/acronym></td><td><p>(.*?) <a class="attribution"', html)[0]
#    subunit_structure = re.findall('Subunit structure<\/acronym></td><td><p>(.*?) <a class="attribution"', html)
    function = re.findall('Function<\/acronym></td><td><p>(.*?) *Subunit', html)[0]
    function = re.sub('<.*?>', '', function)
    function = re.sub('Ref\.\d+', '', function)
    subunit_structure = re.findall('Subunit structure<\/acronym></td><td><p>(.*?) <acronym', html)[0]
    subunit_structure = re.sub('<.*?>', '\t', subunit_structure)
    subunit_structure = re.sub('Ref\.\d+', '\t', subunit_structure)

    #molecular_function = re.findall('Molecular_function(.*?)Complete GO annotation...', subunit_structure)
    molecular_function = re.findall('Molecular function(.*?)Technical term', subunit_structure)
    mo_fun_arr = molecular_function[0].split('\t')
    functions = []
    for word in mo_fun_arr:
        if not (word.startswith('&') or word.strip()==''):
            functions.append(word)
    #for word in mo_fun_arr:
    #mo_fun_arr.remove('')
    #mo_fun_arr.remove('&(.*)')
    print "Molecular_function!!: " + str(functions)
    infos = [full_name, hugo_name, uniprot_id, uniprot_name, ensembl_id, function, subunit_structure]

    return functions #returns the list of functions
    #return html, infos


def test_scrape_uniprot():
#    (html, infos) = scrape_uniprot("P12314")
    (html, infos) = scrape_uniprot("Q9NNX6")
    print infos
    return html, infos


def scrape_list_ids(on_screen=False):
    """
    call scrape_uniprot() on a list of Uniprot IDs saved to a file
    """
    id_list_file = "./uniprot_ids.txt"
    id_list = []
    for line in open(id_list_file, 'r'):
        line = line.strip()
        if line:
            id_list.append(line)

    print "Full_Name\tHugo_Name\tUniprot_ID\tUniprot_name\tEnsembl_ID\tFuntion(Uniprot)\tSubunit_Structure(Uniprot)\n"
    for uniprot_id in id_list:
        (html, infos) = scrape_uniprot(uniprot_id)
        if on_screen: # debugging purposes
            print "-"*80 + "\n" + "%s\t"*5 % tuple(infos[:5]) + "\n\nFUNCTION: \n'%s'\n\nSUBUNITS: '%s'\n\n\n\n" % ( infos[5], infos[6])
        else:
            print "%s\t"*5 % tuple(infos[:5]) + "'%s'\t'%s'\n" % ( infos[5], infos[6])
        time.sleep(2)


if __name__ == '__main__':
#    (html, infos) = test_scrape_uniprot()
    scrape_list_ids(True)