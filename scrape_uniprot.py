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
    #print 'UNPROT: ' + str(uniprot_id)
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

    functions = []
    #print 'html: ' + str(html)
    molecular = re.findall('Molecular_function(.*)Complete GO annotation...', html)
    #subunit = re.findall('Subunit structure<\/acronym></td><td><p>(.*?) <acronym', html)
    #print '!!!!!!!!!!!!!!!!sub unit: ' + str(subunit_structure)

    #print 'subunit: ' + str(subunit)
    if len(molecular) is not 0:
        molecular_function = molecular[0] 
        molecular_function = re.sub('<.*?>', '\t', molecular_function)
        molecular_function = re.sub('Ref\.\d+', '\t', molecular_function)


        #molecular_function = re.findall('Molecular_function(.*?)Complete GO annotation...', subunit_structure)
        #print 'molecular_function: ' + str(molecular_function)
        #molecular_function = re.findall('Molecular function(.*?)Technical term', subunit_structure)
        
        mo_fun_arr = []
        if len(molecular_function) is not 0:
            mo_fun_arr = molecular_function.split('\t')
        functions = []
        for word in mo_fun_arr:
            if not (word.startswith('. Source:') or word.strip()=='' or word.startswith('Inferred')
                or word.startswith('PubMed')):
                functions.append(word)
        
    return functions #returns the list of functions
    #return html, infos


def test_scrape_uniprot():
#    (html, infos) = scrape_uniprot("P12314")
    (html, infos) = scrape_uniprot("Q9NNX6")
    print infos
    return html, infos


def scrape_list_ids(on_screen=False):
    function_dict = {}
    all_functions = []
    """
    call scrape_uniprot() on a list of Uniprot IDs saved to a file
    """
    id_list_file = "./uniprot_ids.txt"
    id_list = []
    for line in open(id_list_file, 'r'):
        line = line.strip()
        if line:
            id_list.append(line)

    #print "Full_Name\tHugo_Name\tUniprot_ID\tUniprot_name\tEnsembl_ID\tFuntion(Uniprot)\tSubunit_Structure(Uniprot)\n"
    for uniprot_id in id_list:
        old_id = uniprot_id
        uniprot_id = uniprot_id.split('_')[0]
        new_list = scrape_uniprot(uniprot_id)
        function_dict[old_id] = new_list

        for item in new_list:
            if item not in all_functions:
                all_functions.append(item)
       
    return [function_dict, all_functions]

    

if __name__ == '__main__':
    scrape_list_ids(True)