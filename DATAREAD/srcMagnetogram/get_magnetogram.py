#!/usr/bin/env python3
############ Prototype: make_submission.py by Maksym Petrenko
############ CMERAL modification of make_submission.py
############ 30/11/2022
############ This python script correspond to a part of the make_submission.py
############ that only download the latest GONG magnetogram for realtime
############ simulation in the run directory as fitsfile.fits
###############################################################################################################

from html.parser import HTMLParser
import requests
import re
import os
import shutil
import gzip
import argparse
from datetime import datetime,timedelta

#modify to change the output directory (run_realtime directory used for realtime simulations)
rundir = os.getcwd()

HEADERS = {"User-Agent":"Mozilla/5.0 (Macintosh Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/97.0.4692.71 Safari/537.36"}

class LinkScrape(HTMLParser):
    def reset(self):
        super().reset()
        self.links = []
        self.in_link = False
    def handle_starttag(self, tag, attrs):
        self.in_link = False
        if tag == 'a':
            self.in_link = True
            for (name, value) in attrs:
                if name == 'href':
                    self.links.append({'href': value})
    def handle_endtag(self, attrs):
        self.in_link = False
    def handle_data(self, data):
        if self.in_link: self.links[-1]['text'] = data
    def clean(self):
        self.links = []

def get_path(parsed_datetime, UseGong):
    """
    Returns "path" - URL before the last slash, for the "parsed" datetime
    """
    if UseGong:
        path = 'https://gong.nso.edu/data/magmap/QR/zqs/'+ \
            parsed_datetime.strftime("%Y%m") + \
            '/mrzqs' +  parsed_datetime.strftime("%y%m%d") + '/'
    else:
        path = f'https://iswa.ccmc.gsfc.nasa.gov/iswa_data_tree/' \
            f'observation/solar/gong/mrzqs/'
        path = path.rstrip('/') + '/' + parsed_datetime.strftime("%Y") \
        + '/' + parsed_datetime.strftime("%m") + '/'
    return path

def get_highest(page_url, pattern, datetimein):
    response = requests.get(page_url)
    page_html = response.text
    link_parser = LinkScrape()
    link_parser.feed(page_html)
    links = link_parser.links
    link_parser.clean()

    last_match = ''
    lats_link = ''
    last_text = ''
    for link in links:
        original_url = link['href']
        text = link['text']
        matches = re.search(pattern, text)
        if (matches):
            if (matches.group(1)>last_match and matches.group(1)<=datetimein):
                last_match = matches.group(1)
                lats_link = original_url
                last_text = text
    return [last_match, last_text, lats_link]

def download_file(url, save_path):
    try:
        # Send GET request to the URL
        response = requests.get(url)
        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            # Write the content of the response to a local file
            with open(save_path, 'wb') as file:
                file.write(response.content)
        else:
            print(
                f"Failed to download file. Status code: {response.status_code}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Use python3 get_magnetogram.py <datetime as yymmddthhmm>")
    parser.add_argument('datetimein', help=
                        "DateTime in format yymmdd't'hhmm")
    parser.add_argument('--UseGong', action='store_true', help=
                        'Use GONG magnetogram repository (ISWA is default)')
    parser.add_argument('--SaveUrl', action='store_true', help=
                        'Save URL to the file')
    args = parser.parse_args()
    datetimein = args.datetimein
    UseGong = args.UseGong
    SaveUrl = args.SaveUrl
    parsed_datetime =  datetime.strptime(datetimein, "%y%m%dt%H%M")
    # Find proper folder with the magnetogram link
    path = get_path(parsed_datetime, UseGong)
    [cr, text, link] = get_highest(
        path, r'(\d\d\d\d\d\dt\d\d\d\d)', datetimein)
    if not link:
        #
        ### No magnetogram in the folder for this day. Look in the previous day
        #
        parsed_datetime -= timedelta(hours=23)
        path = get_path(parsed_datetime, UseGong)
        [cr, text, link] = get_highest(
            path,r'(\d\d\d\d\d\dt\d\d\d\d)',datetimein)
    granule_url = path.rstrip('/') + '/' + link
    # Adjust input files
    fits_file = os.path.join(rundir, "fitsfile.fits")
    fits_file_gz = os.path.join(rundir, link)
    download_file(granule_url, fits_file_gz)
    with gzip.open(fits_file_gz, 'rb') as f_in:
        with open(fits_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.unlink(fits_file_gz)
    if SaveUrl:
        url_file = os.path.join(rundir, "URL.gong")
        fid = open(url_file,'w')
        fid.write(str(granule_url))
        fid.close()
