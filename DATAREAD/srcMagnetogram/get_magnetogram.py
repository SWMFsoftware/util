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

ISWA_DATA_URL = 'https://nispdata.nso.edu/ftp/QR/zqs/'
#202609/mrzqs260901/

#modify to change the output directory (run_realtime directory used for realtime simulations)
OUTPUT_BASE_PATH = os.getcwd()

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
            print(f"File downloaded successfully: {save_path}")
        else:
            print(f"Failed to download file. Status code: {response.status_code}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=
                                     "python3 get_magnetogram.py datetimein")
    parser.add_argument('datetimein', help=
                        "Date_Time in the format yymmdd't'hhmm")
    args = parser.parse_args()
    datetimein = str(args.datetimein)
    matches=re.search(r'(\d\d\d\d\d\dt\d\d\d\d)',datetimein)
    year = '20'+matches.group(1)[0:2]
    month = matches.group(1)[2:4]
    day = matches.group(1)[0:6]
    day_url = ISWA_DATA_URL.rstrip('/')+'/'+\
        str(year)+str(month)+'/mrzqs'+str(day)+'/'
    [cr, text, link] = get_highest(
        day_url,r'(\d\d\d\d\d\dt\d\d\d\d)',matches.group(1))
    if not link:
        #
        ### No magnetogram in the folder for this day. Look in the previous day
        #
        ### Convert datetime input string to the python datetime object
        format_pattern = "%y%m%dt%H%M"
        parsed_datetime = datetime.strptime(datetimein, format_pattern)
        #
        ### Reduce by one day
        parsed_datetime = parsed_datetime - timedelta(hours=23)
        #
        ### Convert to a string
        formatted_string = parsed_datetime.strftime(format_pattern)
        year = '20'+formatted_string[0:2]
        month = formatted_string[2:4]
        day = formatted_string[0:6]
        day_url = ISWA_DATA_URL.rstrip('/')+'/'+\
            str(year)+str(month)+'/mrzqs'+str(day)+'/'
        [cr, text, link] = get_highest(
            day_url,r'(\d\d\d\d\d\dt\d\d\d\d)',matches.group(1))
    granule_url = day_url.rstrip('/') + '/' + link
    # Adjust input files
    fits_file = os.path.join(OUTPUT_BASE_PATH, "fitsfile.fits")
    fits_file_gz = os.path.join(OUTPUT_BASE_PATH,str(link))
    download_file(granule_url,fits_file_gz)
    with gzip.open(fits_file_gz, 'rb') as f_in:
        with open(fits_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.unlink(fits_file_gz)
