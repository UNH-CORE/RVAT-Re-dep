#!/usr/bin/env python
"""
This script uploads raw data to figshare.
"""
from __future__ import division, print_function
import requests
from requests_oauthlib import OAuth1
import os
import json
#import urllib

article = 1256369

def load_credentials():
    with open(os.path.join(os.path.expanduser("~"), ".figsharerc")) as f:
        cred = json.load(f)
    return cred
    
def authenticate():
    cred = load_credentials()
    client_key = cred["client_key"]
    client_secret = cred["client_secret"]
    token_key = cred["token_key"]
    token_secret = cred["token_secret"]
    oauth = OAuth1(client_key=client_key, client_secret=client_secret,
            resource_owner_key=token_key, resource_owner_secret=token_secret,
            signature_type='auth_header')
    return oauth
    
def get_article_details():
    client = requests.session()
    url = "http://api.figshare.com/v1/my_data/articles/{}".format(article)
    oauth = authenticate()
    response = client.get(url, auth=oauth)
    return json.loads(response.content.decode())
    
def get_uploaded_files():
    return get_article_details()["items"][0]["files"]
    
def make_local_file_list():
    file_list = []
    raw_dir = os.path.join("Data", "Raw")
    sections = os.listdir(raw_dir)
    for section in sections:
        section_dir = os.path.join(raw_dir, section)
        runs = os.listdir(section_dir)
        for run in runs:
            run_dir = os.path.join(section_dir, run)
            files = os.listdir(run_dir)
            for f in files:
                fpath = os.path.join(run_dir, f)
                file_list.append(fpath)
    return file_list

def test_upload_file():
    fpath = make_local_file_list()[1]
    fname = "_".join(fpath.split("\\")[-3:])
    client = requests.session()
    files = {"filedata" : (fpath, open(fpath, "rb"))}
    url = "http://api.figshare.com/v1/my_data/articles/{}/files".format(article)
    oauth = authenticate()
    response = client.put(url, auth=oauth, files=files)
    results = json.loads(response.content.decode())
    print(results)

if __name__ == "__main__":
    if os.getcwd()[-7:] == "Scripts":
        print("Changing working directory to experiment root directory")
        os.chdir("../")

    test_upload_file()