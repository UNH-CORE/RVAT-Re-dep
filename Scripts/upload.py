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

if __name__ == "__main__":
    if os.getcwd()[-7:] == "Scripts":
        print("Changing working directory to experiment root directory")
        os.chdir("../")

    print(get_uploaded_files())