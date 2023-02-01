#!/usr/bin/env python
# coding: utf-8

# Copyright 2022
# Maximiliano Isi <max.isi@ligo.org>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.

import argparse
import requests
import os
from tqdm import tqdm
import json
import logging
import wget
import paths

DEFURL = "https://www.gw-openscience.org/eventapi/jsonfull/GWTC"
DEFDIR = paths.pe
parser = argparse.ArgumentParser(description="Download PE data using GWOSC.")
parser.add_argument('-o', '--outdir', default=DEFDIR,
                    help=f"Path target directory (def. {DEFDIR}).")
parser.add_argument('--url', help=f"GWOSC URL (def. {DEFURL}", default=DEFURL)
parser.add_argument("--far-max", help="FAR threshold (def. 1)", default=1.,
                    type=float)
parser.add_argument("--m2-min", help="Minimum mass (def. 2)", default=2.,
                    type=float)
parser.add_argument('--no-download', action='store_true')
parser.add_argument('--force-download', action='store_true')
parser.add_argument('-v', '--verbose', action='store_true')

args = parser.parse_args()

outdir = os.path.abspath(args.outdir)
os.makedirs(outdir, exist_ok=True)

if args.verbose:
    logging.getLogger().setLevel(logging.INFO)

# get data from the current GWTC events, as compiled by GWOSC
logging.info("Pulling events from {}".format(args.url))
r = requests.get(args.url)
rj = r.json()

# let's keep only highly significant BBHs
def is_good(k,v, far_max=args.far_max, m2_min=args.m2_min):
    v_check = (v['far'] <= far_max and v['mass_2_source'] and v['mass_2_source'] > m2_min)
    k_check = 'GW19' in k or 'GW20' in k
    return v_check and k_check

logging.info("Collected {} events from GWOSC.".format(len(rj['events'])))

good_events = {k for k, v in rj['events'].items() if is_good(k, v)}
n = len(good_events)
logging.info(f"Keeping {n} events satisfying selection criteria.")

events = {}
for k in tqdm(good_events):
    v = rj['events'][k]
    events[v['commonName']] = {'id': k}
    # pull additional information for this event
    r_event = requests.get(v['jsonurl'])
    rj_event = r_event.json()
    # locate most recent PE file info
    param_dict = rj_event['events'][k]['parameters']
    pe_dates = []
    for k_pe, v_pe in param_dict.items():
        # ignore parameters from pipelines
        if '_pe_' in k_pe:
            date = v_pe['date_added']
            pe_dates.append((date, k_pe))
    if not pe_dates:
        logging.warning(f"did not find PE for {k}")
        break
    chosen_pe_key = sorted(pe_dates)[-1][1]
    if not param_dict[chosen_pe_key]['is_preferred']:
        logging.warning(f"did not pick preferred PE for {k}")
        pref_key = [k for k, v in param_dict.items()
                    if '_pe_' in k and v['is_preferred']][0]
        logging.info(f"Chose {chosen_pe_key} over {pref_key}")
    # we can now get the URL of the PE samples from the chosen PE run
    data_url = param_dict[chosen_pe_key]['data_url']
    fname = os.path.basename(data_url)
    events[v['commonName']]['data_url'] = data_url
    events[v['commonName']]['file_name'] = fname
    outpath = os.path.join(outdir, fname)
    if not args.no_download:
        if args.force_download or not os.path.exists(outpath):
            logging.info(f"Downloading {k}")
            # response = requests.get(data_url, stream=True)
            # with open(outpath, "wb") as handle:
            #     for data in tqdm(response.iter_content()):
            #         handle.write(data)
            wget.download(data_url, outpath)
        elif os.path.exists(outpath):
            logging.info(f"Skipping preexisting file: {outpath}")
    
logging.info("Recording event data")
p = os.path.join(outdir, "gwosc_catalog.json")
with open(p, "w") as f:
    json.dump(rj['events'], f, indent=4)
logging.info("Saved: {}".format(p))

p = os.path.join(outdir, "event_names.json")
with open(p, "w") as f:
    json.dump(events, f, indent=4)
logging.info("Saved: {}".format(p))
