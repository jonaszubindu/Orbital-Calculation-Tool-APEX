#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 14:16:15 2020

@author: jonaszbinden
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-a", nargs = "+", help="append arguments")
args = parser.parse_args()
print(args.a)