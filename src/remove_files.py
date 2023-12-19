#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 15:04:09 2023

@author: Sam Becker, Vincent van Rie, Divyansh Srivastava 
"""

import os

for istep in range(0, 1000, 1):
    filename = f"snewstellargravhydroTEST_{istep}.hdf5"
    pngfilename = f"snewstellargravhydroTEST_{istep}.png"
    if os.path.exists(filename):
            os.remove(filename)
    if os.path.exists(pngfilename):
        os.remove(pngfilename)