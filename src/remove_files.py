#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
remove_files
Created on 08-12-23

@author: Sam Beckers, Vincent van Rie, Divyansh Srivastava

Short script to remove all the files created by the simulation.
"""
import os

for istep in range(0, 1000, 1):
    filename = f"snewstellargravhydroTEST_{istep}.hdf5"
    pngfilename = f"snewstellargravhydroTEST_{istep}.png"
    if os.path.exists(filename):
            os.remove(filename)
    if os.path.exists(pngfilename):
        os.remove(pngfilename)