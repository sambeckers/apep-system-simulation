#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
remove_files
Created on 08-12-23

@author: Sam Beckers, Vincent van Rie, Divyansh Srivastava

Short script to remove all the files created by the simulation.
"""
import os


def Remove_files(filename):
    for istep in range(0, 1000, 1):
        hdf5filename = f"{filename}{istep}.hdf5"
        pngfilename = f"{filename}{istep}.png"
        if os.path.exists(hdf5filename):
            os.remove(hdf5filename)
        if os.path.exists(pngfilename):
            os.remove(pngfilename)


if __name__ == "__main__":
    # Remove_files("apep_rot_")
    # Remove_files("stellargravhydro_supernova_")
    print("Files removed")
