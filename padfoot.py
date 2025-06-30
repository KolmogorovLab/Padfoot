#!/usr/bin/env python3

#(c) 2023 by Authors
#This file is a part of Severus program.
#Released under the BSD license (see LICENSE file)

"""
This script sets up environment paths
and invokes padfoot without installation.
"""

import os
import sys

#BIN_DIR = "bin"

def main():
    #Setting executable paths
    padfoot_root = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, padfoot_root)
    #bin_absolute = os.path.join(severus_root, BIN_DIR)
    #os.environ["PATH"] = bin_absolute + os.pathsep + os.environ["PATH"]

    #Padffot entry point
    from padfoot.main import main
    sys.exit(main())


if __name__ == "__main__":
    main()
