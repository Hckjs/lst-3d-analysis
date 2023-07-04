#!/usr/bin/env python3

"""
Pretty much this:
https://github.com/cta-observatory/cta-lstchain/blob/main/lstchain/scripts/lstchain_merge_hdf5_files.py
but with slight improvements for my case.
"""

import argparse
import os
from glob import glob
import numpy as np
from sklearn.model_selection import train_test_split

from lstchain.io import auto_merge_h5files
from lstchain.io import get_dataset_keys

parser = argparse.ArgumentParser(description='Merge HDF5 files')

# Required arguments
parser.add_argument(
    '-d', '--input-dir',
    help='path to the source directory of files',
    required=True,
)

# Optional arguments
parser.add_argument(
    '--train-size',
    default=0.5,
    type=float,
)

parser.add_argument(
    '--output-train',
    required=True,
)
parser.add_argument(
    '--output-test',
    required=False,
)

parser.add_argument(
    '-p', '--pattern',
    default='*.h5',
    help='Glob pattern to match files',
)

def main():
    args = parser.parse_args()
    train_size = args.train_size
    file_list = glob(os.path.join(args.input_dir, args.pattern))
    print("Dir:", args.input_dir)
    print("pattern:", args.pattern)
    print("files:", file_list)
    if not args.output_test:
        print("Not creating a test file since none is given")
        train_files = file_list
        test_files = None
    else:
        assert train_size < 1.0
        train_files, test_files = train_test_split(file_list, train_size=train_size)


    # We never want images
    keys = get_dataset_keys(file_list[0])
    keys = {k for k in keys if 'image' not in k}

    auto_merge_h5files(
        train_files,
        args.output_train,
        nodes_keys=keys,
    )
    if test_files:
        auto_merge_h5files(
            test_files,
            args.output_test,
            nodes_keys=keys,
        )



if __name__ == '__main__':
    main()
