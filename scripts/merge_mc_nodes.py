#!/usr/bin/env python3

"""
Pretty much this:
https://github.com/cta-observatory/cta-lstchain/blob/main/lstchain/scripts/lstchain_merge_hdf5_files.py
but with slight improvements for my case:
    - train/test split
    - no metadata check. Clearly this is not optimal, but I am only using this for nodes of the same production. These should be mergable, but due to floating point issues are sometimes not
"""

import argparse
import os
from glob import glob
from sklearn.model_selection import train_test_split

from lstchain.io import auto_merge_h5files
from lstchain.io import get_dataset_keys
from lstchain.io.io import copy_h5_nodes, add_source_filenames, write_metadata, read_metadata
from tables import open_file
import tables
import logging

log = logging.getLogger(__name__)

HDF5_ZSTD_FILTERS = tables.Filters(
    complevel=5,  # enable compression, 5 is a good tradeoff between compression and speed
    complib='blosc:zstd',  # compression using blosc/zstd
    fletcher32=True,  # attach a checksum to each chunk for error correction
    bitshuffle=False,  # for BLOSC, shuffle bits for better compression
)

def auto_merge_h5files(
        file_list,
        output_filename,
        nodes_keys,
        merge_arrays=False,
        filters=HDF5_ZSTD_FILTERS,
):
    """
    Based on https://github.com/cta-observatory/cta-lstchain/blob/v0.9.13/lstchain/io/io.py#L295
    Removed metadata check and other not strictly necessary code.
    """

    file_list = list(file_list)
    keys = set(nodes_keys)

    with open_file(output_filename, 'w', filters=filters) as merge_file:
        with open_file(file_list[0]) as f1:
            copy_h5_nodes(f1, merge_file, nodes=keys)

        for filename in file_list[1:]:
            common_keys = keys.intersection(get_dataset_keys(filename))
            with open_file(filename) as file:
                for k in common_keys:
                    in_node = file.root[k]
                    out_node = merge_file.root[k]
                    try:
                        if isinstance(in_node, tables.table.Table) or merge_arrays:
                            # doing `.astype(out_node.dtype)` fixes an issue
                            # when dtypes do not exactly match but are convertible
                            # https://github.com/cta-observatory/cta-lstchain/issues/671
                            out_node.append(in_node.read().astype(out_node.dtype))
                    except:
                        log.error("Can't append node {} from file {}".format(k, filename))
                        raise
        add_source_filenames(merge_file, file_list)

    # merge global metadata and store source file names
    metadata0 = read_metadata(file_list[0])
    write_metadata(metadata0, output_filename)


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
    assert args.output_train or args.output_test

    if not args.output_test:
        print("Not creating a test file since none is given")
        train_files = file_list
        test_files = None
    elif not args.output_train:
        print("Not creating a train file since none is given")
        train_files = None
        test_files = file_list
    else:
        print("Creating both train and test output")
        assert train_size < 1.0
        train_files, test_files = train_test_split(file_list, train_size=train_size)


    # We never want images
    keys = get_dataset_keys(file_list[0])
    keys = {k for k in keys if 'image' not in k}


    assert train_files or test_files
    if train_files:
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
