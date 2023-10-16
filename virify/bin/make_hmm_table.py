#!/usr/bin/env python3

import pandas as pd
import argparse
import sys
import os

headers = [
    "target name",
    "target accession",
    "tlen",
    "query name",
    "query accession",
    "qlen",
    "full sequence E-value",
    "full sequence score",
    "full sequence bias",
    "#",
    "of",
    "c-Evalue",
    "i-Evalue",
    "domain score",
    "domain bias",
    "hmm coord from",
    "hmm coord to",
    "ali coord from",
    "ali coord to",
    "env coord from",
    "env coord to",
    "acc",
    "description of target"
]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Format hmmscan domain hits table.")
    parser.add_argument("-t", dest="tables", nargs="+",
                        help="hmmscan domain hits tables (multiple files)")
    parser.add_argument("-o", "--outname", dest="outfile_name", required=True,
                        help="Output table name")

    args = parser.parse_args()

    if args.tables is None or len(args.tables) == 0:
        raise Exception("No tables specified.")

    bucket = []
    # tables = args.tables.split(",")

    for table in args.tables:
        with open(table, "r") as file:
            for line in file:
                if line.startswith("#"):
                    continue
                cols = line.split()
                bucket.append(cols)
                # description = " ".join(cols[22:])
                # bucket.append(cols[:22] + [description])

    # df = pd.DataFrame(bucket, columns=headers)
    df = pd.DataFrame(bucket)
    df.to_csv(args.outfile_name + ".tsv", sep="\t", index=False)
