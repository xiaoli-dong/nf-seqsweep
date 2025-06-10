#!/usr/bin/env python
import argparse
import pandas as pd

def main():
    description = "fastq quality control summary report"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "--input-stats",
        default=None,
        help="Path to input read stats file\n"
    )
    parser.add_argument(
        "--trim-stats",
        default=None,
        help="Path to trimmed read stats file\n"
    )
    parser.add_argument(
        "--dehost-stats",
        default=None,
        help="Path to dehosted read stats file\n"
    )
    parser.add_argument(
        "--output",
        default=None,
        help="merged output file in csv format\n"
    )
    args = parser.parse_args()

    # Read and process the first file
    merged_df = pd.read_csv(args.input_stats, sep='\t')
    # Drop 2nd and 3rd columns
    merged_df = merged_df.drop(columns=merged_df.columns[[1, 2]])

    #Add prefix to all columns except the first
    prefix = 'input_'
    merged_df.columns = [merged_df.columns[0]] + [prefix + col for col in merged_df.columns[1:]]

    if args.trim_stats:
        df = pd.read_csv(args.trim_stats, sep='\t')
        df = df.drop(columns=df.columns[[1, 2]])  # Drop 2nd and 3rd columns
        # Add prefix to all columns except the first
        prefix = 'trim_'
        df.columns = [df.columns[0]] + [prefix + col for col in df.columns[1:]]

        """ merged_df = pd.merge(merged_df, df, on='file',
                             suffixes=('_input', '_trim'), how='outer') """
        merged_df = pd.merge(merged_df, df, on='file', how='outer')

        # Add percentage columns
        merged_df['trim_num_seqs_passed_pct'] = (
            merged_df['trim_num_seqs'] / merged_df['input_num_seqs'] * 100).round(2)
        merged_df['trim_sum_len_passed_pct'] = (
            merged_df['trim_sum_len'] / merged_df['input_sum_len'] * 100).round(2)
        print(merged_df)

    if args.dehost_stats:
        df = pd.read_csv(args.dehost_stats, sep='\t')
        df = df.drop(columns=df.columns[[1, 2]])  # Drop 2nd and 3rd columns
        # Add prefix to all columns except the first
        prefix = 'dehost_'
        df.columns = [df.columns[0]] + [prefix + col for col in df.columns[1:]]
        merged_df = pd.merge(merged_df, df, on='file', how='outer')
        merged_df['dehost_num_seqs_passed_pct'] = (
            merged_df['dehost_num_seqs'] / merged_df['input_num_seqs'] * 100).round(2)
        merged_df['dehost_sum_len_passed_pct'] = (
            merged_df['dehost_sum_len'] / merged_df['input_sum_len'] * 100).round(2)
        print(merged_df)
    # Save the result
    merged_df.columns.values[0] = 'sample'
    merged_df.to_csv(args.output, sep=',', index=False)



if __name__ == "__main__":
    main()
