#!/usr/bin/env python
import pandas as pd
import json
import argparse

# Set up command-line argument parsing


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert CSV to JSON with table column headers and data.",
        epilog="""Example:
  python csv_to_json.py input.csv output.json
  python csv_to_json.py input.tsv output.json --sep '\\t'
        """,
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("input_csv", help="Path to the input CSV file")
    parser.add_argument("output_json", help="Path to the output JSON file")
    parser.add_argument("--sep", default=",", help="CSV separator (default: ',')")
    parser.add_argument("--name", default="data", help="data arrayname (default: 'data')")
    parser.add_argument("--colname", default="tablecol", help="table column array name (default: 'tablecol')")

    return parser.parse_args()
# Main function to process the data


def process_csv_to_json(input_csv, output_json):
    # Load the CSV file
    df = pd.read_csv(input_csv, sep=args.sep, encoding='utf-8')

    # Generate `tablecol`
    tablecol = [
        {"title": col.replace("_", " ").title(), "data": col}
        for col in df.columns
    ]

    # Generate `data`
    data = df.to_dict(orient="records")

    # Prepare the final structure in a Python-friendly format
    tablecol_str = f"{args.colname} = {json.dumps(tablecol, indent=2)}"
    data_str = f"{args.name} = {json.dumps(data, indent=2)}"

    # Save to a .json file
    with open(output_json, "w") as f:
        f.write(f"{tablecol_str}\n\n{data_str}\n")

    print(f"Data has been successfully written to {output_json}")
    print(tablecol_str)
    print(data_str)


# Entry point of the script
if __name__ == "__main__":
    args = parse_args()
    process_csv_to_json(args.input_csv, args.output_json)

