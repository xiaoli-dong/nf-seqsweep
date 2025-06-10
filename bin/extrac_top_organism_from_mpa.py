#!/usr/bin/env python
import csv
import sys

if len(sys.argv) != 4:
    print("Usage: python top_species.py <file.tsv> <output.csv> <sample name>")
    sys.exit(1)

tsv_file = sys.argv[1]
output_csv = sys.argv[2]
sname = sys.argv[3]

species_data = {}
#
with open(tsv_file, "r", encoding="utf-8") as file:
    reader = csv.reader(file, delimiter="\t")
    header = next(reader)
    """
        #Classification sample1_count   sample1_perc
        x__cellular_organisms   425705  99.29
        x__cellular_organisms|x__Bacteria       425696  99.29
    """
    for row in reader:
        if len(row) < 3:
            continue
        taxonomy = row[0]
        try:
            count = float(row[1])           # sample1_count
            percentage = float(row[2])      # sample1_perc
        except ValueError:
            continue
        if "s__" not in taxonomy:
            continue
        for level in taxonomy.split("|"):
            if level.startswith("s__"):
                species = level.replace("s__", "").strip()
                if species in species_data:
                    species_data[species]["count"] += count
                    species_data[species]["perc"] += percentage
                else:
                    species_data[species] = {
                        "count": count, "perc": percentage}
                break

# Sort by count and take top 4
top_species = sorted(species_data.items(),
                     key=lambda x: x[1]["count"], reverse=True)[:4]
header = ["Sample", "Match1", "count1", "perc1", "Match2",
          "count2", "perc2", "Match3", "count3", "perc3", "Match4", "count4", "perc4"]
# Write to CSV as single line
with open(output_csv, "w", newline="", encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)
    flat_row = [sname]
    for species, data in top_species:
        flat_row.extend([species, int(data["count"]), round(data["perc"], 2)])
    writer.writerow(flat_row)

# print(f"Top 4 species with percentage written to: {output_csv}")
