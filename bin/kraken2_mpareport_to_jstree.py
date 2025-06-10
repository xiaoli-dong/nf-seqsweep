#!/usr/bin/env python
import json
import sys

def insert_path(tree, path, values, sample_keys):
    node = tree
    for i, level in enumerate(path):
        children = node.setdefault("children", [])
        found = next((child for child in children if child["text"] == level), None)
        if not found:
            found = {
                "text": level,
                "state": {"opened": i <= 2},  # Open up to level 3 (0-based)
                "children": []
            }
            children.append(found)
        node = found

    node["data"] = dict(zip(sample_keys, map(float, values)))

def parse_file_to_tree(filepath):
    with open(filepath, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]

    header = lines[0].split("\t")
    full_sample_names = header[1:]

    # Build original sample key mapping
    original_sample_keys = [name.split('.')[0] for name in full_sample_names]

    # Sort sample keys alphabetically
    sorted_sample_keys = sorted(original_sample_keys)

    # Build index mapping from sorted to original order
    sort_indices = [original_sample_keys.index(key) for key in sorted_sample_keys]

    root = {"text": "Taxonomy", "state": {"opened": True}, "children": []}

    # Replacement mapping
    replacements = {
        "x__cellular_organisms": "R1__cellular_organisms",
        "x__Eukaryota": "R2__Eukaryota",
        "x__Bacteria": "R2__Bacteria",
        "x__Viruses": "R1__Viruses"
    }

    for line in lines[1:]:
        parts = line.split("\t")
        path = parts[0].split("|")
        updated_path = [replacements.get(item, item) for item in path]

        values = parts[1:]
        # Reorder values to match sorted sample keys
        sorted_values = [values[i] for i in sort_indices]

        insert_path(root, updated_path, sorted_values, sorted_sample_keys)

    # Generate JS output
    tree_js = f"var mydata = {json.dumps(root, indent=2)};"

    treecol = [{"header": "Taxonomy"}] + [
        {"value": key, "header": key, "width": "100px"} for key in sorted_sample_keys
    ]
    treecol_js = f"var treecol = {json.dumps(treecol, indent=2)};"

    return tree_js + "\n\n" + treecol_js

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python kraken2_mpareport_to_jstree.py <input_file.tsv> <output_file.json>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_content = parse_file_to_tree(input_file)

    output_file = sys.argv[2]
    with open(output_file, "w", encoding="utf-8") as out:
        out.write(output_content)

    print(f"jstree json file written to {output_file}")
