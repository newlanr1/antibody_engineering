# src/generate_combinations.py
import csv
import argparse
import itertools

def read_csv_columns(input_file):
    """Reads a two-column CSV and returns two lists, skipping the header."""
    col1, col2 = [], []
    with open(input_file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)  # Skip header
        for row in reader:
            if len(row) > 0 and row[0].strip():
                col1.append(row[0].strip())
            if len(row) > 1 and row[1].strip():
                col2.append(row[1].strip())
    return col1, col2

def generate_combinations(list1, list2):
    """Generates all possible combinations (Cartesian product) of two lists."""
    list1 = list1 if list1 else ['']
    list2 = list2 if list2 else ['']
    return list(itertools.product(list1, list2))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate combinations from a 2-column CSV.")
    parser.add_argument("input_csv", help="Path to input CSV")
    parser.add_argument("output_csv", help="Path to output CSV")
    args = parser.parse_args()

    col1, col2 = read_csv_columns(args.input_csv)
    combinations = generate_combinations(col1, col2)
    
    with open(args.output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Sequence1', 'Sequence2'])
        writer.writerows(combinations)
        
    print(f"Generated {len(combinations)} combinations and saved to {args.output_csv}")