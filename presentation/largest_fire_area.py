import csv

# Path to the CSV file
csv_file = '/home/goldma34/fire_insect_co-occurence/data/outputs/on/on_co-occurrences_1986-2012.csv'

# Initialize variables to track the maximum fire area and corresponding Fire_ID
max_area = 0
max_fire_id = None

# Open and read the CSV file
with open(csv_file, mode='r') as file:
    csv_reader = csv.DictReader(file)
    for row in csv_reader:
        fire_area = float(row['Fire_Area'])
        if fire_area > max_area:
            max_area = fire_area
            max_fire_id = row['Fire_ID']

# Print the result
print(f"The Fire_ID with the largest Fire_Area is {max_fire_id} with an area of {max_area}")