
import argparse
import re

txt = "The rain in Spain"
x = re.search("^The.*Spain$", txt)
 
# Initialize parser
parser = argparse.ArgumentParser()

parser.add_argument("-f", "--file" )
parser.add_argument("-o", "--output" )

args = parser.parse_args()

file_input = args.file

with open(file_input, "r") as f:
    data = f.read()
    print(data)
    
process_times = re.findall("\d+.\d+", data)

median = sum([float(process_time) for process_time in process_times]) / len(process_times)


print(median)