import subprocess
import sys
import os


name = sys.argv[1]

files_to_delete = [name + ext for ext in ['.aux', '.bbl', '.log', '.blg', '.out', '.toc']]



for file in files_to_delete:
    try:
        os.remove(file)
    except FileNotFoundError:
        print(f"{file} not found")
        continue