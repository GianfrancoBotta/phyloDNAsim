import gzip
import pandas as pd
import os

def count_reads(fastq_gz_path, gzipped=True):
    """Count the number of reads in a fastq.gz file."""    
    count = 0
    if gzipped:
        with gzip.open(fastq_gz_path, 'rt') as f:
            for i, _ in enumerate(f):
                if i % 4 == 0:  # Every 4th line starting from 0 is a read header
                    count += 1
    else:
        with open(fastq_gz_path, 'rt') as f:
            for i, _ in enumerate(f):
                if i % 4 == 0:  # Every 4th line starting from 0 is a read header
                    count += 1
    
    return count

def get_time(filepath: str) -> float:
    with open(filepath, 'r') as f:
        next(f)  # Skip header
        first_row = next(f)
        return float(first_row.split('\t')[0])
        
def get_folder_size(path):
    total = 0
    for dirpath, dirnames, filenames in os.walk(path):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            if not os.path.islink(filepath):  # skip symlinks
                total += os.path.getsize(filepath)
    return total

def format_size(size_bytes):
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024