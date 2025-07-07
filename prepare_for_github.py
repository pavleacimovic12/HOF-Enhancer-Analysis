#!/usr/bin/env python3
"""
Prepare files for GitHub deployment by moving CSV chunks to main directory
This script moves all CSV chunks from data_chunks/ to the main directory
and the metadata feather file, making them ready for GitHub upload.
"""

import os
import shutil
import glob

def prepare_for_github():
    """Move all data files to main directory for GitHub deployment"""
    print("Preparing files for GitHub deployment...")
    
    # Source directory
    data_chunks_dir = "data_chunks"
    
    if not os.path.exists(data_chunks_dir):
        print(f"Error: {data_chunks_dir} directory not found")
        return False
    
    # Move all CSV files
    csv_files = glob.glob(os.path.join(data_chunks_dir, "*.csv"))
    print(f"Found {len(csv_files)} CSV files to move")
    
    for csv_file in csv_files:
        filename = os.path.basename(csv_file)
        destination = filename
        
        # Check if file already exists in main directory
        if os.path.exists(destination):
            print(f"  {filename} already exists in main directory, skipping")
            continue
            
        # Move file
        shutil.move(csv_file, destination)
        print(f"  Moved {filename}")
    
    # Move metadata feather file
    feather_files = glob.glob(os.path.join(data_chunks_dir, "*.feather"))
    print(f"Found {len(feather_files)} feather files to move")
    
    for feather_file in feather_files:
        filename = os.path.basename(feather_file)
        destination = filename
        
        # Check if file already exists
        if os.path.exists(destination):
            print(f"  {filename} already exists in main directory, skipping")
            continue
            
        # Move file
        shutil.move(feather_file, destination)
        print(f"  Moved {filename}")
    
    print("\nFiles prepared for GitHub deployment!")
    print("You can now upload all files to GitHub repository")
    
    # List files in main directory for verification
    print("\nFiles in main directory:")
    main_files = sorted([f for f in os.listdir(".") if f.endswith((".csv", ".feather"))])
    for file in main_files:
        size_mb = os.path.getsize(file) / (1024 * 1024)
        print(f"  {file} ({size_mb:.1f} MB)")
    
    return True

if __name__ == "__main__":
    success = prepare_for_github()
    if success:
        print("\n✅ Ready for GitHub upload!")
    else:
        print("\n❌ Preparation failed!")