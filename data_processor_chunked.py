"""
Data Processor for GitHub Deployment with Chunked CSV Files
Handles loading and combining chunked CSV files automatically
"""

import pandas as pd
import os
import glob
from pathlib import Path

class DataProcessor:
    """Process and load genomic data from chunked CSV files and metadata"""
    
    def __init__(self):
        # Get absolute path to current directory for Posit Cloud compatibility
        self.base_path = os.path.dirname(os.path.abspath(__file__))
        self.data_dir = self.base_path  # Use absolute path for all operations
        print(f"DataProcessor initialized with base_path: {self.base_path}")
        print(f"Available files: {[f for f in os.listdir(self.base_path) if f.endswith(('.csv', '.feather'))][:5]}...")
        
    def load_all_data(self):
        """Load and process all data files including chunked CSV files"""
        # Load peak data from chunked files
        peak_data = self.load_peak_data()
        
        if peak_data is None or peak_data.empty:
            return None, None, None
        
        # Load enhancer metadata
        enhancer_metadata = self.load_metadata()
        
        # Extract Hall of Fame enhancers
        hof_enhancers = self.extract_hof_enhancers(enhancer_metadata, peak_data)
        
        print(f"=== FINAL VALIDATION ===")
        print(f"Peak data: {len(peak_data) if peak_data is not None else 'None'}")
        print(f"Metadata: {len(enhancer_metadata) if enhancer_metadata is not None else 'None'}")
        print(f"HOF enhancers: {len(hof_enhancers) if hof_enhancers is not None else 'None'}")
        
        # Validate data integrity
        self.validate_data_integrity(peak_data, enhancer_metadata)
        
        # CRITICAL: Return order must match app.py expectations
        return peak_data, enhancer_metadata, hof_enhancers
    
    def load_peak_data(self):
        """Load and combine chunked CSV files"""
        print(f"Loading peak data from directory: {self.data_dir}")
        try:
            # Pattern to match all chunk files in current directory
            chunk_patterns = [
                "part1*chunk*.csv",
                "part2*chunk*.csv", 
                "part3*chunk*.csv",
                "part4*chunk*.csv"
            ]
            
            all_chunks = []
            for pattern in chunk_patterns:
                search_pattern = os.path.join(self.data_dir, pattern)
                chunk_files = glob.glob(search_pattern)
                chunk_files.sort()  # Ensure proper order
                print(f"Pattern {pattern}: found {len(chunk_files)} files")
                
                if chunk_files:
                    print(f"Loading {len(chunk_files)} chunks for pattern {pattern}")
                    
                    for chunk_file in chunk_files:
                        df = pd.read_csv(chunk_file)
                        all_chunks.append(df)
                        print(f"  Loaded {os.path.basename(chunk_file)}: {len(df):,} rows")
            
            if not all_chunks:
                print("No chunk files found in current directory")
                print("Expected files: part1*chunk*.csv, part2*chunk*.csv, part3*chunk*.csv, part4*chunk*.csv")
                print("Current directory contents:")
                for file in os.listdir("."):
                    if file.endswith(".csv"):
                        print(f"  {file}")
                return None
            
            # Combine all chunks
            combined_df = pd.concat(all_chunks, ignore_index=True)
            print(f"Combined all chunks: {len(combined_df):,} total rows")
            
            # Remove any duplicate rows
            original_length = len(combined_df)
            combined_df = combined_df.drop_duplicates()
            if len(combined_df) < original_length:
                print(f"Removed {original_length - len(combined_df):,} duplicate rows")
            
            return combined_df
            
        except Exception as e:
            print(f"Error loading peak data: {str(e)}")
            return None
    
    def load_metadata(self):
        """Load enhancer metadata from feather file"""
        try:
            # Use absolute path for Posit Cloud compatibility
            metadata_path = os.path.join(self.base_path, "Enhancer_and_experiment_metadata_1751579195077.feather")
            print(f"Looking for metadata at: {metadata_path}")
            print(f"File exists: {os.path.exists(metadata_path)}")
            if os.path.exists(metadata_path):
                metadata = pd.read_feather(metadata_path)
                print(f"Loaded metadata: {len(metadata)} records")
                
                # Fix column names to match expected format
                column_mapping = {
                    'Enhancer_ID': 'enhancer_id',
                    'Cargo': 'cargo',
                    'Experiment_Type': 'experiment',
                    'Proximal_Gene': 'proximal_gene',
                    'Image_link': 'image_link',
                    'Neuroglancer 1': 'neuroglancer_1',
                    'Neuroglancer 3': 'neuroglancer_3',
                    'Viewer Link': 'viewer_link',
                    'Coronal_MIP': 'coronal_mip',
                    'Sagittal_MIP': 'sagittal_mip'
                }
                
                # Rename columns that exist
                existing_columns = {old_col: new_col for old_col, new_col in column_mapping.items() if old_col in metadata.columns}
                metadata = metadata.rename(columns=existing_columns)
                
                print(f"Renamed columns: {existing_columns}")
                print(f"Enhanced metadata columns: {list(metadata.columns)}")
                
                return metadata
            else:
                print(f"Metadata file not found at {metadata_path}")
                print("Expected file: Enhancer_and_experiment_metadata_1751579195077.feather")
                print("Current directory contents:")
                for file in os.listdir("."):
                    if file.endswith(".feather"):
                        print(f"  {file}")
                return None
        except Exception as e:
            print(f"Error loading metadata: {str(e)}")
            return None
    
    def extract_hof_enhancers(self, metadata_df, peak_data):
        """Extract the Hall of Fame enhancers and merge with metadata"""
        if peak_data is None or peak_data.empty:
            return pd.DataFrame()
        
        # CRITICAL FIX: Get Hall of Fame enhancers from metadata first
        if metadata_df is None or metadata_df.empty:
            print("No metadata available - cannot identify Hall of Fame enhancers")
            return pd.DataFrame()
        
        # Filter metadata for Hall of Fame enhancers (string 'TRUE', not boolean True)
        hof_metadata = metadata_df[metadata_df['Hall_of_fame'] == 'TRUE']
        print(f"Found {len(hof_metadata)} Hall of Fame enhancer records in metadata")
        
        if hof_metadata.empty:
            print("No Hall of Fame enhancers found in metadata")
            return pd.DataFrame()
        
        # Get unique HOF enhancer IDs from metadata
        hof_enhancer_ids = hof_metadata['enhancer_id'].unique()
        print(f"Found {len(hof_enhancer_ids)} unique Hall of Fame enhancers: {hof_enhancer_ids[:5]}...")
        
        # Create base enhancer information more efficiently
        print(f"Creating base enhancer records for {len(hof_enhancer_ids)} enhancers...")
        
        # Get first occurrence of each HOF enhancer from peak data
        hof_peak_data = peak_data[peak_data['enhancer_id'].isin(hof_enhancer_ids)]
        first_occurrences = hof_peak_data.groupby('enhancer_id').first().reset_index()
        
        # Create base enhancer records
        hof_enhancers = []
        for _, row in first_occurrences.iterrows():
            enhancer_record = {
                'enhancer_id': row['enhancer_id'],
                'chr': row.get('chr', ''),
                'start': row.get('start', 0),
                'end': row.get('end', 0)
            }
            hof_enhancers.append(enhancer_record)
        
        hof_df = pd.DataFrame(hof_enhancers)
        print(f"Created {len(hof_df)} base enhancer records")
        
        # Efficiently merge metadata using pandas merge
        if metadata_df is not None and not metadata_df.empty and not hof_df.empty:
            print(f"Merging metadata for {len(hof_df)} enhancers...")
            
            # Create metadata summary for merge
            metadata_for_merge = metadata_df.groupby('enhancer_id').first().reset_index()
            metadata_cols = ['enhancer_id', 'cargo', 'experiment', 'proximal_gene', 'GC delivered']
            available_cols = [col for col in metadata_cols if col in metadata_for_merge.columns]
            
            # Merge using pandas - much more efficient than individual lookups
            hof_df = hof_df.merge(metadata_for_merge[available_cols], on='enhancer_id', how='left')
            
            print(f"Merged {len(hof_df)} enhancers with metadata")
            
            # Quick validation
            if 'cargo' in hof_df.columns:
                cargo_success = hof_df['cargo'].notna().sum()
                print(f"Cargo merge success: {cargo_success}/{len(hof_df)} enhancers")
        
        print(f"Successfully extracted {len(hof_df)} Hall of Fame enhancers with integrated metadata")
        return hof_df
    
    def get_enhancer_summary(self, peak_data):
        """Generate comprehensive summary statistics for enhancers"""
        if peak_data is None or peak_data.empty:
            return {}
        
        summary = {
            'total_enhancers': peak_data['enhancer_id'].nunique(),
            'total_measurements': len(peak_data),
            'cell_types': peak_data['cell_type'].nunique() if 'cell_type' in peak_data.columns else 0,
            'chromosomes': peak_data['chr'].nunique() if 'chr' in peak_data.columns else 0,
            'max_accessibility': peak_data['accessibility'].max() if 'accessibility' in peak_data.columns else 0,
            'mean_accessibility': peak_data['accessibility'].mean() if 'accessibility' in peak_data.columns else 0
        }
        
        return summary
    
    def validate_data_integrity(self, peak_data, metadata):
        """Validate data integrity and consistency"""
        print("\n=== Data Validation ===")
        
        if peak_data is not None and not peak_data.empty:
            print(f"✓ Peak data loaded: {len(peak_data):,} rows")
            print(f"✓ Unique enhancers: {peak_data['enhancer_id'].nunique()}")
            if 'cell_type' in peak_data.columns:
                print(f"✓ Cell types: {peak_data['cell_type'].nunique()}")
        else:
            print("✗ Peak data missing or empty")
        
        if metadata is not None and not metadata.empty:
            print(f"✓ Metadata loaded: {len(metadata):,} records")
            print(f"✓ Unique enhancers in metadata: {metadata['enhancer_id'].nunique()}")
        else:
            print("✗ Metadata missing or empty")
        
        # Check for common enhancers
        if peak_data is not None and metadata is not None:
            peak_enhancers = set(peak_data['enhancer_id'].unique())
            meta_enhancers = set(metadata['enhancer_id'].unique())
            common_enhancers = peak_enhancers.intersection(meta_enhancers)
            print(f"✓ Common enhancers between datasets: {len(common_enhancers)}")
        
        print("======================\n")