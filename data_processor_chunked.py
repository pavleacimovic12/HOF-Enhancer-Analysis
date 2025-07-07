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
        self.data_dir = "data_chunks"
        
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
        
        # Validate data integrity
        self.validate_data_integrity(peak_data, enhancer_metadata)
        
        return peak_data, hof_enhancers, enhancer_metadata
    
    def load_peak_data(self):
        """Load and combine chunked CSV files"""
        try:
            # Pattern to match all chunk files
            chunk_patterns = [
                "part1*chunk*.csv",
                "part2*chunk*.csv", 
                "part3*chunk*.csv",
                "part4*chunk*.csv"
            ]
            
            all_chunks = []
            for pattern in chunk_patterns:
                chunk_files = glob.glob(os.path.join(self.data_dir, pattern))
                chunk_files.sort()  # Ensure proper order
                
                if chunk_files:
                    print(f"Loading {len(chunk_files)} chunks for pattern {pattern}")
                    
                    for chunk_file in chunk_files:
                        df = pd.read_csv(chunk_file)
                        all_chunks.append(df)
                        print(f"  Loaded {os.path.basename(chunk_file)}: {len(df):,} rows")
            
            if not all_chunks:
                print("No chunk files found")
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
            metadata_path = os.path.join(self.data_dir, "Enhancer_and_experiment_metadata_1751579195077.feather")
            if os.path.exists(metadata_path):
                metadata = pd.read_feather(metadata_path)
                print(f"Loaded metadata: {len(metadata)} records")
                return metadata
            else:
                print(f"Metadata file not found at {metadata_path}")
                return None
        except Exception as e:
            print(f"Error loading metadata: {str(e)}")
            return None
    
    def extract_hof_enhancers(self, metadata_df, peak_data):
        """Extract the Hall of Fame enhancers and merge with metadata"""
        if peak_data is None or peak_data.empty:
            return pd.DataFrame()
        
        # Get unique enhancer IDs from peak data
        hof_enhancer_ids = peak_data['enhancer_id'].unique()
        print(f"Found {len(hof_enhancer_ids)} unique Hall of Fame enhancers")
        
        # Create base enhancer information from peak data
        hof_enhancers = []
        for enhancer_id in hof_enhancer_ids:
            enhancer_data = peak_data[peak_data['enhancer_id'] == enhancer_id]
            if not enhancer_data.empty:
                first_row = enhancer_data.iloc[0]
                hof_enhancers.append({
                    'enhancer_id': enhancer_id,
                    'chr': first_row.get('chr', ''),
                    'start': first_row.get('start', 0),
                    'end': first_row.get('end', 0)
                })
        
        hof_df = pd.DataFrame(hof_enhancers)
        
        # Add metadata if available
        if metadata_df is not None and not metadata_df.empty:
            # Get one representative row per enhancer from metadata
            metadata_summary = metadata_df.groupby('enhancer_id').first().reset_index()
            hof_df = hof_df.merge(metadata_summary[['enhancer_id', 'cargo', 'experiment', 'proximal_gene']], 
                                 on='enhancer_id', how='left')
        
        print(f"Extracted {len(hof_df)} Hall of Fame enhancers with metadata")
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