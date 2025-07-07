"""
Hall of Fame Enhancers Analysis - GitHub Version for Posit Cloud
Interactive genomic data visualization platform for analyzing enhancer accessibility profiles
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
# Import the chunked data processor for GitHub deployment
try:
    from data_processor_chunked import DataProcessor
except ImportError:
    from data_processor import DataProcessor
from visualization import VisualizationGenerator
import os

# Configure Streamlit page
st.set_page_config(
    page_title="Cheat Sheet for Enhancer Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for genomic theme
st.markdown("""
<style>
    .stApp {
        background-color: #f8f9fa;
    }
    .main-header {
        background: linear-gradient(90deg, #2E86AB, #A23B72);
        color: white;
        padding: 1rem;
        border-radius: 10px;
        margin-bottom: 2rem;
    }
    .filter-section {
        background: white;
        padding: 1rem;
        border-radius: 8px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        margin-bottom: 1rem;
    }
    .stSelectbox > div > div {
        background-color: white;
    }
    .enhancer-card {
        background: white;
        padding: 1.5rem;
        border-radius: 8px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        margin-bottom: 1rem;
    }
    .accessibility-plot {
        border: 2px solid #e1e5e9;
        border-radius: 8px;
        padding: 1rem;
        background: white;
    }
    .sidebar .sidebar-content {
        background: linear-gradient(180deg, #f8f9fa 0%, #e9ecef 100%);
    }
</style>
""", unsafe_allow_html=True)

class GitHubDataProcessor(DataProcessor):
    """Modified DataProcessor for GitHub deployment with different data paths"""
    
    def __init__(self):
        super().__init__()
        # Override data paths for GitHub structure
        self.data_dir = "data"
        
    def load_metadata(self):
        """Load enhancer metadata from feather file in data directory"""
        try:
            metadata_path = os.path.join(self.data_dir, "Enhancer_and_experiment_metadata_1751579195077.feather")
            if os.path.exists(metadata_path):
                metadata = pd.read_feather(metadata_path)
                return metadata
            else:
                st.error(f"Metadata file not found at {metadata_path}")
                return None
        except Exception as e:
            st.error(f"Error loading metadata: {str(e)}")
            return None
    
    def load_peak_data(self):
        """Load and process the peak accessibility data from complete high-resolution files"""
        try:
            # File paths for GitHub data structure
            file_paths = [
                os.path.join(self.data_dir, "part1 (1)_1751576434359.csv"),
                os.path.join(self.data_dir, "part2 (1)_1751576437893.csv"),
                os.path.join(self.data_dir, "part3_1751576441401.csv"),
                os.path.join(self.data_dir, "part4_1751576447956.csv")
            ]
            
            # Check if all files exist
            missing_files = [f for f in file_paths if not os.path.exists(f)]
            if missing_files:
                st.error(f"Missing peak data files: {missing_files}")
                return None
            
            # Load and combine all parts
            dfs = []
            for file_path in file_paths:
                if os.path.exists(file_path):
                    df = pd.read_csv(file_path)
                    dfs.append(df)
                    st.info(f"Loaded {len(df)} records from {os.path.basename(file_path)}")
            
            if not dfs:
                st.error("No peak data files found")
                return None
            
            # Combine all dataframes
            combined_df = pd.concat(dfs, ignore_index=True)
            st.success(f"Combined peak data: {len(combined_df)} total records")
            
            return combined_df
            
        except Exception as e:
            st.error(f"Error loading peak data: {str(e)}")
            return None

# Initialize processor and visualization
@st.cache_resource
def load_data():
    """Load and process all data files"""
    processor = GitHubDataProcessor()
    
    # Check if data directory exists
    if not os.path.exists("data_chunks"):
        st.error("Data directory not found. Please ensure the 'data_chunks' folder is present.")
        return None, None, None
    
    # Check if required files exist
    required_files = [
        "data/Enhancer_and_experiment_metadata_1751579195077.feather",
        "data/part1 (1)_1751576434359.csv",
        "data/part2 (1)_1751576437893.csv",
        "data/part3_1751576441401.csv",
        "data/part4_1751576447956.csv"
    ]
    
    missing_files = [f for f in required_files if not os.path.exists(f)]
    if missing_files:
        st.error(f"Missing data files: {', '.join(missing_files)}")
        st.info("Please ensure all data files are uploaded to the 'data' folder")
        return None, None, None
    
    try:
        peak_data, hof_enhancers, enhancer_metadata = processor.load_all_data()
        return peak_data, hof_enhancers, enhancer_metadata
    except Exception as e:
        st.error(f"Error loading data: {str(e)}")
        return None, None, None

def extract_cell_type_number(cell_type):
    """Extract numeric value from cell type for sorting"""
    if pd.isna(cell_type):
        return float('inf')
    
    cell_type_str = str(cell_type)
    for i, char in enumerate(cell_type_str):
        if char.isdigit():
            num_str = ""
            for j in range(i, len(cell_type_str)):
                if cell_type_str[j].isdigit():
                    num_str += cell_type_str[j]
                else:
                    break
            return int(num_str) if num_str else float('inf')
    return float('inf')

def get_filtered_options(selected_filters, base_metadata):
    """Get available options for each filter based on current selections - cell type remains independent"""
    if base_metadata is None or base_metadata.empty:
        return {
            'enhancers': [],
            'cargos': [],
            'experiments': [],
            'genes': [],
            'gc_delivered': [],
            'cell_types': []
        }
    
    def get_options_for_filter(exclude_filter):
        """Get filtered metadata excluding one filter to show available options"""
        temp_metadata = base_metadata.copy()
        
        # Apply all filters except the excluded one
        for filter_name, filter_value in selected_filters.items():
            if filter_name != exclude_filter and filter_value != "All":
                if filter_name == 'enhancer':
                    temp_metadata = temp_metadata[temp_metadata['enhancer_id'] == filter_value]
                elif filter_name == 'cargo':
                    temp_metadata = temp_metadata[temp_metadata['cargo'] == filter_value]
                elif filter_name == 'experiment':
                    temp_metadata = temp_metadata[temp_metadata['experiment'] == filter_value]
                elif filter_name == 'gene':
                    temp_metadata = temp_metadata[temp_metadata['proximal_gene'] == filter_value]
                elif filter_name == 'gc_delivered':
                    temp_metadata = temp_metadata[temp_metadata['GC delivered'] == filter_value]
        
        return temp_metadata
    
    enhancer_metadata = get_options_for_filter('enhancer')
    available_enhancers = sorted([x for x in enhancer_metadata['enhancer_id'].dropna().unique() if x != '']) if not enhancer_metadata.empty else []
    
    cargo_metadata = get_options_for_filter('cargo')
    available_cargos = sorted([x for x in cargo_metadata['cargo'].dropna().unique() if x != '']) if not cargo_metadata.empty else []
    
    experiment_metadata = get_options_for_filter('experiment')
    available_experiments = sorted([x for x in experiment_metadata['experiment'].dropna().unique() if x != '']) if not experiment_metadata.empty else []
    
    gene_metadata = get_options_for_filter('gene')
    available_genes = sorted([x for x in gene_metadata['proximal_gene'].dropna().unique() if x != '']) if not gene_metadata.empty else []
    
    gc_metadata = get_options_for_filter('gc_delivered')
    available_gc_delivered = sorted([x for x in gc_metadata['GC delivered'].dropna().unique() if x != '']) if not gc_metadata.empty else []
    
    # Cell type stays independent - show all available cell types
    base_cell_types = sorted([f"cell_type_{i}" for i in range(1, 35)], key=lambda x: int(x.split('_')[2]))
    available_cell_types = base_cell_types
    
    return {
        'enhancers': available_enhancers,
        'cargos': available_cargos,
        'experiments': available_experiments,
        'genes': available_genes,
        'gc_delivered': available_gc_delivered,
        'cell_types': available_cell_types
    }

# Application header
st.markdown('<div class="main-header"><h1>🧬 CHEAT SHEET FOR ENHANCER ANALYSIS</h1><p>GitHub → Posit Cloud Deployment</p></div>', unsafe_allow_html=True)

# Load data
peak_data, hof_enhancers, enhancer_metadata = load_data()

if peak_data is None:
    st.stop()

# Initialize visualization generator
viz_generator = VisualizationGenerator()

# Initialize session state for filters
if 'filter_state' not in st.session_state:
    st.session_state.filter_state = {
        'enhancer': 'All',
        'cargo': 'All',
        'experiment': 'All',
        'gene': 'All',
        'gc_delivered': 'All',
        'cell_type': 'All'
    }

# Get current filter options based on selections
current_options = get_filtered_options(st.session_state.filter_state, enhancer_metadata)

# Sidebar filters
st.sidebar.markdown("## 🔍 Filter Options")
st.sidebar.markdown("*Deployed from GitHub Repository*")

selected_enhancer = st.sidebar.selectbox(
    "Select Enhancer",
    options=["All"] + current_options['enhancers'],
    index=0 if st.session_state.filter_state['enhancer'] == 'All' or st.session_state.filter_state['enhancer'] not in current_options['enhancers'] else current_options['enhancers'].index(st.session_state.filter_state['enhancer']) + 1,
    help="Choose a specific Hall of Fame enhancer"
)

selected_cargo = st.sidebar.selectbox(
    "Filter by Cargo",
    options=["All"] + current_options['cargos'],
    index=0 if st.session_state.filter_state['cargo'] == 'All' or st.session_state.filter_state['cargo'] not in current_options['cargos'] else current_options['cargos'].index(st.session_state.filter_state['cargo']) + 1,
    help="Filter by delivered genetic cargo"
)

selected_experiment = st.sidebar.selectbox(
    "Filter by Experiment",
    options=["All"] + current_options['experiments'],
    index=0 if st.session_state.filter_state['experiment'] == 'All' or st.session_state.filter_state['experiment'] not in current_options['experiments'] else current_options['experiments'].index(st.session_state.filter_state['experiment']) + 1,
    help="Filter by experiment identifier"
)

selected_gene = st.sidebar.selectbox(
    "Filter by Proximal Gene",
    options=["All"] + current_options['genes'],
    index=0 if st.session_state.filter_state['gene'] == 'All' or st.session_state.filter_state['gene'] not in current_options['genes'] else current_options['genes'].index(st.session_state.filter_state['gene']) + 1,
    help="Filter by nearest gene"
)

selected_gc_delivered = st.sidebar.selectbox(
    "Filter by GC Delivered",
    options=["All"] + current_options['gc_delivered'],
    index=0 if st.session_state.filter_state['gc_delivered'] == 'All' or st.session_state.filter_state['gc_delivered'] not in current_options['gc_delivered'] else current_options['gc_delivered'].index(st.session_state.filter_state['gc_delivered']) + 1,
    help="Filter by genome copies delivered"
)

selected_cell_type = st.sidebar.selectbox(
    "Filter by Cell Type",
    options=["All"] + current_options['cell_types'],
    index=0 if st.session_state.filter_state['cell_type'] == 'All' or st.session_state.filter_state['cell_type'] not in current_options['cell_types'] else current_options['cell_types'].index(st.session_state.filter_state['cell_type']) + 1,
    help="Focus on specific cell type for visualization"
)

# Update session state with new selections
st.session_state.filter_state = {
    'enhancer': selected_enhancer,
    'cargo': selected_cargo,
    'experiment': selected_experiment,
    'gene': selected_gene,
    'gc_delivered': selected_gc_delivered,
    'cell_type': selected_cell_type
}

# Apply metadata filters
enhancer_ids_to_include = set(hof_enhancers['enhancer_id'].unique())

if enhancer_metadata is not None and not enhancer_metadata.empty:
    if any(filter_val != "All" for filter_val in [selected_cargo, selected_experiment, selected_gene, selected_gc_delivered]):
        filtered_metadata = enhancer_metadata.copy()
        
        if selected_cargo != "All":
            filtered_metadata = filtered_metadata[filtered_metadata['cargo'] == selected_cargo]
        
        if selected_experiment != "All":
            filtered_metadata = filtered_metadata[filtered_metadata['experiment'] == selected_experiment]
            
        if selected_gene != "All":
            filtered_metadata = filtered_metadata[filtered_metadata['proximal_gene'] == selected_gene]
            
        if selected_gc_delivered != "All":
            filtered_metadata = filtered_metadata[filtered_metadata['GC delivered'] == selected_gc_delivered]
        
        # Get the enhancer IDs that match the metadata filters
        matching_enhancer_ids = set(filtered_metadata['enhancer_id'].unique())
        enhancer_ids_to_include = enhancer_ids_to_include.intersection(matching_enhancer_ids)
    
    # Filter HOF enhancers to only include those that pass metadata filters - ONE record per enhancer
    filtered_enhancers = hof_enhancers[hof_enhancers['enhancer_id'].isin(list(enhancer_ids_to_include))].drop_duplicates(subset=['enhancer_id'], keep='first').copy()
    relevant_metadata = enhancer_metadata[enhancer_metadata['enhancer_id'].isin(list(enhancer_ids_to_include))] if enhancer_metadata is not None and not enhancer_metadata.empty else pd.DataFrame()

# Display results
if filtered_enhancers.empty:
    st.warning("⚠️ No enhancers match the selected filters. Please adjust your filter criteria.")
elif selected_enhancer == "All":
    st.info(f"📊 Select a specific enhancer from the dropdown above to view detailed analysis")
    st.markdown("### Available Enhancers")
    
    # Show summary table of all enhancers
    if not filtered_enhancers.empty:
        # Get basic info for display
        display_data = []
        for _, row in filtered_enhancers.iterrows():
            enhancer_id = row['enhancer_id']
            chr_info = row['chr']
            start_pos = int(row['start'])
            end_pos = int(row['end'])
            length = end_pos - start_pos
            
            # Get metadata if available
            enhancer_meta = relevant_metadata[relevant_metadata['enhancer_id'] == enhancer_id]
            if not enhancer_meta.empty:
                cargo = enhancer_meta.iloc[0]['cargo']
                experiment = enhancer_meta.iloc[0]['experiment']
                gene = enhancer_meta.iloc[0]['proximal_gene']
                gc_delivered = enhancer_meta.iloc[0]['GC delivered']
                unique_experiments = enhancer_meta['experiment'].nunique()
            else:
                cargo = "N/A"
                experiment = "N/A"
                gene = "N/A"
                gc_delivered = "N/A"
                unique_experiments = 0
            
            display_data.append({
                'Enhancer': enhancer_id,
                'Location': f"{chr_info}:{start_pos:,}-{end_pos:,}",
                'Length (bp)': f"{length:,}",
                'Cargo': cargo,
                'Experiment': experiment,
                'Gene': gene,
                'GC Delivered': gc_delivered,
                'Experiments': unique_experiments
            })
        
        summary_df = pd.DataFrame(display_data)
        st.dataframe(summary_df, use_container_width=True)
        
        st.markdown(f"**Total Enhancers: {len(filtered_enhancers)}**")
        
        # Show overall statistics
        if not relevant_metadata.empty:
            st.markdown("### Summary Statistics")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Unique Cargos", relevant_metadata['cargo'].nunique())
            with col2:
                st.metric("Unique Experiments", relevant_metadata['experiment'].nunique())
            with col3:
                st.metric("Unique Genes", relevant_metadata['proximal_gene'].nunique())

else:
    # Show detailed view for selected enhancer
    selected_enhancer_data = filtered_enhancers[filtered_enhancers['enhancer_id'] == selected_enhancer]
    
    if not selected_enhancer_data.empty:
        enhancer_row = selected_enhancer_data.iloc[0]
        enhancer_id = enhancer_row['enhancer_id']
        chr_info = enhancer_row['chr']
        start_pos = int(enhancer_row['start'])
        end_pos = int(enhancer_row['end'])
        length = end_pos - start_pos
        
        # Get metadata for this enhancer
        enhancer_meta = relevant_metadata[relevant_metadata['enhancer_id'] == enhancer_id]
        
        if not enhancer_meta.empty:
            cargo = enhancer_meta.iloc[0]['cargo']
            experiment = enhancer_meta.iloc[0]['experiment']
            gene = enhancer_meta.iloc[0]['proximal_gene']
            
            # Enhanced header with metadata
            st.markdown(f'<div class="enhancer-card">', unsafe_allow_html=True)
            st.markdown(f"# {enhancer_id}")
            st.markdown(f"**{cargo} • {experiment} • {gene} • {chr_info}:{start_pos:,}-{end_pos:,} ({length:,}bp)**")
            st.markdown('</div>', unsafe_allow_html=True)
            
            # Get all imaging data for this enhancer from the metadata dataframe
            if not relevant_metadata.empty:
                # Get the metadata row for this specific enhancer, experiment type, and GC delivered
                enhancer_meta_row = relevant_metadata[
                    (relevant_metadata['enhancer_id'] == enhancer_id) & 
                    (relevant_metadata['experiment'] == selected_experiment)
                ]
                
                # If GC delivered filter is set, apply it too
                if selected_gc_delivered != "All":
                    enhancer_meta_row = enhancer_meta_row[
                        enhancer_meta_row['GC delivered'] == selected_gc_delivered
                    ]
                
                # If no exact match, use fallback (any experiment, any GC delivered)
                if enhancer_meta_row.empty:
                    enhancer_meta_row = relevant_metadata[relevant_metadata['enhancer_id'] == enhancer_id]
                
                if not enhancer_meta_row.empty:
                    # Use the first matching record
                    meta_row = enhancer_meta_row.iloc[0]
                    
                    # Display imaging if available
                    st.markdown("### 📸 Imaging Data")
                    
                    # Get imaging URLs with priority logic
                    image_link = str(meta_row.get('image_link', ''))
                    neuroglancer_1 = str(meta_row.get('neuroglancer_1', ''))
                    neuroglancer_3 = str(meta_row.get('neuroglancer_3', ''))
                    viewer_link = str(meta_row.get('viewer_link', ''))
                    coronal_mip = str(meta_row.get('coronal_mip', ''))
                    sagittal_mip = str(meta_row.get('sagittal_mip', ''))
                    
                    # Priority logic based on experiment type
                    experiment_type = meta_row.get('experiment', '')
                    
                    if 'lightsheet' in experiment_type.lower():
                        # Lightsheet: prioritize Neuroglancer and MIP projections
                        if neuroglancer_1 and neuroglancer_1 != 'nan':
                            st.markdown("#### Neuroglancer Viewer")
                            st.markdown(f'<iframe src="{neuroglancer_1}" width="100%" height="600"></iframe>', unsafe_allow_html=True)
                        
                        if coronal_mip and coronal_mip != 'nan' and coronal_mip.lower() != 'false':
                            st.markdown("#### Coronal MIP Projection")
                            st.image(coronal_mip, use_column_width=True)
                        
                        if sagittal_mip and sagittal_mip != 'nan' and sagittal_mip.lower() != 'false':
                            st.markdown("#### Sagittal MIP Projection")
                            st.image(sagittal_mip, use_column_width=True)
                    
                    else:
                        # EPI experiments: prioritize contact sheets
                        if image_link and image_link != 'nan':
                            urls = [url.strip() for url in image_link.split(',')]
                            contact_sheet_urls = [url for url in urls if 'contact_sheet' in url.lower()]
                            
                            if contact_sheet_urls:
                                st.markdown("#### Contact Sheet")
                                st.markdown(f'<iframe src="{contact_sheet_urls[0]}" width="100%" height="800"></iframe>', unsafe_allow_html=True)
                        
                        if neuroglancer_1 and neuroglancer_1 != 'nan':
                            st.markdown("#### Neuroglancer Viewer")
                            st.markdown(f'<iframe src="{neuroglancer_1}" width="100%" height="600"></iframe>', unsafe_allow_html=True)
            
            # Create and display accessibility visualization
            st.markdown("### 📊 Peak Accessibility Profile")
            
            # Get peak data for this enhancer
            enhancer_peak_data = peak_data[peak_data['enhancer_id'] == enhancer_id]
            
            if not enhancer_peak_data.empty:
                # Create visualization
                fig = viz_generator.create_peak_visualization(enhancer_peak_data, enhancer_id)
                
                # Apply cell type filter if specified
                if selected_cell_type != "All":
                    # Filter the plot to show only selected cell type
                    cell_type_data = enhancer_peak_data[enhancer_peak_data['cell_type'] == selected_cell_type]
                    if not cell_type_data.empty:
                        fig = viz_generator.create_peak_visualization(cell_type_data, enhancer_id)
                        st.markdown(f"**Showing accessibility for {selected_cell_type} only**")
                
                st.markdown('<div class="accessibility-plot">', unsafe_allow_html=True)
                st.plotly_chart(fig, use_container_width=True)
                st.markdown('</div>', unsafe_allow_html=True)
                
                # Show summary statistics
                st.markdown("### 📈 Summary Statistics")
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.metric("Cell Types", enhancer_peak_data['cell_type'].nunique())
                with col2:
                    st.metric("Peak Positions", len(enhancer_peak_data))
                with col3:
                    max_accessibility = enhancer_peak_data.groupby('cell_type')['accessibility'].max().sort_values(ascending=False)
                    if not max_accessibility.empty:
                        top_cell_type = max_accessibility.index[0]
                        st.metric("Top Cell Type", top_cell_type)
                
                # Show detailed data table
                with st.expander("📋 View Raw Data"):
                    st.dataframe(enhancer_peak_data, use_container_width=True)
                    
                    # Download button for data
                    csv = enhancer_peak_data.to_csv(index=False)
                    st.download_button(
                        label="Download Data as CSV",
                        data=csv,
                        file_name=f"{enhancer_id}_accessibility_data.csv",
                        mime="text/csv"
                    )
            else:
                st.warning("No peak accessibility data available for this enhancer")
        else:
            st.warning("No metadata available for this enhancer")
    else:
        st.error("Selected enhancer not found in filtered results")

# Footer
st.markdown("---")
st.markdown("**Hall of Fame Enhancers Analysis** | GitHub → Posit Cloud Deployment | Built with Streamlit")