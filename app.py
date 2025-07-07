import streamlit as st
import pandas as pd
import numpy as np
import pyarrow.feather as feather
from data_processor_chunked import DataProcessor
from visualization import VisualizationGenerator
import os

# Configure page
st.set_page_config(
    page_title="Hall of Fame Enhancers Analysis",
    page_icon="🧬",
    layout="wide"
)

# Initialize data processor
@st.cache_data
def load_data():
    """Load and process all data files"""
    processor = DataProcessor()
    return processor.load_all_data()

# Load data
try:
    peak_data, hof_enhancers, enhancer_metadata = load_data()
    st.success("App is loaded!")
except Exception as e:
    st.error(f"❌ Error loading data: {str(e)}")
    st.stop()

# Compact title
st.markdown("### Cheat Sheet for Enhancer Analysis")

# Sidebar for filters
st.sidebar.header("🔍 Filters & Selection")
st.sidebar.markdown("Choose filters to focus your analysis")

# Smart filtering logic - get base data first
base_enhancers = sorted(hof_enhancers['enhancer_id'].unique()) if hof_enhancers is not None and 'enhancer_id' in hof_enhancers.columns else []
base_metadata = enhancer_metadata[enhancer_metadata['enhancer_id'].isin(base_enhancers)] if enhancer_metadata is not None and not enhancer_metadata.empty else pd.DataFrame()

# Sort cell types numerically by their leading numbers (1-34)
import re
def extract_cell_type_number(cell_type):
    match = re.match(r'^(\d+)', str(cell_type))
    return int(match.group(1)) if match else 999

base_cell_types = sorted(peak_data['cell_type'].unique(), key=extract_cell_type_number) if peak_data is not None and not peak_data.empty else []

# Initialize session state for filters if not exists
if 'filter_state' not in st.session_state:
    st.session_state.filter_state = {
        'enhancer': 'All',
        'cargo': 'All', 
        'experiment': 'All',
        'gene': 'All',
        'gc_delivered': 'All',
        'cell_type': 'All'
    }

# Smart filter function - updates available options based on current selections
def get_filtered_options(selected_filters, base_metadata):
    """Get available options for each filter based on current selections - cell type remains independent"""
    
    # For each filter, determine what should be available based on OTHER selected filters
    def get_options_for_filter(exclude_filter):
        current_metadata = base_metadata.copy()
        
        # Apply all filters EXCEPT the one we're calculating options for
        for filter_name, filter_value in selected_filters.items():
            if filter_name == exclude_filter or filter_name == 'cell_type':  # Skip cell type and target filter
                continue
            if filter_value != 'All':
                if filter_name == 'enhancer':
                    current_metadata = current_metadata[current_metadata['enhancer_id'] == filter_value]
                elif filter_name == 'cargo':
                    current_metadata = current_metadata[current_metadata['cargo'] == filter_value]
                elif filter_name == 'experiment':
                    current_metadata = current_metadata[current_metadata['experiment'] == filter_value]
                elif filter_name == 'gene':
                    current_metadata = current_metadata[current_metadata['proximal_gene'] == filter_value]
                elif filter_name == 'gc_delivered':
                    current_metadata = current_metadata[current_metadata['GC delivered'] == filter_value]
        
        return current_metadata
    
    # Get available options for each filter
    enhancer_metadata = get_options_for_filter('enhancer')
    available_enhancers = sorted(enhancer_metadata['enhancer_id'].unique()) if not enhancer_metadata.empty else []
    
    cargo_metadata = get_options_for_filter('cargo')
    available_cargos = sorted([x for x in cargo_metadata['cargo'].dropna().unique() if x != '']) if not cargo_metadata.empty else []
    
    experiment_metadata = get_options_for_filter('experiment')
    available_experiments = sorted([x for x in experiment_metadata['experiment'].dropna().unique() if x != '']) if not experiment_metadata.empty else []
    
    gene_metadata = get_options_for_filter('gene')
    available_genes = sorted([x for x in gene_metadata['proximal_gene'].dropna().unique() if x != '']) if not gene_metadata.empty else []
    
    gc_metadata = get_options_for_filter('gc_delivered')
    available_gc_delivered = sorted([x for x in gc_metadata['GC delivered'].dropna().unique() if x != '']) if not gc_metadata.empty else []
    
    # Cell type stays independent - always show all cell types
    available_cell_types = base_cell_types
    
    return {
        'enhancers': available_enhancers,
        'cargos': available_cargos,
        'experiments': available_experiments, 
        'genes': available_genes,
        'gc_delivered': available_gc_delivered,
        'cell_types': available_cell_types
    }

# Get current filter options
current_options = get_filtered_options(st.session_state.filter_state, base_metadata)

# Sidebar filter controls with smart filtering
selected_enhancer = st.sidebar.selectbox(
    "Select Enhancer",
    options=["All"] + current_options['enhancers'],
    index=0 if st.session_state.filter_state['enhancer'] == 'All' or st.session_state.filter_state['enhancer'] not in current_options['enhancers'] else current_options['enhancers'].index(st.session_state.filter_state['enhancer']) + 1,
    help="Choose a specific enhancer to analyze"
)

selected_cargo = st.sidebar.selectbox(
    "Filter by Cargo",
    options=["All"] + current_options['cargos'],
    index=0 if st.session_state.filter_state['cargo'] == 'All' or st.session_state.filter_state['cargo'] not in current_options['cargos'] else current_options['cargos'].index(st.session_state.filter_state['cargo']) + 1,
    help="Filter by experimental cargo type"
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
    help="Filter by cell type for accessibility tracks"
)

# Update session state
st.session_state.filter_state = {
    'enhancer': selected_enhancer,
    'cargo': selected_cargo,
    'experiment': selected_experiment,
    'gene': selected_gene,
    'gc_delivered': selected_gc_delivered,
    'cell_type': selected_cell_type
}

# Initialize filtered_enhancers and relevant_metadata
filtered_enhancers = pd.DataFrame()
relevant_metadata = pd.DataFrame()

# Apply metadata filters only if data is available
if hof_enhancers is not None and not hof_enhancers.empty:
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
        relevant_metadata = enhancer_metadata[enhancer_metadata['enhancer_id'].isin(list(enhancer_ids_to_include))]
    else:
        # No metadata available, use all HOF enhancers
        filtered_enhancers = hof_enhancers.drop_duplicates(subset=['enhancer_id'], keep='first').copy()
else:
    # No HOF enhancers available
    st.error("No Hall of Fame enhancers data available. Please check data files.")

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
            enhancer_meta = relevant_metadata[relevant_metadata['enhancer_id'] == enhancer_id] if not relevant_metadata.empty and 'enhancer_id' in relevant_metadata.columns else pd.DataFrame()
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
                'Location': f"{chr_info}:{start_pos}-{end_pos}",
                'Length (bp)': length,
                'Cargo': cargo,
                'Experiment': experiment,
                'Gene': gene,
                'GC Delivered': gc_delivered,
                'Experiments': unique_experiments
            })
        
        # Display as DataFrame
        summary_df = pd.DataFrame(display_data)
        st.dataframe(summary_df, use_container_width=True)
        
        st.markdown(f"**Total enhancers:** {len(filtered_enhancers)}")
        
else:
    # Show detailed view for selected enhancer
    selected_row = filtered_enhancers[filtered_enhancers['enhancer_id'] == selected_enhancer]
    
    if not selected_row.empty:
        enhancer_row = selected_row.iloc[0]
        enhancer_id = enhancer_row['enhancer_id']
        chr_info = enhancer_row['chr']
        start_pos = int(enhancer_row['start'])
        end_pos = int(enhancer_row['end'])
        length = end_pos - start_pos
        
        # Get metadata for this enhancer
        enhancer_meta = relevant_metadata[relevant_metadata['enhancer_id'] == enhancer_id] if not relevant_metadata.empty and 'enhancer_id' in relevant_metadata.columns else pd.DataFrame()
        
        if not enhancer_meta.empty:
            cargo = enhancer_meta.iloc[0]['cargo']
            experiment = enhancer_meta.iloc[0]['experiment']
            gene = enhancer_meta.iloc[0]['proximal_gene']
            
            # Enhanced header with metadata
            st.markdown(f'<div style="background: linear-gradient(90deg, #667eea 0%, #764ba2 100%); padding: 20px; border-radius: 10px; margin-bottom: 20px;">', unsafe_allow_html=True)
            st.markdown(f'<h2 style="color: white; margin: 0; font-size: 24px;">{enhancer_id}</h2>', unsafe_allow_html=True)
            st.markdown(f'<p style="color: white; margin: 5px 0; font-size: 16px;"><strong>Location:</strong> {chr_info}:{start_pos:,}-{end_pos:,} ({length:,} bp)</p>', unsafe_allow_html=True)
            st.markdown(f'<p style="color: white; margin: 5px 0; font-size: 16px;"><strong>Cargo:</strong> {cargo} | <strong>Experiment:</strong> {experiment} | <strong>Gene:</strong> {gene}</p>', unsafe_allow_html=True)
            st.markdown('</div>', unsafe_allow_html=True)
            
            # Get all imaging data for this enhancer from the metadata dataframe
            if not relevant_metadata.empty:
                # Get the metadata row for this specific enhancer, experiment type, and GC delivered
                if 'enhancer_id' in relevant_metadata.columns and 'experiment' in relevant_metadata.columns:
                    enhancer_meta_row = relevant_metadata[
                        (relevant_metadata['enhancer_id'] == enhancer_id) & 
                        (relevant_metadata['experiment'] == selected_experiment)
                    ]
                    
                    # If GC delivered filter is set, apply it too
                    if selected_gc_delivered != "All" and 'GC delivered' in relevant_metadata.columns:
                        enhancer_meta_row = enhancer_meta_row[
                            enhancer_meta_row['GC delivered'] == selected_gc_delivered
                        ]
                    
                    # If no exact match, use fallback (any experiment, any GC delivered)
                    if enhancer_meta_row.empty:
                        enhancer_meta_row = relevant_metadata[relevant_metadata['enhancer_id'] == enhancer_id]
                else:
                    enhancer_meta_row = pd.DataFrame()
                
                if not enhancer_meta_row.empty:
                    # Use the first matching record
                    meta_row = enhancer_meta_row.iloc[0]
                    
                    # Display imaging if available
                    st.markdown("### 📸 Imaging Data")
                    
                    # Get the imaging URLs from the metadata
                    image_link = meta_row.get('image_link', '')
                    neuroglancer_1 = meta_row.get('neuroglancer_1', '')
                    neuroglancer_3 = meta_row.get('neuroglancer_3', '')
                    viewer_link = meta_row.get('viewer_link', '')
                    coronal_mip = meta_row.get('coronal_mip', '')
                    sagittal_mip = meta_row.get('sagittal_mip', '')
                    
                    # Initialize viewer components
                    viewer_components = []
                    contact_sheets = []
                    
                    # Process image_link (contact sheets)
                    if image_link and image_link != 'FALSE' and pd.notna(image_link):
                        if ',' in str(image_link):
                            contact_sheets.extend([url.strip() for url in str(image_link).split(',')])
                        else:
                            contact_sheets.append(str(image_link).strip())
                    
                    # Process neuroglancer viewers
                    if neuroglancer_1 and neuroglancer_1 != 'FALSE' and pd.notna(neuroglancer_1):
                        if ',' in str(neuroglancer_1):
                            viewer_components.extend([url.strip() for url in str(neuroglancer_1).split(',')])
                        else:
                            viewer_components.append(str(neuroglancer_1).strip())
                    
                    if neuroglancer_3 and neuroglancer_3 != 'FALSE' and pd.notna(neuroglancer_3):
                        if ',' in str(neuroglancer_3):
                            viewer_components.extend([url.strip() for url in str(neuroglancer_3).split(',')])
                        else:
                            viewer_components.append(str(neuroglancer_3).strip())
                    
                    # Process viewer_link
                    if viewer_link and viewer_link != 'FALSE' and pd.notna(viewer_link):
                        if ',' in str(viewer_link):
                            viewer_components.extend([url.strip() for url in str(viewer_link).split(',')])
                        else:
                            viewer_components.append(str(viewer_link).strip())
                    
                    # Process MIP projections
                    mip_projections = []
                    if coronal_mip and coronal_mip != 'FALSE' and pd.notna(coronal_mip):
                        if ',' in str(coronal_mip):
                            mip_projections.extend([url.strip() for url in str(coronal_mip).split(',')])
                        else:
                            mip_projections.append(str(coronal_mip).strip())
                    
                    if sagittal_mip and sagittal_mip != 'FALSE' and pd.notna(sagittal_mip):
                        if ',' in str(sagittal_mip):
                            mip_projections.extend([url.strip() for url in str(sagittal_mip).split(',')])
                        else:
                            mip_projections.append(str(sagittal_mip).strip())
                    
                    # Display imaging based on experiment type
                    experiment_type = meta_row.get('experiment', '').upper()
                    
                    if experiment_type == 'LIGHTSHEET':
                        # For lightsheet: prioritize Neuroglancer viewers and MIP projections
                        if viewer_components:
                            st.markdown("#### 🧠 Neuroglancer Viewers")
                            for i, viewer_url in enumerate(viewer_components):
                                if viewer_url:
                                    st.markdown(f"**Viewer {i+1}:**")
                                    st.components.v1.iframe(viewer_url, height=600)
                        
                        if mip_projections:
                            st.markdown("#### 📊 MIP Projections")
                            for i, mip_url in enumerate(mip_projections):
                                if mip_url:
                                    st.markdown(f"**MIP Projection {i+1}:**")
                                    st.image(mip_url, use_column_width=True)
                        
                        if contact_sheets:
                            st.markdown("#### 📷 Contact Sheets")
                            for i, sheet_url in enumerate(contact_sheets):
                                if sheet_url:
                                    st.markdown(f"**Contact Sheet {i+1}:**")
                                    st.image(sheet_url, use_column_width=True)
                                    
                    elif experiment_type == 'EPI':
                        # For EPI: prioritize contact sheets first, then Neuroglancer
                        if contact_sheets:
                            st.markdown("#### 📷 Contact Sheets")
                            for i, sheet_url in enumerate(contact_sheets):
                                if sheet_url:
                                    st.markdown(f"**Contact Sheet {i+1}:**")
                                    st.image(sheet_url, use_column_width=True)
                        
                        if viewer_components:
                            st.markdown("#### 🧠 Neuroglancer Viewers")
                            for i, viewer_url in enumerate(viewer_components):
                                if viewer_url:
                                    st.markdown(f"**Viewer {i+1}:**")
                                    st.components.v1.iframe(viewer_url, height=600)
                        
                        if mip_projections:
                            st.markdown("#### 📊 MIP Projections")
                            for i, mip_url in enumerate(mip_projections):
                                if mip_url:
                                    st.markdown(f"**MIP Projection {i+1}:**")
                                    st.image(mip_url, use_column_width=True)
                    
                    else:
                        # Default: show all available imaging
                        if contact_sheets:
                            st.markdown("#### 📷 Contact Sheets")
                            for i, sheet_url in enumerate(contact_sheets):
                                if sheet_url:
                                    st.markdown(f"**Contact Sheet {i+1}:**")
                                    st.image(sheet_url, use_column_width=True)
                        
                        if viewer_components:
                            st.markdown("#### 🧠 Neuroglancer Viewers")
                            for i, viewer_url in enumerate(viewer_components):
                                if viewer_url:
                                    st.markdown(f"**Viewer {i+1}:**")
                                    st.components.v1.iframe(viewer_url, height=600)
                        
                        if mip_projections:
                            st.markdown("#### 📊 MIP Projections")
                            for i, mip_url in enumerate(mip_projections):
                                if mip_url:
                                    st.markdown(f"**MIP Projection {i+1}:**")
                                    st.image(mip_url, use_column_width=True)
                    
                    # Show message if no imaging data available
                    if not any([contact_sheets, viewer_components, mip_projections]):
                        st.info("No imaging data available for this enhancer.")
                    
                else:
                    st.info("No imaging metadata available for this enhancer.")
        
        else:
            st.warning("No metadata available for this enhancer.")
        
        # Get peak data for selected enhancer
        if peak_data is not None and not peak_data.empty:
            enhancer_peaks = peak_data[peak_data['enhancer_id'] == enhancer_id]
            
            if not enhancer_peaks.empty:
                st.markdown("### 📊 Peak Accessibility Analysis")
                
                # Apply cell type filter if specified
                if selected_cell_type != "All":
                    enhancer_peaks = enhancer_peaks[enhancer_peaks['cell_type'] == selected_cell_type]
                
                if not enhancer_peaks.empty:
                    # Generate visualization
                    viz_generator = VisualizationGenerator()
                    
                    if selected_cell_type == "All":
                        # Show all cell types
                        fig = viz_generator.create_peak_visualization(enhancer_peaks, enhancer_id)
                    else:
                        # Show specific cell type
                        fig = viz_generator.create_cell_type_specific_view(enhancer_peaks, selected_cell_type)
                    
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.info(f"No accessibility data available for {selected_cell_type} cell type.")
            else:
                st.info("No peak accessibility data available for this enhancer.")
        else:
            st.info("Peak data not available.")
    else:
        st.error("Selected enhancer not found in the data.")