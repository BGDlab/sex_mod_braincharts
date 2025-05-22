#written by claude
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from collections import defaultdict
import glob
from PyPDF2 import PdfReader, PdfWriter

def main():
    # Base directories to search
    base_dirs = ['cv_sample_B_test', 'cv_sample_A_test']
    
    # Dictionary to store files by their prefix (first part of filename)
    files_by_prefix = defaultdict(lambda: {'csv': [], 'png': []})
    
    # Walk through directories and collect files
    for base_dir in base_dirs:
        if not os.path.exists(base_dir):
            print(f"Warning: Directory {base_dir} does not exist. Skipping.")
            continue
            
        for root, dirs, files in os.walk(base_dir):
            for file in files:
                if file.endswith('.csv') or file.endswith('.png'):
                    # Extract the first part of the filename (e.g., "CSF", "GMV")
                    # assuming format is PREFIX_*
                    match = re.match(r'^([^_]+)', file)
                    if match:
                        prefix = match.group(1)
                        file_path = os.path.join(root, file)
                        
                        if file.endswith('.csv'):
                            files_by_prefix[prefix]['csv'].append(file_path)
                        else:  # PNG
                            files_by_prefix[prefix]['png'].append(file_path)
    
    # First, create a temporary PDF with the table of contents
    with PdfPages('toc.pdf') as pdf:
        # Create a table of contents as the first page
        plt.figure(figsize=(11, 8.5))
        plt.suptitle("Table of Contents", fontsize=16)
        
        content_text = []
        for prefix in sorted(files_by_prefix.keys()):
            content_text.append(f"{prefix} Group:")
            for csv_file in sorted(files_by_prefix[prefix]['csv']):
                content_text.append(f"  CSV: {os.path.basename(csv_file)}")
            for png_file in sorted(files_by_prefix[prefix]['png']):
                content_text.append(f"  PNG: {os.path.basename(png_file)}")
            content_text.append("")  # Add empty line between groups
        
        plt.text(0.1, 0.9, "\n".join(content_text), 
                 va='top', fontsize=10, 
                 transform=plt.gca().transAxes)
        
        plt.axis('off')
        pdf.savefig()
        plt.close()
        
    # Then create the content PDF
    with PdfPages('content.pdf') as pdf:
        # Process each prefix group
        for prefix, files in sorted(files_by_prefix.items()):
            print(f"Processing group: {prefix}")
            
            # Create a summary page for this prefix
            plt.figure(figsize=(11, 8.5))
            plt.suptitle(f"{prefix} Group Summary", fontsize=16)
            plt.axis('off')
            pdf.savefig()
            plt.close()
            
            # Display PNG files
            for i, png_file in enumerate(sorted(files['png'])):
                try:
                    img = plt.imread(png_file)
                    plt.figure(figsize=(11, 8.5))
                    plt.imshow(img)
                    plt.title(f"{os.path.basename(png_file)}")
                    plt.axis('off')
                    pdf.savefig()
                    plt.close()
                    print(f"  Added PNG: {png_file}")
                except Exception as e:
                    print(f"  Error processing {png_file}: {e}")
            
            # Display CSV data
            for i, csv_file in enumerate(sorted(files['csv'])):
                try:
                    # Read CSV with pandas
                    df = pd.read_csv(csv_file)
                    
                    # Display first few rows in a table
                    fig, ax = plt.subplots(figsize=(11, 8.5))
                    ax.axis('off')
                    
                    # Create a table from the DataFrame
                    table_rows = min(10, len(df))  # Display at most 10 rows
                    
                    # Add table to the plot
                    table = ax.table(
                        cellText=df.head(table_rows).values,
                        colLabels=df.columns,
                        loc='center',
                        cellLoc='center',
                        colColours=['lightblue'] * len(df.columns)
                    )
                    
                    # Set font size for table
                    table.auto_set_font_size(False)
                    table.set_fontsize(9)
                    
                    # Adjust layout
                    plt.title(f"CSV Preview: {os.path.basename(csv_file)}")
                    plt.tight_layout()
                    pdf.savefig()
                    plt.close()
                    print(f"  Added CSV: {csv_file}")
                except Exception as e:
                    print(f"  Error processing {csv_file}: {e}")
    
    # Now use PyPDF2 to merge the PDFs with TOC first
    try:
        pdf_writer = PdfWriter()
        
        # Add TOC first
        pdf_reader = PdfReader('toc.pdf')
        for page in pdf_reader.pages:
            pdf_writer.add_page(page)
            
        # Then add content
        pdf_reader = PdfReader('content.pdf')
        for page in pdf_reader.pages:
            pdf_writer.add_page(page)
            
        # Write the combined PDF
        with open('combined_results.pdf', 'wb') as f:
            pdf_writer.write(f)
            
        # Clean up temporary files
        os.remove('toc.pdf')
        os.remove('content.pdf')
        
        print(f"PDF created: combined_results.pdf")
    except Exception as e:
        print(f"Error merging PDFs: {e}")
        print("The individual PDFs 'toc.pdf' and 'content.pdf' have been created instead.")

if __name__ == "__main__":
    main()
