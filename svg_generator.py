#!/usr/bin/env python3
"""
Standalone SVG generator for QGIS processing results.
This script uses system Python to avoid OSGeo4W matplotlib issues.
"""

import sys
import os
import json
import argparse
from pathlib import Path

def generate_svg_from_data(data_file, output_file):
    """
    Generate SVG from processing data using system Python matplotlib.
    """
    try:
        # Import matplotlib using system Python
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
        from matplotlib.ticker import MultipleLocator
        
        # Load processing data
        with open(data_file, 'r') as f:
            data = json.load(f)
        
        # Create figure
        fig, ax = plt.subplots(1, 1, figsize=(8, 5))
        
        # Plot data (simplified - you'd need to adapt this to your actual data structure)
        if 'plot_s' in data and 'plot_dem' in data:
            ax.plot(data['plot_s'], data['plot_dem'], 'b-', linewidth=1, label='DEM')
        
        # Add title and labels
        ax.set_title(f"Profile at Station {data.get('station', 'Unknown')}")
        ax.set_xlabel("Distance (m)")
        ax.set_ylabel("Elevation (m)")
        
        # Set equal aspect ratio and grid
        ax.set_aspect("equal", adjustable="datalim")
        ax.xaxis.set_major_locator(MultipleLocator(1.0))
        ax.yaxis.set_major_locator(MultipleLocator(1.0))
        ax.grid(True, which="major", axis="both", alpha=0.3)
        
        # Save SVG
        fig.savefig(output_file, format="svg", bbox_inches="tight", pad_inches=0.1)
        plt.close(fig)
        
        print(f"SVG generated: {output_file}")
        return True
        
    except Exception as e:
        print(f"SVG generation failed: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Generate SVG from QGIS processing data')
    parser.add_argument('data_file', help='JSON file with processing data')
    parser.add_argument('output_file', help='Output SVG file')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.data_file):
        print(f"Data file not found: {args.data_file}")
        return 1
    
    success = generate_svg_from_data(args.data_file, args.output_file)
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
