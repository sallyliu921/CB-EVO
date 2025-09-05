#!/usr/bin/env python3
"""
Extract optimization results from CBTune log files.
This script parses log files and extracts LUT area, level and optimization time data.
"""

import pandas as pd
import argparse
from dataclasses import dataclass
from typing import List, Optional

@dataclass
class OptimizationResult:
    """Data class to store optimization results."""
    lut_area: int
    lut_level: int
    time: float

class LogParser:
    """Parser for CBTune log files."""
    
    def __init__(self, log_path: str):
        self.log_path = log_path
        self.results: List[OptimizationResult] = []

    def parse(self) -> None:
        """Parse the log file and extract optimization results."""
        try:
            with open(self.log_path, 'r') as file:
                lines = file.readlines()
                
            for i in range(1, len(lines)):
                if "Total decision and optimization time :" in lines[i]:
                    # Extract time
                    time_str = lines[i].split(':')[-1].strip('s').strip()
                    
                    # Extract LUT area and level
                    prev_line = lines[i-1].strip()
                    parts = prev_line.split()
                    lut_area = int(parts[2])
                    lut_level = int(parts[5])
                    
                    # Store result
                    result = OptimizationResult(
                        lut_area=lut_area,
                        lut_level=lut_level,
                        time=float(time_str)
                    )
                    self.results.append(result)
        except FileNotFoundError:
            print(f"Error: Log file not found at {self.log_path}")
            return

    def save_to_excel(self, output_path: str) -> None:
        """Save parsed results to Excel file."""
        if not self.results:
            print("No results to save")
            return

        # Convert results to DataFrame
        data = [[r.lut_area, r.lut_level, r.time] for r in self.results]
        df = pd.DataFrame(data, columns=['LUT Area', 'LUT Level', 'Time (s)'])
        df.to_excel(output_path, index=False)
        print(f"Results saved to '{output_path}'")

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Extract optimization results from CBTune log files',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        'input_file',
        help='Path to the CBTune log file'
    )
    parser.add_argument(
        '-o', '--output',
        default='output.xlsx',
        help='Path to save the output Excel file'
    )
    return parser.parse_args()

def main() -> None:
    """Main entry point."""
    args = parse_args()
    
    # Parse log file
    parser = LogParser(args.input_file)
    parser.parse()
    
    # Save results
    parser.save_to_excel(args.output)

if __name__ == '__main__':
    main()