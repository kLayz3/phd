#!/usr/bin/env python3

# Dependencies:
# python-pandas python-openpyxl python-odfpy

import pandas as pd
import json
import os

from pathlib import Path

script_dir = Path(__file__).resolve().parent

df = pd.read_excel(
    script_dir / '../../FRS_run_sheet_2024_carbon.ods',
    sheet_name='S111_2024Feb_EXPERT',
    header=None,
    engine='odf',
)

key_row       = 1 # Excel row  2 -> pandas index  1
starting_row = 25 # Excel row 26 -> pandas index 25
final_row    = 69 # Excel row 70 -> pandas index 69

data = {}

for i in range(starting_row, final_row+1):
	record = {
		str(key): None if pd.isna(value) else value
		for key, value in zip(df.iloc[1], df.iloc[i])
		if pd.notna(key)
	}

	data[str(i + 1)] = record

with open(script_dir / '../params/runsheet.json', 'w') as f:
	json.dump(data, f, indent=2, default=str, allow_nan=False)
