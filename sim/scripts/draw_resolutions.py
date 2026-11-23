#!/usr/bin/env python3

import argparse
import json
import subprocess
from math import sqrt
from pathlib import Path

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure


def read_params(filename):
	path = Path(__file__).resolve().parents[2] / 'params' / filename
	# FOOT parameter files use C-style comments.
	return json.loads(subprocess.check_output(['cpp', '-P', '-x', 'c', path]))


def main():
	default_output = Path('foot_resolutions.png')
	parser = argparse.ArgumentParser(
		description='Plot the per-ion position resolutions from params/foot_cov.json as a single PNG.',
		epilog=f'Example: {Path(__file__).name} -1 -o resolutions.png',
	)
	parser.add_argument('-1', '--one-index', action='store_true',
	                    help='label the detectors 1–8 instead of 0–7')
	parser.add_argument('-o', '--output', type=Path, default=default_output,
	                    help=f'output PNG path (default: {default_output})')
	args = parser.parse_args()

	pairs = read_params('foot_cov.json')['values']
	layers = read_params('foot_setup.json')['DSSSD']['composition']
	silicon = next(layer for layer in layers if layer['material'] in ('Si', 'G4_Si'))
	pitch = silicon['pitch']
	effective_pitch = silicon['effective pitch']
	styles = {
		'12C': ('#0072B2', 'o', r'$^{12}$C'),
		'9C' : ('#D55E00', 's', r'$^{9}$C'),
		'4He': ('#009E73', '^', r'$^{4}$He'),
		'1H' : ('#8C509C', 'D', r'$^{1}$H'),
	}
	resolutions = {ion: [] for ion in styles}
	for pair in range(1, 5):
		entries = {entry['ion']: entry for entry in pairs[f'pair {pair}']['position']}
		for ion, widths in resolutions.items():
			resolution = entries[ion]['resolution']
			# Physical detector order within each FOOT pair is v, then u.
			widths.extend((resolution['v'], resolution['u']))

	first = 1 if args.one_index else 0
	indices = range(first, first + 8)
	figure = Figure(figsize=(8, 4.8), layout='constrained')
	FigureCanvasAgg(figure)
	ax = figure.subplots()
	for ion, (color, marker, label) in styles.items():
		ax.plot(indices, resolutions[ion], color=color, marker=marker, linestyle='-', label=label)
	ax.axhline(pitch / sqrt(12), color='firebrick', linestyle='--', zorder=1,
	           label=r'$\mathrm{Pitch}/\sqrt{12}$')
	ax.axhline(effective_pitch / sqrt(12), color='tomato', linestyle='--', zorder=1,
	           label=rf'${effective_pitch / pitch:g}\times\mathrm{{Pitch}}/\sqrt{{12}}$')
	ax.set_xticks(indices, labels=[f'FOOT{i}' for i in indices])
	ax.set(title='Detector Resolutions per Ion Species', ylabel='Intrinsic position resolution [µm]')
	ax.spines[['top', 'right']].set_visible(False)
	ax.set_axisbelow(True)
	ax.grid(axis='y', alpha=0.3)
	ax.legend()
	figure.savefig(args.output, format='png', dpi=180)


if __name__ == '__main__':
	main()
