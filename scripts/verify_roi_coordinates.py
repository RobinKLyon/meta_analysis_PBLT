"""Check MNI-coordinate mapping from the project root without loading TFR data.

Usage:
    python scripts/verify_roi_coordinates.py --roi-file path/to/PBLT_ROIs_all_sheets_insula.csv
"""

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from roi_coordinates_old import build_roi_coordinate_lookup, normalize_patient, normalize_contact


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--roi-file', type=Path, required=True)
    parser.add_argument('--patient', default='GRE_2014_AUDe')
    parser.add_argument('--contacts', nargs='*', default=['X02','X03','X07'])
    args = parser.parse_args()

    lookup = build_roi_coordinate_lookup(args.roi_file)
    missing = []
    for contact in args.contacts:
        key = (normalize_patient(args.patient), normalize_contact(contact))
        coords = lookup.get(key)
        if coords is None:
            print(f'[MISSING] {args.patient}/{contact}: not found in ROI CSV')
            missing.append(contact)
        else:
            print(f'[FOUND] {args.patient}/{contact}: ' + ', '.join(f'{axis}={coords[axis]}' for axis in ('X','Y','Z')))
    if missing:
        raise SystemExit(f'{len(missing)} contact(s) not found. Confirm --roi-file points to the expected CSV.')


if __name__ == '__main__':
    main()
