import os
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(ROOT_DIR, "DATA")
path_to_PBLT = os.path.join(DATA_DIR, 'PBLT')
Source_dir = r'D:\sEEG_data\PBLT' # the source directory where the PBLT data is stored
path_to_julich_output = os.path.join(ROOT_DIR, "..", "JULICH_ATLAS_PARTIALISATION")
PBLT_ROI_csv_name = "PBLT_ROIs.csv"