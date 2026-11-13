
"""
Copy PBLT patient data and create config.json with sampling frequencies.
"""
from pathlib import Path
from datetime import datetime
import shutil
import json
import re
from scipy.io import loadmat
from src.paths_data import path_to_PBLT, Source_dir

# ============================================================
# CONFIGURATION
# ============================================================

SOURCE_DIR = Path(Source_dir)
DEST_DIR = Path(path_to_PBLT)
OVERWRITE = False
CONFIG_PATH = DEST_DIR / "config.json"

def extract_sampling_frequency(patient_dir, patient_id):
    mat_path = patient_dir / f"{patient_id}_PBLT.mat"

    if not mat_path.is_file():
        return {
            "sampling_frequency_hz": None,
            "metadata_file": str(mat_path),
            "status": "file_not_found",
        }

    D = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    fsample = D["D"].Fsample if hasattr(D["D"], "Fsample") else "Fsample not found in .mat file"  # Assuming 'fsample' is a field in the MATLAB struct
    # Convertir la valeur MATLAB en nombre Python
    frequency = float(fsample)

    return {
        "sampling_frequency_hz": frequency,
        "metadata_file": str(mat_path),
        "status": "extracted",
    }


def copy_patient_data(source_dir=SOURCE_DIR, dest_dir=DEST_DIR):
    source_dir = Path(source_dir)
    dest_dir = Path(dest_dir)

    if not source_dir.is_dir():
        raise FileNotFoundError(f"Dossier source introuvable : {source_dir}")

    dest_dir.mkdir(parents=True, exist_ok=True)

    mat_files = list(source_dir.rglob("*_PBLT_TFR.mat"))
    patients = {}

    for mat in mat_files:
        patient_name = mat.name.removesuffix("_TFR.mat")
        patient_id = patient_name.removesuffix("_PBLT")

        key = (str(mat.parent), patient_name)
        patients[key] = {
            "patient_name": patient_id,
            "mat": mat,
            "wya": mat.with_suffix(".wya"),
            "pos": mat.parent / f"{patient_id}_PBLT.pos",
            "source_patient_dir": mat.parent,
        }

    config = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "source_dir": str(source_dir.resolve()),
        "destination_dir": str(dest_dir.resolve()),
        "patients": {},
    }

    copied = skipped = errors = 0

    for (_, patient_name), files in patients.items():
        patient_id = files["patient_name"]

        # One destination folder per patient
        patient_dir = dest_dir / patient_id
        patient_dir.mkdir(parents=True, exist_ok=True)

        records = []

        # Only these three files are copied; original names are preserved.
        for ext in ("mat", "wya", "pos"):
            src = files[ext]
            dst = patient_dir / src.name

            record = {
                "type": f".{ext}",
                "source": str(src),
                "destination": str(dst),
            }

            if not src.is_file():
                record["status"] = "missing"
                print(f"[MANQUANT] {src}")
            else:
                record["size_bytes"] = src.stat().st_size

                try:
                    if dst.exists() and not OVERWRITE:
                        record["status"] = "already_present"
                        skipped += 1
                        print(f"[EXISTE] {dst}")
                    else:
                        shutil.copy2(src, dst)
                        record["status"] = "copied"
                        copied += 1
                        print(f"[COPIÉ] {src.name} -> {dst}")
                except OSError as error:
                    record["status"] = "error"
                    record["error"] = str(error)
                    errors += 1
                    print(f"[ERREUR] {src}: {error}")

            records.append(record)

        # Extract sampling frequency from source metadata.
        sampling_info = extract_sampling_frequency(
            files["source_patient_dir"],
            patient_id,
        )

        config["patients"][patient_id] = {
            "patient_name": patient_name,
            "sampling_frequency": sampling_info,
            "files": records,
            "missing_extensions": [
                r["type"] for r in records
                if r["status"] == "missing"
            ],
        }

        print(
            f"[FRÉQUENCE] {patient_id}: "
            f"{sampling_info['sampling_frequency_hz']} Hz "
            f"({sampling_info['status']})"
        )

    config["summary"] = {
        "number_of_patients": len(patients),
        "copied": copied,
        "already_present": skipped,
        "errors": errors,
        "sampling_frequency_extracted": sum(
            1 for patient in config["patients"].values()
            if patient["sampling_frequency"]["sampling_frequency_hz"]
            is not None
        ),
    }

    with open(CONFIG_PATH, "w", encoding="utf-8") as f:
        json.dump(config, f, indent=4, ensure_ascii=False)

    print("\n========== RÉSUMÉ ==========")
    print(f"Patients : {len(patients)}")
    print(f"Copiés : {copied} | Déjà présents : {skipped} | Erreurs : {errors}")
    print(
        "Fréquences extraites : "
        f"{config['summary']['sampling_frequency_extracted']}"
    )
    print(f"Config : {CONFIG_PATH}")

    return config


config = copy_patient_data()
