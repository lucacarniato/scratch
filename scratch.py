import os
import shutil
import zipfile

from pathlib import Path

def collect_vfgp_files(source_dir):
    return list(Path(source_dir).rglob("*.vfgp"))

def copy_files(files, dest_dir):
    os.makedirs(dest_dir, exist_ok=True)
    copied_files = []
    for file in files:
        dest_path = Path(dest_dir) / file.name
        shutil.copy(file, dest_path)
        copied_files.append(dest_path)
    return copied_files

def unzip_files(zip_files, extract_to):
    os.makedirs(extract_to, exist_ok=True)
    for zip_file in zip_files:
        try:
            with zipfile.ZipFile(zip_file, 'r') as zip_ref:
                zip_ref.extractall(extract_to / zip_file.stem)
        except zipfile.BadZipFile:
            print(f"Warning: {zip_file} is not a valid zip file")

def collect_and_copy_vdb_files(search_dir, dest_dir):
    os.makedirs(dest_dir, exist_ok=True)
    for vdb_file in Path(search_dir).rglob("*.vdb"):
        shutil.copy(vdb_file, Path(dest_dir) / vdb_file.name)

def main(source_dir, temp_dir, vdb_output_dir):
    # Step 1: Collect .vfgp files
    vfgp_files = collect_vfgp_files(source_dir)
    print(f"Found {len(vfgp_files)} .vfgp files.")

    # Step 2: Copy to temp directory
    copied_vfgp_files = copy_files(vfgp_files, temp_dir)
    
    # Step 3: Unzip files in temp directory
    unzip_dir = Path(temp_dir) / "unzipped"
    unzip_files(copied_vfgp_files, unzip_dir)

    # Step 4: Collect and copy .vdb files to final output directory
    collect_and_copy_vdb_files(unzip_dir, vdb_output_dir)

    # Step 5: Delete the temp directory
    shutil.rmtree(temp_dir)
    print(f"Temporary directory '{temp_dir}' deleted.")

if __name__ == "__main__":
    source_directory = "path/to/source"
    temporary_directory = "path/to/temp"
    vdb_output_directory = "path/to/output_vdbs"
    
    main(source_directory, temporary_directory, vdb_output_directory)
