import os
import shutil
from pathlib import Path

def collect_and_copy_vdb_files(source_dir, dest_dir):
    os.makedirs(dest_dir, exist_ok=True)
    vdb_files = list(Path(source_dir).rglob("*.vdb"))

    for vdb_file in vdb_files:
        dest_file = Path(dest_dir) / vdb_file.name
        shutil.copy(vdb_file, dest_file)
        print(f"Copied: {vdb_file} -> {dest_file}")

    print(f"Total .vdb files copied: {len(vdb_files)}")

def delete_directory(dir_path):
    shutil.rmtree(dir_path)
    print(f"Deleted directory: {dir_path}")

def main(unzipped_source_dir, vdb_output_dir, delete_temp=False):
    collect_and_copy_vdb_files(unzipped_source_dir, vdb_output_dir)

    if delete_temp:
        delete_directory(unzipped_source_dir)

if __name__ == "__main__":
    unzipped_folder = "path/to/unzipped"
    vdb_output_folder = "path/to/vdb_output"
    delete_unzipped_when_done = False  # Set to True if you want to delete the source

    main(unzipped_folder, vdb_output_folder, delete_unzipped_when_done)

