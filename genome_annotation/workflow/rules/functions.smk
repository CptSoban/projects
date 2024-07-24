import os
import glob

def get_final_output():
    final_output = expand(
        "{work_dir}/workflow/results/braker/{run_id}/braker.gtf", work_dir = config["work_dir"], run_id = config["Run ID"]
    )
    return final_output


def get_basenames(directory_path):
    # Get a list of all files in the given directory
    files = glob.glob(os.path.join(directory_path, "*"))

    # Initialize an empty list to store the base IDs
    base_ids = []

    # Iterate over each file path in the list of files
    for file_path in files:
        # Check if the current path is a file (not a directory)
        if os.path.isfile(file_path):
            # Get the file name from the file path
            file_name = os.path.basename(file_path)
            # Split the file name into the base name and extension
            basename, _ = os.path.splitext(file_name)
            # Extract the base ID by splitting the base name at the first underscore
            base_id = basename.split('_')[0]
            # Add the base ID to the list of base IDs
            base_ids.append(base_id)

    # Join the list of base IDs into a single string, separated by commas
    base_ids_string = ','.join(base_ids)
    return base_ids_string
