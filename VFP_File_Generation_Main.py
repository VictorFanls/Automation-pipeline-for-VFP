import VFP_File_Generation_Utils as u

if __name__ == '__main__':
    # Generate AoA .dat files
    u.run_aoa_generation()

    u.copy_generated_files_to_main_dir()