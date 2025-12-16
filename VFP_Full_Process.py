import subprocess

# 1. Run AoA generation script (interactive)
print("=== Step 1: Generating AoA .dat files ===")
subprocess.run(["python", "VFP_File_Generation_Main.py"])

# 2. Run VFP solver batch execution (automated loop)
print("\n=== Step 2: Running VFP .bat automation ===")
subprocess.run(["python", "VFP_Run_Main_V3.py"])

# 3. Run data extraction from generated .forces and wavedrg files
print("\n=== Step 3: Extracting data and saving CSV ===")
subprocess.run(["python", "VFP_Data_Extraction_Main.py"])
