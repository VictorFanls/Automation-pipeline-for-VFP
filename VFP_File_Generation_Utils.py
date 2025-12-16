import re
import os
import shutil

#Change lines in file and change index from 0 -> 1 for continuatino run
def process_file(original_file, mach_str, aoa_str, apply_filter=False): #Filter makes sure that the given file will not change
    pattern = re.compile(rf"{re.escape(mach_str)}\s+{re.escape(aoa_str)}")

    with open(original_file, 'r') as f:
        lines = f.readlines()

    if not apply_filter:
        return lines
            
    indices = [i for i, line in enumerate(lines) if pattern.search(line)]
    if len(indices) >= 3:
        start, end = indices[0], indices[2]
        
        # Change index from 0 to 1 for continuation run

        parts = lines[end].split(' ')
        tmp = list(parts[4])
        index = tmp[0]
        if index == '0' : tmp[0] = '1'
        parts[4] = ''.join(tmp)
        lines[end] = ' '.join(parts)

        filtered = [l for i, l in enumerate(lines) if i < start or i >= end]
    else:
        filtered = lines

    return filtered

#Create different AoA's 
def iterate_AoA_modifications(original_file, mach_str, initial_aoa_str, output_base_name):
    
    aoa_low = float(input("Enter the value for the min angle of attack: "))
    aoa_high = float(input("Enter the value for max angle of attack: "))
    step = float(input("Enter the value for the step between iterations: "))
    aoa_range = aoa_high - aoa_low
    itterations = aoa_range / step
    aoa = float(initial_aoa_str)
    base_pattern = re.compile(rf"({re.escape(mach_str)})\s+({re.escape(initial_aoa_str)})")

    # Use script directory, not current working directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    flow_dir = os.path.join(script_dir, "Flow_Conditions")
    os.makedirs(flow_dir, exist_ok=True)

    for i in range(int(itterations) + 1):
        new_aoa = aoa + i * step
        new_aoa_str = f"{new_aoa:.4f}"

        lines = process_file(original_file, mach_str, initial_aoa_str, apply_filter=(i != 0))
        #print(lines)
        sign = "-" if new_aoa < 0 else "+"
        aoa_for_filename = f"{abs(new_aoa):.2f}".replace('.', 'p')
        output_filename = f"{output_base_name}{sign}{aoa_for_filename}.dat"

        output_file_path = os.path.join(flow_dir, output_filename)

        updated_lines = []
        updated = False

        for line in lines:
            if not updated and base_pattern.search(line):
                line = base_pattern.sub(lambda m: f"{m.group(1)}   {new_aoa_str}", line)
                updated = True
            updated_lines.append(line)

        with open(output_file_path, 'w') as f:
            f.writelines(updated_lines)

        

# Run full AoA generation workflow
def run_aoa_generation():
    original_file = input("Enter the original .dat filename: ").strip()

    match = re.search(r"M(\d{3})[^-+]*([-+]?\d+p\d+)\.dat$", original_file)
    if not match:
        print("Filename format not recognized. Expected format like 'M085Re19p8ma-1p00.dat'")
        return

    mach_str = f"{int(match.group(1)) / 100:.4f}"
    aoa_str = match.group(2).replace('p', '.')
    aoa_str = f"{float(aoa_str):.4f}"
    print(mach_str,aoa_str)
    if not aoa_str.startswith('-') and not aoa_str.startswith('+'):
        if not aoa_str.startswith('-') and not aoa_str.startswith('+'):
            aoa_float = float(aoa_str)
        if aoa_float == 0.0:
            aoa_str = f"{aoa_float:.4f}"  
        elif aoa_float > 0:
            aoa_str = f"+{aoa_float:.4f}"
        else:
            aoa_str = f"{aoa_float:.4f}"


    output_base_name = os.path.splitext(original_file)[0]
    output_base_name = re.sub(r"[-+]\d+p\d+$", "", output_base_name)

    try:
        iterate_AoA_modifications(original_file, mach_str, aoa_str, output_base_name)
    except FileNotFoundError:
        print(f"File '{original_file}' not found.")

def copy_generated_files_to_main_dir():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    flow_dir = os.path.join(script_dir, "Flow_Conditions")

    if not os.path.exists(flow_dir):
        print("Flow_Conditions directory does not exist.")
        return

    dat_files = [f for f in os.listdir(flow_dir) if f.endswith(".dat")]

    for fname in dat_files:
        src_path = os.path.join(flow_dir, fname)
        dst_path = os.path.join(script_dir, fname)
        shutil.copy2(src_path, dst_path)

    print(f"Copied {len(dat_files)} .dat files to: {script_dir}")