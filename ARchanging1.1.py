#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
import sys
from pathlib import Path
from typing import List, Tuple

def parse_numbers(s: str) -> List[float]:
    """Parse all floating point numbers from a string."""
    return [float(x) for x in re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", s)]

def format_like(orig: str, values: List[float]) -> str:
    """Keep original spacing/column style where possible.
    Fallback to space-separated with 6 decimals."""
    parts = re.split(r"([-+]?\d*\.\d+|[-+]?\d+)", orig.rstrip("\n"))
    nums = [p for p in parts if re.fullmatch(r"[-+]?\d*\.\d+|[-+]?\d+", p)]
    if len(nums) == len(values):
        vi = 0
        for i, p in enumerate(parts):
            if re.fullmatch(r"[-+]?\d*\.\d+|[-+]?\d+", p):
                # Preserve original style (e.g., column width) but force 6 decimals for precision
                parts[i] = f"{values[vi]:.6f}"
                vi += 1
        return "".join(parts) + "\n"
    return ("  ".join(f"{v:.6f}" for v in values) + "\n")

def read_geo_blocks(path: Path):
    """
    Reads wing stations from the .GEO file, handling header and tail blocks.
    The tail is defined as the content starting from the first line with only one number.
    """
    lines = path.read_text(errors="ignore").splitlines(True)
    i = 0
    stations = []
    head_before_first_hdr = []
    
    # --- State machine for parsing ---
    while i < len(lines):
        line = lines[i]
        nums = parse_numbers(line)
        
        # 1. Check for Tail section (starts with a line containing only one number)
        if len(nums) == 1 and stations:
            # Once a station is found, a single number line marks the start of the tail
            tail = lines[i:]
            return stations, tail, head_before_first_hdr
        
        # 2. Check for Standard Station Header (4 numbers: y, zoff, chord, xle)
        if len(nums) == 4:
            y, zoff, chord, xle = nums
            hdr_line = line
            i += 1
            if i >= len(lines): break
            
            cnt_nums = parse_numbers(lines[i])
            if len(cnt_nums) < 3:
                raise RuntimeError(f"Unexpected counts at line {i+1}: {lines[i]}")
            
            cnt_line = lines[i]
            # Assumes nu and nl are at index 1 and 2
            nu, nl = int(cnt_nums[1]), int(cnt_nums[2]) 
            coords = lines[i+1 : i+1+nu+nl]
            
            stations.append(dict(
                y=y, zoff=zoff, chord=chord, xle=xle,
                hdr_line=hdr_line, cnt_line=cnt_line, coords=coords
            ))
            i += 1 + nu + nl
            
        # 3. Check for Default/Root Station (5 numbers: 0, nu, nl, ..., ...)
        elif len(nums) >= 3 and int(nums[0]) == 0:
            # Treat this as the root station definition
            cnt_line = line
            nu, nl = int(nums[1]), int(nums[2])
            coords = lines[i+1 : i+1+nu+nl]
            
            # The GEO format often omits the 4-number header for the root,
            # using default values of 0/1 for y/zoff/chord/xle.
            stations.append(dict(
                y=0.0, zoff=0.0, chord=1.0, xle=0.0,
                hdr_line=None, cnt_line=cnt_line, coords=coords
            ))
            i += 1 + nu + nl
            
        # 4. Content before first station is considered header
        else:
            if not stations:
                head_before_first_hdr.append(line)
            i += 1
            
    # If no tail was explicitly found (e.g., file ends after stations)
    return stations, [], head_before_first_hdr

def compute_span_area(stations: List[dict]) -> Tuple[float, float]:
    """Computes half-span (b) and total area (S) using the trapezoidal rule."""
    ys = [s["y"] for s in stations]
    chords = [s["chord"] for s in stations]
    
    # Calculate half-wing area (S_half) using trapezoidal rule
    S_half = 0.0
    for i in range(len(stations)-1):
        dy = ys[i+1] - ys[i]
        # Sum of chords * dy for each panel
        S_half += 0.5 * (chords[i] + chords[i+1]) * dy
        
    b = 2.0 * max(ys)
    S = 2.0 * S_half
    return b, S

def scale_stations(stations: List[dict], ky: float, kc: float, sweep_mode: str, keep_angles: bool=False):
    """Scales the station geometry based on ky and kc."""
    scaled = []
    for s in stations:
        y = s["y"] * ky
        chord = s["chord"] * kc
        
        # FIX: To preserve the sweep angle, Xle must scale with Y (ky), not chord (kc).
        if sweep_mode == "scale":
            xle = s["xle"] * ky
        else:
            xle = s["xle"] # Fixed (sweep_mode == "fixed")

        # Zoff scaling: Scales with ky to keep dihedral angle (if keep_angles=True), otherwise with kc.
        zoff = s["zoff"] * (ky if keep_angles else kc)
        
        # Re-format header line if it exists
        if s["hdr_line"] is not None:
            # Format using y, zoff, chord, xle (note: order matters for format_like)
            new_hdr = format_like(s["hdr_line"], [y, zoff, chord, xle])
        else:
            new_hdr = None
            
        scaled.append(dict(
            y=y, zoff=zoff, chord=chord, xle=xle,
            hdr_line=new_hdr if new_hdr else s["hdr_line"],
            cnt_line=s["cnt_line"], coords=s["coords"]
        ))
    return scaled

def write_geo(path_out: Path, head_before_first_hdr, stations: List[dict], tail: List[str]):
    """Writes the scaled GEO file."""
    out_lines = []
    out_lines.extend(head_before_first_hdr)
    for s in stations:
        if s["hdr_line"] is not None:
            out_lines.append(s["hdr_line"])
        out_lines.append(s["cnt_line"])
        out_lines.extend(s["coords"])
    out_lines.extend(tail)
    path_out.write_text("".join(out_lines), encoding="utf-8")
    return path_out

def update_map(path_in: Path, path_out: Path, ky: float, kc: float):
    """Scales Y/X/Z coordinates in the .map file based on nearby labels."""
    lines = path_in.read_text(errors="ignore").splitlines(True)
    out = []
    # Basic state to know whether previous line had XF1 labels etc.
    prev = ""
    for ln in lines:
        stripped = ln.strip()
        nums = parse_numbers(stripped)
        
        # Y1, Y2, X3, X4, Z3, Y3 line
        if len(nums) == 6 and ("Y1" in prev and "Y2" in prev and "X3" in prev):
            Y1, Y2, X3, X4, Z3, Y3 = nums
            Y1 *= ky; Y2 *= ky; Y3 *= ky # Y coords (span direction) scale with ky
            X3 *= kc; X4 *= kc; Z3 *= kc # X/Z coords (chord/thickness direction) scale with kc
            out.append(format_like(ln, [Y1, Y2, X3, X4, Z3, Y3]))
            
        # XF1, XF2, YF1, YF2, XF3, YCL line
        elif len(nums) == 6 and ("XF1" in prev and "XF2" in prev):
            XF1, XF2, YF1, YF2, XF3, YCL = nums
            XF1 *= kc; XF2 *= kc; XF3 *= kc # X coords scale with kc
            YF1 *= ky; YF2 *= ky; YCL *= ky # Y coords scale with ky
            out.append(format_like(ln, [XF1, XF2, YF1, YF2, XF3, YCL]))
            
        else:
            out.append(ln)
            
        prev = ln
    path_out.write_text("".join(out), encoding="utf-8")
    return path_out

def main():
    if len(sys.argv) < 5:
        print("""python wing_ar_param.py <input.GEO> <input.map> <target_AR> <out_prefix> [--preserve-area|--preserve-span|--preserve-chord] [--sweep scale|fixed]
Defaults: --preserve-area --sweep scale

Meaning:
  --preserve-area: keep wing area S constant (recommended). Sets ky = sqrt(AR_target/AR0), kc = 1/ky.
  --preserve-span: keep span b constant. ky = 1, kc = AR0/AR_target.
  --preserve-chord: keep chord scale. kc = 1, ky = AR_target/AR0.

  --sweep scale: scale x_le with ky to preserve sweep angle (non-dimensional sweep). (FIXED: now scales with ky)
  --sweep fixed: keep x_le unchanged in absolute coordinates.

  --keep-angles: (Optional) If provided, Zoff scales with ky to preserve dihedral/incidence angles. Default is to scale Zoff with kc.
""")
        sys.exit(1)

    geo_in = Path(sys.argv[1])
    map_in = Path(sys.argv[2])
    AR_target = float(sys.argv[3])
    out_prefix = sys.argv[4]
    
    # --- Parse Optional Arguments ---
    preserve = "area"
    sweep_mode = "scale"
    keep_angles = False
    for arg in sys.argv[5:]:
        if arg == "--preserve-area": preserve = "area"
        elif arg == "--preserve-span": preserve = "span"
        elif arg == "--preserve-chord": preserve = "chord"
        elif arg == "--sweep": pass # Handled by the next argument
        elif arg in ("scale", "fixed"): sweep_mode = arg
        elif arg == "--keep-angles": keep_angles = True # Option for zoff scaling

    # --- 1. Read Data ---
    try:
        stations, tail, head = read_geo_blocks(geo_in)
    except RuntimeError as e:
        raise SystemExit(f"Error reading GEO file: {e}")
        
    if not stations:
        raise SystemExit("No wing stations parsed from GEO.")

    # --- 2. Calculate Scaling Factors ---
    b0, S0 = compute_span_area(stations)
    AR0 = (b0*b0) / S0

    if preserve == "area":
        # S1 = S0 => ky = sqrt(AR1/AR0), kc = 1/ky
        ky = (AR_target/AR0) ** 0.5
        kc = 1.0 / ky
    elif preserve == "span":
        # b1 = b0 => ky = 1, kc = AR0/AR1
        ky = 1.0
        kc = AR0 / AR_target
    else: # preserve == "chord"
        # c_avg1 = c_avg0 => kc = 1, ky = AR1/AR0
        kc = 1.0
        ky = AR_target / AR0

    # --- 3. Scale and Write Files ---
    stations_scaled = scale_stations(stations, ky, kc, sweep_mode=sweep_mode, keep_angles=keep_angles)

    geo_out = Path(out_prefix + ".GEO")
    write_geo(geo_out, head, stations_scaled, tail)

    map_out = None
    if map_in.exists():
        map_out = Path(out_prefix + ".map")
        update_map(map_in, map_out, ky=ky, kc=kc)

    # --- 4. Report Results ---
    b1, S1 = compute_span_area(stations_scaled)
    AR1 = (b1*b1)/S1
    print(f"Original: b={b0:.6f}, S={S0:.6f}, AR={AR0:.6f}")
    print(f"Target:   AR={AR_target:.6f} -> Applied ky={ky:.6f}, kc={kc:.6f} (sweep={sweep_mode})")
    print(f"New:      b={b1:.6f}, S={S1:.6f}, AR={AR1:.6f} (Mismatch due to trapezoidal rule/precision if AR1 != AR_target)")
    if map_out:
        print(f"Wrote: {geo_out} and {map_out}")
    else:
        print(f"Wrote: {geo_out}")

if __name__ == "__main__":
    main()