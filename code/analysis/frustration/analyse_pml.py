import re
from collections import Counter
import pandas as pd

# ===== Get PML file path =====
pml_file = "tertiary_frustration.pml"

# ===== Read file content =====
try:
    with open(pml_file, "r", encoding="utf-8") as f:
        lines = f.readlines()
except FileNotFoundError:
    print(f"Error: File {pml_file} not found.")
    exit(1)
except UnicodeDecodeError:
    print(f"Error: File {pml_file} encoding incompatible. Please check encoding.")
    exit(1)

# ===== Regular expression to match any atom name and uppercase Chain =====
pattern = re.compile(
    r"draw_links\s+resi\s+(\d+)\s+and\s+name\s+\w+\s+and\s+Chain\s+([A-Z]),\s+resi\s+(\d+)\s+and\s+name\s+\w+\s+and\s+Chain\s+([A-Z]),\s+color=(\w+),\s*color2=\w+"
)

# ===== Initialize counter =====
counter = Counter()

# ===== Iterate over matched lines and count =====
for line in lines:
    match = pattern.search(line)
    if match:
        resi1, chain1, resi2, chain2, color = match.groups()

        # Count A-A, B-B, C-C, A-B, A-C, B-C types
        if chain1 == chain2:  # Same-chain combination
            if chain1 in ["A", "B", "C"]:
                key = (chain1, chain2, color)
                counter[key] += 1
        else:  # Different-chain combination
            chains = tuple(sorted([chain1, chain2]))  # Ensure A-B and B-A are treated as the same
            if set(chains).issubset({"A", "B", "C"}):
                key = (chains[0], chains[1], color)
                counter[key] += 1

# ===== Construct DataFrame =====
data = []
for (chain1, chain2, color), count in counter.items():
    data.append({"Group": f"{chain1}-{chain2}", "Color": color, "Count": count})

df = pd.DataFrame(data)
if not df.empty:
    df = df.sort_values(by=["Group", "Color"])
else:
    print("Warning: No valid data to process.")

# ===== Print and save to TXT file =====
output_txt = "draw_links_summary.txt_

