import re
from collections import Counter
import pandas as pd

# ===== 获取 PML 文件路径 =====
pml_file = "tertiary_frustration.pml"

# ===== 读取文件内容 =====
try:
    with open(pml_file, "r", encoding="utf-8") as f:
        lines = f.readlines()
except FileNotFoundError:
    print(f"错误：文件 {pml_file} 未找到。")
    exit(1)
except UnicodeDecodeError:
    print(f"错误：文件 {pml_file} 编码不兼容，请检查编码。")
    exit(1)

# ===== 正则表达式，匹配任意原子名、大小写 Chain =====
pattern = re.compile(
    r"draw_links\s+resi\s+(\d+)\s+and\s+name\s+\w+\s+and\s+Chain\s+([A-Z]),\s+resi\s+(\d+)\s+and\s+name\s+\w+\s+and\s+Chain\s+([A-Z]),\s+color=(\w+),\s*color2=\w+"
)

# ===== 初始化计数器 =====
counter = Counter()

# ===== 遍历匹配行并计数 =====
for line in lines:
    match = pattern.search(line)
    if match:
        resi1, chain1, resi2, chain2, color = match.groups()

        # 统计 A-A, B-B, C-C, A-B, A-C, B-C 类型
        if chain1 == chain2:  # 同链组合
            if chain1 in ["A", "B", "C"]:
                key = (chain1, chain2, color)
                counter[key] += 1
        else:  # 异链组合
            chains = tuple(sorted([chain1, chain2]))  # 确保 A-B 和 B-A 视为相同
            if set(chains).issubset({"A", "B", "C"}):
                key = (chains[0], chains[1], color)
                counter[key] += 1

# ===== 构造 DataFrame =====
data = []
for (chain1, chain2, color), count in counter.items():
    data.append({"Group": f"{chain1}-{chain2}", "Color": color, "Count": count})

df = pd.DataFrame(data)
if not df.empty:
    df = df.sort_values(by=["Group", "Color"])
else:
    print("警告：无有效数据可处理。")

# ===== 打印并保存到 TXT 文件 =====
output_txt = "draw_links_summary.txt"
try:
    with open(output_txt, "w", encoding="utf-8") as out:
        out.write("# Draw Links Summary\n")
        out.write("# Group\tColor\tCount\n")
        out.write(df.to_string(index=False) if not df.empty else "无数据")
    print(f"结果已保存到：{output_txt}")
except IOError as e:
    print(f"错误：保存文件 {output_txt} 失败 - {e}")