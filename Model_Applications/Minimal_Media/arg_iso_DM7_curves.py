########### Generate Growth curves from a large dataset exploring arginine and isoleucine inclusion in the LGG minimal medium ###########

import matplotlib.pyplot as plt
import statistics
from matplotlib.lines import Line2D
import pandas as pd

excel_path = "arg_iso_cleaned.xlsx"

df = pd.read_excel(excel_path, sheet_name=sheet_name)


df.rename(columns={df.columns[0]: "Duration"}, inplace=True)

dm7_cols = [c for c in df.columns if "DM7" in c and "iso" not in c and "arg" not in c]
iso_cols = [c for c in df.columns if "iso" in c]
arg_cols = [c for c in df.columns if "arg" in c]

raw_data = {
    "DM7": {row["Duration"]: [row[c] for c in dm7_cols if pd.notna(row[c])] for _, row in df.iterrows()},
    "DM7_iso": {row["Duration"]: [row[c] for c in iso_cols if pd.notna(row[c])] for _, row in df.iterrows()},
    "DM7_arg": {row["Duration"]: [row[c] for c in arg_cols if pd.notna(row[c])] for _, row in df.iterrows()},
}
processed_data={}
for media, data in raw_data.items():
    time_df={}
    for time, replicates in data.items():
        average= sum(replicates)/len(replicates)
        st_dev = statistics.stdev(replicates)
        time_df[time]=average, st_dev
    processed_data[media]=time_df

colors = {'DM7': '#FF0000', 'DM7_iso': '#0000FF', 'DM7_arg': '#FFA500'}
plt.figure(figsize=(10,6))

legend_handles = []

for media, points in processed_data.items():
    if not isinstance(points, dict):
        print(f"Warning: {media} points is not a dictionary, skipping")
        continue

    times = sorted(points.keys())
    means = []
    sds = []

    for t in times:
        val = points[t]
        if isinstance(val, (list, tuple)) and len(val) == 2:
            mean, sd = val
        else:
            raise ValueError(f"Unexpected format at {media} time {t}: {val}")
        means.append(mean)
        sds.append(abs(sd))

    means = [float(m) for m in means]
    sds = [float(s) for s in sds]

    line, = plt.plot(times, means, color=colors.get(media, 'black'), label=media)
    legend_handles.append(line)

    lower = [m - s for m, s in zip(means, sds)]
    upper = [m + s for m, s in zip(means, sds)]
    plt.fill_between(times, lower, upper, color=colors.get(media, 'black'), alpha=0.2)


plt.xlabel('Duration (hr)', fontsize=18)
plt.ylabel('OD600', fontsize=18)
plt.xticks(fontsize=18)  
plt.yticks(fontsize=18)  
plt.legend(handles=legend_handles, fontsize=20)
plt.tight_layout()
plt.savefig('arg_iso_DM7_curve.svg', format='svg') 
