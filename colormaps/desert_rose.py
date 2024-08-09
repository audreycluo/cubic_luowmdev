
import json
from matplotlib import colors as mcolors
from matplotlib.colors import LinearSegmentedColormap


hex_colors = [
   "#8b0e0e", "#a42d39", "#d2644f", "#ee806b", # desert roses
   "#efc0a5", "#f5d9b5", "#d1f5f4", "#c8ede5", # transitional 
   "#a8e4d7",  "#72bdc0", "#009DA8FF", "#0c7390", "#075e77"] # blues
 
 
hex_colors = hex_colors[::-1]

# convert hex color codes to a list of RGB tuples
rgb_colors = [mcolors.hex2color(hex_color) for hex_color in hex_colors]

# Create a color map from these colors
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", rgb_colors, N=256)


# Extract RGB colors from the colormap
rgb_colors = [custom_cmap.colors[i] for i in range(custom_cmap.N)]

# Save to JSON file
with open('/cbica/projects/luo_wm_dev/code/tract_profiles/colormaps/desert_rose.json', 'w') as f:
    json.dump(hex_colors, f)

