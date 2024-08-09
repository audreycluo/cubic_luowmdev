import json
from matplotlib import colors as mcolors
from matplotlib.colors import LinearSegmentedColormap



hex_colors = [ "#052632", "#084559","#045F79"  "#017C9C", 
              "#1C9CB0", "#3CB0B7", "#5AC1BB", "#77CFBE", 
              "#92DBC1", "#ABE4C4", "#C1EBC7", "#D4EFCC", 
              "#E7EA9D", "#F5F1A6", "#FBFCC6"]


# convert hex color codes to a list of RGB tuples
rgb_colors = [mcolors.hex2color(hex_color) for hex_color in hex_colors]

# Create a color map from these colors
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", rgb_colors, N=256)
 
# Save to JSON file
with open('/cbica/projects/luo_wm_dev/code/tract_profiles/colormaps/aquamarine.json', 'w') as f:
    json.dump(hex_colors, f)

