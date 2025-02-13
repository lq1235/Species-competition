import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Define data
data = np.array([
    [0.3333,0.1000,0.2000],           
    [0.2,0.3333,0.1000],
    [0.1,0.2000,0.3333]
])

# Create a heat map
plt.figure(figsize=(8, 6))  # Resize the image to better show larger fonts
heatmap = sns.heatmap(data, annot=True, cmap='viridis', cbar=True, linewidths=2, annot_kws={"size": 27})

# Resize the image to better show larger words Add labels and titles
plt.title('Resource Quantity Demanded', fontsize=22)
plt.xlabel('Resource', fontsize=27)
plt.ylabel('Species', fontsize=27)

# Resize the image to better show larger words Set the scale labels on the x and y axes to start at 1 and make the font larger
heatmap.set_xticklabels([1, 2, 3], fontsize=27)
heatmap.set_yticklabels([1, 2, 3], fontsize=27)

# Gets the color bar object and resizes the label font
cbar = heatmap.collections[0].colorbar
cbar.ax.tick_params(labelsize=27)

# Adjust the font size of the scale label
heatmap.tick_params(axis='both', which='major', labelsize=27)

# Save the image in PDF format
plt.savefig('ResourceQuantityDemanded.pdf')

# Display image
plt.show()
