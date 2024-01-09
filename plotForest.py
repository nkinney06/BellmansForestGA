import pandas as pd
import matplotlib.pyplot as plt
import math
import matplotlib.gridspec as gridspec

# Your dataframe
df = pd.read_csv("forestData.txt",header=None,sep="\t")
df.columns =['x', 'y', 'type']

# Separate data based on type
circle_data = df[df['type'] == 'circle'].reset_index()
system_info = df[df['type'] == 'system'].reset_index()
forest_data = df[df['type'] == 'forest'].reset_index()
points_data = df[df['type'] == 'points'].reset_index()
strats_data = df[df['type'] == 'strats'].reset_index()
scores_info = df[df['type'] == 'generation'].reset_index()
solved_data = df[df['type'] == 'solved'].reset_index()
degrees_of_freedom = int(system_info['y'][0])

# Create subplots with 1 row and 3 columns
fig, axes = plt.subplots(2, 3, figsize=(18, 9))

# Plot polygon
axes[0,0].plot(forest_data['x'].tolist() + [forest_data['x'].iloc[0]],
             forest_data['y'].tolist() + [forest_data['y'].iloc[0]],
             color='blue', label='forest', marker='o')
axes[0,0].scatter(points_data['x'], points_data['y'], color='green', label='points')
axes[0,0].set_title('forest and {} Grid points'.format(len(points_data)))
axes[0,0].set_xlabel('')
axes[0,0].set_ylabel('')
axes[0,0].grid(True)
axes[0,0].set_aspect('equal',adjustable='datalim')

# Plot path
axes[0,1].plot(strats_data['x'], strats_data['y'], color='red', label='strats', marker='o')
axes[0,1].set_title('best escape strategy with {} segments'.format(degrees_of_freedom))
axes[0,1].set_xlabel('')
axes[0,1].set_ylabel('')
axes[0,1].grid(True)
axes[0,1].set_aspect('equal',adjustable='datalim')

# plot circle data
axes[0,2].scatter(circle_data['x'], circle_data['y'], color='green', label='circle')
axes[0,2].set_title('{} test angles'.format(len(circle_data)/(degrees_of_freedom+1)))
axes[0,2].set_xlabel('')
axes[0,2].set_ylabel('')
axes[0,2].plot(circle_data['x'][0:degrees_of_freedom+1], circle_data['y'][0:degrees_of_freedom+1], color='red', label='Solution', marker='o')
axes[0,2].plot(circle_data['x'][degrees_of_freedom+1:(degrees_of_freedom+1)*2], circle_data['y'][degrees_of_freedom+1:(degrees_of_freedom+1)*2], color='blue', label='Solution', marker='o')
axes[0,2].grid(True)
axes[0,2].set_aspect('equal',adjustable='datalim')

# plot fitness
axes[1,0].plot(scores_info['x'], scores_info['y'], color='black', label='fitness', marker='o')
axes[1,0].set_title('fitness minimization')
axes[1,0].set_xlabel('generation')
axes[1,0].set_ylabel('escape distance')
axes[1,0].grid(True)

# Plot solution
axes[1,1].plot(solved_data['x'], solved_data['y'], color='red', label='solved', marker='o')
axes[1,1].plot(forest_data['x'].tolist() + [forest_data['x'].iloc[0]],
             forest_data['y'].tolist() + [forest_data['y'].iloc[0]],
             color='blue', label='Polygon', marker='o')
axes[1,1].set_title('escape distance {}'.format(system_info['x'][0]))
axes[1,1].grid(True)
axes[1,1].set_aspect('equal',adjustable='datalim')

# plot the final path
axes[1,2].plot(solved_data['x'], solved_data['y'], color='red', label='solved', marker='o')
axes[1,2].set_title('escape distance {}'.format(system_info['x'][0]))
axes[1,2].grid(True)
axes[1,2].set_aspect('equal',adjustable='datalim')

# Save the plot as PNG
plt.savefig('example_plot.png')
plt.close()
exit()