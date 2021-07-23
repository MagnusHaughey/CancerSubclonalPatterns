

import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt 
import sys
import argparse



# Parse command line argument
parser = argparse.ArgumentParser()
parser.add_argument('--path', help='File path to cells.csv file to be plotted', type=str)
args = parser.parse_args()


# Import raw data
x , y , mutations = np.loadtxt(args.path , delimiter="," , unpack=True)


# matplotlib.pyplot.pcolor() function requires the input data to be square. So first find if x or y is the larger dimension.
minX = int(np.min(x))
minY = int(np.min(y))
maxX = int(np.max(x))
maxY = int(np.max(y))

largest_dim = np.max([maxX-minX , maxY-minY])






#~~~~~~~~~~~~~~~~~~~~	Plotting    ~~~~~~~~~~~~~~~~~~~~~#


# Initiate input array
toPlot = np.zeros((int(largest_dim)+1 , int(largest_dim)+1))

# Populate array
for i in range(len(x)):
	toPlot[int(x[i])-minX][int(y[i])-minY] = mutations[i] + 1



x = np.linspace(0 , int(largest_dim)+1 , int(largest_dim)+1)
y = np.linspace(0 , int(largest_dim)+1 , int(largest_dim)+1)
X , Y = np.meshgrid(x,y)


# Set up colour map
cmap = mpl.colors.ListedColormap(['white' , (250/255.0, 236/255.0, 45/255.0) , (237/255.0, 86/255.0, 40/255.0)])
norm = mpl.colors.Normalize(vmin = 0, vmax = 2)



plt.figure(figsize=( 8 , 8 ))

#plt.xlim([ 0 , maxX-minX ])
#plt.ylim([ 0 , maxY-minY ])

plt.gca().axis('off')

plt.pcolor(X , Y , toPlot , cmap = cmap , shading = 'auto')



#plt.show()
plt.savefig(args.path + ".png" , dpi = 300 , format = 'png')






