import numpy as np
import matplotlib.pyplot as plt

from matplotlib.animation import FuncAnimation # matplotlib's movie maker

# if(len(sys.argv) != 2):
	# sys.exit("Error! Usage: ./plotArrivalTime fname")

# fname = sys.argv[1]

# dt = np.dtype([('rLenFront', np.uint32, (1)),
               # ('energy', np.float64, (1)),
               # ('theta', np.float64, (1)),
               # ('time', np.float64, (1)),
               # ('eperp', np.float64, (1)),
               # ('x', np.float64, (1)),
               # ('y', np.float64, (1)),
               # ('z', np.float64, (1)),
               # ('zoff', np.float64, (1)),
               # ('nhit', np.int32, (1)),
               # ('nhitBot', np.int32, (1)),
               # ('nhitTop', np.int32, (1)),
               # ('rLenBack', np.uint32, (1))])

# data = np.fromfile(fname, dtype=dt)

X = np.arange(-1.5,1.5,0.1)
Y = np.arange(-0.5,1.0,0.1)
X, Y = np.meshgrid(X, Y)
Z = zeros(np.size(X)) # Copy sizing
for x in X:
	if x > 0:
		Z = 1





# Stealing from the internet! (Josh Borrow @ Durham)
# Some global variables to define the whole run
# total_number_of_frames = 100
# all_data = [np.random.rand(512, 512) for x in range(100)]

# def animate(frame):
	# """
	# Animation function. Takes the current frame number (to select the 
	# potion of data to plot) and a line object to update. 
	# """
	
	# # Not strictly necessary, just so we know we are stealing these from
	# # the global scope
	# global all_data, image
	
	# # We want up-to and _including_ the frame'th element
	# image.set_array(all_data[frame])
	
	# return image
	
# # Now we can do the plotting!
# fig, ax = plt.subplots(1, figsize=(1,1))
# # Remove a bunch of stuff to make sure we only 'see' the actual imshow
# # Stretch to fit the whole plane
# fig.subplots_adjust(0, 0, 1, 1)
# # Remove bounding line
# ax.axis("off")

# # Initialise our plot. Make sure you set vmin and vmax!
# image = ax.imshow(all_data[0], vmin=0, vmax=1)

# animation = FuncAnimation( 
	# # Your Matplotlib Figure object
	# fig,
	# # The function that does the updating of the figure
	# animate,
	# # Frame information (here just frame number)
	# np.arange(total_number_of_frames),
	# # Extra arguments to the animate function
	# fargs=[],
	# # Frame-time in ms; i.e. for a given frame-rate x, 1000/x
	# interval=1000 / 25
# )

# # Try to set the DPI to the actual number of pixels you're plotting
# animation.save("out_2dgrid.mp4", dpi=512)

	
	
