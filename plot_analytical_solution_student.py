from matplotlib import pyplot as plt
import numpy as np
import AnalyticalSolutions as ansol
import matplotlib.animation as Animation
import os

# Data
plot_video = False
font_size = 14
text_color = 'r'
hpos,vpos = 0.45,0.9 # horizontal and vertical positions of texts
title = "Initial wave function - m=10 , n=10"
video_name = 'Ex7_7-2-b-m10_n10-analytical-solution-2periods'
Nx = 317
Ny = 1000
xL = 0.
xR = 10.
yL = 0.
yU = 10.
kx = 10
ky = 10
F0 = 1.
u = 4.
time = 2*0.353553390593274
Nt = 100

# initialise
x_grid = np.linspace(xL,xR,Nx)
y_grid = np.linspace(yL,yU,Ny)
time_grid = np.linspace(0.,time,Nt)
solution_time = []

# compute solution
for t in time_grid:
  print("".join(["doing time: ",str(t)]))
  solution = np.zeros((Ny,Nx))
  for iDy,y in enumerate(y_grid):
    for iDx,x in enumerate(x_grid):
      solution[iDy][iDx] = 0.
    
# plot
im = []
colorbar_min = np.amin(solution_time)
colorbar_max = np.amax(solution_time)
plt.rcParams.update({'font.size': font_size})
fig = plt.figure(1)
for iD,frame in enumerate(solution_time):
  ax = plt.gca()
  plot = ax.imshow(frame,interpolation=None,extent=[xL,xR,yL,yU],animated=True)
  time_text = ax.text(hpos,vpos,' t: '+str(time_grid[iD]),\
                horizontalalignment='center',verticalalignment='bottom',\
                transform=ax.transAxes,color=text_color)
  plot.set_clim(colorbar_min, colorbar_max)              
  ax.set_title(title)
  ax.set_xlabel("x")
  ax.set_ylabel("y")
  im.append([plot,time_text])

# generate video
animation = Animation.ArtistAnimation(fig,im,repeat=False,blit=True) 

# save video in mp4
metadata = dict(title=title,artisti='EPFLPhysNum') # initialise movie metadata
Writer = Animation.writers['ffmpeg'] # select the kind of movie writer
writer = Writer(fps=10, metadata=metadata) # initialising the writer
video_path = ''.join([video_name,'.mp4']) # path to the video
if os.path.isfile(video_path):
  os.remove(video_path) #if it does exist, remove it
animation.save(video_path,writer=writer) # save the new video
  
# Play video
if(plot_video):
  plt.draw()  
  plt.show()
