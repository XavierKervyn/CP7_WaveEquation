# This tool is useful for playing 2D videos
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as Animation

# Data
plot_video = False
font_size = 10
text_color = 'r'
hpos,vpos = 0.45,0.9 # horizontal and vertical positions of texts
title = "Onde type 1, Nx=300, Ny=300, CFL=0.95"
video_name = 'Ex7_7-3-a-onde-type-1-Nx300-Ny300-CFL9dot5en1-T500s'
delimiter = " "
meshfilename = "output_mesh_7-3-a.out"
filename = "output_f_7-3-a.out"

# Read data
with open(meshfilename,'r') as mesh_file:
  x_mesh = np.fromstring(mesh_file.readline(),dtype=np.float64,sep=delimiter)
  y_mesh = np.fromstring(mesh_file.readline(),dtype=np.float64,sep=delimiter)
gloabl_array = np.loadtxt(filename,delimiter=delimiter);
time = np.unique(gloabl_array[:,0])
ntime = time.size
video = []
for t in time:
  print("".join(["Reading time: ",str(t)," of: ",str(time[-1])]))
  frame = gloabl_array[gloabl_array[:,0]==t,1:]
  video.append(frame)  
  
# Plot
im = []
colorbar_min = np.amin(video)
colorbar_max = np.amax(video)
plt.rcParams.update({'font.size': font_size})
fig = plt.figure(1)
for iD,frame in enumerate(video):
  print("".join(["Doing frame: ",str(iD)," of: ",str(ntime)]))
  ax = plt.gca()
  plot = ax.imshow(frame,interpolation=None,extent=[np.amin(x_mesh),np.amax(x_mesh),\
  np.amin(y_mesh),np.amax(y_mesh)],aspect='equal',animated=True)
  time_text = ax.text(hpos,vpos,' t: '+str(time[iD]),\
                horizontalalignment='center',verticalalignment='bottom',\
                transform=ax.transAxes,color=text_color)
  plot.set_clim(colorbar_min, colorbar_max)              
  ax.set_title(title)
  ax.set_xlabel("x")
  ax.set_ylabel("y")
  im.append([plot,time_text])
fig.colorbar(plot)  
print("".join(["Maximum wave high [max(f)]: ",str(colorbar_max)]))
 
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
