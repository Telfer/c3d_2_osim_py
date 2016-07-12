# Convert c3d file to files required for import to OpenSim
# Author: Scott Telfer
# Last Update: 2016_07_11

# NOTES:
# Biomechanical ToolKit needs to be installed in your module library. Get this
# from https://pypi.python.org/pypi/btk#downloads (version used 0.3)

# =============================================================================

# Import required modules
import easygui
import btk
import numpy as np
import csv

# =============================================================================

# Import c3d file
c3d_file = easygui.fileopenbox()
reader = btk.btkAcquisitionFileReader()	
reader.SetFilename(str(c3d_file))	
reader.Update()	
acq = reader.GetOutput()


# =============================================================================

## Extract file name
file_name = str(c3d_file)
file_name = file_name.split("\\")
file_name = file_name[-1][:-4]


# =============================================================================

## Extract data for marker co-ordinates
# Make list of marker names
marker_names = [] # Create empty list to store marker names
for i in btk.Iterate(acq.GetPoints()):
    marker_names.append(i.GetLabel())

# Make array of all coordinates 
marker_coords = [] # Create empty list to store marker coordinates
for i in range(0, len(marker_names)):
    coords = acq.GetPoint(marker_names[i]).GetValues()
    marker_coords.append(coords)
marker_coords = np.concatenate(marker_coords, axis = 1)

# Add Frame# column
frame_nos = range(1, (len(marker_coords) + 1), 1)

# Add Time column
time_frame = len(marker_coords) / acq.GetPointFrequency()
time_frames = np.arange(0, time_frame, 1 / acq.GetPointFrequency())

# Add frame_nos and time_frames to array
marker_coords = np.c_[frame_nos, time_frames, marker_coords]


# =============================================================================

## Additional details for .trc file
# Marker collection data
NumMarkers = len(marker_names)
DataRate = acq.GetPointFrequency()
CameraRate = DataRate
NumFrames = len(marker_coords)
Units = 'mm'
OrigDataRate = DataRate
OrigDataStartFrame = 1
OrigNumFrames = NumFrames

# Header row 1
output_filename = file_name + ".trc"
trc_1 = ['PathFileType', '4', '(X/Y/Z)'] + [output_filename]

# Header row 2
trc_2 = ['DataRate', 'CameraRate', 'NumFrames', 'NumMarkers', 
         'Units', 'OrigDataRate', 'OrigDataStartFrame', 'OrigNumFrames']

#Header row 3
trc_3 = [DataRate, CameraRate, NumFrames, NumMarkers, Units, OrigDataRate, 
         OrigDataStartFrame, OrigNumFrames]

# Header row 4
marker_names2 = []
for e in marker_names:
    marker_names2.append(e)
    marker_names2.extend(['', ''])
trc_4 = ['Frame#', 'Time'] + marker_names2

#Header row 5
xyz_strs = ['X', 'Y', 'Z'] * NumMarkers
xyz_nos = ['1', '2', '3'] * NumMarkers
xyz = [x + y for x, y in zip(xyz_strs, xyz_nos)]
trc_5 = ['', ''] + xyz

# Header row 6
trc_6 = ['', ''] + (['', '', ''] * NumMarkers)


# =============================================================================

## Extract force data


# =============================================================================

## Extract EMG data


# =============================================================================

## Write .trc file (marker coordinate data)
with open('C:/Users/telfe/Dropbox/My_Projects/VA_Dose_Response/OpenSim/Training/c3d examples/mydata2.trc', 'wb') as mycsvfile:
    thedatawriter = csv.writer(mycsvfile, delimiter = '\t')
    thedatawriter.writerow(trc_1)
    thedatawriter.writerow(trc_2)
    thedatawriter.writerow(trc_3)
    thedatawriter.writerow(trc_4)
    thedatawriter.writerow(trc_5)
    thedatawriter.writerow(trc_6)    
    for row in marker_coords:
        thedatawriter.writerow(row)


# =============================================================================
