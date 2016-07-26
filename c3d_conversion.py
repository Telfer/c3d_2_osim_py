# Convert c3d file to files required for import to OpenSim
# Author: Scott Telfer
# Last Update: 2016_07_21

# NOTES:
# Biomechanical ToolKit needs to be installed in your module library. Get this
# from https://pypi.python.org/pypi/btk#downloads (version used: 0.3)


# =============================================================================

def c3d_to_opensim(c3d_filepath, trc_filepath, mot_filepath, threshold):

    # =========================================================================
    
    # Import required modules
    import btk
    import numpy as np
    import csv
    import math
    
    # =========================================================================    
    
    ## Helper functions
    def format_list(list_of_strings): 
        d = map(str, list_of_strings)
        for x in xrange(0, len(d)):
            length = len(d[x])
            if length > 10:
                d[x] = d[x][: - (length - 10)]
            if length < 10 and '.' in d[x]:
                zeropad = '0' * (10 - length)        
                d[x] = d[x] + zeropad
            if length < 10 and '.' not in d[x]:
                zeropad = '.' + '0' * (9 - length)
                d[x] = d[x] + zeropad
        return d
    
    
    # =========================================================================
    
    # Import c3d file    
    reader = btk.btkAcquisitionFileReader()	
    reader.SetFilename(str(c3d_filepath))	
    reader.Update()	
    acq = reader.GetOutput()
    
    
    # =========================================================================
    
    ## Extract file name
    file_name = str(c3d_filepath)
    file_name = file_name.split("\\")
    file_name = file_name[-1][:-4]
    
    
    # =========================================================================
    
    ## Extract data for marker co-ordinates
    # Make list of marker names
    marker_names = [] # Create empty list to store marker names
    for i in btk.Iterate(acq.GetPoints()):
        marker_names.append(i.GetLabel())
    
    # Make array of all coordinates 
    marker_coords = [] # Create empty list to store marker coordinates
    axis_conversion = True
    for i in range(0, len(marker_names)):
        coords = acq.GetPoint(marker_names[i]).GetValues()
        if axis_conversion == True:
            coords_x = coords[:, 1]     
            coords_y = coords[:, 2]
            coords_z = coords[:, 0]
            coords_x = format_list(coords_x)
            coords_y = format_list(coords_y)
            coords_z = format_list(coords_z)
            coords_str = np.column_stack([coords_x, coords_y, coords_z])
        marker_coords.append(coords_str)
    marker_coords = np.concatenate(marker_coords, axis = 1)
    
    # Add Frame# column
    frame_nos = range(1, (len(marker_coords) + 1), 1)
    frame_nos = map(str, frame_nos)
    
    # Add Time column
    time_frame = len(marker_coords) / acq.GetPointFrequency()
    time_frames = np.arange(0, time_frame, 1 / acq.GetPointFrequency())
    time_frames = time_frames[0:len(marker_coords)]
    time_frames = format_list(time_frames)
    
    # Add frame_nos and time_frames to array
    marker_coords = np.c_[frame_nos, time_frames, marker_coords]
    
    
    # =========================================================================
    
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
             'Units', 'OrigDataRate', 'OrigDataStartFrame', 
             'OrigNumFrames']
    
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
    xyz_nos = map(str, np.repeat(list(range(1, NumMarkers + 1)), 3))
    xyz = [x + y for x, y in zip(xyz_strs, xyz_nos)]
    trc_5 = ['', ''] + xyz
    
    # Header row 6
    trc_6 = ['', ''] + (['', '', ''] * NumMarkers)
    
    
    # =========================================================================
    
    ## Extract force data
    pfe = btk.btkForcePlatformsExtractor()
    pfe.SetInput(acq)
    pfe.Update()
    pfc = pfe.GetOutput()
    fp_present = pfc.IsEmpty()
    fp_data = []    
    
    if fp_present == False: 
        fp_types = []
        
        for i in btk.Iterate(pfe.GetOutput()):
            fp_type = i.GetType()
            fp_types.append(fp_type)
        
        # Find number of force plates in collection
        NumForcePlates = len(fp_types)
        
        # Find 
        analog_freq = acq.GetAnalogFrequency()
        no_analog_frames = acq.GetAnalogFrameNumber()
        time_frame_a = no_analog_frames / analog_freq
        time_frames_a = np.arange(0, time_frame_a, 1 / analog_freq)
        time_frames_a = format_list(time_frames_a)
        
        # Loop across all force plates
        fp_data1 = []        
        for j in range(0, NumForcePlates):
            fp = pfc.GetItem(j)
            
            # Extract and calcualte required variables from plate data 
            # (translate to opensim coordinate frame)
            if fp_types[j] == 2:
                # Plate thickness        
                dz = fp.GetOrigin()
                dz = (dz[2] / 1000) * -1
                
                # Create rotational matrix from corners
                fp_cors = fp.GetCorners()
                cor_43 = fp_cors[:, 3] - fp_cors[:, 2]
                cor_43_mag = cor_43[0]**2 + cor_43[1]**2 + cor_43[2]**2
                cor_43_mag = math.sqrt(cor_43_mag)
                x_vec = cor_43 / cor_43_mag
                cor_23 = fp_cors[:, 1] - fp_cors[:, 2]
                cor_23_mag = cor_23[0]**2 + cor_23[1]**2 + cor_23[2]**2
                cor_23_mag = math.sqrt(cor_23_mag)
                y_vec = cor_23 / cor_23_mag
                z_vec = [(x_vec[1] * y_vec[2]) - (x_vec[2] * y_vec[1]),
                         (x_vec[2] * y_vec[0]) - (x_vec[0] * y_vec[2]),
                         (x_vec[0] * y_vec[1]) - (x_vec[1] * y_vec[0])]  
                rot_mat = np.column_stack([x_vec, y_vec, z_vec])
                
                # Find force plate center in lab coordinate system
                fp_x_cen = ((fp_cors[0, 0] + fp_cors[0, 2]) / 2) / 1000
                fp_y_cen = ((fp_cors[1, 0] + fp_cors[1, 3]) / 2) / 1000

                # Force and torque vectors
                channels = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
                for k in range(0, 6):
                    channels[k] = fp.GetChannel(k)
                    channels[k] = channels[k].GetValues()
                    channels[k] = channels[k][:, 0]
                Fx_o = channels[0] * -1
                Fy_o = channels[1]
                Fz_o = channels[2] * -1
                Mx_o = channels[3] / 1000
                My_o = channels[4] / 1000
                Mz_o = channels[5] / 1000

                # threshold force and torque data
                thres_indexes = [i for i, v in enumerate(Fz_o) if v < 20.0]
                for h in thres_indexes:
                    Fx_o[h] = 0
                    Fy_o[h] = 0
                    Fz_o[h] = 0
                    Mx_o[h] = 0
                    My_o[h] = 0
                    Mz_o[h] = 0

                # Calculate COP coordinates
                COPx1 = (((Fx_o * dz) - My_o) / Fz_o) * -1
                COPx1 = np.nan_to_num(COPx1)
                COPy1 = (((Fy_o * dz) + Mx_o) / Fz_o) * -1
                COPy1 = np.nan_to_num(COPy1)
                COPz1 = [0.0] * len(COPy1)
                COP_mat = np.row_stack([COPx1, COPy1, COPz1])
                COP_mat = np.dot(rot_mat, COP_mat)
                COPx1 = COP_mat[0, ]
                COPy1 = COP_mat[1, ]
                COPz1 = COP_mat[2, ]
                
                # Apply offset to lab coordinates
                COPx1 = COPx1 + fp_x_cen
                COPy1 = COPy1 + fp_y_cen

                # Make array to be output (convert to OpenSim Ref frame here)
                Fx = format_list(Fy_o)
                Fy = format_list(Fz_o)
                Fz = format_list(Fx_o)
                Mx = format_list(My_o)
                My = format_list(Mz_o)
                Mz = format_list(Mx_o)
                COPx = format_list(COPy1)
                COPy = ['0.00000000'] * len(Fz)
                COPz = format_list(COPx1)

                forceplate_data = np.column_stack([Fx, Fy, Fz, COPx, COPy, 
                                                   COPz, Mx, My, Mz])
                fp_data1.append(forceplate_data)

            # Format and add time column
            #fp_data = np.column_stack(fp_data1)
            #fp_data = np.column_stack([time_frames_a, fp_data])

            if fp_types[j] == 3:
                # Origin coords
                d = fp.GetOrigin()
                dx = d[0] / 1000
                dy = d[1] / 1000
                dz = d[2] / 1000

                # Create rotational matrix from corners
                fp_cors = fp.GetCorners()
                cor_43 = fp_cors[:, 3] - fp_cors[:, 2]
                cor_43_mag = cor_43[0]**2 + cor_43[1]**2 + cor_43[2]**2
                cor_43_mag = math.sqrt(cor_43_mag)
                x_vec = cor_43 / cor_43_mag
                cor_23 = fp_cors[:, 1] - fp_cors[:, 2]
                cor_23_mag = cor_23[0]**2 + cor_23[1]**2 + cor_23[2]**2
                cor_23_mag = math.sqrt(cor_23_mag)
                y_vec = cor_23 / cor_23_mag
                z_vec = [(x_vec[1] * y_vec[2]) - (x_vec[2] * y_vec[1]),
                         (x_vec[2] * y_vec[0]) - (x_vec[0] * y_vec[2]),
                         (x_vec[0] * y_vec[1]) - (x_vec[1] * y_vec[0])]  
                rot_mat = np.column_stack([x_vec, y_vec, z_vec])
                
                # Find force plate center in lab coordinate system
                fp_x_cen = abs(fp_cors[0, 0] - fp_cors[0, 2]) / 2000 
                fp_y_cen = abs(fp_cors[1, 0] - fp_cors[1, 2]) / 2000

                # Output channel vectors
                channels = ['F12x', 'F34x', 'F14y', 'F23y', 
                            'F1z', 'F2z', 'F3z', 'F4z']
                for i in range(0, 8):
                    channels[i] = fp.GetChannel(i)
                    channels[i] = channels[i].GetValues()
                    channels[i] = channels[i][:, 0]
                F12x = channels[0]
                F34x = channels[1]
                F14y = channels[2]
                F23y = channels[3]
                F1z = channels[4] * -1
                F2z = channels[5] * -1
                F3z = channels[6] * -1
                F4z = channels[7] * -1

                # Force vectors
                Fx_o = F12x + F34x
                Fy_o = F14y + F23y
                Fz_o = F1z + F2z + F3z + F4z

                # Torque vectors
                Mx_o = (dy * (F1z + F2z - F3z - F4z)) 
                My_o = (dx * (-F1z + F2z + F3z - F4z))
                Mz_o = ((dy * (-F12x + F34x)) + (dx * (F14y - F23y)))

                # Threshold force and moment vectors                
                thres_indexes = [i for i, v in enumerate(Fz_o) if v < 20.0]
                for h in thres_indexes:
                    Fx_o[h] = 0
                    Fy_o[h] = 0
                    Fz_o[h] = 0
                    Mx_o[h] = 0
                    My_o[h] = 0
                    Mz_o[h] = 0

                # Center of pressure
                COPx1 = ((Fx_o * dy) - My_o) / Fz_o
                COPx1 = np.nan_to_num(COPx1)                
                COPy1 = ((Fy_o * dy) + Mx_o) / Fz_o
                COPy1 = np.nan_to_num(COPy1)
                COPz1 = [0.0] * len(COPy1)
                COP_mat = np.row_stack([COPx1, COPy1, COPz1])
                COP_mat = np.dot(rot_mat, COP_mat)
                COPx1 = COP_mat[0, ]
                COPy1 = COP_mat[1, ]
                COPz1 = COP_mat[2, ]

                # Apply offset to lab coordinates
                COPx1 = COPx1 + fp_x_cen
                COPy1 = COPy1 + fp_y_cen

                # Make array to be output
                Fx = format_list(Fy_o)
                Fy = format_list(Fz_o)
                Fz = format_list(Fx_o)
                Mx = format_list(My_o)
                My = format_list(Mz_o)
                Mz = format_list(Mx_o)
                COPx = format_list(COPy1)
                COPy = ['0.00000000'] * len(Fz)
                COPz = format_list(COPx1)

                forceplate_data = np.column_stack([Fx, Fy, Fz, COPx, COPy, 
                                                   COPz, Mx, My, Mz])
                fp_data1.append(forceplate_data)

            # Format and add time column
            fp_data = np.column_stack(fp_data1)
            fp_data = np.column_stack([time_frames_a, fp_data])


    # =========================================================================
    
    ## Additional details for .mot file
    # Header lines
    mot_1 = [file_name + '.mot']
    mot_2 = ['version=1']
    mot_3 = ['nRows=' + str(len(fp_data))]
    mot_4 = ['nColumns=' + str((9 * NumForcePlates) + 1)] 
    mot_5 = ['inDegrees=yes']
    mot_6 = ['endheader']
    
    # Column headings
    suffixes = []
    for ii in range(0, NumForcePlates):
        suffix = 'fp' + str(ii + 1) + '_' 
        suffixes.append(suffix)
    
    suffixes = np.repeat(suffixes, 9)
    
    col_heads = ['ground_force_vx', 'ground_force_vy', 
                 'ground_force_vz', 'ground_force_px', 'ground_force_py', 
                 'ground_force_pz','ground_torque_x', 'ground_torque_y', 
                 'ground_torque_z'] * NumForcePlates
    
    col_headers = [m + n for m, n in zip(suffixes, col_heads)]  
    col_headers = ['time'] + col_headers
    
    
    # =========================================================================

    ## Extract EMG data
    
    
    # =========================================================================
    
    ## Write .trc file (marker coordinate data)
    with open(trc_filepath, 'wb') as mycsvfile:
        thedatawriter = csv.writer(mycsvfile, delimiter = '\t')
        thedatawriter.writerow(trc_1)
        thedatawriter.writerow(trc_2)
        thedatawriter.writerow(trc_3)
        thedatawriter.writerow(trc_4)
        thedatawriter.writerow(trc_5)
        thedatawriter.writerow(trc_6)
        for row in marker_coords:
            thedatawriter.writerow(row)
    
    
    # =========================================================================
    
    ## Write .mot file (force plate and analogue data)
    with open(mot_filepath, 'wb') as mycsvfile:
        thedatawriter = csv.writer(mycsvfile, delimiter = '\t')    
        thedatawriter.writerow(mot_1)
        thedatawriter.writerow(mot_2)
        thedatawriter.writerow(mot_3)
        thedatawriter.writerow(mot_4)
        thedatawriter.writerow(mot_5)
        thedatawriter.writerow(mot_6)
        thedatawriter.writerow(col_headers)
        for row in fp_data:
            thedatawriter.writerow(row)


# =============================================================================

# Example
import easygui
trc_filepath = ('C:/Users/telfe/Dropbox/My_Projects/VA_Dose_Response/Open'
                'Sim/Pilot_data/ST_2016_07_26/ST_dynamic.trc')
mot_filepath = ('C:/Users/telfe/Dropbox/My_Projects/VA_Dose_Response/Open'
                'Sim/Pilot_data/ST_2016_07_26/ST_dynamic.mot')
c3d_file = easygui.fileopenbox()
c3d_to_opensim(c3d_file, trc_filepath, mot_filepath, threshold = -20)


