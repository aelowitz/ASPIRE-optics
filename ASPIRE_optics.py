

'''
This optics code was written by Amy Lowitz (lowitz@arizona.edu).  
The key underlying Gaussian optics equations are derived in many textbooks and 
papers and are standard to the field.  In most cases, however, I point to specific 
equation numbers from Paul Goldsmith's textbook, Quasioptical Systems, in inline comments, 
as a guidepost to any user needing an entry-point for learning more about the relevant
physical underpinnings.  


Author: Amy Lowitz - Spring 2025

Updates:



Summary:
This module calculates and plots various quantities useful in the design of a Gaussian optical
system.  This code was developed pecifically with the ASPIRE receiver (Advanced 
South Pole Integrated Receiver for EHT) in mind, but could be adapted to more general 
systems.  


Usage:
    1) edit hardcoded mirror parameters in generate_parameter_dict()
    2) run either main_forward() or main_reverse()





#TODO: Make agnostic to number of lenses
       Make agnostic to direction forward/reverse
       General cleanup


'''



import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt


def optimize_ff():
    '''
    
    '''
    
    #generate the ranges to vary over
    #d_5h_b7_list = np.linspace(90, 92, 25)
    #f_M5_b7_list = np.linspace(80.8, 81, 10)
    f_M3_list = np.linspace(100, 140, 10)
    f_M4_list = np.linspace(100,140, 10)
    
    
    #insert varying parameters into the beams dictionary
    beams = generate_parameter_dict()    
   
    #calculate the non-changing beam
    #beams['B6'] = beam_propagation(beams['B6'], edge_taper_dB = 30, Nz = 10000, truncate = 2000)
    
    #run through the varying params
    min_sum = 100000000
    min_f3 = 0
    min_f4 = 0        
    for f4 in f_M4_list:
        #print('next f')
        for f3 in f_M3_list:
        
            #insert varying parameters into the beams dictionary
            beams = generate_parameter_dict()
            beams['B7']['focal_lengths'][2] = f4
            beams['B6']['focal_lengths'][2] = f4
            beams['B7']['focal_lengths'][1] = f3
            beams['B6']['focal_lengths'][1] = f3
            
            
            #calculate the beam
            beams['B7'] = beam_propagation(beams['B7'], edge_taper_dB = 30, Nz = 10000, truncate = 2000)
            beams['B6'] = beam_propagation(beams['B6'], edge_taper_dB = 30, Nz = 10000, truncate = 2000)
            
            #calculate figure of merit
            M34_waist_diff = np.abs((beams['B6']['z_w0'][2] - beams['B6']['z_elements'][-2])-(beams['B7']['z_w0'][2] - beams['B7']['z_elements'][-2]))
            M23_waist_diff = np.abs((beams['B6']['z_w0'][3] - beams['B6']['z_elements'][-2])-(beams['B7']['z_w0'][3] - beams['B7']['z_elements'][-2]))            
            if min_sum > M34_waist_diff + M23_waist_diff:
                min_sum = M34_waist_diff + M23_waist_diff
                min_f3 = f3
                min_f4 = f4
                print(np.round(min_sum,3), np.round(min_f3,3), np.round(min_f4,3))
    
    return(min_f3, min_f4, min_sum)
    
    
def optimize_dd():
    '''
    
    '''
    
    #generate the ranges to vary over
    #d_5h_b7_list = np.linspace(90, 92, 25)
    #f_M5_b7_list = np.linspace(80.8, 81, 10)
    d_h5_B6_list = np.linspace(80,110, 25)
    d_45_B6_list = np.linspace(990, 1010, 25)
    
    #start with the base beams dictionary
    beams = generate_parameter_dict()
    
   
    #calculate the non-changing beam
    beams['B7'] = beam_propagation(beams['B7'], edge_taper_dB = 30, Nz = 10000, truncate = 2000)
    
    #run through the varying params
    min_sum = 100000000
    min_d1 = 0
    min_d2 = 0        
    for d1 in d_h5_B6_list:
        #print('next f')
        for d2 in d_45_B6_list:
        
            #beams['B7']['element_spacings'][2] = f4
            beams['B6']['element_spacings'][0] = d1
            #beams['B7']['element_spacings'][1] = f3
            beams['B6']['element_spacings'][1] = d2
            
            
            #calculate the beam
            #beams['B7'] = beam_propagation(beams['B7'], edge_taper_dB = 30, Nz = 10000, truncate = 2000)
            beams['B6'] = beam_propagation(beams['B6'], edge_taper_dB = 30, Nz = 10000, truncate = 2000)
            
            #calculate figure of merit
            M34_waist_diff = np.abs((beams['B6']['z_w0'][2] - beams['B6']['z_elements'][-2])-(beams['B7']['z_w0'][2] - beams['B7']['z_elements'][-2]))
            M23_waist_diff = np.abs((beams['B6']['z_w0'][3] - beams['B6']['z_elements'][-2])-(beams['B7']['z_w0'][3] - beams['B7']['z_elements'][-2]))            
            if min_sum > M34_waist_diff + M23_waist_diff:
                min_sum = M34_waist_diff + M23_waist_diff
                min_d1 = d1
                min_d2 = d2
                print(np.round(min_sum,3), np.round(min_d1,3), np.round(min_d2,3))
    
    return(min_d1, min_d2, min_sum)
    
def optimize_dfd():
    '''
    
    '''
    
    #generate the ranges to vary over
    d_h5_B6_list = np.linspace(80,105, 26)
    f5_B6_list = [80] #np.linspace(85, 95, 11)
    d_h5_B7_list = np.linspace(80, 105, 26)
    
    
    #start with the base beams dictionary
    beams = generate_parameter_dict()  
   
    #calculate the non-changing beam
    #beams['B7'] = beam_propagation(beams['B7'], edge_taper_dB = 30, Nz = 10000, truncate = 2000)
    
    #run through the varying params
    min_sum = 100000000
    min_d1 = 0
    min_f1 = 0
    min_d2 = 0        
    for d1 in d_h5_B6_list:
        for d2 in d_h5_B7_list:
            for f1 in f5_B6_list:

                beams['B6']['element_spacings'][0] = d1
                #beams['B6']['focal_lengths'][0] = f1
                beams['B7']['element_spacings'][0] = d2
                
                
                #calculate the beam
                beams['B7'] = beam_propagation(beams['B7'], edge_taper_dB = 30, Nz = 10000, truncate = 2000)
                beams['B6'] = beam_propagation(beams['B6'], edge_taper_dB = 30, Nz = 10000, truncate = 2000)
                
                #calculate figure of merit
                M34_waist_diff = np.abs((beams['B6']['z_w0'][2] - beams['B6']['z_elements'][-2])-(beams['B7']['z_w0'][2] - beams['B7']['z_elements'][-2]))
                M23_waist_diff = np.abs((beams['B6']['z_w0'][3] - beams['B6']['z_elements'][-2])-(beams['B7']['z_w0'][3] - beams['B7']['z_elements'][-2]))            
                if min_sum > M34_waist_diff + M23_waist_diff:
                    min_sum = M34_waist_diff + M23_waist_diff
                    min_d1 = d1
                    min_d2 = d2
                    min_f1 = f1
                    print(np.round(min_sum,3), np.round(min_d1,3), np.round(min_f1, 3), np.round(min_d2,3))
    
    return(min_d1, min_f1, min_d2, min_sum)    


def generate_parameter_dict():
    '''
        ##ALMA band 6: 211-275 GHz.  Band 7: 275-373.  Band 2: 67-116.
    
    '''
    
    #physics parameters and unchangable telescope and horn parameters
    c = 299792458 #m/s, speed of light in a vacuum
    f1 = 8003.21 #mm, focal length of the primary, from Padin 2008
    f2 = 818.2808 #mm, effective (thin lens equiv) focal length of the secondary, from Stark memo 2018.  f1 = 1310.0 mm, f2 = 2180.0 mm
    d_12 = 8003.21 + 1310 # mm, distance between the primary and secondary, from Stark memo. Focal length + distance from prime focus to M2
    dg = 2180 # mm, distance from secondary to gregorian focus.  From Stark memo.
    a_h_B6 = 3.54 #mm, horn aperture radius, from Greer memo
    a_h_B7 = 3.00225 # horn aperture radius, from ALMA Band 7 Cartridge Preliminary design report p.29
    R_h_B6 = 46.672 #mm, slant length of horn, from Greer memo
    R_h_B7 = 45.7817 # slant length of horn, from ALMA Band 7 cartridge Preliminary design report, p.29


    ######################################################################################
    ############# user-provided system parameters  #######################################
    

    #mirror parameters that are common between bands (but not common with SPT, so we can change them)
    f_M3 = 103#115#111.25  # focal length of M3  #135 in Gene's design
    f_M4 = 130#130# 133.75  # focal length of M4 #151 in Gene's
    d_23 = dg+135 #160   # distance between M2 and M3
    d_34 =  600 #570      # distance between M3 and M4 #600 in Gene's design
    
    
    #mirror parameters that are not common between bands
    freq_b6 = 227.1e9    # band 6 frequency, Hz
    f_M5_b6 = 3055.25 #72.6244        # thin lens equivalent focal length of M5, band 6 (cold mirror), mm
    d_45_b6 = 900 #995       # distance between M4 and M5 for band 6, mm.  918.903 in gene's
    d_5h_b6 = 62.122 #60       # distance between M5 and horn aperture for band 6, mm.  
    
    freq_b7 = 348.6e9    # band 7 frequency, Hz
    f_M5_b7 = 34.5      # thin lens equivalent focal length of M5, band 7 (cold mirror), mm
    d_45_b7 = 900      # distance between M4 and M5 for band 7, mm
    d_5h_b7 = 37       # distance between M5 and horn aperture for band 7, mm
    
    
    
    ######################################################################################
    ######################################################################################
    
    
    
    
    ############ put all the parameters into a dict  #################
    beams = {'B6':{}, 'B7':{}}

    beams['B6']['focal_lengths'] = [f1, f2, f_M3, f_M4, f_M5_b6] #  [f_M5_b6, f_M4, f_M3, f2, f1]
    beams['B6']['element_spacings'] = [d_12, d_23, d_34, d_45_b6, d_5h_b6]   #[d_5h_b6, d_45_b6, d_34, d_23, d12]

    beams['B7']['focal_lengths'] = [f1, f2, f_M3, f_M4, f_M5_b7]    #[f_M5_b7, f_M4, f_M3, f2, f1]
    beams['B7']['element_spacings'] = [d_12, d_23, d_34, d_45_b7, d_5h_b7]   #[d_5h_b7, d_45_b7, d_34, d_23, d12]
    
    
    # calculated system parameters  #
    
    # band 6
    beams['B6']['frequency'] = freq_b6
    beams['B6']['greg_foc'] = dg
    beams['B6']['lambda'] = (c/beams['B6']['frequency'])*1000 #wavelength in mm
    beams['B6']['a_h'] = a_h_B6
    beams['B6']['R_h'] = R_h_B6
    beams['B6']['w_h'] = beams['B6']['a_h']*0.644 # beam radius at the aperture of the horn.  # Goldsmith eq 7.41.  The .644 comes from fig 7.6.  Assumes corrugated circular horn.
    
    # band 7
    beams['B7']['frequency'] = freq_b7
    beams['B7']['greg_foc'] = dg
    beams['B7']['lambda'] = (c/beams['B7']['frequency'])*1000 #wavelength in mm
    beams['B7']['a_h'] = a_h_B7 # horn aperture radius, from ALMA Band 7 Cartridge Preliminary design report p.29
    beams['B7']['R_h'] = R_h_B7
    beams['B7']['w_h'] = beams['B7']['a_h']*0.644 #beam radius at the aperture of the horn.  # Goldsmith eq 7.41.  The .644 comes from fig 7.6.  Assumes corrugated circular horn.        
    
    return beams
    
    
def main_forward(plot = True):
    '''
    From primary, moving towards horn
    '''

    beams = generate_parameter_dict()
    
    beams['B6'] = forward_beam_propagation(beams['B6'], Nz = 10000)
    beams['B7'] = forward_beam_propagation(beams['B7'], Nz = 10000)
    
    if plot:
        thin_lens_plot_forward(beams)
    
    
    return beams
    


    

def main_reverse(plot = True):
    '''
    from horn moving towards primary
    
    '''
    
    beams = generate_parameter_dict()
    
    beams['B6'] = reverse_beam_propagation(beams['B6'], edge_taper_dB = 30, Nz = 10000)
    beams['B7'] = reverse_beam_propagation(beams['B7'], edge_taper_dB = 30, Nz = 10000)
    
    M34_waist_diff = (beams['B6']['z_w0'][2] - beams['B6']['z_elements'][-2])-(beams['B7']['z_w0'][2] - beams['B7']['z_elements'][-2])
    M23_waist_diff = (beams['B6']['z_w0'][3] - beams['B6']['z_elements'][-2])-(beams['B7']['z_w0'][3] - beams['B7']['z_elements'][-2])
    
    
    #beam waist offsets between bands
    print('Waist 3-4: {}'.format(np.round(M34_waist_diff, 2)))
    print('Waist 2-3: {}'.format(np.round(M23_waist_diff, 2)))
    
    #beam widths at fixed mirrors
    secondary_index_b6 = beams['B6']['z_list'].index(beams['B6']['z_elements'][4])
    secondary_index_b7 = beams['B7']['z_list'].index(beams['B7']['z_elements'][4])
    w_secondary_b6 = beams['B6']['w'][secondary_index_b6]
    w_secondary_b7 = beams['B7']['w'][secondary_index_b7]
    print('beam radius at secondary: {}, {}'.format( np.round(w_secondary_b6, 2) , np.round(w_secondary_b7, 2) ))
    print('beam radius at primary: {}, {}'.format( np.round(beams['B6']['w'][-1], 2) , np.round(beams['B7']['w'][-1],2)  ))
    
    if plot:
        thin_lens_plot_reverse(beams)
    
    return beams
    
    
    

def thin_lens_plot_reverse(beams, plot_phi = False):
    '''
    
    
    
    '''
    
    plt.ion()
    plt.figure()
    
    plt.plot(np.array(beams['B6']['z_list']) - beams['B6']['z_elements'][-2], beams['B6']['w'], '-r', linewidth = 2, label = 'B6 beam width')
    plt.plot(np.array(beams['B7']['z_list']) - beams['B7']['z_elements'][-2], beams['B7']['w'], '-b', linewidth = 2, label = 'B7 beam width')
    
    plt.text(-2424, 65, "M3", fontsize = 14, color = 'black')
    plt.text(-3050, 35, "M4", fontsize = 14, color = 'black')
    
    plt.vlines(beams['B6']['z_w0'][2] - beams['B6']['z_elements'][-2], ymin = 0, ymax = 10, color = 'magenta', label = 'B6 waist')
    plt.vlines(beams['B7']['z_w0'][2] - beams['B7']['z_elements'][-2], ymin = 0, ymax = 10, color = 'cyan', label = 'B7 waist')
    
    #beam waists
    plt.vlines(beams['B6']['z_w0'][3] - beams['B6']['z_elements'][-2], ymin = 0, ymax = 10, color = 'magenta')
    plt.vlines(beams['B7']['z_w0'][3] - beams['B7']['z_elements'][-2], ymin = 0, ymax = 10, color = 'cyan')
    
    #mirrors
    plt.vlines(8310, ymin = 0, ymax = 5000, color = 'black') #primary
    plt.hlines(4443, xmin = 8210, xmax = 8410, color='black') # primary 11dB
    plt.vlines(0, ymin=0, ymax = 875, color = 'black') #secondary
    plt.hlines(437.5, xmin = -100, xmax = +100, color='black')
    plt.vlines(beams['B6']['z_elements'][-3] - beams['B6']['z_elements'][-2], ymin = 0, ymax = 35.5, color = 'black') # tertiary
    plt.vlines(beams['B6']['z_elements'][-4] - beams['B6']['z_elements'][-2], ymin = 0, ymax = 20.3, color = 'black') # quaternary
    plt.vlines(beams['B6']['z_elements'][-5] - beams['B6']['z_elements'][-2], ymin = 0, ymax = 7.62, color = 'black') # quinary
    
    
    plt.grid()
    plt.xlabel('distance from secondary [mm]')
    plt.ylabel('beam radius [mm]')
    plt.xlim([-4200,-1800])
    plt.ylim([-10, 150])
    plt.legend()
    
    if plot_phi:
        plt.figure()
        plt.plot(np.array(beams['B6']['z_list']) - beams['B6']['z_elements'][-2], beams['B6']['phi'], '-r', linewidth = 2, label = 'B6 phi')
        plt.plot(np.array(beams['B7']['z_list']) - beams['B7']['z_elements'][-2], beams['B7']['phi'], '-b', linewidth = 2, label = 'B7 phi')
        
        plt.text(-2424, 0, "M3", fontsize = 14, color = 'black')
        plt.text(-3000, 0, "M4", fontsize = 14, color = 'black')
        
        
        plt.xlabel('distance from secondary [mm]')
        plt.ylabel('phi [radians]')
        plt.legend()
    

def thin_lens_plot_forward(beams):
    '''


    '''

    plt.ion()
    plt.figure()
    
  
    #plot beamwidths    
    plt.plot(beams['B6']['z_list'], beams['B6']['w'], '-r', linewidth = 2, label = 'B6 beamwidth')
    plt.plot(beams['B7']['z_list'], beams['B7']['w'], '-b', linewidth = 2, label = 'B7 beamwidth')
    
    #beam waists
    plt.vlines(beams['B6']['z_w0'][2], ymin = 0, ymax = 10, color = 'magenta', label = 'B6 waist')
    plt.vlines(beams['B6']['z_w0'][3], ymin = 0, ymax = 10, color = 'magenta')
    plt.vlines(beams['B7']['z_w0'][2], ymin = 0, ymax = 10, color = 'cyan', label = 'B7 waist')
    plt.vlines(beams['B7']['z_w0'][3], ymin = 0, ymax = 10, color = 'cyan')
    
    #mirrors
    plt.vlines(0, ymin = 0, ymax = 5000, color = 'black') #primary
    plt.hlines(4443, xmin = -20, xmax = 20, color='black') #primary 11dB
    plt.vlines(beams['B6']['z_elements'][1], ymin=0, ymax = 875, color = 'black') #secondary
    plt.vlines(beams['B6']['z_elements'][2], ymin = 0, ymax = 35.5, color = 'black') #tertiary 5w
    plt.hlines(44.375, xmin = beams['B6']['z_elements'][2]-10, xmax = beams['B6']['z_elements'][2]+10, color='black') #tert 4w
    plt.vlines(beams['B6']['z_elements'][3], ymin = 0, ymax = 20.3, color = 'black') #quaternary
    plt.hlines(25.375, xmin = beams['B6']['z_elements'][3]-10, xmax = beams['B6']['z_elements'][3]+10, color='black') #quat 4w
    plt.vlines(beams['B6']['z_elements'][4], ymin = 0, ymax = 7.62, color = 'black') #quinary
    plt.hlines(9.525, xmin = beams['B6']['z_elements'][4]-10, xmax = beams['B6']['z_elements'][4]+10, color='black') #quin 4w
    
    #mirror labels
    plt.text(beams['B6']['z_elements'][2]-20, -10, 'M3', fontsize = 12, color = 'black')
    plt.text(beams['B6']['z_elements'][3]-20, -10, 'M4', fontsize = 12, color = 'black')
    plt.text(beams['B6']['z_elements'][4]-20, -10, 'M5', fontsize = 12, color = 'black')
    plt.text(beams['B6']['z_w0'][-1]-20, -10, 'horn', fontsize = 12, color = 'black')
    
    #tidy plot
    plt.xlim([11450, 13500])
    plt.ylim([-32, 175])
    plt.grid()
    plt.xlabel('distance from primary [mm]')
    plt.ylabel('beam radius [mm]')
    plt.legend()
    plt.tight_layout()



  
     
def forward_beam_propagation(single_freq_beam, Nz = 100, primary_truncate_dB = 11, truncate = 0, calculate_m5 = True):

    focal_lengths = single_freq_beam['focal_lengths']
    element_spacings = single_freq_beam['element_spacings']
    lambd = single_freq_beam['lambda']
    R_h = single_freq_beam['R_h']
    w_h = single_freq_beam['w_h']
    dg = single_freq_beam['greg_foc']


    z_elements = [0] #primary is at z=0
    for distance in element_spacings:
        new_z  = z_elements[-1] + distance # add the next distance onto the last z-position
        z_elements.append(new_z)
        
    #calculate waist radius and z-location of beam waists between each optical element
    z_w0 = [focal_lengths[0], z_elements[1] + dg] #first two waists are at prime and gregorian foci
    w_primary = 5000/(0.3393 * (primary_truncate_dB**.5))  #goldsmith 2.35b, beamwidth at primary
    w0_primary = min(w0_calc(lambd, w_primary, focal_lengths[0])) #primary beam waist radius
    w0_secondary = lens(w0_primary, z_elements[1]-focal_lengths[0], focal_lengths[1], lambd)[0]
    w0 = [w0_primary, w0_secondary]
    focal_length_counter = 2
    for element_z in z_elements:
        if element_z > z_elements[1] and element_z < z_elements[-1]: # don't do anything with the first two or the last one except store it
            d_in = element_z - previous_z_w0 #location of previous beam waist
            w0_out, d_out = lens(previous_w0, d_in, focal_lengths[focal_length_counter], lambd)
            w0.append(w0_out)
            z_w0.append(d_out + element_z)
            focal_length_counter += 1
        previous_z_w0 = z_w0[-1]
        previous_w0 = w0[-1]
    
    print('w0: {}'.format(np.round(w0, 2)))
    print('z_w0: {}'.format(np.round(z_w0, 2)))
     
    #generate points along z-axis
    z_list = [0] #primary position
    previous_element_z = 0
    for element_z in z_elements[1:-1]: #from after the primary up to the last mirror before the horn
        if element_z - previous_element_z <= 500: # for small element spacings, use Nz
            n_points = Nz+1
        elif element_z - previous_element_z <= 2000: #for medium element spacings, use slightly denser points
            n_points = (Nz*2)+1
        else: #for big element spacings, use more points
            n_points = (Nz*6)+1
        z_steps = np.linspace(previous_element_z, element_z, n_points).tolist()[1:] #chop off the first one because it's already in z
        z_list.extend(z_steps)
        previous_element_z = element_z
    
    #add points between horn aperture and horn waist
    z_steps = np.linspace(z_list[-1], z_w0[-1], 200)[1:]
    z_list.extend(z_steps)

    
    if truncate > 0:
        z_list = [x for x in z_list if x <= truncate]
    
    # for every point in z, calculate the beam radius (w), radius of wavefront curvature (R), phase slippage (phi), and edge taper radius (edge_taper_w)
    
    w = []
    R = []
    phi = []
    
    #edge_taper_w = []
    #edge_taper_factor = 0.3393*(edge_taper_dB**.5) #Goldsmith 2.35b
    
    element_number = 1
    for z in z_list: #step through each z-position in the system
        if z > z_elements[element_number] and element_number < len(z_elements)-1: #move on to the next element
            element_number += 1
        w_tmp, R_tmp = get_wR(lambd, z-z_w0[element_number-1], w0[element_number-1])
        #phi_tmp = phi_slippage(lambd, z-z_w0[element_number-1], w0[element_number-1])
        w.append(w_tmp)
        R.append(R_tmp)
        #phi.append(phi_tmp)
        #edge_taper_w.append(w_tmp*edge_taper_factor)
        
    single_freq_beam['w0'] = w0
    single_freq_beam['z_w0'] = z_w0
    single_freq_beam['z_elements'] = z_elements
    single_freq_beam['z_list'] = z_list
    single_freq_beam['w'] = w
    single_freq_beam['R'] = R
    #single_freq_beam['phi'] = phi
    #single_freq_beam['edge_taper_w'] = edge_taper_w
    
    if calculate_m5:
        
        #calculate horn parameters
        w0_in, z_offset_horn = w0z0_calc(lambd, R_h, w_h) #z_offset is the distance inside the horn (measured from the aperture) where the beam waist falls
        w0_out = w0[3]# waist between M4 and M5
        print('Horn waist for {} mm: {} mm'.format(np.round(lambd,2), np.round(w0_in,2)))
        d_in = z_offset_horn + element_spacings[-1] #distance from horn waist to last mirror
        d_out = z_elements[-2]-z_w0[-2]
        f_M5 = reverse_lens(wo_in, d_in, d_out, w0_out, lambd)  #returns a tuple, both roots
        
        print('Suggested M5 focal lengths for {} mm:  {}'.format(np.round(lambd,2), np.round(f_M5,4)))
        
    return single_freq_beam     
        
        
def reverse_beam_propagation(single_freq_beam, edge_taper_dB = 30, Nz = 100, truncate = 0):
    '''
    
    
    Parameters:
    -----------
    single_freq_beam: dict containing the following keys:
        focal_lengths: list of floats, focal lengths of each optical element [primary, secondary, ...] in mm
        element_spacings: list of floats, distances between each optical element [primary to secondary, secondary to tertiary, ...]
        lambda: float, wavelength in mm
        R_h: float, slant length of horn in mm, equal to the radius of curvature at the horn aperture for conical corrugated horns
        w_h: float, beam radius at the aperture of the horn in mm
        Nz: int, how many points along the z-axis to calculate between each optical element
            e.g. if there are 3 optical elements, then (2*Nz)+1 points will be calculated.  
    edge_taper_dB: int or float, how many dB out to go for the edge tapered radius
    Nz: int, how many points to calculate between elements
    truncate: int or float, z-position (measured from horn) after which to truncate the 
        calculations of w, R, etc. (makes things run faster for optimization stuff, when you don't
        care about the secondary and primary).  If truncate <= 0, then no truncation is done.
    
    
    
    
    Returns:
    --------
    The original dict with the following keys added:
        w0: list of floats, beam waist radii, typically in mm, starting from the horn waist
        z_w0: list of floats, z-positions of each beam waist, typically in mm, starting 
            from the horn waist
        z_elements: list of floats, z-positions of each element, typically in mm, starting 
            from the horn aperture
        z_list: list of floats, every z-position we will calculate w, R, and phi for.  
            0 is at the horn aperture. Typically in mm.  
        w: list of floats, the beam radius at every z-position in z-list.  Typically in mm
        R: list of floats, the wavefront radius of curvatureat every z-position in 
            z-list.  Typically in mm.
        phi: list of floats, the phase shift or "slippage" at every z-position in z-list

    
    '''
    
    focal_lengths = single_freq_beam['focal_lengths']
    focal_lengths.reverse()
    element_spacings = single_freq_beam['element_spacings']
    element_spacings.reverse()
    lambd = single_freq_beam['lambda']
    R_h = single_freq_beam['R_h']
    w_h = single_freq_beam['w_h']
    
    
    #calculate the location of each optical element along the z-axis (optical axis)
    z_elements = [0] #horn aperture is at z=0  #called z_e in Zeng code
    for distance in element_spacings:
        new_z  = z_elements[-1] + distance # add the next distance onto the last z-position
        z_elements.append(new_z)
    

    #calculate the waist radius and z-location of the beam waists between each optical element
    w0_horn, z_offset_horn = w0z0_calc(lambd, R_h, w_h) #z_offset is the distance inside the horn (measured from the aperture) where the beam waist falls
    w0 = [w0_horn]
    z_w0 = [-z_offset_horn] #negative because the waist will be located inside the horn and we're calling the horn aperture z=0  #list of the z-position of each waist
    focal_length_counter = 0
    for element_z in z_elements:
    
        if element_z > 0: # don't do anything with the first one except store it
            d_in = element_z - previous_z_w0 # where is the previous beam waist  #TODO# re-check this for minus signs.  might need to be + instead of -.  
            w0_i, d_out = lens(previous_w0, d_in, focal_lengths[focal_length_counter], lambd)
            w0.append(w0_i)
            z_w0.append(d_out + element_z)
            focal_length_counter += 1            
        previous_z_w0 = z_w0[-1]
        previous_w0 = w0[-1]


    #generate points along z-axis
    for element_z in z_elements:
        if element_z == 0: #for the first one, deal with horn waist location at negative z
            z_list = np.linspace(-z_offset_horn, 0, Nz+1).tolist()
        else:
            if element_z - previous_element_z <= 1000: # for small element spacings, use Nz
                z_steps = np.linspace(previous_element_z, element_z, Nz+1).tolist()[1:] #chop off the first one because it's already in z
            elif element_z - previous_element_z <= 2000: #for medium element spacings, use slightly denser points
                z_steps = np.linspace(previous_element_z, element_z, (Nz*2)+1).tolist()[1:] #chop off the first one because it's already in z
            else: #for big element spacings, use much denser points
                z_steps = np.linspace(previous_element_z, element_z, (Nz*6)+1).tolist()[1:] #chop off the first one because it's already in z
            z_list.extend(z_steps)
        previous_element_z = element_z
    if truncate > 0:
        z_list = [x for x in z_list if x <= truncate]
        
    

    
    #for every point in z, calculate the beam radius (w), radius of curvature (R), phase slippage (phi), and edge taper radius (edge_taper_w)
    w = []
    R = []
    phi = []
    edge_taper_w = []
    
    edge_taper_factor = 0.3393*(edge_taper_dB**.5) #Goldsmith 2.35b
    
    element_number = 1
    for z in z_list:  #step through each z-position in the system
        if z > z_elements[element_number]: #move on to the next element
            element_number += 1
            #print('element number: {}.  z = {}'.format(element_number, z))
        #print(z_w0[element_number], element_number)
        w_tmp, R_tmp = get_wR(lambd, z-z_w0[element_number-1], w0[element_number-1])   
        phi_tmp = phi_slippage(lambd, z-z_w0[element_number-1], w0[element_number-1])
        w.append(w_tmp)
        edge_taper_w.append(w_tmp*edge_taper_factor)
        R.append(R_tmp)
        phi.append(phi_tmp)
    
    
    single_freq_beam['w0'] = w0
    single_freq_beam['z_w0'] = z_w0
    single_freq_beam['z_elements'] = z_elements
    single_freq_beam['z_list'] = z_list
    single_freq_beam['w'] = w
    single_freq_beam['R'] = R
    single_freq_beam['phi'] = phi
    single_freq_beam['edge_taper_w'] = edge_taper_w
    
    return single_freq_beam
   
   
    

def w0z0_calc(lambd, R, w):
    '''
    Calculate the beam waist radius and location when you know the radius of 
    curvature (R) at some other particular beam radius (w). 
    
    Parameters:
    -----------
    lambd: float, wavelength in mm
    R: float or int, the radius of curvature in mm at the location with beam radius w
    w: float or int, the beam radius in mm at the location with radius of curvature R
    
    Returns:
    --------
    w0 float, the beam diameter in mm at the beam waist
    z0: float, the offset along the z axis in mm between the beam waist and the 
        location with beam radius w
    
    
    '''
    
    w0 = w/((1+(np.pi*w**2/(lambd*R))**2)**.5)  #goldsmith table 2.3 line 6
    z0 = R/(1+(lambd*R/(np.pi*w**2))**2)        # goldsmith table 2.3 line 6.  Note I have taken out a minus sign compared to the Lingzhen Zeng code, to match Goldsmith and treat this as a non-directed offset
    
    return w0, z0


def w0_calc(lambd, w, z):
    '''
    Calculate the beam waist radius, w0, given the beam radius at a given distance, z, 
    away from the beam waist.  Not currently (5/25) called elsewhere in the code, but
    useful for "by hand" calculations.
    
    Parameters:
    ----------
    lambd: float, wavelength typically in mm
    w: float, 1/e beam width radius at a distance z from the beam waist, typically in mm
    z: float, distance from the beam waist to where w is measured, typically in mm
    
    '''
    
    w0_pos = ((w**2/2) * (1 + ((1 - ((2*lambd*z/(3.1415*w**2))**2))**.5)))**.5
    w0_neg = ((w**2/2) * (1 - ((1 - ((2*lambd*z/(3.1415*w**2))**2))**.5)))**.5
    
    return w0_pos, w0_neg


def lens(w0_in, d_in, f, lambd):
    '''
    Calculates the beam transformation through a paraxial thin lens (or equivalent mirror).
    Note this function by itself works with any consistent units, but in typical usage 
    within this module all of the parameters are in mm, so that's why the Parameters
    and Returns descriptions specify those units. 
    
    Parameters:
    ----------
    w0_in: float, the beam waist radius in mm, before the lens/mirror.  See 
        Goldsmith fig 3.6
    d_in: float, the distance between the lens/mirror and w0_in, along the z axis, in mm.
        See Goldsmith fig 3.6
    f: float, the focal length of the lens/mirror in mm
    lambd: float, the wavelength in mm
    
    
    Returns:
    --------
    w0_out: float, the beam waist radius in mm, after the mirror/lens.  See 
        Goldsmith fig 3.6
    d_out: float, the distance between the lens/mirror and w0_out in mm.  See
        Goldsmith fig 3.6 
        
    '''
    
    z = 3.14159*(w0_in**2)/lambd
    c = d_in/f
    
    d_out = (1 + ((c-1) / (((c-1)**2) + ((z**2)/(f**2))))) * f      # Goldsmith 3.31a
    w0_out = w0_in / ((((c-1)**2) + ((z**2)/(f**2)))**.5)           # Goldsmith 3.31b
    
    if d_out < 0:
        print('WARNING: d_out < 0: d_out = {} for f = {}'.format(d_out, f))
    
    return w0_out, d_out
    
    

def reverse_lens(w0_in, d_in, d_out, w0_out, lambd):

    zc = 3.14159*(w0_in**2)/lambd
    
    #using the d_out equation (given d_in, d_out, w0_in.  W0_out and f free)
    A = d_out + d_in
    B = (-2*d_in*d_out) - (d_in**2) - (zc**2)
    C = ((d_in**2) * d_out) + (d_out * (zc**2))  # plus or minus sign here is questionable
    
    
    f_pos_d = (-B + ((B**2 - (4*A*C))**.5)) / (2*A)
    f_neg_d = (-B - ((B**2 - (4*A*C))**.5)) / (2*A)
    
    
    #using the w0_out equation.  Given w0_in, w0_out, d_in.  d_out, and f free.  
    A = -1 - ((w0_in/w0_out)**2)
    B = -2*d_in
    C = (zc**2) + (d_in**2)
    
    f_pos_w = (-B + ((B**2 - (4*A*C))**.5)) / (2*A)
    f_neg_w = (-B - ((B**2 - (4*A*C))**.5)) / (2*A)
    
    [w0_out_neg, d_out] = lens(w0_in, d_in, f_neg_w, lambd)
    [w0_out_pos, d_out] = lens(w0_in, d_in, f_pos_w, lambd)
    print('M5 to horn waist dist for waist matching: {}, {}'.format(np.round(w0_out_pos, 2), np.round(w0_out_neg, 2)))
       
    return f_pos_d, f_neg_d, f_pos_w, f_neg_w
    
def lens_explorer(w0_in, d_in, f_list, lambd):
    '''
    Makes a plot of focal lengths vs d_outs for a given d_in and w0_in
    
    '''
    plt.ion()
    d_out = []
    w0_out = []
    for f in f_list:
        w0_out_tmp, d_out_tmp = lens(w0_in, d_in, f, lambd)
        d_out.append(d_out_tmp)
        w0_out.append(w0_out_tmp*10)
    
    plt.figure()
    plt.plot(f_list, d_out, label = 'd_out')
    plt.grid()
    plt.xlabel('focal length [mm]')
    plt.tight_layout()
    plt.plot(f_list, w0_out, label = 'w0')
    plt.legend()
    plt.ylabel('w0_out*10 or d_out [mm]')
    plt.tight_layout()
    

        
def M5_focal_vs_mirror_size(lambd, w0_in, d_out):
    '''
    Makes a plot showing the M5 focal length required for the mirror closest to the horn
    vs the required 5w mirror diameter for a given d_out (which should be roughly the distance
    from the mirror to the window)
    
    '''           
    plt.figure()
    pos = []
    neg = []
    diameter = []
    d_in_list = np.linspace(.1,150,3000)
    for d_in in d_in_list:
        w_tmp, R_tmp = get_wR(lambd, d_in, w0_in)
        diameter.append(w_tmp*2*5/10) #cm, dimeter of 5w mirror
        a,b,c,d = reverse_lens(w0_in, d_in, d_out, 1, lambd)
        pos.append(a)
        neg.append(b)
    plt.plot(diameter, pos, label = 'positive root')
    plt.plot(diameter, neg, label = 'negative root')
    plt.grid()
    plt.legend()
    plt.xlabel('5w mirror diameter [cm]')
    plt.ylabel('focal length of M5 [mm]')
    plt.title('Requirement for d_out of {}'.format(d_out))
    
    plt.figure()
    plt.plot(d_in_list, pos, label = 'positive root')
    plt.plot(d_in_list, neg, label = 'negative root')
    plt.grid()
    plt.legend()
    plt.xlabel('d_in [mm]')
    plt.ylabel('focal length of M5 [mm]')
    plt.title('Requirement for d_out of {}'.format(d_out))
    
    
    indices = np.where(~np.isnan(neg))
    min_5w = diameter[indices[0][0]]
    min_5w_inches = min_5w/2.54
    min_focal_neg = neg[indices[0][0]]
    min_focal_pos = pos[indices[0][0]]
    d_in_min = d_in_list[indices[0][0]]
    
    print('Minimum 5w diameter: {} cm, {} in.'.format(np.round(min_5w, 2), np.round(min_5w_inches, 2)))
    print('focal lengths: {} mm (pos), {} mm (neg).'.format(np.round(min_focal_pos,2), np.round(min_focal_neg, 2)))
    print('d_in {} mm'.format(np.round(d_in_min, 2)))

    #return(diameter, neg)
    
             
    
    
def constraint(vars, X, z_c, w_oi, w_oo):
    '''
    
    '''
    d_i, d_o, f = vars
    
    eq_spacing = d_i + d_o - X
    
    denominator = ((d_i/f)**2) + ((z_c/f)**2)
    
    eq_d_out = 1 + (((d_i/f)-1) / denominator) - (d_o/f)
    
    eq_w0_out = (w_oi / (denominator**.5)) - w_oo
    
    return [eq_spacing, eq_d_out, eq_w0_out]


def solve_constraint(X, z_c, w_oi, w_oo, initial_guess = (100, 1000, 100)):
    '''
    
    
    '''
    d_i, d_o, f = opt.fsolve(constraint, x0 = initial_guess, args = (X, z_c, w_oi, w_oo)) 
    
    return d_i, d_o, f
    
    
    

    
    

def get_wR(lambd, z, w_0): 
    '''
    Get the beam radius (w) and radius of curvature (R) 
    when you have the beam waist (w_0), propagation distance (z), and wavelength (lambd).
    Typically use mm for the lengths, but this function on its own works for any consistent unit. 

    Parameters:
    ----------
    lambd: int or float, wavelength (lambd instead of lambda because lambda is a special word in python)
    z: int or float, distance along the optical axis from the beam waist
    w_0: int or float, the radius of the beam waist

    Returns:
    --------
    w: int or float, the 1/e Gaussian beam radius at location z
    R: int or float, the radius of curvature of the wavefront at location z
    
    '''
    w = w_0*np.sqrt(1+(lambd*z/(3.14159*w_0**2))**2)    #goldsmith eq 2.21b and 2.26c
    try:
        R = z*(1+(3.14159*w_0**2/(lambd*z))**2)             #goldsmith eq 2.21a and 2.26b
    except:
        print('exception1')
        if z == 0:
            R = float('nan')
        else:
            raise
    
    return w, R
    
    
def phi_slippage(lambd, z, w0):
    '''
    Calculate the phase slippage when you know the beam waist radius, z offset from the waist, and wavelength
    
    
    Parameters:
    -----------
    lambd: float or int, wavelength typically in mm
    z: float or int, distance along the optical axis from the beam waist, typically in mm
    w0: float or int, radius of the beam waist, typically in mm
    
    Returns:
    -------
    phi: float, Gaussian beam phase shift in radians
    
    
    '''
    
    phi = np.arctan((lambd*z)/(np.pi*w0**2))    # Goldsmith 2.26d
    
    return phi
    
    
    
    
    
    
    
    
    

    
    
