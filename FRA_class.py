import pandas as pd
import seaborn as sns
import numpy as np
from numpy import load
import scipy as sp
import matplotlib.pyplot as plt

from scipy import signal
from scipy.signal import welch


class FRA_irregularities():
    def __init__(self,L_min=1.524,L_max=304.8,N=3000,k=0.25,signal=None,dt=None,signal_type='vert'):
        
        self.L_min=L_min # The default values of L_min and L_max are the the maximum permissible range of wavelengths for FRA standards
        self.L_max=L_max
        self.irreg_type = None
        self.N_harmonics = N
        self.k = k
        
        self.omega_max = None
        self.omega_min = None
        self.d_omega = None
        self.omega = None
        
        self.wave=None
        
        self.s_vert=None
        self.s_lat=None
        self.vert_irreg = None
        self.lat_irreg = None
        
        self.class_list = None
        
        self.signal = signal
        self.dt = dt
        self.signal_type = signal_type
        
    def PSD(self,type_irreg=6):
        
        self.omega_max = 2*np.pi/self.L_min # Maximum angular frequency (spatial wavenumber) in rad/m
        self.omega_min = 2*np.pi/self.L_max # Minimum angular frequency (spatial wavenumber) in rad/m
    
        self.d_omega = (self.omega_max-self.omega_min)/self.N_harmonics      # Frequency increment (rad/m)
        n = np.arange(1,self.N_harmonics+1,1)                      # index vector
        self.omega = self.omega_min + (n-0.5)*self.d_omega    # discrete angular frequency (rad/m)
    
        # Creating wavelength domain vector
        self.wave = 2*np.pi/self.omega
    
        if type_irreg == 6:
            Av = 0.0339*10**-4  # m^2 * (rad/m)
            Aa = 0.0339*10**-4  # m^2 * (rad/m)
            omega_c = 0.8245    # rad/m
            omega_s = 0.438     # rad/m
       
        elif type_irreg == 5:
            Av = 0.2095*10**-4  
            Aa = 0.0762*10**-4 
            omega_c = 0.8245   
            omega_s = 0.8209    
       
        elif type_irreg == 4:
            Av = 0.5376*10**-4  
            Aa = 0.3027*10**-4 
            omega_c = 0.8245   
            omega_s = 1.1312    

        elif type_irreg == 3:
            Av = 0.6816*10**-4  
            Aa = 0.4128*10**-4 
            omega_c = 0.8245   
            omega_s = 0.852    
    
        
        elif type_irreg == 2:
            Av = 1.0181*10**-4  
            Aa = 1.2107*10**-4 
            omega_c = 0.8245   
            omega_s = 0.9308    
    
        
        elif type_irreg == 1:
            Av = 1.2107*10**-4  
            Aa = 3.3634*10**-4 
            omega_c = 0.8245   
            omega_s = 0.6046
            
        else:
            print('Provide a FRA classe between 6 and 1')    
            return None
    
        self.s_vert = 2*np.pi*(self.k*Av*omega_c**2)/((self.omega**2)*(self.omega**2+omega_c**2)) # m^2/(1/m)
        self.s_lat = 2*np.pi*(self.k*Aa*omega_c**2)/((self.omega**2)*(self.omega**2+omega_c**2)) # m^2/(1/m)
        
        return self.wave,self.omega,self.s_vert,self.s_lat
        
    def _create_PSD(self,class_list=[6,5,4]):
        
        self.class_list=class_list
        self.vert_irreg = []
        self.lat_irreg = []
        
        for item in self.class_list:
            _,_,vert, lat = FRA_irregularities.PSD(self,type_irreg=item)
            self.vert_irreg.append(vert)
            self.lat_irreg.append(lat)
        
        print('Classes {} were created'.format(class_list))
        
        return self.wave,self.omega,self.vert_irreg,self.lat_irreg
    
    # def plot_PSD(self,x_axis='wavelength',v_line=True,scale='log',plot_signal=False):
        
    #     if self.vert_irreg==None:
    #         FRA_irregularities._create_PSD(self,class_list=[6,5,4,3]) # Create PSDs
        
    #     fig, ax = plt.subplots(2,1,figsize=(30,20))
    #     ax[0].set_ylim(10**(-14),10**0)
    #     ax[1].set_ylim(10**(-14),10**0);  
        
    #     ax[0].set_title('PSD standard curves for vertical irregularities')
    #     ax[0].set_ylabel('PSD $(m^2/(1/m)$')
    #     ax[0].set_xticks
    #     ax[1].set_title('PSD standard curves for lateral irregularities')
    #     ax[1].set_ylabel('PSD $(m^2/(1/m)$')                 
        
        
    #     if x_axis == 'spatial_angular_frequency':
            
    #         for idx,item in enumerate(self.class_list):
    #             ax[0].plot(self.omega,self.vert_irreg[idx],label='FRA class {}'.format(self.class_list[idx]))
    #             ax[1].plot(self.omega,self.lat_irreg[idx],label='FRA class {}'.format(self.class_list[idx]))
            
    #         if v_line == True:
    #             ax[0].axvline(self.omega.min(), label='Lower spatial frequency: {} rad/m'.format(np.round(self.omega_min,3)),color='m',linestyle='--')
    #             ax[0].axvline(self.omega.max(), label='Upper spatial frequency: {} rad/m'.format(np.round(self.omega_max,3)),color='m',linestyle='--')
    #             ax[1].axvline(self.omega.min(), label='Lower spatial frequency: {} rad/m'.format(np.round(self.omega_min,3)),color='m',linestyle='--')
    #             ax[1].axvline(self.omega.max(), label='Upper spatial frequency: {} rad/m'.format(np.round(self.omega_max,3)),color='m',linestyle='--')            

            
    #         ax[0].set_xlabel('Spatial angular frequency (rad/m)')
    #         ax[0].legend() 
            
    #         ax[1].set_xlabel('Spatial angular frequency (rad/m)')    
    #         ax[1].legend()
        
    #     if x_axis == 'wavelength':
            
    #         for idx,item in enumerate(self.class_list):
    #             ax[0].plot(self.wave,self.vert_irreg[idx],label='FRA class {}'.format(self.class_list[idx]))
    #             ax[1].plot(self.wave,self.lat_irreg[idx],label='FRA class {}'.format(self.class_list[idx]))
                
    #         if v_line == True:
    #             ax[0].axvline(self.wave.min(), label='Lower FRA wavelength: {} m'.format(self.L_min),color='m',linestyle='--')
    #             ax[0].axvline(self.wave.max(), label='Upper FRA wavelength: {} m'.format(self.L_max),color='m',linestyle='--')
    #             ax[1].axvline(self.wave.min(), label='Lower FRA wavelength: {} m'.format(self.L_min),color='m',linestyle='--')
    #             ax[1].axvline(self.wave.max(), label='Upper FRA wavelength: {} m'.format(self.L_max),color='m',linestyle='--') 
                
    #         ax[0].set_xlabel('Wavelength (m)')
    #         ax[1].set_xlabel('Wavelength (m)')      
        
        
    #     if scale=='log':
    #         ax[0].set_xscale('log')
    #         ax[0].set_yscale('log')  
    #         ax[1].set_xscale('log')
    #         ax[1].set_yscale('log')
            
    #     ax[0].legend()         
    #     ax[1].legend()   
          
    #     if plot_signal==True:
                    
    #         f,welch_coef = FRA_irregularities.Welch_PSD_signal(self)
            
    #         if self.signal_type =='vert':
    #             ax[0].plot(1/f,welch_coef,label = 'Signal PSD')
    #         elif self.signal_type=='lat':
    #             ax[1].plot(1/f,welch_coef,label = 'Signal PSD')
    #         else:
    #             print('Only vertical and lateral irregularities can be plotted')
        
    #     ax[0].legend()         
    #     ax[1].legend() 
        
        
    def Welch_PSD_signal(self,window_size_frac=0.2,overlap_frac=0.5):
               
        segment_size = np.int32(window_size_frac*len(self.signal)) 
        fft_size = 2 ** (int(np.log2(segment_size)) + 1) # round up to next highest power of 2 - used for zero padding
    
        overlap_size = overlap_frac*segment_size

        f_signal, welch_coef_signal = welch(self.signal, 
                            1/self.dt,
                            nperseg=segment_size, 
                            noverlap=overlap_size,
                            nfft=fft_size,
                            return_onesided=True,
                            scaling='density',
                            detrend='constant',
                            window='hann',
                            average='mean')
        
        return f_signal, welch_coef_signal    