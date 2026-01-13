SyntheticHaystack.m is an idealized model developed in MATLAB to simulate a sub-bottom profiler operating at a constant height, passing over a collection of resonators positioned on a flat seafloor. To represent an irregular collection of resonating objects, a specified number of resonators were randomly assigned a resonance frequency, Q-factor, resonance amplitude coefficient, and spatial position within specified ranges. The model places the SBP along a sequence of discrete positions along the range direction, corresponding to successive pings acquired at a constant platform speed and ping rate. At each position, the SBP transmits a downward-looking chirp with a user-specified beam width, centered at normal-incidence. The water column is assumed to be homogeneous, with a constant sound speed. The seafloor is represented as a flat, planar interface characterized by a uniform reflection coefficient. A set of discrete lithic resonators are positioned on the seafloor over a finite horizontal extent.  Each resonator produces an angular-dependent acoustic response when excited by the incident chirp. Resonator responses are strongest at high incidence angles (e.g., phi = ~55 deg) and weaker at near-normal angles (phi = 0 deg), reflecting directional resonance behavior. As the SBP passes over the resonator field, individual resonators enter and exit the SBP beam footprint, contributing time-delayed resonant energy to the received signal depending on SBP position.  For each step along the SBP path, the following calculations are made.

1. The outgoing pulse, a chirp beginning at time t=0, is generated. The chirp is defined as a linear frequency sweep between two specified frequencies with a Gaussian envelope and specified pulse duration. The SBP beam pattern is a function of beam width sbp_bw and polar angle phi. The beam pattern is defined as a parabola with maximum amplitude of 1 at phi = 0 deg and 0.5 at phi = sbp_bw/2, where sbp_bw = the SBP beam width.
2. The normal-incidence backscattered wave from the seafloor is calculated by setting it equal to the outgoing chirp pulse scaled by the reflection coefficient (specified as 0.9) and spherical spreading. The backscattered wave is received starting at time t = 2h_sbp/c_sw, where h_sbp = 10.6 m is the SBP height above seafloor and c_sw = 1500 m/s. 
3. The resonator responses are calculated for each resonator in the model. 
        3.1 Time is calculated as the time that the resonance frequency from the chirp would reach the resonator position if the pulse were heading straight from the SBP to the resonator. 
        3.2 The resonance amplitude is equal to the SBP beam pattern at the direct-incidence angle squared (transmitting and receiving beam pattern are identical due to the principal of reciprocity) times the specified resonator amplitude coefficient times the resonator beam pattern. The resonator beam pattern is a predefined function of polar angle phi. Spherical spreading is also taken into account on the incident and return trips. 
        3.3 The resonance response is created by multiplying the outgoing chirp pulse by a gaussian window of the calculated amplitude centered around the time at which the resonance frequency from the chirp would reach the resonator position and get back to the SBP. The length of the gaussian is a function of resonator Q.
4. The received signal is set equal to the scattered wave plus the sum of the waves from each resonator. 
5. The final signal} is  calculated from the outgoing pulse cross-correlated with the received signal (backscattered wave + resonator responses).

Problem Setup: the user may specify the following

save_run        = 1 to save | 0 no save. saves settings, resonator characteristics and seismic plot; 
save_path       = path for saving output; 
f_s             = sampling frequency; 
chirp_freq1     = start frequency of linear chirp sweep; 
chirp_freq2     = end frequency of linear chirp sweep; 
chirp_length    = duration of linear chirp sweep; 
chirp_alpha     = alpha of chirp gaussian created by gausswin; 
chirp_t         = time array of outgoing chirp; 
receive_length  = duration of sbp receiving; 
receive_t       = time array of received signal receive_t=0 at beginning of outgoing pulse; 
sbp_dx          = step size for sbp; 
sbp_x           = array defining horizontal path of sbp; 
sbp_height      = sbp height above seafloor; 
sbp_beamwidth   = angle of half-max beam pattern from the sbp; 
n_resonators    = number of resonators; 
res_freq_min    = minimum resonance frequency; 
res_freq_max    = maximum resonance frequency; 
res_amps        = resonator amplitudes relative to the transmit pulse; 
res_depths      = resonator depths below seafloor; 
res_xs          = x positions of resonators from start of chirp path; 
water_c         = speed of sound in water; 
sediment_reflection_coeff = sediment reflection coefficient; 

Output:
text file containing resonance characteristics, 
text file containing model settings, 
2D sbp profile, 

