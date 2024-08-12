# lowfield_seq_development
This project aims at developing 
    (1) a testing suite for pypulseq sequence development including 
        (a) building  different sequences, 
        (b) analyzing sequences (analyze_seq.py) with regard to gradient moments and slew rate, PNS, SAR and mechanical resonances, 
        (c) simulating sequences using MRzero package (simulate_seq.py), 
        (d) plot results (plot_sim.py),
        (e) exporting sequences to Virtual Machine Siemenes IDEA for VE11C sequence simulation using poet
    (2) further develop 3D TSE sequence for Ultra Lowfield MRI at PTB
    (3) compare 3D TSE performance between 3T and ULF
    (4) develop QA procedure for ULF in comparison to 3T

Major sequence adaptions
    - added crusher gradients around 180 refocussing pulses
    - added crusher after last 180 ref pulse of echo train
    - changed phase of 90 and 180 pulse to pi/2 and 0
    - added oversampling factor for ADC readout and replaced full sampling of 64k ADC samples
    - added linear readout in phase encoding direction 1 and sequential (slice-wise) readout in phase encoding 2
    - added calculation of minimal ESP, max ETL and recommended ETL
