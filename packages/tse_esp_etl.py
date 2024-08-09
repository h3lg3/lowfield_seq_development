# Print min_esp (echo spacing) and max_etl (echo train length) and recommended etl for given PE steps
import numpy as np
import math 

def get_esp_etl(
    tau_1: float = 3e-3,
    tau_2: float = 4e-3,
    tau_3: float = 5e-3,
    echo_time: float = 25e-3,
    T2: int = 100,              # T2 value of main tissue of interest
    n_enc_pe1: int = 64,
) -> dict:
    
    # minimum echo spacing to accomodate gradients
    min_esp = -(min(tau_1, tau_2, tau_3) - echo_time/2)*2   
    min_esp = np.ceil(min_esp*1e3)*1e-3    # round to 1 ms -> new echo time
    
    # sampling duration [ms] till signal drops to 20%
    max_sampling = -math.log(0.2)*T2*1e-3  
               
    # maximum numbers of 180° echoes fitting in sampling duration     
    max_etl = int(np.floor(max_sampling/min_esp))
    
    # ETL that is multiple of n_enc_pe1 and closest to max_etl
    cc = [(i, n_enc_pe1%i, np.abs(max_etl-i)) for i in np.arange(1, 2*max_etl)]
    # Filter the list to find tuples where the second element is zero
    filtered_cc = [item for item in cc if item[1] == 0]
    # Closest ETL value that is multiple of n_enc_pe1
    rec_etl = int(min(filtered_cc, key=lambda x: x[2])[0])     

    return dict({'min_esp':min_esp, 'max_etl':max_etl, 'rec_etl':rec_etl})