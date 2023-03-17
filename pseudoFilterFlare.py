# multi-band flare i.d. design

# sources: 
# (1) NOAA database process: https://www.mdpi.com/1996-1073/9/1/14 
# (2) Original VNF: https://www.mdpi.com/2072-4292/5/9/4423 

# m10 raw VIIRS detection -> (1) if multiband dect then planck fit -> find T and RH (2) else: no temp, use 10 km distance sort -> 15 arc second persistence -> 10km performance on remaining
# (1) T > 1400K assume true flare (2) < 1400k industrial sites (3) only M10 and M10 not calcd

def flare_processing_main(raw_VIIRS_SDR, database_path):
    
    # calculate temperature and radiant heat per VIIRS detection
    detections_multib, detections_single = VNF_temperature_processing(raw_VIIRS_SDR)
    # filter results to gas flares by T and time persistence
    filter_vnf_results(detections_multib, detections_single)
    # append results to database 
    
    return database_path


def VNF_temperature_processing(raw_VIIRS_SDR):
    # apply VNF ID using M7, M8, M10, M12 and M13
    # requires peak radiant emission ID
    
    for detection in raw_VIIRS_SDR
        
        num_bands_tempt = numbands(detection)
        if num_bands_tempt == 1:
            # cannot calc single band with plank, skip
            detections_single.append()
        else:
            # M10 is strongest channel detection
            
            # Candidate hot pixels are identified as those with digital numbers exceeding the mean plus four standard deviations.
            
            # Planck Curve Fitting
            # TODO
            
            # calculate T and RH using Planck fitting
            
            # Wien’s Displacement Law 
            # displacement constant
            b = 2897.8 # K microm
            # "The M10 spectral band records the peak radiant emissions from the typical gas flare"
            lambda_max = # wavelength of peak radiant emissions
            # temperature (K)
            T = b / lambda_max

            # ε = emission scaling factor 
            # ratio btw observed radiances vs radiances for an object at that temperature filling the entire field of view. 
            # S = ε * size of pixel footprint
            # sigma = Stefan–Boltzmann constant
            RH = sigma*S*T^4
            
            # append with calculated values
            detections_multib.append(detection(T,RH))
        

    
    return detections_multib, detections_single



def filter_vnf_results(detections_multib, detections_single)
    
    # < 1400 falls into biomass range
    if temperature_per_viirs_dect < 1400: # in Kelvin
        temperature_per_viirs_dect.remove()
        
    elif temperature_per_viirs_dect is None:
        # M10 no calculation case
        no_temp.append(temperature_per_viirs_dect)
    
    # 1500 K to 2000 K ideal T range
    # any insecurity can be verified w/ persistence
    else:
        binned = bin_by_date_time(temperature_per_viirs_dect)
        persisted_flares, unpersisted_flares = persistance_det(binned)
        
        # for unper in unpersisted_flares:
            # filter additional outliers...
    return persisted_flares, unpersisted_flares, no_temp
        
def persistance_det(sorted_potential_flares):
    # global 15 arc second grid 
    # tallying the number of times an M10 detection occurred during the year.
    
    for detection in sorted_potential_flares:
        # is a same location flagged in a different date-time bin?
        if detection.coordinates() in sorted_potential_flares.bins():
            persisted_flares.append(detection)
        else:
            unpersisted_flares.append(detection)
    
    
    return persisted_flares, unpersisted_flares
        
