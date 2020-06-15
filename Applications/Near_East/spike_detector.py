def find_spikes(age_input, intensity_input, window_size=10, intensity_threshold=5, dFdt_threshold=0.6, duration_dFdt_threshold=0.12):
    from scipy.signal import find_peaks
    import pandas as pd
    import numpy as np
    
    # use a rolling-averaged timeseries to find the spikes.
    
    # intensity:
    intensity = pd.Series(intensity_input).rolling(window_size).mean().to_numpy()
    # remove first window_size-1 elements (as Nan)
    intensity = intensity[window_size-1:]
    
    # but we also need to make sure that the ages are of the correct length. We can take the same average as all ages are equally spaced...
    age = pd.Series(age_input).rolling(window_size).mean().to_numpy()
    # remove first window_size-1 elements (as Nan)
    age = age[window_size-1:]
    
    grad = np.gradient(intensity,age)
    # find the peaks
    maxima,_ =  find_peaks(intensity)
    minima,_ =  find_peaks(-intensity)
    
    # initialise output:
    start_age, end_age = [age_input[0]]*len(maxima),[age_input[-1]]*len(maxima)
    spike_number = -1
    
    spikes = np.zeros_like(age, dtype=bool); spikes[:] = False
    
    
    # To find a spike, we consider each local maximum in turn, and scan to the left and right up to the nearest local minima. To the left, for each age (index i), we consider whether it has:
    # (a) a difference of "intensity_threshold" (e.g. 5 micro T) between intensity[i] and that of the local maximum
    # (b) whether dF/dt[i] > "duration_dFdt_threshold" AND dF/dt[i-1] < "duration_dFdt_threshold".
    #     If such a point is found satisfying (a) and (b), then set start_spike_index = i
    
    # Search to the right, looking for:
    # (a) a difference of "intensity_threshold" (e.g. 5 micro T) between intensity[i] and that of the local maximum
    # (b) whether -dF/dt[i] > "duration_dFdt_threshold" AND -dF/dt[i+1] < "duration_dFdt_threshold"
    #     If such a point is found satisfying (a) and (b), then set end_spike_index = i
    
    
    # If no point that satisfies the above is found, then no spike.
    # If between start_spike_index and end_spike_index |dF/dt| > dFdt_threshold, then we have a spike
    
    # only test if there are some local minima which can be used to test the intensity jump:
    if minima.size > 0:
        #
        for thispeak in maxima:
    #   check to see if there is a minimum to the left or right. If not, try the next peak.
            if thispeak > minima.max() or thispeak < minima.min():
                continue

            start_spike_index, end_spike_index = -1,-1
            # find index of local minima closest to each peak on either side.
            left = minima[minima < thispeak].max()
            right = minima[minima > thispeak].min()
            
            for i in range(thispeak, left-1,-1):   #begin at spike and work backwards in time; include the end points.
                if intensity[thispeak] - intensity[i] > intensity_threshold and grad[i] > duration_dFdt_threshold and grad[i-1] < duration_dFdt_threshold:
                    start_spike_index = i
                    break
        
            for i in range(thispeak, right+1,1):   #begin at spike and work forwards in time.
                if intensity[thispeak] - intensity[i] > intensity_threshold and -grad[i] > duration_dFdt_threshold and -grad[i+1] < duration_dFdt_threshold:
                    end_spike_index = i
                    break

            # test to see if spike gradient exceeded:
            if start_spike_index >-1 and end_spike_index > -1:
                if max(np.abs( grad[start_spike_index:end_spike_index+1])) > dFdt_threshold:
                    spike_number += 1
                    start_age[spike_number] = age[start_spike_index]
                    end_age[spike_number]=age[end_spike_index]
                    spikes[thispeak] = True

       
   
    return spike_number+1, start_age[0:spike_number+1], end_age[0:spike_number+1], age[spikes]
