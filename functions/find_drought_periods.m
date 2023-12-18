function [duration, intensity, frequency, drought_start_end] = find_drought_periods(spi, nmonths_end_drought, min_drought_duration, min_intensity)

    % Initialize variables
    start_dates = [];
    end_dates = [];
    intensity = [];

    drought = 0;
    p_count = 0;
    n_count = 0;
    % Iterate through SPI time series
    for i = 1:length(spi)
        if i == 500
            i;
        end


        if drought == 0 % it's not an ongoing drought
            if spi(i) < 0 % SSI is negative
                if n_count == 0
                    start = i;
                end
                n_count = n_count + 1;
                p_count = max(0, p_count-1);
                if n_count >= min_drought_duration % Start of drought period
                    start_dates(end+1) = start;
                    drought = 1;
                    p_count = 0;
                    n_count = 0;
                end
            else % SSI is positive
                p_count = p_count+1;
                n_count = max(0, n_count-1); 
                if p_count > 6
                    n_count = 0;
                    p_count = 0;
                end
            end
        else % ongoing drought
            if p_count == 0
                endd = i;
            end
            if spi(i) >= 0
                p_count = p_count + 1;
                
                if p_count >= nmonths_end_drought % End of drought period
                    end_dates(end+1) = endd-1;
                    drought = 0; 
                    p_count = 0;
                    n_count = 0;
                end
            else
                p_count = max(0, p_count-1); 
            end
        end
    end


    

    if length(start_dates) > length(end_dates)
        end_dates(end+1) = length(spi);
    end

    for i = 1:length(start_dates)
        intensity(i) = mean(spi(start_dates(i):end_dates(i)) );
    end


    % check on intensity
    if isempty(start_dates)
        frequency = 0;
        duration  = 0;
        intensity = 0;
        drought_start_end = [];
    else

        check_on_intensity = intensity<min_intensity;
    
        drought_start_end = [start_dates', end_dates'];
    
        drought_start_end = drought_start_end(check_on_intensity, :);
        frequency = sum(check_on_intensity);
        duration  = drought_start_end(:,2) - drought_start_end(:,1);
        intensity = intensity(check_on_intensity);
        
        if isempty(intensity)
            frequency = 0;
            duration  = 0;
            intensity = 0;
            drought_start_end = [];
    
        end
    end
