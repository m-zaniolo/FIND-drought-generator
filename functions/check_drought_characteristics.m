function check_drought_characteristics(obs, param)
% this function performs a number of checks on data and user-set parameter
% quality.

if param.target_intensity > param.min_drought_intensity
    warning( ['Target drought intensity is lower than minimum drought intensity. '])
end

if param.target_duration<param.min_drought_duration
    warning( ['Target drought duration is lower than minimum drought duration. '])
end

if param.target_duration*param.target_ndroughts > 0.8*param.nyear_togenerate*12
    warning( ['For better results, the desired number of dry months ' ...
        '(desired_duration*desired_n_droughts) should amount to less than 80% ' ...
        'of length of the time horizon of the synthetic scenario defined in ' ...
        'the param structure. '] )
end

if param.target_intensity>3.0
     warning( ['A desired intensity larger than 3.0 is in most cases not likely to achieve' ...
         'good results.'])
end

if length(obs)<30*12
    warning( ['For SSI calculation, 30 years of data is usually requried but 50 is recommended.'] )
end

