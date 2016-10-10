function trajectory = calculate_det_traj(time_span,nn_init,input_filename)

pp = read_input_parameters_2(input_filename);

[tt nn] = ode45(@calculate_det_traj_core,time_span,nn_init);
trajectory = [tt nn];

    function dn = calculate_det_traj_core(tt,nn)
        
        number_species = length(nn);
        dn = zeros(number_species,1); % must be a column vector
        
        for species = 1:number_species
            growth_rate = (-1+2*pp.rr)*pp.ss;
            dn(species) = ...
                nn(species)*growth_rate;  % auto-catalytic growth
%                 - nn(species)*(1+growth_rate)*pp.uu(species); % mutational outflux
            if species > 1
                dn(species) = dn(species) + nn(species-1)*pp.uu(species-1); % mutational influx
            end
        end        
    end
end
