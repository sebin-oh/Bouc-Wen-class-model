function fitness = GA_fitness(pop, disp, force, up_down_limit)
    
    fitness = zeros([size(pop,1),1]);  
    
    
    for ii = 1:size(pop,1)
        % fprintf('%d\n',ii);
        
        x = pop(ii,:);
        
        % apply upper and lower limits of the populations
        for jj =1:size(x,2)
            if x(jj) < up_down_limit(jj,2)
                x(jj) = up_down_limit(jj,2);
            elseif x(jj)>up_down_limit(jj,1)
                x(jj) = up_down_limit(jj,1);
            end
        end
        
        results_iBWBN = BoucWen(x, disp);

        % Fitness function: mean squred error
        fitness(ii) = norm(force-results_iBWBN(:,1));
    end
    
end