% Developed by Taeyong Kim at SNU
% Mar 3, 2020

% Performing muation using random numbers

% offspring: the chromosome of offspring


function crossover = GA_mutation(offspring, learning_rate, jj, up_down_limit)
    
    % Mutation changes a single gene in each offspring randomly
    % The random value to be added to the gene.
    
    offspring2 = [];
    
    for ii = 1:size(offspring,1)
        
        if rand<0.3
            
            for kk = 1:size(offspring,2)
                
                x(1,kk) = offspring(ii,kk) + 2*(rand-0.5)* ...
                          min([up_down_limit(kk,1)-offspring(ii,kk),offspring(ii,kk)-up_down_limit(kk,2)]);
            end

        else
        
            x = offspring(ii,:)+2*(rand(1,size(offspring,2))-0.5).*offspring(ii,:) ...
                                *exp(-(jj-1)*learning_rate);
        end
                
        % apply upper and lower limits of the populations
        for jj =1:size(x,2)
            if x(jj) < up_down_limit(jj,2) % period
                x(jj) = up_down_limit(jj,2);
            elseif x(jj)>up_down_limit(jj,1)
                x(jj) = up_down_limit(jj,1);
            end
        end
        
        offspring2 = [offspring2; x];
    end
    
    crossover = offspring2;
    
end
