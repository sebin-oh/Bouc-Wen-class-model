% Developed by Taeyong Kim at SNU
% Jan 31, 2020

% Performing crossover using the parents chromosome

% parents: parent chromosome
% offspring_size: size of the off spring


function offspring = GA_crossover(parents, offspring_size)
    
    offspring = zeros([offspring_size, size(parents,2)]);
    % The point at which crossover takes place between two parents. 
    % Usually it is at the center.
    crossover_point = round(size(offspring,2)/2);
    
    for ii = 1:offspring_size
        
        % Index of the first parent to mate.
        parent1_idx = ceil(rand(1)*(size(parents,1)-1));%rem(ii,size(parents,1));
        % Index of the second parent to mate.
        parent2_idx = ceil(rand(1)*(size(parents,1)-1));%rand();rem(ii+1,size(parents,1));
        
        if parent1_idx == 0
            parent1_idx = size(parents,1);
        end
        if parent2_idx == 0
            parent2_idx = size(parents,1);
        end
        
        % The new offspring will have its first half of its genes taken from the first parent
        offspring(ii, 1:crossover_point) = parents(parent1_idx, 1:crossover_point);
        % The new offspring will have its second half of its genes taken from the second parent.
        offspring(ii, crossover_point+1:end) = parents(parent2_idx, crossover_point+1:end);
        
    end
    

end

