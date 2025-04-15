% Developed by Taeyong Kim at SNU
% Jan 31, 2020

% Selecting the best individuals in the current generation as parents for 
% producing the offspring of the next generation.

function parents = GA_mating(pop, fitness, num_parents, up_down_limit)
       
    best_fitness = min(fitness);
    
    pop(:,14) = 1./pop(:,14);
    % Estiamte ratio of parameters
    idx = [6;7;9;10;11;12;13;14];
    fitness2 = zeros([length(fitness),1]);
    for taeyong = 1:size(pop,1)
        tmp = pop(taeyong,:);
        tmp_sum = 0;
        for kim = 1:length(idx)
            tmp_idx = idx(kim);
            tmp_sum = tmp_sum + (tmp(tmp_idx)-up_down_limit(tmp_idx,2))...
                       /(up_down_limit(tmp_idx,1)-up_down_limit(tmp_idx,2));
        end
        fitness2(taeyong,1) = tmp_sum;
    end
    pop(:,14) = 1./pop(:,14);
    
    % FInd the parent based on both fitness1 and fitness2
    fitness11 = fitness;
    tmp_parents = [];
    for taeyong = 1:size(pop,1)
        max_fitness2_idx = find(fitness2 == min(fitness2));
        if length(max_fitness2_idx)>=2
            max_fitness2_idx = max_fitness2_idx(1);
        end
        
        if abs(fitness11(max_fitness2_idx)-best_fitness)/best_fitness < 0.2
            tmp_parents = [tmp_parents; pop(max_fitness2_idx,:)];
        end
        fitness2(max_fitness2_idx) = 99999999999;
        fitness11(max_fitness2_idx) = 99999999999;
    end
    
    if size(tmp_parents,1)>=ceil(num_parents/2)
        tmp_parents = tmp_parents(1:ceil(num_parents/2),:);
        tmp_parents2 = [];
        for parent_num = 1:num_parents-size(tmp_parents,1)

            max_fitness_idx = find(fitness == min(fitness));
            if length(max_fitness_idx)>=2
                tmp_parents2 = [tmp_parents2; pop(max_fitness_idx(1),:)];
            else
                tmp_parents2 = [tmp_parents2; pop(max_fitness_idx,:)];
            end
            fitness(max_fitness_idx) = 99999999999;

        end
        
        parents = [tmp_parents; tmp_parents2];        
    else
        tmp_parents2 = [];
        for parent_num = 1:num_parents-size(tmp_parents,1)

            max_fitness_idx = find(fitness == min(fitness));
            if length(max_fitness_idx)>=2
                tmp_parents2 = [tmp_parents2; pop(max_fitness_idx(1),:)];
            else
                tmp_parents2 = [tmp_parents2; pop(max_fitness_idx,:)];
            end
            fitness(max_fitness_idx) = 99999999999;

        end
        
        parents = [tmp_parents; tmp_parents2];
    end
    
end
