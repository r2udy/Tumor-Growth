clear all;

% Get data's file names
cd data\'IF Files'\
fields = dir();
list_filename_intermd = [fields.name];
list_filename = regexp(list_filename_intermd, '[d][a-z]+\d+?_\d+?_\d+? \d+?_\d+?_\d+?.txt', 'match');

%% Paramters 
% Create video
%0 -> no video
%1 -> video
%2 -> video of orignal data
makevideo = 2;
% Save data
%0 -> no saving
%1 -> classic saving
%2 -> h5py compatible saving
data_saving = 2;
% Filter
filter_activ = 0;

% for k = 1:length(list_filename)/3
k = 1;
    %% Get the data
    data = fileread(string(list_filename{k}));
    matches = regexp(data, '[-+]?\d+(\.\d+)?([eE][-+]?\d+)?', 'match');
    data_new = str2double(matches);
    
    %% Data features
    n = length(data_new);
    data = reshape(data_new, 3, n/3)';
    nb_occurence = histc(data(:,1), unique(data(:,1)));
    Nt = size(nb_occurence,1);
    
    %% Extrapolation
    [max_value, i_max] = max(nb_occurence);
    data_i_max = data(sum(nb_occurence(1:i_max-1))+1 : sum(nb_occurence(1:i_max-1)) + nb_occurence(i_max), 2:3 );
    [x_i_max, ind_sorted] = sort(data_i_max(:, 1));
    n = size(x_i_max,1);
    nb_occurence_loop = [0; nb_occurence];
    
    data_3d = data_i_max(ind_sorted, :);
    
    for i = 2:(Nt+1)
        
        data_j = data(sum(nb_occurence_loop(1:i-1))+1 : sum(nb_occurence_loop(1:i-1)) + nb_occurence_loop(i), 2:3 ); % j = i-1
        [x_j, index, ~] = unique(data_j(:, 1));
        y_j = data_j(index, 2);

        % Butterworth filter
        if filter_activ
            %Spectre of data_j
            n_j = length(x_j);
            freq = (0: n_j-1)/n_j;
            magn_y_j = abs(fft(y_j));
                        
            %Butterworth Filter
            fc = 0.15*n_j;
            fs = n_j;
            [b, a] = butter(9, fc/(fs/2));
            y_j_butter = filter(b, a, y_j);
            
            y_j = y_j_butter;
        end

        % Chebyshev Interpolation
        m = size(y_j,1);
        ks = -1+1/m:2/m:1-1/m; xk = 0.5*(x_i_max(end) - x_i_max(1))*cos(pi*(ks+1)/2) + 0.5*(x_i_max(end) + x_i_max(1));

        %Barycentric interpolation
        c = (-1).^(0:m-1).*sin(pi*(k+1)/2);
        numer = zeros(size(x_i_max));
        denom = zeros(size(x_i_max));
        exact = zeros(size(x_i_max));
        for j = 1:m
            xdiff = x_i_max - xk(j);
            temp = c(j)./xdiff;
            numer = numer + temp*y_j(j);
            denom = denom + temp;
            exact(xdiff==0) = j;
        end

        y_j_cheb = numer./denom; jj = find(exact); y_j_cheb(jj) = y_j(exact(jj));
        
        data_3d = cat(3, data_3d, [x_i_max,y_j_cheb]);
    end
    
    data3d = data_3d(:, :, 2:end); % data3d(:, :, 1) was just for implementing the technique
    heights = squeeze(data3d(:, 2, :)); % isolate the heights on one variable

    %% cleaning data to remove 'nan' values
    if isnan(heights)
        for l = 1:size(heights, 2)
            index_i_nan = find(isnan(heights(:,l)));
            mean_heights = mean(heights(setdiff(1:size(heights,1), index_i_nan),l));
            heights(index_i_nan,l) = mean_heights;
        end
    end
    
    %% removing duplication
    x_interm = data_3d(:, 1, 2:end);
    % space variable
    [x_heights, ix, ~] = unique(squeeze(x_interm(:, :, 1))');
    % height variable 
    heights = heights(ix,:);
    % time variable
    t_heights = 0:10:700;
    
    if data_saving >0 
        if data_saving == 1
            %% save data
            save(sprintf('Surface_growth_' + "filtered" + '_(x,t).mat'), 'heights', 't_heights', 'x_heights');
        elseif data_saving == 2
            mfw = dsp.MatFileWriter(sprintf("Verif_algo" + '_h_t.mat'), 'heights_t');
            mfw(diff(heights'));
            release(mfw)
            mfw = dsp.MatFileWriter(sprintf("Verif_algo" + '_h.mat'), 'heights');
            mfw(heights(:, 1:end-1)');
            release(mfw)
            mfw = dsp.MatFileWriter(sprintf("Verif_algo" + '_x.mat'), 'x_heights');
            mfw(x_heights);
            release(mfw)
            mfw = dsp.MatFileWriter(sprintf("Verif_algo" + '_t.mat'), 't_heights');
            mfw(t_heights(1:end-1)');
            release(mfw)
        end
    end

    %% make movie with all the fields 
    if makevideo>0    
        %First plot
        height_max = max(max(heights)) + 0.1*max(max(heights));
        height_min = min(min(heights));
        
        figure();
        x = x_heights;
        y = heights;
        ylim([height_min, height_max]);
        frame = getframe(gcf);
        %Video maker
        fps = 10;
        %video_name = 'Surface_growth_' + "filtered";
        video_name = string(list_filename{k}) + '_video_clean'
        v = VideoWriter(video_name, 'MPEG-4'); 
        v.FrameRate = fps;
        open(v);               
        if ~makevideo>2
                for i=1:Nt
                    plot(x, y(:, i))
                    ylim([height_min, height_max]);
                    frame = getframe(gcf);
                    writeVideo(v,frame);
                end
        else
            height_max = max(max(data(:,3))) + 0.1*max(max(data(:,3)));
            height_min = min(min(data(:,3)));
            for i=2:Nt+1
                    data_j = data(sum(nb_occurence_loop(1:i-1))+1 : sum(nb_occurence_loop(1:i-1)) + nb_occurence_loop(i), 2:3 ); % j = i-1
                    [x_j, index, ~] = unique(data_j(:, 1));
                    y_j = data_j(index, 2);
                    plot(x_j, y_j);
                    ylim([height_min, height_max]);
                    ylabel('Heights of invasive chains');
                    xlabel('1-D field, outest layer');
                    frame = getframe(gcf);
                    writeVideo(v,frame);
            end
        end
        
                close(v);
    end
% end


