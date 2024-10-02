clear all;

%% make movie with all the fields 

%Get the data
data = fileread('data/IF FILES/defectors06_07_2023 13_31_125.txt');
matches = regexp(data, '[-+]?\d+(\.\d+)?([eE][-+]?\d+)?', 'match');
data_new = str2double(matches);

%Data features
n = length(data_new);
data = reshape(data_new, 3, n/3)';
nb_occurence = histc(data(:,1), unique(data(:,1)));
Nt = size(nb_occurence,1);
makevideo = false;

if makevideo
    %First plot
    height_max = max(data(:,3));
    height_min = min(data(:,3));
    
    figure();
    x = data(1:nb_occurence(1), 2);
    y = data(1:nb_occurence(1), 3);
    p1=plot(x, y);
    caxis([height_min, height_max])
    xlim([min(data(:,2)), max(data(:,2))]);
    ylim([min(data(:,3)), max(data(:,3) + 10)]);
    %Video maker
    fps = 10;
    video_name = 'data/IF FILES/defectors06_07_2023 13_31_125_video';
    v = VideoWriter(video_name, 'MPEG-4'); 
    v.FrameRate = fps;
    frame = getframe(gcf);
    open(v);               
    
            for i=1:Nt-1
               i
               data_i = data(sum(nb_occurence(1:i))+1 : sum(nb_occurence(1:i)) + nb_occurence(i+1), 2:3 );
               [x, index] = sort(data_i(:,1));
               y = data_i(:,2);
               plot(x,y(index)); 
               
    
               frame = getframe(gcf);
               writeVideo(v,frame);
            end
    
            close(v);
end
