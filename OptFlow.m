close all;
clear;
clc;
%%%%% Opt-Flow but based on Images

%% Pre-processing

images_path = 'C:/Users/Jacopo/Desktop/Weizmann/Movies/Long One';
output_path = 'C:/Users/Jacopo/Desktop/Weizmann/Code/Optical Flow/MatlabOut/';
lt = dir(images_path);
lt = lt(3 : end, :);
lt = struct2table(lt);
img_idx = [];

for i = 1 : height(lt)
    if endsWith(lt.name(i), '.tif')
        img_idx = [img_idx, i];
    end
end

lt = lt(img_idx, :);

%%

flowFarn = opticalFlowFarneback('NumPyramidLevels', 4);

h = figure;

j = 0;
nFrames = height(lt);
loop = true;

while loop
    j = j + 1;

    img = append(images_path, '/', char(lt.name(j)));
    frRGB = read(Tiff(img, 'r'));
    frGray = im2gray(frRGB);  
    flow = estimateFlow(flowFarn, frGray);
    
    imshow(frRGB)
    hold on
    plot(flow, DecimationFactor = [21, 21], ScaleFactor = 25)
    f = get(gca, 'Children');
    set(f(1), 'Color', 'y')
    hold off
    
    name_save = append(output_path, sprintf('%04d', j), '.png');
    ax = gca;
    exportgraphics(ax, name_save ,'Resolution',300)
 
    if j == nFrames
        loop = false;
    end
end