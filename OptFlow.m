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
clc;
%%

flowFarn = opticalFlowFarneback('NumPyramidLevels', 4);

h = figure;

j = 0;
nFrames = height(lt);
loop = 1;

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
    
    if j == nFrames
        loop = 0;
    end

    pause(0.1)
end