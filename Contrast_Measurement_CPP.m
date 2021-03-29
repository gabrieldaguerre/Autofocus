function Contrast = Contrast_Measurement_CPP(Image)
% image gray
% contrast per pixel https://uk.mathworks.com/matlabcentral/answers/81133-how-to-calculate-contrast-per-pixel-cpp-of-an-image#:~:text=To%20measure%20the%20quality%20of,pixel%20and%20its%20adjacent%20pixel.


        kernel = [-1, -1, -1, -1, 8,-1, -1, -1]/8;
        difffocus_image = conv2(double(Image),kernel,'same');
        Contrast = mean2(difffocus_image);
end