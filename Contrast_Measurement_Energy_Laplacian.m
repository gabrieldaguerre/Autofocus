function Contrast = Contrast_Measurement_Energy_Laplacian(Image)
        LAP= fspecial('laplacian');
        C=imfilter(Image,LAP,'replicate','conv');
        Contrast = mean2(C.^2);
end