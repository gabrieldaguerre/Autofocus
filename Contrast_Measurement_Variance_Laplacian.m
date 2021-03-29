function Contrast = Contrast_Measurement_Variance_Laplacian(Image)
        LAP= fspecial('laplacian');
        ILAP=imfilter(Image,LAP,'replicate','conv');
        Contrast = std2(ILAP)^2;
end