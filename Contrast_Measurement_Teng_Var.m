function Contrast = Contrast_Measurement_Teng_Var(Image)

        Sx = fspecial ('sobel');
        Gx = imfilter(double(Image), Sx, 'replicate', 'conv');
        Gy = imfilter(double(Image), Sx, 'replicate', 'conv');
        C = Gx.^2 + Gy.^2;
        Contrast = std2(C)^2;
end