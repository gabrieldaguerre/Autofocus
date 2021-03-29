function Contrast = Contrast_Measurement_Teng(Image)

        Sx = fspecial ('sobel');
        Gx = imfilter(double(Image), Sx, 'replicate', 'conv');
        Gy = imfilter(double(Image), Sx, 'replicate', 'conv');
        C = Gx.^2 + Gy.^2;
        Contrast = mean2(C);
end