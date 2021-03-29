function Contrast = Contrast_Measurement_NormedVar(Image)
        A= size(Image);
        W = A(1,2);
        H = A(1,1);
        b=0;
        mean = mean2< (Image);
        for j=1:W
            for i=1:H
                p = double(Image (i,j));
                t = (p-mean)*(p-mean);
                b = b+ t;
            end
        end
        Contrast = b/(H*W*mean);
end