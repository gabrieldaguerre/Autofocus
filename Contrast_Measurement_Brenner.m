function Contrast = Contrast_Measurement_Brenner(Image)
    [M,N] = size(Image);
    DH = zeros(M,N);
    DV = zeros(M,N);
    DV(1:M-2,:) = Image(3:end,:)-Image(1:end-2,:);
    DH(:,1:N-2) = Image(:,3:end)-Image(:,1:end-2);
    FM=max(DH,DV);
    FM=FM.^2;
    Contrast = mean2(FM);
end
    
