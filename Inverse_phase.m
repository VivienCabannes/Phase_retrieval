name1 = 'images/Ek_Balam.jpg';
name2 = 'images/Valladolid.jpg';

% Read 2 images
X1 = imread(name1);
X2 = imread(name2);
imwrite(X1, 'image1.jpeg')
imwrite(X2, 'image2.jpeg')

% Perform their fourier transform
F1 = fft2(X1);
F2 = fft2(X2);

% Inverse their phases
phase1 = F1 ./ abs(F1);
phase2 = F2 ./ abs(F2);
Ff1 = abs(F1).*phase2;
Ff2 = abs(F2).*phase1;

% Inverse Fourier
iX1 = ifft2(Ff1);
iX2 = ifft2(Ff2);
imwrite(uint8(iX1), 'image1_phase2.jpeg')
imwrite(uint8(iX2), 'image2_phase1.jpeg')