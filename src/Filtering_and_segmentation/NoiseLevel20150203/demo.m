img = double(imread('lena.png'));
mskflg = 0; % if you want to produce the mask, put one for mskflg

img = rgb2gray(img);

level = [5,10,20,40];
for i=1:size(level,2);
 noise = img + randn(size(img)) * level(i);
 tic;
 [nlevel th] = NoiseLevel(noise);
 t=toc;
 fprintf('True: %5.2f  R:%5.2f G:%5.2f B:%5.2f\n', level(i), nlevel(1), nlevel(1), nlevel(1) );
 fprintf('Calculation time: %5.2f [sec]\n\n', t );
 
 if( mskflg )
   msk = WeakTextureMask( noise, th );
   imwrite(uint8(msk*255), sprintf('msk%02d.png', level(i)));
 end
end
