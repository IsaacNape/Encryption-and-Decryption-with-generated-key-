
%% images encryption for QKD demo
% generate keys from BVs code... or  load previosly generated codes.
% by simply commenting QKD() and loading QKD2020Data.mat -> specificly generated for the wits image.
% The QKD procedure was written by  Bienvenue Ndagano. Here we use it to generate encryotion keys 
clear
clc
Im = double( rgb2gray(imread('Wits image.png')) ); %% load message 'image'
Im =uint8(Im./max(Im(:)).*249); %% load message 'image'
ypixels = size(Im, 1);
xpixels = size(Im, 2);

% GetKeys
%[SentMat, RecievedMat ] = QKD( ypixels, xpixels); %% generate keys from BVs code
load QKD2020Data; %% for wits image or wait longer by using the QKD code

% convert both image and encrypt... decrypt codes to uint8
ImageImageEncrypted = uint8(zeros(size(Im))) ;
Decoded = uint8(zeros(size(Im))) ;
EncCode = uint8( SentMat./max(SentMat(:)) .*255); %% Encryption Key
DecCode = uint8( RecievedMat./max(RecievedMat(:)).*255); %% Decryption Key

close all
InputImage(:,:) = imresize(Im(:,:),size(EncCode)); % make image same size as encryption code
ImageEncrypted(:,:) = bitxor(InputImage (:,:), EncCode ); %% bitXoRencryption

fig1=figure(5)
imagesc(InputImage);
axis off
axis image
%set(gcf, 'OuterPosition', [100, 100, 900, 400], 'color', 'w')
title('Input image')
print(fig1, 'InputImage.png', '-dpng' ,'-r600' )
imwrite(InputImage,'InputImage.png');
caxis([0, 255])

fig1=figure(1)
imagesc(EncCode);
axis off
axis image
%set(gcf, 'OuterPosition', [100, 100, 900, 400], 'color', 'w')
title('Alices key')
% print(fig1, 'EncCode.png', '-dpng' ,'-r600' )
imwrite(EncCode ,'EncCode.png');
caxis([0, 257])


fig2=figure(2)
imagesc(DecCode);
axis off
axis image
%set(gcf, 'OuterPosition', [100, 100, 900, 400], 'color', 'w')
title('Bob key')
% print(fig2, 'DecCode.png', '-dpng' ,'-r600' )
imwrite(DecCode ,'DecCode.png');
caxis([0, 255])


fig3=figure(3)
imagesc(ImageEncrypted);
axis off
title('Encrypted image')
axis image
%set(gcf, 'OuterPosition', [100, 100, 900, 400], 'color', 'w')
% print(fig3, 'ImageEncrypted.png', '-dpng' ,'-r600' )
imwrite(ImageEncrypted ,'ImageEncrypted.png');
caxis([0, 255])


% image decryption
ImageDecoded(:,:) = bitxor(ImageEncrypted(:,:), DecCode );
ImageDecrypted(:,:) = bitxor(ImageEncrypted(:,:), DecCode);

fig4=figure(4)
imagesc(ImageDecoded);
axis off
axis image
caxis([0, 255])
title('Decrypted image')
%set(gcf, 'OuterPosition', [100, 100, 900, 400], 'color', 'w')
print(fig4, 'ImageDecoded.png', '-dpng' ,'-r600' )
imwrite(ImageDecoded ,'ImageDecoded.png');
