function [RecievedMat, SentMat] = QKD( ypixel, xpixel)
% Written by Bienvecue Ndagano and refined by Isaac Nape 
% Creates an  ypixel by xpixel matrix for encrytion (SentdMat)  and decrytion (RecievedMat)

%slice is the subset of data from the crosstalk matrix 
%path = 'C:\Users\NMashaba\Google Drive\PHD\Written works\2018\Bessel Self healing paper\QKD_withScalarModes\Matrices For HD selfhealing experiment Eileen\19072018\19072018\'; % file path
file = 'A1All2'; % Crosstalk Matrix file
crosstalk = load([file '.mat']); % load the crosstalk data
slice=1;
% Select a given crosstalk matrix in the file by uncommenting the
% section below
%slice = 2;
data = crosstalk.(file);
rawcrosstalk = data(:,:,slice);
%rawcrosstalk(rawcrosstalk<0)=0;
% OR average over all the measurements uncomment line below
% rawcrosstalk = mean(crosstalk.(file),3);

s = size(rawcrosstalk); % Find the size of the crosstalk array

% Split the crosstalk array into detection is different bases for Alice and
% Bob
normcrosstalk1 = rawcrosstalk(1:s(1)/2,1:s(1)/2);
normcrosstalk2 = rawcrosstalk(1:s(1)/2,s(1)/2+1:s);
normcrosstalk3 = rawcrosstalk(s(1)/2+1:s,1:s(1)/2);
normcrosstalk4 = rawcrosstalk(s(1)/2+1:s,s(1)/2+1:s);

% Example of ideal crosstalk data
% normcrosstalk1 = eye(s(1)/2);
% normcrosstalk2 = ones(s(1)/2,s(1)/2)./(s(1)/2);
% normcrosstalk3 = ones(s(1)/2,s(1)/2)./(s(1)/2);
% normcrosstalk4 = eye(s(1)/2);

% Normalize the crosstalk into probabilities.
for lauf = 1:s(1)/2
    normcrosstalk1(lauf,:) = normcrosstalk1(lauf,:)./sum(normcrosstalk1(lauf,:));
    normcrosstalk2(lauf,:) = normcrosstalk2(lauf,:)./sum(normcrosstalk2(lauf,:));
    normcrosstalk3(lauf,:) = normcrosstalk3(lauf,:)./sum(normcrosstalk3(lauf,:));
    normcrosstalk4(lauf,:) = normcrosstalk4(lauf,:)./sum(normcrosstalk4(lauf,:));
end
t1 = trace(normcrosstalk1)./sum(normcrosstalk1(:));
t2 = trace(normcrosstalk4)./sum(normcrosstalk4(:));
t = mean([t1 t2]);

if t < 0.89
    warning(['Please note that your average measurement fidelity of ' num2str(round(t,3)) ' is below the 0.89 threshold for BB84'])
else
end
% Plot the normalized crosstalk
%figure; imagesc([normcrosstalk1 normcrosstalk2; normcrosstalk3 normcrosstalk4])
%caxis([0 1]); colorbar;


bases = [1 2]; % Label the bases
states = (1:21);% Label the states within a basis.

% Combine the crosstalk data into one array M(a,b,c,d) where a is the state
% sent, b the state measured, c and d are Alice and Bob's bases,
% respectively.
normcrosstalk = zeros([size(normcrosstalk1) length(bases) length(bases)]);
normcrosstalk(:,:,1,1) = normcrosstalk1;
normcrosstalk(:,:,1,2) = normcrosstalk2;
normcrosstalk(:,:,2,1) = normcrosstalk3;
normcrosstalk(:,:,2,2) = normcrosstalk4;


bitsent = xpixel*ypixel*3; % Number of bits sent
PbasisAlice = ones(size(bases))/length(bases); % Probability of Alice picking a given basis
PstateAlice = ones(size(states))/length(states); % Probability of Alice picking a given state within a basis

AliceBasis = randsample(bases,bitsent,true,PbasisAlice); % Alice picks a basis
Alicestate = randsample(states,bitsent,true,PstateAlice); % Alice picks a state within a basis

BobBasis = randsample(bases,bitsent,true,PbasisAlice); % Bob picks a basis

% Bob measures the state lauf according to crosstalk matrix, state sent by Alice and basis chosen by Bob
for lauf = 1:bitsent
    Bobmeasurement(lauf) = randsample(states,1,true,normcrosstalk(Alicestate(lauf),:,AliceBasis(lauf),BobBasis(lauf)));
end

% Proceeed with sifting by comparing the basis choices
siftbases = mod(AliceBasis+BobBasis,2);
counter = 0;
for lauf = 1:bitsent
    if siftbases(lauf) == 0 % prepare and measure bases match.
        counter = counter+1;
        Sentbit(counter) = Alicestate(lauf); % retain Alice's prepared state
        Receivedbit(counter) = Bobmeasurement(lauf); % retain Bob's measured state
    else
    end
end
% reshape as matrices for keys
figure(3)
RecievedMat = reshape(Receivedbit(1:xpixel*ypixel), [ypixel,xpixel]);
figure(4)
SentMat = reshape(Sentbit(1:xpixel*ypixel), [ypixel,xpixel]);
% bar([Sentbit;Receivedbit]') % Plot the states numbers for visual comparison
disp([num2str(counter) ' states were sifted out of the ' num2str(round(bitsent,3)) ' that were sent.'])
subplot(1,2,1)
imagesc(RecievedMat); colormap(hot);
colorbar
title('Alice key')
axis image
subplot(1,2,2)
imagesc(SentMat); colormap(hot);
colorbar
axis image
title('Bob key')
% Compute the error between Alice's and Bob's keys.
error = 0;
for lauf = 1:length(Receivedbit)
    if Receivedbit(lauf) ~= Sentbit(lauf) % Key values are different
        error = error +1;
    else
    end
end
ErrorFraction = error./length(Receivedbit); % Probability of error in sifted key
disp(['The probability of errors in the sifted key is: ' num2str(round(ErrorFraction,3))])
save('QKD2020Data2.mat', 'RecievedMat', 'Sentbit', 'SentMat', 'Receivedbit', 'ErrorFraction' );