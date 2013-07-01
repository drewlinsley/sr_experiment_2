function noise_array = generate_noise(params,stims)
%generate array of noise images... turn these into 1/f noise eventually
noise_array = cell(numel(stims.db_seq),1);
phase = 1i*random('Uniform',0,2*pi,params.pic_size(1),params.pic_size(2));  % Random uniform distribution 0 - 2pi
x = ones(params.pic_size(1),1) * (-params.pic_size(2)/2 : (params.pic_size(2)/2 - 1)); % x = x/(cols/2);
y = (-params.pic_size(1)/2 : (params.pic_size(1)/2 - 1))' * ones(1,params.pic_size(2));% y = y/(rows/2);

radius = sqrt(x.^2 + y.^2);         % Matrix values contain radius from centre.
radius(params.pic_size(1)/2+1,params.pic_size(2)/2+1) = 1;      % .. avoid division by zero.

amp = 1./(radius.^1);          % Construct the amplitude spectrum
amp = fftshift(amp);

phasemod = round(fftshift(radius.^.1 + 1));

phasechange = 2*pi*((random('unid',4+1,params.pic_size(1),params.pic_size(2)) -1 - 4/2) .* phasemod );

dphase = 1i*phasechange/(100-1);   % premultiply by i to save time in th eloop

for m = 1:numel(noise_array),
    %noise_array{m} = ((picSize(1)*picSize(2))*randn(picSize(1), picSize(2)) + 255);
    %noise_array{m} = abs(spatialPattern(picSize, -1)).*255; %make some pink noise...
    phase = phase + dphase;
    newfft =  amp .* exp(phase);
    noise_array{m} = mat2gray(real(ifft2(newfft)));
    noise_array{m} = gray2ind(noise_array{m},256);
end
