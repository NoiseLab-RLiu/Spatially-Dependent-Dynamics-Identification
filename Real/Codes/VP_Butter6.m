%%  Please load the dataset from "Data" folder before running this code
%cfreq = 3;
[b,a] = butter(6,[cfreq*10000-1000 cfreq*10000+1000]/150000,'bandpass');
dx = 1e-3;
dy = dx;
dt = 1/3e5;
% Rebuild VP in time domain
% load("Xsus.mat");
X = zeros(100,100,3001);
for i=1:100
    for j=1:100
        X(i,j,:) = ifft(Xsus(:,i,j));
    end
end


%% ADD NOISE
rng(0);
noise = sqrt(100)*randn(size(X));
Xn = X+noise;

% 
U_filt = zeros(size(X,1),size(X,2),size(X,3));
for i=1:size(X,1)
    for j=1:size(X,2)
        dataOut = filter(b,a,squeeze(Xn(i,j,:)));
        U_filt(i,j,:) = dataOut;
    end
end

