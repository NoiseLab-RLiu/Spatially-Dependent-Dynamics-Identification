%Uused = H(:,1:1:eind); %Heat eq.
Uused = Waves; %% Burgers
[N_x,Mused] = size(Uused);
method = 'FD';

Ux=zeros(N_x,Mused);
Uxx=zeros(N_x,Mused);

for k = 1:Mused
    Ux(:,k)=numder(Uused(:,k), dx, 1,method);
    Uxx(:,k)=numder(Uused(:,k), dx, 2,method); 
end

% time derivatives
Ut=zeros(N_x,Mused);
Utt=zeros(N_x,Mused);
for i = 1:N_x
        Ut(i,:) = numder(Uused(i,:),dt,1,method);
        Utt(i,:) = numder(Uused(i,:),dt,2,method);
end

% time-spatial derivatives
Utx=zeros(N_x,Mused);
for k = 1:Mused
        Utx(:,k) = numder(Ut(:,k),dx,1,method);
end

Sinu=zeros(N_x,Mused);
for k = 1:Mused
        Sinu(:,k) = sin(Uused(:,k)*dx);
end

U1 = ones(N_x,Mused);

%% ROI for Heat eq.
% Ix = 3:98;
% Iy = 1;
% It=2:1999;
%% ROI for Burgers eq.
Ix = 41:91;
Iy = 1;
It=2:199;
%% Build dictionary
m = length(It);
N = length(Ix);
Phi_tensor = zeros(m,6,N);
for i=1:length(Ix)
        Phi_tensor(:,1,i) = squeeze(Ut(Ix(i),It));
        Phi_tensor(:,2,i) = squeeze(Utt(Ix(i),It));
        Phi_tensor(:,3,i) = squeeze(Uused(Ix(i),It).*Ux(Ix(i),It));
        Phi_tensor(:,4,i) = squeeze(Uxx(Ix(i),It));
        Phi_tensor(:,5,i) = squeeze(Utx(Ix(i),It));
        Phi_tensor(:,6,i) = squeeze(Sinu(Ix(i),It));
end