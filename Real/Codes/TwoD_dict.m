%%
ind  = 1:1000;
Uused = U_filt(:,:,ind);

%%
% Uused = 10*Waves(11:110,11:110,:);
[N_x,N_y,Mused] = size(Uused);
method = 'FD';

UxUy=zeros(N_x,N_y,Mused);
UUxUUy=zeros(N_x,N_y,Mused);
UxxUyy=zeros(N_x,N_y,Mused);
UtxUty = zeros(N_x,N_y,Mused);
dy=dx;
for i=2:N_x-1
    for j=2:N_y-1
        UxxUyy(i,j,:) = (Uused(i-1,j-1,:)+Uused(i+1,j+1,:)-2*Uused(i,j,:))/(2*dx^2)+(Uused(i-1,j,:)+Uused(i+1,j,:)-2*Uused(i,j,:))/dx^2+(Uused(i-1,j+1,:)+Uused(i+1,j-1,:)-2*Uused(i,j,:))/(2*dx^2)+(Uused(i,j-1,:)+Uused(i,j+1,:)-2*Uused(i,j,:))/dx^2;
        UxUy(i,j,:) = (Uused(i-1,j,:)-Uused(i+1,j,:))/(2*dx)+(Uused(i,j-1,:)-Uused(i,j+1,:))/(2*dx);
        UUxUUy(i,j,:) = Uused(i,j,:).*UxUy(i,j,:);
    end
end


Ux=zeros(N_x,N_y,Mused);
Uy=zeros(N_x,N_y,Mused);
Uxx=zeros(N_x,N_y,Mused);
Uyy=zeros(N_x,N_y,Mused);
for k = 1:Mused
    Utmp=Uused(:,:,k);
    for i = 1:N_y
        Ux(:,i,k)=numder(Utmp(:,i), dx, 1,method);
        Uxx(:,i,k)=numder(Utmp(:,i), dx, 2,method);
    end
    for i = 1:N_x
        Uy(i,:,k)=numder(Utmp(i,:), dy, 1,method);
        Uyy(i,:,k)=numder(Utmp(i,:), dy, 2,method);
    end
    Utmp=Ux(:,:,k);
    for i = 1:size(Ux,1)
        uxy=numder(Utmp(i,:), dy, 1,method);
        Uxy(i,:,k)=uxy;
    end    
end

% time derivatives
Ut=zeros(N_x,N_y,Mused);
Utt=zeros(N_x,N_y,Mused);
for i = 1:N_x
    for j=1:N_y
        Ut(i,j,:) = numder(Uused(i,j,:),dt,1,method);
        Utt(i,j,:) = numder(Uused(i,j,:),dt,2,method);
       % UtxUty(i,j,:) = numder(UxUy(i,j,:),dt,1,method);
    end
end

%time-spatial derivatives
Utx=zeros(N_x,N_y,Mused);
Uttxx=zeros(N_x,N_y,Mused);
for k = 1:Mused
    for j=1:N_y
        Utx(:,j,k) = numder(Ut(:,j,k),dx,1,method);
        Uttxx(:,j,k) = numder(Utt(:,j,k),dx,2,method);
    end
end

Uty=zeros(N_x,N_y,Mused);
Uttyy=zeros(N_x,N_y,Mused);
for k = 1:Mused
    for i=1:N_x
        Uty(i,:,k) = numder(Ut(i,:,k),dy,1,method);
        Uttyy(i,:,k) = numder(Utt(i,:,k),dy,2,method);
    end
end

SinU = sin(Uused);

Ix = 6:size(Uused,1)-5;
Iy = 6:size(Uused,2)-5;
It = 6:size(Uused,3)-5;

%U1 = ones(30,30,200);
D = 9;
Mdic = length(It);
Ndic = length(Ix)*length(Iy);
Phi_tensor = zeros(Mdic,D,Ndic);
for i=1:length(Ix)
    for j=1:length(Iy)
        Phi_tensor(:,1,length(Iy)*(i-1)+j) = squeeze(Ut(Ix(i),Iy(j),It));
        Phi_tensor(:,2,length(Iy)*(i-1)+j) = squeeze(Utt(Ix(i),Iy(j),It));
        Phi_tensor(:,3,length(Iy)*(i-1)+j) = squeeze(Uused(Ix(i),Iy(j),It).*Ux(Ix(i),Iy(j),It));
        Phi_tensor(:,4,length(Iy)*(i-1)+j) = squeeze(Uused(Ix(i),Iy(j),It).*Uy(Ix(i),Iy(j),It));
        Phi_tensor(:,5,length(Iy)*(i-1)+j) = squeeze(Uxx(Ix(i),Iy(j),It));
        Phi_tensor(:,6,length(Iy)*(i-1)+j) = squeeze(Uyy(Ix(i),Iy(j),It));
        Phi_tensor(:,7,length(Iy)*(i-1)+j) = squeeze(Utx(Ix(i),Iy(j),It));
        Phi_tensor(:,8,length(Iy)*(i-1)+j) = squeeze(Uty(Ix(i),Iy(j),It));
        Phi_tensor(:,9,length(Iy)*(i-1)+j) = squeeze(SinU(Ix(i),Iy(j),It));
    end
end

