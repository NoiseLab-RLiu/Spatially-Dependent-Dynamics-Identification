%% Wave Eq.
Waves=zeros(25,200);
x_tmp=0:12;
k=2*pi/12;
Waves(1+6:13+6,1)=-5*(sin(k*x_tmp+pi/2)-1);
figure
plot(Waves(:,1))
Waves(:,2)=Waves(:,1);
alpha = zeros(25,1);
c = 1*ones(25,1);
c(1:5)=1.5;
c(6:10) = 1.4:-0.1:1;
c(16:20) = 1.2:0.2:2;
c(21:end)=2;

c = 1*c;
alpha(16:25)=0.01:0.01:0.1;%50:1/50:1;

dt=1;
dx=2;
for t=3:199
    tmp=2*Waves(2:end-1,t-1)/dt^2-Waves(2:end-1,t-2)/dt^2+alpha(2:end-1).*Waves(2:end-1,t-2)/(2*dt)+c(2:end-1).^2.*(Waves(3:end,t-1)-2*Waves(2:end-1,t-1)+Waves(1:end-2,t-1))/(2*dx);
    Waves(2:end-1,t) = tmp./(ones(23,1)./dt^2+alpha(2:end-1)./(2*dt));
end

rng(0);
noise = sqrt(var(Waves(:))*0.001)*randn(size(Waves)); %sqrt(0.38*0.01)
Xn = Waves+noise;

figure
imagesc(Xn)
colorbar
ylabel('x (\Delta x)')
xlabel('t (\Delta t)')
ax=gca
ax.FontSize=14;
%%exportgraphics(gcf,'U_att.png')
%% Heat Eq
H=zeros(101,2001);
x = -2.5:0.05:2.5;
alpha = 0.5*normpdf(x);
alpha = alpha';
figure
plot(alpha)

H(1,:) = 1000;
H(101,:) = 1000;
H(26,1:2) = 1000;
H(51,1:2) = 1000;
H(76,1:2) = 1000;

dx = 0.1;
dt = 5e-5;

for t=2:1001
    H(2:end-1,t+1) = 2*dt/(dx^2).*alpha(2:end-1).*(H(3:end,t)-2*H(2:end-1,t)+H(1:end-2,t))+H(2:end-1,t-1);
end

figure
imagesc(log10(H(:,1:1001)))
%colorbar
ax=gca
cbh = colorbar('XTickLabel',{'$\leq -6$','-3','0','3'}, ...
               'XTick', -6:3:3,...
              'TickLabelInterpreter','latex')
caxis([-6,3])
set(gca,'YDir','normal') 
xticks([1,200,400,600,800,1000])
xticklabels({'0','10','20','30','40','50'})
xlabel('$t~(\rm ms)$','interpreter','latex')
ylabel('$x~(\rm m)$','interpreter','latex')
yticks([1,26,51,76,101])
yticklabels({'0','2.5','5','7.5','10'})
title('$log_{10}(\mathbf U)$','interpreter','latex')
ax.FontSize = 22
ax.TickLabelInterpreter = 'latex';
exportgraphics(ax,'heat_data_gauss.png')
%% Burgers Eq
Waves = zeros(101,201);
u0 = 10*normpdf(x);
Waves(:,1) = u0';
Waves(:,2) = u0';
dx = 1;
dt = 0.05;
nu = 0.1;
for t=2:200
    Waves(2:end-1,t+1) = nu/(dx^2)*2*dt*(Waves(3:end,t)-2*Waves(2:end-1,t)+Waves(1:end-2,t))-Waves(2:end-1,t).*(Waves(3:end,t)-Waves(1:end-2,t))./dx.*dt+Waves(2:end-1,t-1);
end
figure
imagesc(Waves(:,1:151))
caxis([0 4.5])
hold on
plot(1:152,21*ones(152,1),'r','LineWidth',2)
hold on
plot(1:152,91*ones(152,1),'r','LineWidth',2)
cbh = colorbar('XTick', 0:1:4,...
              'TickLabelInterpreter','latex')
set(gca,'YDir','normal') 
xticks([1,51,101,151])
xticklabels({'0','2.5','5','7.5'})
yticks([1,21,41,61,81,101])
yticklabels({'0','20','40','60','80','100'})
xlabel('$t~(\rm s)$','interpreter','latex')
ylabel('$x~(\rm m)$','interpreter','latex')
title('$\mathbf U$','interpreter','latex')
ax = gca
ax.FontSize = 22
ax.TickLabelInterpreter = 'latex';

exportgraphics(ax,'Burgers_data_new.png')
% figure
% plot(Waves(:,end)-Waves(:,1))
%%
Ix = 3:23;
alphas = [1e-7,1e-6,1e-5,1e-4,1e-3,1e-2];
alphas = [alphas,0.1:0.1:1];
diff = 1e10;
for a=10:10
    alpha = alphas(a);
u_slicex = Xn(3:23,5);
u_true = Waves(3:23,5);
ux = fdtvr(u_slicex,dx,alpha,1e-6,1000);
uxx = fdtvr(ux,dx,alpha,1e-6,1000);
ux_true = numder(Waves(1:25,5),dx,1,'FD');
uxx_true = numder(Waves(1:25,5),dx,2,'FD');
if(norm(ux_true(3:23)-ux)<diff)
    diff = norm(ux_true(3:23)-ux);
    sel_x = a;
end
end
ax = figure
subplot(3,1,1)
plot(Ix,u_slicex)
hold on
plot(Ix,u_true)
legend({'noisy','true'},'FontSize',12)
xlim([3,23])
ylabel('u')
subplot(3,1,2)
plot(Ix,ux)
hold on
plot(Ix,ux_true(3:23))
legend({'FD-TVR','true'},'FontSize',12)
xlim([3,23])
ylabel('u_x')
subplot(3,1,3)
plot(Ix,uxx)
hold on
plot(Ix,uxx_true(3:23))
legend({'FD-TVR','true'},'FontSize',12)
xlim([3,23])
ylabel('u_{xx}')
xlabel('$x ({\rm d}x)$','Interpreter','latex')
%ax = gca
exportgraphics(ax,'TVR_ux_alpha04.png')
%%
alphas = [1e-7,1e-6,1e-5,1e-4,1e-3,1e-2];
alphas = [alphas,0.1:0.1:1];
diff = 1e10;
for b=1:16
    alpha = alphas(b);
u_slicet = Xn(3,5:195);
u_true = Waves(3:23,5);
ut = fdtvr(u_slicet',dt,alpha,1e-6,1000);
utt = fdtvr(ut,dt,alpha,1e-6,1000);
ut_true = numder(Waves(3,:),dt,1,'FD');
utt_true = numder(Waves(3,:),dt,2,'FD');
if(norm(ut_true(5:195)-ut)<diff)
    diff = norm(ut_true(5:195)-ut);
    sel_t = b;
end
end
figure
plot(u_slicet)
hold on
plot(Waves(3,5:195))
%%
%%
Waves=zeros(100,150);
x_tmp=0:57;
k=2*pi/19;
amp = 5;
Waves(1:58,1)=amp*sin(k*x_tmp+pi/2)-amp;
figure
plot(Waves(:,1))
Waves(:,2)=Waves(:,1);
alpha = zeros(100,1);
c = 1*ones(100,1);
c(1:30)=1.5;
c(30:35) = 1.5:-0.1:1;
c(55:60) = 1:0.2:2;
c(60:end)=2;
alpha(45:50)=0:0.01:0.05;
alpha(51:end) = 0.05;

figure
subplot(1,2,1)
plot(c)
subplot(1,2,2)
plot(alpha)

dt=1;
dx=2;
for t=3:149
    tmp=2*Waves(2:99,t-1)/dt^2-Waves(2:99,t-2)/dt^2+alpha(2:99).*Waves(2:99,t-2)/(2*dt)+c(2:99).^2.*(Waves(3:100,t-1)-2*Waves(2:99,t-1)+Waves(1:98,t-1))/(2*dx);
    Waves(2:99,t) = tmp./(ones(98,1)./dt^2+alpha(2:99)./(2*dt));
end

figure
imagesc(Waves)
colorbar
ylabel('x (\Delta x)')
xlabel('t (\Delta t)')
ax=gca
ax.FontSize=14;
%%
loc=0;
mindiff=1e10;
for i=1:epoch
    if(mse(LAMBDA1(i,1:2:end),alpha(3:23))<mindiff)
        loc = i;
        mindiff = mse(LAMBDA1(i,1:2:end),alpha(3:23));
    end
end
%%
%% Spline inpainting
load('wave1D_calpha_NN.mat', 'x');
xq = 0:1:198;
WavesInt = zeros(199,150);
for i=1:150
    v = Xn(:,i);
    WavesInt(:,i) = interp1(x,v,xq,'spline');
end
figure
imagesc(WavesInt)
title('Spline Interpolation');
%% Coef domain
Coef = zeros(100,7);
Coef(1:30,5) = -1.5^2;
Coef(60:end,5) = -2^2;
Coef(31:59,5) = -1^2;
Coef(51:end,2) = 1/100;
Coef(:,3) = 1;
figure
imagesc(Coef)
colorbar
fc  = fftshift(fft2(Coef));
figure
imagesc(abs(fc))
colorbar
%% PINN result
PINN_RES=[[ 1.2966331e-05]
 [ 9.8960977e-03]
 [ 1.9425319e-01]
 [ 2.8519400e-03]
 [-7.8916609e-01]
 [-6.7222449e-03]
 [ 1.5093613e-04]];
%%
PINN_RES0 = [[-7.2086274e-05]
 [-4.0381824e-04]
 [ 2.9703006e-01]
 [-1.2721907e-03]
 [-6.8207133e-01]
 [-1.1969993e-02]
 [ 4.5715256e-05]];
PINN_RES1 = [[-6.7528621e-05]
 [-1.7391150e-04]
 [ 2.7364215e-01]
 [-1.8668197e-03]
 [-5.4130405e-01]
 [-9.8567478e-02]
 [ 4.6448589e-05]];
PINN_RES2 = [[-1.9707459e-05]
 [-9.2646998e-04]
 [ 3.7714267e-01]
 [-8.7637309e-04]
 [-4.8886976e-01]
 [-8.0407284e-02]
 [ 1.2458164e-04]];
PINN_RES3 = [[-2.9422250e-04]
 [ 7.5461231e-03]
 [-6.5358495e-03]
 [ 1.5838008e-02]
 [-4.1236731e-01]
 [-3.9530143e-01]
 [-1.5973023e-03]];
PINN_RES4 = [[-9.2512555e-06]
 [ 9.0210605e-03]
 [-2.5284439e-02]
 [ 2.1516183e-02]
 [-4.4008544e-01]
 [-4.8965529e-01]
 [-5.9831371e-05]];
PINN_RES5 = [[-3.5587222e-05]
 [-7.5885968e-04]
 [-3.0162880e-01]
 [ 7.2904542e-02]
 [-3.9240018e-01]
 [-8.2147968e-01]
 [-3.7805026e-04]];
PINN_RES6 = [[ 6.7373847e-05]
 [ 1.5451058e-02]
 [-1.8976757e-02]
 [-7.9309992e-02]
 [-6.7462373e-01]
 [-4.0931901e-01]
 [ 8.6937926e-04]];
PINN_RES7 =  [[ 5.1366878e-06]
 [-7.0691537e-03]
 [ 6.1091565e-02]
 [ 1.4289336e-01]
 [-5.9147328e-01]
 [-2.1188062e-01]
 [-9.1862091e-04]];

PINN_RES = [PINN_RES0,PINN_RES1,PINN_RES2,PINN_RES3,PINN_RES4,PINN_RES5,PINN_RES6,PINN_RES7];
figure
imagesc(abs(PINN_RES))
colorbar
ylabel('dict index')
xlabel('location interval index')
ax = gca;
ax.FontSize = 15;
norms = [];
for i=1:7
    norms = [norms,norm(PINN_RES(i,:))];
end
norms
%% Build dict- FD-TVR
alpha=0.3;
clear Waves
Waves = Xn(:,:);
%method = 'FD';
N_x = size(Waves,1);
Mused = size(Waves,2);
Ux=zeros(N_x,Mused);
Uxx=zeros(N_x,Mused);
for k = 1:size(Waves,2)
    Utmp=Waves(:,k);
    ux = fdtvr(Utmp,dx,alpha,1e-6,1000);%alphas(sel_x)
    Ux(:,k)=ux;
    Uxx(:,k)=fdtvr(ux,dx,alpha,1e-6,1000);
end

% time derivatives
Ut=zeros(N_x,Mused);
Utt=zeros(N_x,Mused);
for i = 1:N_x
    ut = fdtvr(Waves(i,:)',dt,alpha,1e-6,1000);
    sum(ut)
    Ut(i,:) = ut;
    Utt(i,:) = fdtvr(ut,dt,alpha,1e-6,1000);
end
% 
% Utx=zeros(N_x,Mused);
% %Usin=zeros(N_x,Mused);
% for k = 1:size(Waves,2)
%     Utmp=Ut(:,k);
%     Utx(:,k)=numder(Utmp, dx, 1,method);    
% end
% Usin = sin(Waves);
%% Build dict 
clear Waves
Waves = Xn(:,:);
method = 'FD';
N_x = size(Waves,1);
Mused = size(Waves,2);
Ux=zeros(N_x,Mused);
Uxx=zeros(N_x,Mused);
for k = 1:size(Waves,2)
    Utmp=Waves(:,k);
    Ux(:,k)=numder(Utmp, dx, 1,method);
    Uxx(:,k)=numder(Utmp, dx, 2,method);
end

% time derivatives
Ut=zeros(N_x,Mused);
Utt=zeros(N_x,Mused);
for i = 1:N_x
    Ut(i,:) = numder(Waves(i,:),dt,1,method);
    Utt(i,:) = numder(Waves(i,:),dt,2,method);
end

Utx=zeros(N_x,Mused);
%Usin=zeros(N_x,Mused);
for k = 1:size(Waves,2)
    Utmp=Ut(:,k);
    Utx(:,k)=numder(Utmp, dx, 1,method);    
end
Usin = sin(Waves);
%% GD_waves
tsample = 50:1:140;
Ds = cell(length(tsample),96);
for i=1:96
    for t=1:length(tsample)
        dict = [Ut(2+i,tsample(t)), Uxx(2+i,tsample(t))];
        Ds{t,i} = dict;
    end
end
Ts = cell(length(tsample),96);
for i=1:96
    for t=1:length(tsample)
        Ts{t,i} = [Utt(2+i,tsample(t))];
    end
end
C = zeros(96,2);
delta = 0.5;
Loss = zeros(500000,1);
%Loss(1:200000) = Loss200000;
tic
for iter =1:500000
    dC = zeros(96,2);
    if(mod(iter,100)==0)
        iter
        loss
        toc
    end
     if(mod(iter,2500)==0)
        delta = delta/2;
        %toc
    end
    loss = 0;
    for p=1:length(tsample)
        for r = 1:96
            for c=1:2            
                dC(r,c) = dC(r,c) + (C(r,1)*Ds{p,r}(1)+C(r,2)*Ds{p,r}(2)-Ts{p,r})*Ds{p,r}(c);
            end
        end
        loss = loss + norm(C(r,1)*Ds{p,r}(1)+C(r,2)*Ds{p,r}(2)-Ts{p,r})^2;
    end
    Loss(iter,1) = loss;
    C = C - delta*dC;
end
figure
plot(Loss)
%title('Loss')
ylabel('Loss')
xlabel('Iterations')
ax=gca;
ax.FontSize = 18;

figure
imagesc(C)
colorbar

Ctrue = zeros(96,2);
Ctrue(:,1) = -alpha(3:98);
Ctrue(:,2) = c(3:98).^2;
figure
imagesc(C)
colorbar
caxis([-0.05,4])
%title('Coefficients at iter=3*10^5')%Coefficients at iter=10^5
xticks([1,2])
yticks([8:10:88]);
yticklabels({'10','20','30','40','50','60','70','80','90'})
xlabel('Atom indices')
ylabel('Location indices')
ax=gca;
ax.FontSize=18;

%norm(Ctrue-C100000,'F')
%%
C=[1,2;3,4];
D1=[5,6;7,8]+nd(:,:,1);
D2=[4,3;2,1]+nd(:,:,2);
D3=[1,3;5,7]+nd(:,:,3);
T1 = sum(C.*D1,2)+nt(:,:,1);
T2 = sum(C.*D2,2)+nt(:,:,2);
T3 = sum(C.*D3,2)+nt(:,:,3);
Cini = zeros(2,2);
delta = 0.0002;
Cpath = zeros(8001,4);
for i=1:16000
dC11 = 2*(sum(Cini(1,:).*D1(1,:))-T1(1))*D1(1,1)+2*(sum(Cini(1,:).*D2(1,:))-T2(1))*D2(1,1)+2*(sum(Cini(1,:).*D3(1,:))-T3(1))*D3(1,1);
dC12 = 2*(sum(Cini(1,:).*D1(1,:))-T1(1))*D1(1,2)+2*(sum(Cini(1,:).*D2(1,:))-T2(1))*D2(1,2)+2*(sum(Cini(1,:).*D3(1,:))-T3(1))*D3(1,2);
dC21 = 2*(sum(Cini(2,:).*D1(2,:))-T1(2))*D1(2,1)+2*(sum(Cini(2,:).*D2(2,:))-T2(2))*D2(2,1)+2*(sum(Cini(2,:).*D3(2,:))-T3(2))*D3(2,1);
dC22 = 2*(sum(Cini(2,:).*D1(2,:))-T1(2))*D1(2,2)+2*(sum(Cini(2,:).*D2(2,:))-T2(2))*D2(2,2)+2*(sum(Cini(2,:).*D3(2,:))-T3(2))*D3(2,2);
DER = [dC11,dC12;dC21,dC22];
Cini = Cini - delta*DER;
Cpath(i+1,1) = Cini(1,1);
Cpath(i+1,2) = Cini(1,2);
Cpath(i+1,3) = Cini(2,1);
Cpath(i+1,4) = Cini(2,2);
end
%%
% figure
% plot(1:100,alpha)
% A=[1,2,3;4,5,6]
% B=reshape(A,[],1)%1 4 2 5 3 6
% Utt = Utt';
% Ut = Ut';
% Uxx = Uxx';
% largeInd = length(Utt);
% rang = 3:largeInd-2;
uttv = reshape(Utt(3:end-2,3:end-2),[],1);
utv = reshape(Ut(3:end-2,3:end-2),[],1);
uxxv = reshape(Uxx(3:end-2,3:end-2),[],1);
alpha1 = repmat(alpha',[150,1]);
alpha2 = 0.05-alpha1;
c1 = repmat((c'.^2),[150,1]);
c2 = 4-c1;
nums = length(utv);
% alpha1 = zeros(nums,1);
% alpha1((50-2)*(150-4)+1:end)=0.05;
% alpha2 = 0.05*ones(nums,1)-alpha1;
% %debug
% % res = uttv-c1.*uxxv;
% % q=res./utv;
% %
% c1 = ones(nums,1);
% c1(1:(30-2)*(150-4)) = 1.5^2;
% c1((60-2)*(150-4)+1:end) = 2^2;
% c2 = ones(nums,1);
% c2(25*(150-4)+1:75*(150-4)) = 0.5;
alpha1 = reshape(alpha1(3:end-2,3:end-2),[],1);
alpha2 = reshape(alpha2(3:end-2,3:end-2),[],1);
c1 = reshape(c1(3:end-2,3:end-2),[],1);
c2 = reshape(c2(3:end-2,3:end-2),[],1);

dict_r = [-alpha1.*utv, -alpha2.*utv, c1.*uxxv, c2.*uxxv];
norms = [];
for i=1:4
    norms = [norms;norm(dict_r(:,i))];
end
dict_r = normc(dict_r);

lam=0.01*max(abs(dict_r'*uttv))/length(dict_r);
% lasso
[a_ext_tmp, stru] = lasso(dict_r,uttv,'Lambda',lam,'Intercept',false,'RelTol',1e-8,'MaxIter',10^6);
a_ext_tmp
a_ext_tmp./norms
% recalc
[-alpha1.*utv, c1.*uxxv]\uttv
%% de-block
s_utt = sum(Utt(3:end-2,3:end-2));
s_uxx = sum(Uxx(3:end-2,3:end-2));
s_ut = sum(Ut(3:end-2,3:end-2));
dict_sr0 = [(c(3:end-2)).^2.*s_uxx',(2-c(3:end-2)).^2.*s_uxx',(alpha(3:end-2)).*s_ut',(0.025-alpha(3:end-2)).*s_ut'];
norms = [];
for i=1:4
    norms = [norms;norm(dict_sr0(:,i))];
end
dict_sr = normc(dict_sr0);
lam=0.1*max(abs(dict_sr'*s_utt'))/length(dict_sr);
% lasso
[a_ext_tmp, stru] = lasso(dict_sr,s_utt','Lambda',lam,'Intercept',false,'RelTol',1e-8,'MaxIter',10^6);
a_ext_tmp
a_ext_tmp./norms
% recalc
[dict_sr0(:,1), dict_sr0(:,3)]\s_utt'
%% recover coefs
step=5;
dicts = cell(21,2);
cnt = 1;
for i=3:23
    Uxx_slice = Uxx(i, 5:195);
    Ut_slice = Ut(i, 5:195);
    dict = [Ut_slice', Uxx_slice'];
    dicts{cnt,1} = dict;
    Utt_slice = Utt(i, 5:195);
    dicts{cnt,2} = Utt_slice';
%     Ut_slice = Ut(i,5:195);
%     dicts{cnt,2} = Ut_slice';
    cnt = cnt+1;
end
coef_fd = zeros(21,2);
for i=1:21
    coef_fd(i,:) = dicts{i,1}\dicts{i,2};
end
figure
subplot(2,1,1)
plot(3:23,coef_fd(:,1))
subplot(2,1,2)
plot(3:23,coef_fd(:,2))
%% Plotting
ax = figure
subplot(2,1,1)
plot(3:23,-alpha(3:23))
hold on
plot(3:23,coef_fd(:,1))
hold on
plot(3:23, coef_tvr(:,1))
hold on
plot(3:23,-coef_sdpinn1)
xlim([3,23])
%ylim([-4.1,-0.9])
legend({'True','Recovered by LSQ','Recovered by LSQ with TVR','Recovered by SD-PINN'},'Interpreter','latex','NumColumns',2,'FontSize',16)
%xlabel('x (\rm dx)','Interpreter','latex')
ylabel('$-\alpha$','Interpreter','latex')
xticks([])
ax2=gca
ax2.FontSize=14
ax2.TickLabelInterpreter='latex'
subplot(2,1,2)
plot(3:23,c(3:23).^2)
hold on
plot(3:23,coef_fd(:,2))
hold on
plot(3:23, coef_tvr(:,2))
hold on
plot(3:23,-coef_sdpinn2)
xlim([3,23])
%ylim([-4.1,-0.9])
%legend({'True','Recovered by LSQ','Recovered by LSQ with TVR','Recovered by SD-PINN'},'Interpreter','latex')
xlabel('$x (\Delta  x)$','Interpreter','latex')
ylabel('$c^2$','Interpreter','latex')
xticks([3,5,10,15,20,23])
ax2=gca
ax2.FontSize=14
ax2.TickLabelInterpreter='latex'


% 
%ax.TickLabelInterpreter='latex'
exportgraphics(ax,'recovered-paras-noisymeas_2atten_c.png')
%%
figure
plot(Xn(:,100))
figure
plot(Uxx(:,100))

CO = zeros(96,2);
CO(:,1) = -alpha(3:end-2);
CO(:,2) = c(3:end-2).^2;
COR = randn(96,2)
% error
Error = 0;
error = zeros(10,1);
LIS = [];
for i=1:10
    lis = (dicts{i,2} - CO(:,2).*dicts{i,1}(:,2))./dicts{i,1}(:,1);
    LIS = [LIS,lis];
    error(i) = norm(CO(:,1).*dicts{i,1}(:,1)+CO(:,2).*dicts{i,1}(:,2)-dicts{i,2});
    Error = sum(error);
end
Error
%% GD test
rng(0)
nd = .5*randn(2,2,3);
nt = .5*randn(2,1,3);
C=[1,2;3,4];
D1=[5,6;7,8]+nd(:,:,1);
D2=[4,3;2,1]+nd(:,:,2);
D3=[1,3;5,7]+nd(:,:,3);
T1 = sum(C.*D1,2)+nt(:,:,1);
T2 = sum(C.*D2,2)+nt(:,:,2);
T3 = sum(C.*D3,2)+nt(:,:,3);
Cini = zeros(2,2);
delta = 0.0002;
Cpath = zeros(8001,4);
for i=1:16000
dC11 = 2*(sum(Cini(1,:).*D1(1,:))-T1(1))*D1(1,1)+2*(sum(Cini(1,:).*D2(1,:))-T2(1))*D2(1,1)+2*(sum(Cini(1,:).*D3(1,:))-T3(1))*D3(1,1);
dC12 = 2*(sum(Cini(1,:).*D1(1,:))-T1(1))*D1(1,2)+2*(sum(Cini(1,:).*D2(1,:))-T2(1))*D2(1,2)+2*(sum(Cini(1,:).*D3(1,:))-T3(1))*D3(1,2);
dC21 = 2*(sum(Cini(2,:).*D1(2,:))-T1(2))*D1(2,1)+2*(sum(Cini(2,:).*D2(2,:))-T2(2))*D2(2,1)+2*(sum(Cini(2,:).*D3(2,:))-T3(2))*D3(2,1);
dC22 = 2*(sum(Cini(2,:).*D1(2,:))-T1(2))*D1(2,2)+2*(sum(Cini(2,:).*D2(2,:))-T2(2))*D2(2,2)+2*(sum(Cini(2,:).*D3(2,:))-T3(2))*D3(2,2);
DER = [dC11,dC12;dC21,dC22];
Cini = Cini - delta*DER;
Cpath(i+1,1) = Cini(1,1);
Cpath(i+1,2) = Cini(1,2);
Cpath(i+1,3) = Cini(2,1);
Cpath(i+1,4) = Cini(2,2);
end
figure
for i=1:4
    plot(Cpath(:,i));
    hold on
end
ylim([0,5])
legend({'1','2','3','4'})
% matrix inverse
T = [T1(1), T2(1), T1(2), T2(2)];
%D1=[5,6;7,8];D2=[4,3;2,1];
D = [D1(1,1),D2(1,1),0,0;D1(1,2),D2(1,2),0,0;0,0,D1(2,1),D2(2,1);0,0,D1(2,2),D2(2,2)];
recC = T*inv(D);
recC = reshape(recC',2,2)'
%%
Ix = 10:19;
It = 3:148;
lhs = reshape(Utt(Ix,It),[],1);
utv = reshape(Ut(Ix,It),[],1);
uxxv = reshape(Uxx(Ix,It),[],1);
rhs = [utv, uxxv];
rhs\lhs
%% Spatially Dependent Direct LSQ
It = 3:148;
coefs = zeros(51,2);
for i=1:51
    utv = Ut(i+19,It);
    uxxv = Uxx(i+19,It);
    uttv = Utt(i+19,It);
    rhs = [utv', uxxv'];
    coefs(i,:) = rhs\uttv';
end
figure
subplot(2,1,1)
plot(GT1)
hold on
plot(-coefs(:,1))
legend({'True','Recovered'})
xlabel('spatial index')
ylabel('\alpha')
title('Least Squares')
ax=gca
ax.FontSize=14;
subplot(2,1,2)
plot(GT2)
hold on
plot(-coefs(:,2))
ylim([-4.1,-0.9])
legend({'True','Recovered'})
xlabel('spatial index')
ylabel('-c^2')
ax=gca
ax.FontSize=14;
%%
dict = [ones((size(Waves(3:end-2,3:end-2),1))*(size(Waves(3:end-2,3:end-2),2)),1),reshape(Ut(3:end-2,3:end-2),[],1),reshape(Utt(3:end-2,3:end-2),[],1),reshape(Waves(3:end-2,3:end-2),[],1).*reshape(Ux(3:end-2,3:end-2),[],1),reshape(Uxx(3:end-2,3:end-2),[],1),reshape(Utx(3:end-2,3:end-2),[],1), reshape(Usin(3:end-2,3:end-2),[],1)];
s = [1,1,1,1,-1,-1,1];

indvalid = [2,3,5];

dicts = [dict;s];
dicts = dicts(:,indvalid);

dicnorm = [];
for i=1:3
    dicnorm = [dicnorm, norm(dict(:,i))];
end


Phi_sn=zeros(size(dict,1),size(dict,2));
for d=1:size(dict,2)
    Phi_sn(:,d) = dict(:,d)./norm(dict(:,d)); % introduce nan column if the column is all 0.
end
Phi_ext = [Phi_sn;s];
tar = zeros(size(dicts,1),1);
tar(end) = 1;

del_c = [];%1,4,6,7
Phi_ext(:,del_c)=[];
% set lambda
lam=0.2*max(abs(Phi_ext'*tar))/length(tar);
% lasso
[a_ext_tmp, stru] = lasso(Phi_ext,tar,'Lambda',lam,'Intercept',false,'RelTol',1e-8,'MaxIter',10^6);
a_ext_tmp
% figure
% imagesc(Uxx)
% colorbar
dicnorm(del_c) = [];
a_ext_tmp./dicnorm'

%%
It=10:140;
Ix=3:98;
m = length(It);
N = length(Ix);
%veclen = length(Ix)*length(Iy)*length(It);
utm = squeeze(reshape(Ut(Ix,It),N,m));
uttm = squeeze(reshape(Utt(Ix,It),N,m));
uxm = squeeze(reshape(Ux(Ix,It),N,m));
uxxm = squeeze(reshape(Uxx(Ix,It),N,m));

Theta_e_tens = zeros(size(utm,2),3,size(utm,1));
Theta_e_tens(:,1,:)=utm';
Theta_e_tens(:,2,:)=uttm';
Theta_e_tens(:,3,:)=uxxm';
%% Recover PDEs
a_raw=zeros(size(Theta_e_tens,2),length(Ix));
a_res=zeros(size(Theta_e_tens,2),length(Ix));
vec_constraint = [1,1,-1];
for i = 1:length(Ix)
    i
    %out(i) = doCvxStuff(i);
    [a_raw(:,i),a_res(:,i)] = CVX_L1min(squeeze(Theta_e_tens(:,:,i)),vec_constraint);
end
%% LASSO
a_raw=zeros(size(Theta_e_tens,2),length(Ix)*length(Iy));
a_res=zeros(size(Theta_e_tens,2),length(Ix)*length(Iy));
vec_constraint = [1,1,1,1,1,-1,-1,-1,-1,-1,-1];
tic
for i=1:length(Ix)*length(Iy)
    %i
    Theta_s = squeeze(Theta_e_tens(:,:,i));
    n = size(Theta_s,2);
%vec1row = ones(1,n);
    Theta_sn=zeros(size(Theta_s,1),size(Theta_s,2));
    for d=1:size(Theta_s,2)
        Theta_sn(:,d) = Theta_s(:,d)./norm(Theta_s(:,d));
    end
    Theta_ext = [Theta_sn;vec_constraint];
    tar = zeros(size(Theta_s,1)+1,1);
    tar(end) = 1;
    epsilon=100;
    lam=1e-3;
    a_ext_tmp = lasso(Theta_ext,tar,'Lambda',lam,'Intercept',false);
    rescale_factor = vecnorm(Theta_s)./vecnorm(Theta_sn);
    a_eres_tmp=a_ext_tmp./(rescale_factor');
    a_raw(:,i) = a_ext_tmp;
    a_res(:,i) = a_eres_tmp;
end
toc

indicator = zeros(11,length(Ix)*length(Iy));
for i=1:length(Ix)*length(Iy)
    [maxv,maxi] = max(abs(a_raw(:,i)));
    ind=find(abs(a_raw(:,i))>=0.01*maxv);
    indicator(ind,i) = 1;
end
figure
imagesc(indicator)

ind_perc = zeros(11,1);
for i=1:11
    ind_perc(i)=sum(indicator(i,:))/(length(Ix)*length(Iy));
end
figure
stem(ind_perc,'filled')

ind_waveeq=zeros(length(Ix)*length(Iy),1);
label = [0,1,1,0,0,1,1,0,0,0,0]';
for i=1:length(Ix)*length(Iy)
    if(sum(indicator(:,i)&label)==4)
        ind_waveeq(i,1) = 1;
    end
end
    

figure
for i=1:11
    plot(a_res(i,:));
    hold on
end




speeds_rec = zeros(size(a_res,2),1);
for i=1:size(a_res,2)
    speeds_rec(i) = sqrt(-a_res(3,i)/a_res(2,i));
end

figure
plot(c,'LineWidth',1)
hold on
plot(Ix,speeds_rec,'rx','LineWidth',1)
hold on
plot(alpha,'r-','LineWidth',1)
hold on
plot(Ix,a_res(1,:)./a_res(2,:),'bx','LineWidth',1)
ylim([0 2.5])
legend({'true c','estimated c','true \alpha','estimated \alpha'})
xlabel('Spatial locations (\Delta x)')
ylabel('Phase speeds (m/s) and attenuating factor')
ax=gca
ax.FontSize=15;
exportgraphics(gcf,'Speeds_att.png')
