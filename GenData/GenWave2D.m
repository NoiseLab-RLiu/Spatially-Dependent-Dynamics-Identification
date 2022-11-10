Waves=zeros(120,120,1000);

alpha = zeros(120,120);
c = 2*ones(120,120);
for i=1:5
    c(25+i:end,25+i:end) = 2+0.2*i;
end
c(31:end,31:end)=3;
alpha(60:end,60:end)=0.5;

figure
imagesc(c(23:32,23:32))
axis square
colorbar

figure
imagesc(alpha)
axis square
colorbar
%% Plotting
axx = figure
subplot(1,2,1)
imagesc(c(2:31,2:31))
caxis([2 4])
cbh = colorbar('XTick', 2:1:4,...
              'TickLabelInterpreter','latex')
set(gca,'YDir','normal') 
xticks([1,10,20,30])
xticklabels({'0.1','1','2','3'})
yticks([1,10,20,30])
yticklabels({'0.1','1','2','3'})
xlabel('$x~(\rm m)$','interpreter','latex')
ylabel('$y~(\rm m)$','interpreter','latex')
title('True $c$ (m/s)','interpreter','latex')
axis square
ax = gca
ax.FontSize = 22
ax.TickLabelInterpreter = 'latex';

subplot(1,2,2)
imagesc(alpha(2:31,2:31))
caxis([0 0.2])
cbh = colorbar('XTick', 0:0.1:0.2,...
              'TickLabelInterpreter','latex')
set(gca,'YDir','normal') 
xticks([1,10,20,30])
xticklabels({'0.1','1','2','3'})
yticks([ ])
% yticklabels({'0.1','1','2','3'})
xlabel('$x~(\rm m)$','interpreter','latex')
%ylabel('$x~(\rm m)$','interpreter','latex')
title('True $\alpha$','interpreter','latex')
axis square
ax = gca
ax.FontSize = 22
ax.TickLabelInterpreter = 'latex';

%exportgraphics(axx,'ac_gt_new.png')
%%

% for i=0:99
%     for j=0:99
%         c(i+1,j+1) = 2.5+.5*sin(2/99*pi*i)*sin(2/99*pi*j);
%         alpha(i+1,j+1) = .125*(1+cos(2/99*pi*i)*cos(2/99*pi*j));
%     end
% end
% 
% figure
% imagesc(alpha)
% colorbar



dt=0.2;
dx=1;
for t=3:999
    tmp=2*Waves(2:119,2:119,t-1)/dt^2-Waves(2:119,2:119,t-2)/dt^2+alpha(2:119,2:119).*Waves(2:119,2:119,t-2)/(2*dt)...
        +c(2:119,2:119).^2.*(Waves(3:120,2:119,t-1)+Waves(1:118,2:119,t-1)-4*Waves(2:119,2:119,t-1)...
        +Waves(2:119,3:120,t-1)+Waves(2:119,1:118,t-1))/(dx^2);
    Waves(2:119,2:119,t) = tmp./(ones(118,118)./dt^2+alpha(2:119,2:119)./(2*dt));
    Waves(19,5,t) = 10*sin(.002*pi*t*dt*t*dt);
    Waves(81,115,t) = 20*sin(.001*pi*t*dt*t*dt);
end

figure
for i=400:599
    imagesc(Waves(21:35,21:35,i));
    %axis square
    set(gca, 'YDir','normal')
    colorbar
    caxis([-4, 4]);
    title(num2str(i))
    pause(.1)
end
%%
axx = figure
%subplot(1,3,1)
subplot('position', [0.06 0.55 0.3 0.3])
imagesc(Waves(:,:,1))
caxis([0 3])
cbh = colorbar('XTick', 0:1:3,...
              'TickLabelInterpreter','latex')
set(gca,'YDir','normal') 
xticks([1,11,21,32])
xticklabels({'0','1','2','3.1'})
yticks([1,11,21,32])
yticklabels({'0','1','2','3.1'})
xlabel('$x~(\rm m)$','interpreter','latex')
ylabel('$y~(\rm m)$','interpreter','latex')
title('$t=0,i_t=0$','interpreter','latex')
axis square
ax = gca
ax.FontSize = 22
ax.TickLabelInterpreter = 'latex';

%subplot(1,3,2)
subplot('position', [0.38 0.55 0.3 0.3])
imagesc(Waves(:,:,31))
cbh = colorbar('XTick', -1:1:1,...
              'TickLabelInterpreter','latex')
axis square
caxis([-1 1])
set(gca,'YDir','normal') 
xticks([1,11,21,32])
xticklabels({'0','1','2','3.1'})
yticks([ ])
%yticklabels({'0.1','1','2','3'})
xlabel('$x~(\rm m)$','interpreter','latex')
%ylabel('$y~(\rm m)$','interpreter','latex')
title('$t=0.3~{\rm s},i_t=30$','interpreter','latex')
axis square
ax = gca
ax.FontSize = 22
ax.TickLabelInterpreter = 'latex';

%subplot(1,3,3)
subplot('position', [0.7 0.55 0.3 0.3])
imagesc(Waves(:,:,51))
cbh = colorbar('XTick', -1:1:1,...
              'TickLabelInterpreter','latex')
axis square
caxis([-1 1])
set(gca,'YDir','normal') 
xticks([1,11,21,32])
xticklabels({'0','1','2','3.1'})
yticks([ ])
%yticklabels({'0.1','1','2','3'})
xlabel('$x~(\rm m)$','interpreter','latex')
%ylabel('$y~(\rm m)$','interpreter','latex')
title('$t=0.5~{\rm s},i_t=50$','interpreter','latex')
axis square
ax = gca
ax.FontSize = 22
ax.TickLabelInterpreter = 'latex';

exportgraphics(axx,'Waves_32_all_new.png')
%%
Var = zeros(1000,1);
for i=1:1000
    tmp = Waves(:,:,i);
    Var(i) = var(tmp(:));
end
figure
plot(Var)

rng(0);
noise = sqrt(2)*randn(size(Waves));
%Xn = X+noise;
Varn = zeros(1000,1);
for i=1:1000
    tmp = noise(:,:,i);
    Varn(i) = var(tmp(:));
end
figure
plot(Varn)


figure
imagesc(alpha)
colorbar
%%


dt=0.2;
dx=1;
for t=3:999
    tmp=2*waves(2:99,2:99,t-1)/dt^2-waves(2:99,2:99,t-2)/dt^2+alpha(2:99,2:99).*waves(2:99,2:99,t-2)/(2*dt)...
        +c(2:99,2:99).^2.*(waves(3:100,2:99,t-1)+waves(1:98,2:99,t-1)-4*waves(2:99,2:99,t-1)...
        +waves(2:99,3:100,t-1)+waves(2:99,1:98,t-1))/(dx^2);
    waves(2:99,2:99,t) = tmp./(ones(98,98)./dt^2+alpha(2:99,2:99)./(2*dt));
    waves(63,71,t) = 10*sin(.0025*pi*t*dt*t*dt);
    waves(32,19,t) = 10*sin(.00125*pi*t*dt*t*dt);
end

figure
for i=999:999
imagesc(waves(:,:,i))
colorbar
title(num2str(i))
pause(.1)
end

var = zeros(1000,1);
for i=1:1000
    tmp = Waves(:,:,i);
    VAR(i) = var(tmp(:));
end
figure
plot(VAR)





%% Build dict
Uused  = Waves(:,:,400:599);
[N_x, N_y, Mused] = size(Uused);
method = 'FD';
%N_x = size(Uused,1);
%Mused = size(Uused,2);
Ux=zeros(N_x,N_y,Mused);
Uxx=zeros(N_x,N_y,Mused);
for t=1:size(Uused,3)
    Utmp = Uused(:,:,t);
    for k = 1:size(Utmp,2)
        Utmp2=Utmp(:,k);
        Ux(:,k,t)=numder(Utmp2, dx, 1,method);
        Uxx(:,k,t)=numder(Utmp2, dx, 2,method);
    end
end

Uy=zeros(N_x,N_y,Mused);
Uyy=zeros(N_x,N_y,Mused);
for t=1:size(Uused,3)
    Utmp = Uused(:,:,t);
    for k = 1:size(Utmp,1)
        Utmp2=Utmp(k,:);
        Uy(k,:,t)=numder(Utmp2, dx, 1,method);
        Uyy(k,:,t)=numder(Utmp2, dx, 2,method);
    end
end


% time derivatives
Ut=zeros(N_x,N_y,Mused);
Utt=zeros(N_x,N_y,Mused);
for ix=1:size(Uused,1)
    for iy = 1:size(Uused,2)
        Utmp=Uused(ix,iy,:);
        Ut(ix,iy,:)=numder(Utmp, dt, 1,method);
        Utt(ix,iy,:)=numder(Utmp, dt, 2,method);
    end
end

It = 400:599;
coefhat = zeros(10,10,2);
chat = zeros(10,10);
biasx = 22;
biasy = 22;
for ix = 23:32
    for iy = 23:32
        uxx_vec = reshape(squeeze(Uxx(ix,iy,:)), length(It),1);
        uyy_vec = reshape(squeeze(Uyy(ix,iy,:)), length(It),1);
        dict = [uxx_vec, uyy_vec];
        utt_vec = reshape(squeeze(Utt(ix,iy,:)), length(It),1);
        coefhat(ix-biasx,iy-biasy, :) = dict\utt_vec;
        chat(ix-biasx,iy-biasy) = sqrt((coefhat(ix-biasx,iy-biasy,1)+coefhat(ix-biasx,iy-biasy,2))/2);
    end
end
figure
imagesc(chat)
colorbar
%%
It=6:495;
Ix=6:95;
Iy=Ix;
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

%% test relationship
%W = (((x-x_cen)/(5*dx)).^2-1).^3.*(((y-y_cen)/(5*dy)).^2-1).^3;
x_cen = dx*Ix(1); %x_center
x_inf =  x_cen - half_interval*dx;
x_sup = x_cen+half_interval*dx;
y_cen = dy*Iy(1); %x_center
y_inf =  y_cen - half_interval*dy;
y_sup = y_cen+half_interval*dy;
t_cen = dt*It(1);
t_inf = t_cen - half_interval*dt;
t_sup = t_cen + half_interval*dt;
x = x_inf:dx:x_sup;
x = x';%x must be a column vector to be consistant with U 
y = y_inf:dy:y_sup;
t = t_inf:dt:t_sup;
W_pow=3;
W = zeros(length(x),length(y),length(t));
for i=1:length(x)
    for j=1:length(y)
        for k=1:length(t)
            W(i,j,k) = (((x(i)-x_cen)/(half_interval*dx))^2-1)^W_pow* (((y(j)-y_cen)/(half_interval*dy))^2-1)^W_pow* (((t(k)-t_cen)/(half_interval*dt))^2-1)^W_pow;
        end
    end
end

Wt = zeros(11,11,11);
for i=1:11
    for j=1:11
        Wt(i,j,:) = numder(W(i,j,:), dt, 1,'FD');
    end
end

i_x=50;
i_y=50;
i_t=500;
half_interval = 5;
interval=11;
relation = zeros(interval,interval,interval);
relationW = zeros(interval,interval,interval);
for i=-half_interval:half_interval
    for j=-half_interval:half_interval
        for k=-half_interval:half_interval
            relation(i+half_interval+1,j+half_interval+1,k+half_interval+1) =Utt(i_x+i,i_y+j,i_t+k)/(Uxx(i_x+i,i_y+j,i_t+k)+Uyy(i_x+i,i_y+j,i_t+k));
            relationW(i+half_interval+1,j+half_interval+1,k+half_interval+1)=Utt(i_x+i,i_y+j,i_t+k)*W(i+half_interval+1,j+half_interval+1,k+half_interval+1)/(Uxx(i_x+i,i_y+j,i_t+k)*W(i+half_interval+1,j+half_interval+1,k+half_interval+1)+Uyy(i_x+i,i_y+j,i_t+k)*W(i+half_interval+1,j+half_interval+1,k+half_interval+1));
        end
    end
end

U = zeros(100,100,1000);
Ux = zeros(100,100,1000);
Uxx = zeros(100,100,1000);
Uy = zeros(100,100,1000);
Uyy = zeros(100,100,1000);
Ut = zeros(100,100,1000);
Utt = zeros(100,100,1000);
for i=1:100
    for j=1:100
        for k=1:1000
            U(i,j,k) = sin(2*i)*sin(2*j)*sin(2*k);
            Utt(i,j,k)=-4*sin(2*k)*sin(2*i)*sin(2*j);
            Uxx(i,j,k) = -4*sin(2*i)*sin(2*j)*sin(2*k);
            Uyy(i,j,k)=sin(2*i)*sin(2*k)*(-4*sin(2*j));
            Ut(i,j,k) = 2*cos(2*k)*sin(2*i)*sin(2*j);
            Ux(i,j,k) = 2*cos(2*i)*sin(2*j)*sin(2*k);
            Uy(i,j,k) = 2*cos(2*j)*sin(2*i)*sin(2*k);
        end
    end
end

int_utt = trapz(y,trapz(x,trapz(t,Utt(i_x-half_interval:i_x+half_interval,i_y-half_interval:i_y+half_interval,i_t-half_interval:i_t+half_interval),3),2))
int_uxx = trapz(y,trapz(x,trapz(t,Uxx(i_x-half_interval:i_x+half_interval,i_y-half_interval:i_y+half_interval,i_t-half_interval:i_t+half_interval),3),2))
int_uyy =  trapz(y,trapz(x,trapz(t,Uyy(i_x-half_interval:i_x+half_interval,i_y-half_interval:i_y+half_interval,i_t-half_interval:i_t+half_interval),3),2))
rel_int = int_utt/(int_uxx+int_uyy)

% UW = Uused(i_x-half_interval:i_x+half_interval,i_y-half_interval:i_y+half_interval,i_t-half_interval:i_t+half_interval).*W;
% UWt = zeros(interval,interval,interval);
% for i=2:4
%     for j=2:4
%         for k=2:4
%             UWt(i,j,k) = (UW(i,j,k+1)-UW(i,j,k-1))/(2*dt);
%         end
%     end
% end
% int_uttw = trapz(y,trapz(x,trapz(10*t,Utt(i_x-half_interval:i_x+half_interval,i_y-half_interval:i_y+half_interval,i_t-half_interval:i_t+half_interval).*W,3),2))
% int_uxxw = trapz(y,trapz(x,trapz(t,Uxx(i_x-half_interval:i_x+half_interval,i_y-half_interval:i_y+half_interval,i_t-half_interval:i_t+half_interval).*W,3),2))
% int_uyyw = trapz(y,trapz(x,trapz(t,Uyy(i_x-half_interval:i_x+half_interval,i_y-half_interval:i_y+half_interval,i_t-half_interval:i_t+half_interval).*W,3),2))
ker_rel = zeros(8100,1);
raw_rel = zeros(8100,1);

%%
ker_rel2 = zeros(8100,1);
for i_x=6:95
    i_x
    tic
    for i_y=6:95
        kernalt =  zeros(interval,interval);
        kernalx =  zeros(interval,interval);
        kernaly =  zeros(interval,interval);
        approx_kernalt2 =  zeros(interval,interval);
%         kernalt20 =  zeros(interval,interval);
        approx_kernalx2 =  zeros(interval,interval);
        approx_kernaly2 =  zeros(interval,interval);
        intt = zeros(interval,interval);
        intx = zeros(interval,interval);
        inty = zeros(interval,interval);
        for i=1:interval
            for j=1:interval
        %         kernalt(i,j) = W(i,j,5)*Ut(i_x-3+i,i_y-3+j,i_t+half_interval)-W(i,j,1)*Ut(i_x-3+i,i_y-3+j,i_t-half_interval)-trapz(Wt(i,j,:).*Ut(i_x-3+i,i_y-3+j,i_t-half_interval:i_t+half_interval));
        %         kernalx(i,j) = -trapz(Wx(:,i,j).*Ux(i_x-half_interval:i_x+half_interval,i_y-3+i,i_t-3+j));
        %         kernaly(i,j) = -trapz(Wy(i,:,j).*Uy(i_x-3+i,i_y-half_interval:i_y+half_interval,i_t-3+j));
        %         intt(i,j) = trapz(t,Utt(i_x-3+i,i_y-3+j,i_t-half_interval:i_t+half_interval).*W(i,j,:));
        %         intx(i,j) = trapz(x,Uxx(i_x-half_interval:i_x+half_interval,i_y-3+i,i_t-3+j).*W(:,i,j));
        %         inty(i,j) = trapz(y,Uyy(i_x-3+i,i_y-half_interval:i_y+half_interval,i_t-3+j).*W(i,:,j));
                % rectangle
                kernalt(i,j) = -sum(Wt(i,j,:).*Ut(i_x-(half_interval+1)+i,i_y-(half_interval+1)+j,i_t-half_interval:i_t+half_interval));
                kernalx(i,j) = -sum(Wx(:,i,j).*Ux(i_x-half_interval:i_x+half_interval,i_y-(half_interval+1)+i,i_t-(half_interval+1)+j));
                kernaly(i,j) = -sum(Wy(i,:,j).*Uy(i_x-(half_interval+1)+i,i_y-half_interval:i_y+half_interval,i_t-(half_interval+1)+j));
        
                intt(i,j) = sum(Utt(i_x-(half_interval+1)+i,i_y-(half_interval+1)+j,i_t-half_interval:i_t+half_interval).*W(i,j,:));
                intx(i,j) = sum(Uxx(i_x-half_interval:i_x+half_interval,i_y-(half_interval+1)+i,i_t-(half_interval+1)+j).*W(:,i,j));
                inty(i,j) = sum(Uyy(i_x-(half_interval+1)+i,i_y-half_interval:i_y+half_interval,i_t-(half_interval+1)+j).*W(i,:,j));
        
                uslicet = squeeze(U(i_x-(half_interval+1)+i,i_y-(half_interval+1)+j,i_t-half_interval:i_t+half_interval));
                uslicex = squeeze(U(i_x-half_interval:i_x+half_interval,i_y-(half_interval+1)+i,i_t-(half_interval+1)+j));
                uslicey = squeeze(U(i_x-(half_interval+1)+i,i_y-half_interval:i_y+half_interval,i_t-(half_interval+1)+j));
                tp=polyfit(0:1:10,uslicet,5);
                pt = polyder(tp);
                ptt = polyder(pt);
                xp=polyfit(0:1:10,uslicex,5);
                px = polyder(xp);
                pxx = polyder(px);
                yp=polyfit(0:1:10,uslicey,5);
                py = polyder(yp);
                pyy = polyder(py);
                %syms te
                te = 495:.1:505;
                %expr_wutt = (((te-495)./(half_interval))^2-1).^(W_pow).*(ptt(1)*(te-495).^3+ptt(2)*(te-495).^2+ptt(3)*(te-495)+ptt(4));
                %expr_wtut = 2*W_pow*(te-495)./(half_interval).*(((te-495)./(half_interval)).^2-1).^(W_pow-1).*(pt(1)*(te-495).^4+pt(2)*(te-495).^3+pt(3)*(te-495).^2+pt(4)*(te-495)+pt(5));
                expr_wttu1 = 2*(((i-1-half_interval)./(half_interval))^2-1)^W_pow.*(((j-1-half_interval)./(half_interval))^2-1)^W_pow.*W_pow*((2*W_pow-1).*((te-i_t)./(half_interval)).^2-1).*(((te-i_t)./(half_interval)).^2-1).^(W_pow-2).*(tp(1)*(te-(i_t-half_interval)).^5+tp(2)*(te-(i_t-half_interval)).^4+tp(3)*(te-(i_t-half_interval)).^3+tp(4)*(te-(i_t-half_interval)).^2+tp(5)*(te-(i_t-half_interval))+tp(6));
                %syms te
                %expr_wttu = 2*(((i-1-half_interval)./(half_interval))^2-1)^W_pow.*(((j-1-half_interval)./(half_interval))^2-1)^W_pow.*W_pow*((2*W_pow-1).*((te-i_t)./(half_interval)).^2-1).*(((te-i_t)./(half_interval)).^2-1).^(W_pow-2).*(tp(1)*(te-(i_t-half_interval)).^5+tp(2)*(te-(i_t-half_interval)).^4+tp(3)*(te-(i_t-half_interval)).^3+tp(4)*(te-(i_t-half_interval)).^2+tp(5)*(te-(i_t-half_interval))+tp(6));
                te = i_x-half_interval:.1: i_x+half_interval;
                expr_wxxu1 = 2*(((i-1-half_interval)./(half_interval))^2-1)^W_pow.*(((j-1-half_interval)./(half_interval))^2-1)^W_pow.*W_pow*((2*W_pow-1).*((te-i_x)./(half_interval)).^2-1).*(((te-i_x)./(half_interval)).^2-1).^(W_pow-2).*(xp(1)*(te-(i_x-half_interval)).^5+xp(2)*(te-(i_x-half_interval)).^4+xp(3)*(te-(i_x-half_interval)).^3+xp(4)*(te-(i_x-half_interval)).^2+xp(5)*(te-(i_x-half_interval))+xp(6));
                %syms te
                %expr_wxxu = 2*(((i-1-half_interval)./(half_interval))^2-1)^W_pow.*(((j-1-half_interval)./(half_interval))^2-1)^W_pow.*W_pow*((2*W_pow-1).*((te-i_x)./(half_interval)).^2-1).*(((te-i_x)./(half_interval)).^2-1).^(W_pow-2).*(xp(1)*(te-(i_x-half_interval)).^5+xp(2)*(te-(i_x-half_interval)).^4+xp(3)*(te-(i_x-half_interval)).^3+xp(4)*(te-(i_x-half_interval)).^2+xp(5)*(te-(i_x-half_interval))+xp(6));
                te = i_y-half_interval:.1: i_y+half_interval;
                expr_wyyu1 = 2*(((i-1-half_interval)./(half_interval))^2-1)^W_pow.*(((j-1-half_interval)./(half_interval))^2-1)^W_pow.*W_pow*((2*W_pow-1).*((te-i_y)./(half_interval)).^2-1).*(((te-i_y)./(half_interval)).^2-1).^(W_pow-2).*(yp(1)*(te-(i_y-half_interval)).^5+yp(2)*(te-(i_y-half_interval)).^4+yp(3)*(te-(i_y-half_interval)).^3+yp(4)*(te-(i_y-half_interval)).^2+yp(5)*(te-(i_y-half_interval))+yp(6));
                %syms te
                %expr_wyyu = 2*(((i-1-half_interval)./(half_interval))^2-1)^W_pow.*(((j-1-half_interval)./(half_interval))^2-1)^W_pow.*W_pow*((2*W_pow-1).*((te-i_y)./(half_interval)).^2-1).*(((te-i_y)./(half_interval)).^2-1).^(W_pow-2).*(yp(1)*(te-(i_y-half_interval)).^5+yp(2)*(te-(i_y-half_interval)).^4+yp(3)*(te-(i_y-half_interval)).^3+yp(4)*(te-(i_y-half_interval)).^2+yp(5)*(te-(i_y-half_interval))+yp(6));
                % figure
                % plot(expr_wtut)
                %kernalt2(i,j) = int(expr_wttu,[i_t-half_interval i_t+half_interval]);
                approx_kernalt2(i,j) = trapz(i_t-half_interval:.1: i_t+half_interval, expr_wttu1);
                %kernalx2(i,j)  = int(expr_wxxu,[i_x-half_interval i_x+half_interval]);
                approx_kernalx2(i,j) = trapz(i_x-half_interval:.1: i_x+half_interval, expr_wxxu1);
                %kernaly2(i,j)  = int(expr_wyyu,[i_y-half_interval i_y+half_interval]);
                approx_kernaly2(i,j) = trapz(i_y-half_interval:.1: i_y+half_interval, expr_wyyu1);
        
                %ker_rel(90*(i_y-1-5)+i_x-5) = int_kernalt/(int_kernalx+int_kernaly);
        
                
               %  kernalt20(i,j) = sum(Wtt(i,j,:).*U(i_x-(half_interval+1)+i,i_y-(half_interval+1)+j,i_t-half_interval:i_t+half_interval));
        %         kernalx2(i,j) = sum(Wxx(:,i,j).*U(i_x-half_interval:i_x+half_interval,i_y-(half_interval+1)+i,i_t-(half_interval+1)+j));
        %         kernaly2(i,j) = sum(Wyy(i,:,j).*U(i_x-(half_interval+1)+i,i_y-half_interval:i_y+half_interval,i_t-(half_interval+1)+j));
            end
        end
        int_kernalt2 = trapz(y,trapz(x,approx_kernalt2,2),1);
        int_kernalx2 = trapz(y,trapz(t,approx_kernalx2,2),1);
        int_kernaly2 = trapz(x,trapz(t,approx_kernaly2,2),1);
        chat =  sqrt((1/dt)*int_kernalt2/(int_kernalx2+int_kernaly2))
        ker_rel2(i_x-5,i_y-5) = chat;    
    end
     toc
end
figure
imagesc(reshape(real(ker_rel2),90,90))
colorbar
caxis([1.5 3.5])

%%
i=6;
j=6;
figure
plot(squeeze(U(i_x-(half_interval+1)+i,i_y-(half_interval+1)+j,i_t-half_interval:i_t+half_interval)));%Wtt(i,j,:).*
hold on
plot(fitted)
figure
plot(squeeze(Ut(i_x-(half_interval+1)+i,i_y-(half_interval+1)+j,i_t-half_interval:i_t+half_interval)));%Wt(i,j,:).*
hold on
plot(10*fittedt,'LineWidth',5)
figure
plot(squeeze(Utt(i_x-(half_interval+1)+i,i_y-(half_interval+1)+j,i_t-half_interval:i_t+half_interval)));%.*W(i,j,:)
hold on
plot(fittedtt)

syms te
%te=498:502
expr_wutt = (((te-500)./(half_interval))^2-1).^(W_pow).*(ptt(1)*(te-495).^3+ptt(2)*(te-495).^2+ptt(3)*(te-495)+ptt(4));
expr_wtut = 2*W_pow*(te-500)./(half_interval).*(((te-500)./(half_interval)).^2-1).^(W_pow-1).*(pt(1)*(te-495).^4+pt(2)*(te-495).^3+pt(3)*(te-495).^2+pt(4)*(te-495)+pt(5));
expr_wttu = 2*W_pow*((2*W_pow-1)*((te-500)./(half_interval)).^2-1)*(((te-500)./(half_interval)).^2-1).^(W_pow-2).*(p(1)*(te-495).^5+p(2)*(te-495).^4+p(3)*(te-495).^3+p(4)*(te-495).^2+p(5)*(te-495)+p(6));
% figure
% plot(expr_wtut)
Fint = int(expr_wutt,[i_t-half_interval i_t+half_interval])
Fint = int(expr_wtut,[i_t-half_interval i_t+half_interval])
Fint = int(expr_wttu,[i_t-half_interval i_t+half_interval])

te=495:1:505;
figure
plot(495:.1:505,2*W_pow*((2*W_pow-1)*((te-500)./(half_interval)).^2-1).*(((te-500)./(half_interval)).^2-1).^(W_pow-2).*(p(1)*(te-495).^5+p(2)*(te-495).^4+p(3)*(te-495).^3+p(4)*(te-495).^2+p(5)*(te-495)+p(6)))
sum(2*W_pow*((2*W_pow-1)*((te-500)./(half_interval)).^2-1).*(((te-500)./(half_interval)).^2-1).^(W_pow-2).*(p(1)*(te-495).^5+p(2)*(te-495).^4+p(3)*(te-495).^3+p(4)*(te-495).^2+p(5)*(te-495)+p(6)))

%% polynomial fitting for Uslice
uslice = squeeze(U(i_x-(half_interval+1)+i,i_y-(half_interval+1)+j,i_t-half_interval:i_t+half_interval));
p=polyfit(0:1:10,uslice,5)
ts = 0:1:10;
fitted = zeros(11,1);
for i=1:11
    for pind=1:6
        fitted(i) = fitted(i)+p(pind)*(ts(i)^(6-pind));
    end
end
sum(abs(squeeze(U(i_x-(half_interval+1)+i,i_y-(half_interval+1)+j,i_t-half_interval:i_t+half_interval))-fitted))
pt = polyder(p)
ptt = polyder(pt)
fittedt = zeros(11,1);
fittedtt = zeros(11,1);
for i=1:11
    for pind=1:5
        fittedt(i) = fittedt(i)+pt(pind)*(ts(i)^(5-pind));
    end
end
for i=1:11
    for pind=1:4
        fittedtt(i) = fittedtt(i)+ptt(pind)*(ts(i)^(4-pind));
    end
end
% i=3;
% j=3;
% figure
% plot(-squeeze(Utt(i_x-3+i,i_y-3+j,i_t-half_interval:i_t+half_interval).*W(i,j,:)))
% hold on
% plot(squeeze(Wt(i,j,:).*Ut(i_x-3+i,i_y-3+j,i_t-half_interval:i_t+half_interval)))
% -sum(Utt(i_x-(half_interval+1)+i,i_y-(half_interval+1)+j,i_t-half_interval:i_t+half_interval).*W(i,j,:))
% sum(Wt(i,j,:).*Ut(i_x-3+i,i_y-3+j,i_t-half_interval:i_t+half_interval))
% trapz(498:502,Wt(i,j,:).*Ut(i_x-3+i,i_y-3+j,i_t-half_interval:i_t+half_interval))
% %UWt-Wtt.*Uused(i_x-half_interval:i_x+half_interval,i_y-half_interval:i_y+half_interval,i_t-half_interval:i_t+half_interval);
int_kernalt = trapz(y,trapz(x,kernalt,2),1);
int_kernalx = trapz(y,trapz(t,kernalx,2),1);
int_kernaly = trapz(x,trapz(t,kernaly,2),1);
int_kernalt/(int_kernalx+int_kernaly)
ker_rel(90*(i_y-1-5)+i_x-5) = int_kernalt/(int_kernalx+int_kernaly);

int_kernalt2 = trapz(y,trapz(x,kernalt2,2),1);
int_kernalx2 = trapz(y,trapz(t,kernalx2,2),1);
int_kernaly2 = trapz(x,trapz(t,kernaly2,2),1);
int_kernalt2/(int_kernalx2+int_kernaly2)

int_uttw = trapz(y,trapz(x,intt,2),1);
int_uxxw = trapz(y,trapz(t,intx,2),1);
int_uyyw = trapz(x,trapz(t,inty,2),1);
int_uttw/(int_uxxw+int_uyyw)
raw_rel(90*(i_y-1-5)+i_x-5) = int_uttw/(int_uxxw+int_uyyw);


figure
plot(ker_rel)
ylim([4 10])
speed = zeros(90,90);
for i=1:90
    for j=1:90
        speed(i,j) = real(sqrt(ker_rel(90*(j-1)+i)));
    end
end
figure
imagesc(speed)
colorbar
caxis([1.5 3.5])

syms te
%te=498:502
expr_wtut = 2*W_pow*(te-500)./2.*(((te-500)./2).^2-1).^(W_pow-1).*(sin(2*i_x)*sin(2*i_y)*2*cos(2.*te));
% figure
% plot(expr_wtut)
Fint = int(expr_wtut,[i_t-half_interval i_t+half_interval])
te=495:505;
figure
plot(te,2*W_pow*(0^2-1)^W_pow*(0^2-1)^W_pow*(te-500)./2.*(((te-500)./2).^2-1).^(W_pow-1).*(sin(2*i_x)*sin(2*i_y)*2*cos(2.*te)))
hold on
te=498:502;
plot(te,2*W_pow*(0^2-1)^W_pow*(0^2-1)^W_pow*(te-500)./2.*(((te-500)./2).^2-1).^(W_pow-1).*(sin(2*i_x)*sin(2*i_y)*2*cos(2.*te)))


trapz(te,2*W_pow*(te-500)./2.*(((te-500)./2).^2-1).^(W_pow-1).*(sin(2*i_x)*sin(2*i_y)*2*cos(2.*te)))
trapz(te,(2*W_pow*(te-500)./half_interval.*(((te-500)./half_interval).^2-1).^(W_pow-1))'.*squeeze(Ut(i_x,i_y,i_t-half_interval:i_t+half_interval)))
i=6;
j=6;
-trapz(te,Utt(i_x-(half_interval+1)+i,i_y-(half_interval+1)+j,i_t-half_interval:i_t+half_interval).*W(i,j,:))
figure
plot(squeeze(Ut(i_x-(half_interval+1)+i,i_y-(half_interval+1)+j,i_t-half_interval:i_t+half_interval)))
figure
plot((2*W_pow*(te-500)./half_interval.*(((te-500)./half_interval).^2-1).^(W_pow-1)))
figure
plot(squeeze(W(i,j,:)))
sum(Utt(i_x-(half_interval+1)+i,i_y-(half_interval+1)+j,i_t-half_interval:i_t+half_interval).*W(i,j,:))
sum((2*W_pow*(te-500)./half_interval.*(((te-500)./half_interval).^2-1).^(W_pow-1))'.*squeeze(Ut(i_x,i_y,i_t-half_interval:i_t+half_interval)))
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
    if(sum(indicator(:,i)==label)==11)
        ind_waveeq(i,1) = 1;
    end
end

%     for i=1:length(Ix)
%         for j=1:length(Iy)
%             Theta_e_tens(:,3,length(Iy)*(i-1)+j) = squeeze(Utt(Ix(i),Iy(j),It));

%c_rec=zeros(100,100);
a_rec=zeros(100,100); 
A_lsq=[];
tar = zeros(199,1);
tar(end) = 1;
vec_cons_3 = [1,-1,-1];
vec_cons_4 = [1,1,-1,-1];
for i=1:length(Ix)
    for j=1:length(Iy)
        ind_trans = length(Iy)*(i-1)+j;
        if(ind_waveeq(ind_trans)==1)
            %Theta_s=[Theta_e_tens(:,3,ind_trans),Theta_e_tens(:,6,ind_trans),Theta_e_tens(:,7,ind_trans)];
            Theta_s=[Theta_e_tens(:,2,ind_trans),Theta_e_tens(:,3,ind_trans),Theta_e_tens(:,6,ind_trans),Theta_e_tens(:,7,ind_trans)];
            %Theta_s=[Theta_s;vec_cons_3];
            Theta_s=[Theta_s;vec_cons_4];
            a_lsq = Theta_s\tar;
            A_lsq=[A_lsq,a_lsq];
            %c_rec(Ix(i),Iy(j)) = sqrt((abs(a_res(6,ind_trans))+abs(a_res(7,ind_trans)))/2*(abs(a_res(3,ind_trans))));
            %c_rec(Ix(i),Iy(j)) = sqrt((abs(a_lsq(2))+abs(a_lsq(3)))/(2*(abs(a_lsq(1)))));
            c_rec(Ix(i),Iy(j)) = sqrt((abs(a_lsq(3))+abs(a_lsq(4)))/(2*(abs(a_lsq(2)))));
            a_rec(Ix(i),Iy(j)) = a_lsq(1)/a_lsq(2);
        end
    end
end

figure
imagesc(c_rec)%(Ix,Iy)
axis square
colorbar

c_indicat=ones(100,100);
c_indicat=c_indicat&c_rec;
c_valid = sum(sum(c_indicat))/(98^2);
           

figure
imagesc(a_rec)%(Ix,Iy)
axis square
colorbar


%%
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