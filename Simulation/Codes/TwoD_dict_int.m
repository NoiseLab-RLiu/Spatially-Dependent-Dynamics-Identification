U = Waves;
%%
interval = 11;
half_interval = floor(interval/2);
W_pow = 3;
x=1:interval;
y=1:interval;
t=1:interval;
numt = 49;
numx = 90;
numy = 90;
bias=5;
ker_rel2 = zeros(numx,numy,numt);
int_kernalt2 = zeros(numx,numy,numt);
int_kernalt0 = zeros(numx,numy,numt);
int_kernalx2 = zeros(numx,numy,numt);
int_kernaly2 = zeros(numx,numy,numt);
int_kernalt = zeros(numx,numy,numt);
int_kernaltt = zeros(numx,numy,numt);
int_kernalxx2 = zeros(numx,numy,numt);
int_kernalyy2 = zeros(numx,numy,numt);

int_kernaltx = zeros(numx,numy,numt);
int_kernalty = zeros(numx,numy,numt);
int_kernalsin = zeros(numx,numy,numt);
int_kernalxu2 = zeros(numx,numy,numt);
int_kernalyu2 = zeros(numx,numy,numt);
indt = 0;
for i_t =20:20:980
    indt = indt+1;
    i_t
    tic
    for i_x=6:95
        for i_y=6:95
            kernalt =  zeros(interval,interval);
            kernalt2 =  zeros(interval,interval);
            kernalx =  zeros(interval,interval);
            kernaly =  zeros(interval,interval);
            approx_kernalt2 =  zeros(interval,interval);
    %         kernalt20 =  zeros(interval,interval);
            approx_kernalx2 =  zeros(interval,interval);
            approx_kernaly2 =  zeros(interval,interval);
            approx_kernalt =  zeros(interval,interval);
            approx_kernaltx =  zeros(interval,interval);
            approx_kernalty =  zeros(interval,interval);
            approx_kernalsin =  zeros(interval,interval);
            approx_kernalxu2 =  zeros(interval,interval);
            approx_kernalyu2 =  zeros(interval,interval);

            intt = zeros(interval,interval);
            intx = zeros(interval,interval);
            inty = zeros(interval,interval);
            for i=1:interval
                for j=1:interval

            
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
                    te = i_t-half_interval:.1:i_t+half_interval;
                    expr_wtu1 = -(((i-1-half_interval)./(half_interval))^2-1)^W_pow.*(((j-1-half_interval)./(half_interval))^2-1)^W_pow.*2*W_pow.*((te-i_t)./(half_interval)).*(((te-i_t)./(half_interval)).^2-1).^(W_pow-1).*(tp(1)*(te-(i_t-half_interval)).^5+tp(2)*(te-(i_t-half_interval)).^4+tp(3)*(te-(i_t-half_interval)).^3+tp(4)*(te-(i_t-half_interval)).^2+tp(5)*(te-(i_t-half_interval))+tp(6));
                    
                    expr_wttu1 = 2*(((i-1-half_interval)./(half_interval))^2-1)^W_pow.*(((j-1-half_interval)./(half_interval))^2-1)^W_pow.*W_pow*((2*W_pow-1).*((te-i_t)./(half_interval)).^2-1).*(((te-i_t)./(half_interval)).^2-1).^(W_pow-2).*(tp(1)*(te-(i_t-half_interval)).^5+tp(2)*(te-(i_t-half_interval)).^4+tp(3)*(te-(i_t-half_interval)).^3+tp(4)*(te-(i_t-half_interval)).^2+tp(5)*(te-(i_t-half_interval))+tp(6));
                    expr_wtxu1 = 4*W_pow^2.*((i-1-half_interval)./(half_interval)).*(((i-1-half_interval)./(half_interval))^2-1)^(W_pow-1).*(((j-1-half_interval)./(half_interval))^2-1)^W_pow.*((te-i_t)./(half_interval)).*(((te-i_t)./(half_interval)).^2-1).^(W_pow-1).*(tp(1)*(te-(i_t-half_interval)).^5+tp(2)*(te-(i_t-half_interval)).^4+tp(3)*(te-(i_t-half_interval)).^3+tp(4)*(te-(i_t-half_interval)).^2+tp(5)*(te-(i_t-half_interval))+tp(6));
                    expr_wtyu1 = 4*W_pow^2.*((j-1-half_interval)./(half_interval)).*(((j-1-half_interval)./(half_interval))^2-1)^(W_pow-1).*(((i-1-half_interval)./(half_interval))^2-1)^W_pow.*((te-i_t)./(half_interval)).*(((te-i_t)./(half_interval)).^2-1).^(W_pow-1).*(tp(1)*(te-(i_t-half_interval)).^5+tp(2)*(te-(i_t-half_interval)).^4+tp(3)*(te-(i_t-half_interval)).^3+tp(4)*(te-(i_t-half_interval)).^2+tp(5)*(te-(i_t-half_interval))+tp(6));
                    expr_wsinu = (((i-1-half_interval)./(half_interval))^2-1).^W_pow.*(((j-1-half_interval)./(half_interval))^2-1).^W_pow.*(((te-i_t)./(half_interval)).^2-1).^W_pow.*sin(tp(1)*(te-(i_t-half_interval)).^5+tp(2)*(te-(i_t-half_interval)).^4+tp(3)*(te-(i_t-half_interval)).^3+tp(4)*(te-(i_t-half_interval)).^2+tp(5)*(te-(i_t-half_interval))+tp(6));
                    
                    expr_wxu2 = -0.5*2*W_pow*((i-1-half_interval)./(half_interval)).*(((i-1-half_interval)./(half_interval))^2-1)^(W_pow-1).*(((j-1-half_interval)./(half_interval))^2-1)^W_pow.*(((te-i_t)./(half_interval)).^2-1).^W_pow.*((tp(1)*(te-(i_t-half_interval)).^5+tp(2)*(te-(i_t-half_interval)).^4+tp(3)*(te-(i_t-half_interval)).^3+tp(4)*(te-(i_t-half_interval)).^2+tp(5)*(te-(i_t-half_interval))+tp(6)).^2);
                    expr_wyu2 = -0.5*2*W_pow*((j-1-half_interval)./(half_interval)).*(((j-1-half_interval)./(half_interval))^2-1)^(W_pow-1).*(((i-1-half_interval)./(half_interval))^2-1)^W_pow.*(((te-i_t)./(half_interval)).^2-1).^W_pow.*((tp(1)*(te-(i_t-half_interval)).^5+tp(2)*(te-(i_t-half_interval)).^4+tp(3)*(te-(i_t-half_interval)).^3+tp(4)*(te-(i_t-half_interval)).^2+tp(5)*(te-(i_t-half_interval))+tp(6)).^2);
                    te = i_x-half_interval:.1: i_x+half_interval;
                    expr_wxxu1 = 2*(((i-1-half_interval)./(half_interval))^2-1)^W_pow.*(((j-1-half_interval)./(half_interval))^2-1)^W_pow.*W_pow*((2*W_pow-1).*((te-i_x)./(half_interval)).^2-1).*(((te-i_x)./(half_interval)).^2-1).^(W_pow-2).*(xp(1)*(te-(i_x-half_interval)).^5+xp(2)*(te-(i_x-half_interval)).^4+xp(3)*(te-(i_x-half_interval)).^3+xp(4)*(te-(i_x-half_interval)).^2+xp(5)*(te-(i_x-half_interval))+xp(6));
                    te = i_y-half_interval:.1: i_y+half_interval;
                    expr_wyyu1 = 2*(((i-1-half_interval)./(half_interval))^2-1)^W_pow.*(((j-1-half_interval)./(half_interval))^2-1)^W_pow.*W_pow*((2*W_pow-1).*((te-i_y)./(half_interval)).^2-1).*(((te-i_y)./(half_interval)).^2-1).^(W_pow-2).*(yp(1)*(te-(i_y-half_interval)).^5+yp(2)*(te-(i_y-half_interval)).^4+yp(3)*(te-(i_y-half_interval)).^3+yp(4)*(te-(i_y-half_interval)).^2+yp(5)*(te-(i_y-half_interval))+yp(6));
                    approx_kernalt2(i,j) = trapz(i_t-half_interval:.1: i_t+half_interval, expr_wttu1);% NOTE: ...
                    approx_kernalx2(i,j) = trapz(i_x-half_interval:.1: i_x+half_interval, expr_wxxu1);
                    approx_kernaly2(i,j) = trapz(i_y-half_interval:.1: i_y+half_interval, expr_wyyu1);
                    approx_kernalt(i,j) = trapz(i_t-half_interval:.1: i_t+half_interval, expr_wtu1);

                    approx_kernaltx(i,j) = trapz(i_t-half_interval:.1: i_t+half_interval, expr_wtxu1);
                    approx_kernalty(i,j) = trapz(i_t-half_interval:.1: i_t+half_interval, expr_wtyu1);
                    approx_kernalsin(i,j) = trapz(i_t-half_interval:.1: i_t+half_interval, expr_wsinu);
                    approx_kernalxu2(i,j) = trapz(i_t-half_interval:.1: i_t+half_interval, expr_wxu2);
                    approx_kernalyu2(i,j) = trapz(i_t-half_interval:.1: i_t+half_interval, expr_wyu2);        
            
                    
                end
            end
            int_kernalt2(i_x-bias,i_y-bias,indt) = trapz(y,trapz(x,approx_kernalt2,2),1)/(half_interval*dt)^2;
            int_kernalt0(i_x-bias,i_y-bias,indt) = trapz(y,trapz(x,approx_kernalt,2),1)/(half_interval*dt);
            
            int_kernalt(i_x-bias,i_y-bias,indt) = trapz(y,trapz(x,kernalt,2),1);%/dt;
            int_kernaltt(i_x-bias,i_y-bias,indt) = trapz(y,trapz(x,kernalt2,2),1);
            int_kernalxx2(i_x-bias,i_y-bias,indt) = trapz(y,trapz(t,kernalx,2),1);%/dt;
            int_kernalyy2(i_x-bias,i_y-bias,indt) = trapz(x,trapz(t,kernaly,2),1);
            
            int_kernaltx(i_x-bias,i_y-bias,indt) = trapz(y,trapz(x,approx_kernaltx,2),1)/(half_interval*dx);
            int_kernalty(i_x-bias,i_y-bias,indt) = trapz(y,trapz(x,approx_kernalty,2),1)/(half_interval*dy);
            int_kernalsin(i_x-bias,i_y-bias,indt) = trapz(y,trapz(x,approx_kernalsin,2),1);
            int_kernalxu2(i_x-bias,i_y-bias,indt) = trapz(y,trapz(x,approx_kernalxu2,2),1)/(half_interval*dx);
            int_kernalyu2(i_x-bias,i_y-bias,indt) = trapz(y,trapz(x,approx_kernalyu2,2),1)/(half_interval*dy);
            int_kernalx2(i_x-bias,i_y-bias,indt) =  trapz(y,trapz(t,approx_kernalx2,2),1)/(half_interval*dx)^2;
            int_kernaly2(i_x-bias,i_y-bias,indt) =  trapz(x,trapz(t,approx_kernaly2,2),1)/(half_interval*dy)^2;
            chat =  sqrt(int_kernalt2(i_x-bias,i_y-bias,indt)/(int_kernalx2(i_x-bias,i_y-bias,indt)+int_kernaly2(i_x-bias,i_y-bias,indt)));
            ker_rel2(i_x-bias,i_y-bias,indt) = chat;    
        end         
    end
    toc
end


%%
vec_intkert = zeros(numx*numy,numt);
vec_intkertt = zeros(numx*numy,numt);
vec_intkerxx2 = zeros(numx*numy,numt);
vec_intkeryy2 = zeros(numx*numy,numt);


vec_intkert0 = zeros(numx*numy,numt);
vec_intkeruux = zeros(numx*numy,numt);
vec_intkeruuy = zeros(numx*numy,numt);
vec_intkert2 = zeros(numx*numy,numt);
vec_intkerx2 = zeros(numx*numy,numt);
vec_intkery2 = zeros(numx*numy,numt);
vec_intkertx = zeros(numx*numy,numt);
vec_intkerty = zeros(numx*numy,numt);
vec_intkersin = zeros(numx*numy,numt);
for i=1:numx
    for j=1:numy
        vec_intkert(numy*(i-1)+j,:) = int_kernalt(i,j,:);
        vec_intkertt(numy*(i-1)+j,:) = int_kernaltt(i,j,:);
        vec_intkerxx2(numy*(i-1)+j,:) = dt*int_kernalxx2(i,j,:);
        vec_intkeryy2(numy*(i-1)+j,:) = dt*int_kernalyy2(i,j,:);
        
        vec_intkert0(numy*(i-1)+j,:) = int_kernalt0(i,j,:);
        vec_intkeruux(numy*(i-1)+j,:) = int_kernalxu2(i,j,:);
        vec_intkeruuy(numy*(i-1)+j,:) = int_kernalyu2(i,j,:);

        vec_intkert2(numy*(i-1)+j,:) = int_kernalt2(i,j,:);
        vec_intkerx2(numy*(i-1)+j,:) = int_kernalx2(i,j,:);
        vec_intkery2(numy*(i-1)+j,:) = int_kernaly2(i,j,:);

        vec_intkertx(numy*(i-1)+j,:) = int_kernaltx(i,j,:);
        vec_intkerty(numy*(i-1)+j,:) = int_kernalty(i,j,:);
        vec_intkersin(numy*(i-1)+j,:) = int_kernalsin(i,j,:);
    end
end
save_intpt2 = cell(numx*numy,1);
save_intpx2 = cell(numx*numy,1);
save_intpy2 = cell(numx*numy,1);
save_intpt = cell(numx*numy,1);
save_intpuux = cell(numx*numy,1);
save_intpuuy = cell(numx*numy,1);
save_intptx = cell(numx*numy,1);
save_intpty = cell(numx*numy,1);
save_intpsin = cell(numx*numy,1);

Integral_res = zeros(numt,9,numx*numy);

for i=1:numx*numy
    save_intpt{i} = vec_intkert0(i,:)';
    save_intpuux{i} = vec_intkeruux(i,:)';
    save_intpuuy{i} = vec_intkeruuy(i,:)';

    save_intpt2{i} = vec_intkert2(i,:)';
    save_intpx2{i} = vec_intkerx2(i,:)';
    save_intpy2{i} = vec_intkery2(i,:)';

    save_intptx{i} = vec_intkertx(i,:)';
    save_intpty{i} = vec_intkerty(i,:)';
    save_intpsin{i} = vec_intkersin(i,:)';

    Integral_res(:,1,i) = save_intpt{i};
    Integral_res(:,2,i) = save_intpt2{i};
    Integral_res(:,3,i) = save_intpuux{i};
    Integral_res(:,4,i) = save_intpuuy{i};
    Integral_res(:,5,i) = save_intpx2{i};
    Integral_res(:,6,i) = save_intpy2{i};
    Integral_res(:,7,i) = save_intptx{i};
    Integral_res(:,8,i) = save_intpty{i};
    Integral_res(:,9,i) = save_intpsin{i};
end


%% LSQ
% vec_coef_lsq = zeros(numx*numy,2);
% vec_speed  = zeros(numx*numy,1);
% wtsum = zeros(numx*numy,1);
% wttsum = zeros(numx*numy,1);
% for i=1:numx*numy
%     wti = save_intpt{i};
%     wtti  = save_intpt2{i};
%     wtsum(i) = sum(wti);
%     wttsum(i) = sum(wtti);
% end
% 
% for i=1:numx*numy
%     vec_coef_lsq(i,:)  = pinv([save_intpx2{i},save_intpy2{i}])*save_intpt2{i};
%     vec_speed(i) = sqrt((vec_coef_lsq(i,1)+vec_coef_lsq(i,2))/2);
% end
% 
% speed_rec_lsq = zeros(numx,numy);
% for i=1:numx
%     for j=1:numy
%         speed_rec_lsq(i,j) = vec_speed(numy*(i-1)+j);
%     end
% end
% figure
% imagesc(real(speed_rec_lsq))
% colorbar
% caxis([1.5 3.5])
