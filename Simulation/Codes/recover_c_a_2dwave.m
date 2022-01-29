clear c_rec;
clear a_rec;
c1_rec=nan(2,size(Uused,1),size(Uused,2)); % recovered c
c2_rec=nan(2,size(Uused,1),size(Uused,2)); % recovered c
c3_rec=nan(size(Uused,1),size(Uused,2)); % recovered c

a_rec=nan(size(Uused,1),size(Uused,2)); % recovered c
a1 = nan(size(Uused,1),size(Uused,2));
a2 = nan(size(Uused,1),size(Uused,2));

a3_rec = nan(size(Uused,1),size(Uused,2)); 
a3_lsq = nan(3,size(Uused,1),size(Uused,2));

%a_rec=nan(length(Ix),length(Iy)); % recovered alpha

ind_waves = zeros(N,1);
ind_waveeq3 = zeros(N,1);
ind_waveeq4 = zeros(N,1);
for i=1:N
    if(all(indicator(:,i)==[0,1,0,0,1,1,0,0,0]') )
        ind_waveeq3(i) = 1;
    end
    if(all(indicator(:,i)==[1,1,0,0,1,1,0,0,0]'))
        ind_waveeq4(i) = 1;
    end
    if(all(indicator(:,i)==[0,1,0,0,1,1,0,0,0]') || all(indicator(:,i)==[1,1,0,0,1,1,0,0,0]'))% add a 0
        ind_waves(i) = 1;
    end
end
%% number of locations where the wave eq. is successfully recovered
percent = sum(ind_waves)/N;
% sum(ind_waves)/N
% sum(ind_waveeq4)/45^2
Phi_tensor = dict;
%% recover speed c and attenuating factor alpha
%if(percent(cfreq-2)~=0)
    tar = zeros(M+1,1);
    tar(end) = 1;
    vec_cons_3 = [1,-1,-1];
    vec_cons_4 = [1,1,-1,-1];
    %diff56=0;
    for i=1:length(Ix)
        for j=1:length(Iy)
            ind_trans = length(Iy)*(i-1)+j;
           % if(sum(sum(Integral_res(:,:,ind_trans)))~=0)
            if(ind_waveeq3(ind_trans)==1)
                Phi_s=[Phi_tensor(:,2,ind_trans),Phi_tensor(:,5,ind_trans),Phi_tensor(:,6,ind_trans)];
                Phi_s=[Phi_s;vec_cons_3];
                a_lsq = Phi_s\tar;
                a_rec(Ix(i),Iy(j)) = 0;
                c1_rec(:,Ix(i),Iy(j)) = -a_lsq(2:3)./a_lsq(1);%sqrt((abs(a_lsq(2))+abs(a_lsq(3))/(2*(abs(a_lsq(1)))))); %sqrt(1/dt)*
                c2_rec(:,Ix(i),Iy(j))  = pinv([Phi_tensor(:,5,ind_trans),Phi_tensor(:,6,ind_trans)])*Phi_tensor(:,2,ind_trans);%sqrt(mean(pinv([Phi_tensor(:,5,ind_trans),Phi_tensor(:,6,ind_trans)])*Phi_tensor(:,2,ind_trans)));
                c3_rec(Ix(i),Iy(j)) = -(a_lsq(2)+a_lsq(3))/a_lsq(1)/2;%CORRECT
                a1(Ix(i),Iy(j)) = 0;
                a2(Ix(i),Iy(j)) = a_lsq(1);
                a3_rec(Ix(i),Iy(j)) = 0;
            elseif(ind_waveeq4(ind_trans)==1)
                Phi_s=[Phi_tensor(:,1,ind_trans),Phi_tensor(:,2,ind_trans),Phi_tensor(:,5,ind_trans),Phi_tensor(:,6,ind_trans)];
                Phi_s=[Phi_s;vec_cons_4];
                a_lsq = Phi_s\tar;
                a3_lsq(:,Ix(i),Iy(j))  = [(1/dt)*Phi_tensor(:,1,ind_trans),Phi_tensor(:,5,ind_trans),Phi_tensor(:,6,ind_trans)]\Phi_tensor(:,2,ind_trans);
                a_rec(Ix(i),Iy(j)) = a_lsq(1)/a_lsq(2);% UNABLE TO RECOVER CORRECTLY USING 49 TIME POINTS
                c1_rec(:,Ix(i),Iy(j)) = -a_lsq(3:4)./(a_lsq(2));%sqrt((abs(a_lsq(3))+abs(a_lsq(4))/(2*(abs(a_lsq(2))))));      
                c2_rec(:,Ix(i),Iy(j))  =  pinv([Phi_tensor(:,5,ind_trans),Phi_tensor(:,6,ind_trans)])*Phi_tensor(:,2,ind_trans);%sqrt(mean(pinv([Phi_tensor(:,5,ind_trans),Phi_tensor(:,6,ind_trans)])*Phi_tensor(:,2,ind_trans)));
                c3_rec(Ix(i),Iy(j)) = -(a_lsq(3)+a_lsq(4))/a_lsq(2)/2;
                a1(Ix(i),Iy(j)) = a_lsq(1);
                a2(Ix(i),Iy(j)) = a_lsq(2);
                a3_rec(Ix(i),Iy(j)) = a3_lsq(1,Ix(i),Iy(j));
            end
            %diff56 = abs(a_lsq(end)-a_lsq(end-1));
           % end
        end
    end
%end
%%
figure
imagesc(sqrt(c3_rec),'AlphaData',~isnan(c3_rec))
colorbar
%caxis([450 700])
axis square
%set(gca, 'YDir','normal')
ylabel('y (mm)','interpreter','latex','FontSize',22)
xlabel('x (mm)','interpreter','latex','FontSize',22)
%name = strcat(num2str(10*cfreq),'kHz');
%title(name, 'interpreter','latex','FontSize',22)
% xticks([10 30 ])
% xticklabels({'30','40','50','60','70'})
ylim([6,95])
xlim([6,95])
ax=gca
ax.TickLabelInterpreter = 'latex';
ax.FontSize=28;
%%
% figure
% imagesc(a_rec,'AlphaData',~isnan(a_rec))
% colorbar
% %caxis([450 700])
% axis square
% %set(gca, 'YDir','normal')
% ylabel('y (m)','interpreter','latex','FontSize',22)
% xlabel('x (m)','interpreter','latex','FontSize',22)
% title('$50 \rm{kHz}$ ', 'interpreter','latex','FontSize',22)
% % xticks([10 30 ])
% % xticklabels({'30','40','50','60','70'})
% ylim([6,95])
% xlim([6,95])
% ax=gca
% ax.TickLabelInterpreter = 'latex';
% ax.FontSize=28;
%%
% 
% figure
% % subplot(1,2,1)
% imagesc(squeeze(sqrt(c3_rec(:,:))),'AlphaData',~isnan(squeeze(c1_rec(1,:,:))))
% %hold on
% %plot(71,63,'r.','MarkerSize',15);
% %hold on
% %plot(19,32,'r.','MarkerSize',15);
% colorbar
% caxis([450 750])
% axis square
% set(gca, 'YDir','normal')
% ylabel('y (mm)','interpreter','latex','FontSize',22)
% xlabel('x (mm)','interpreter','latex','FontSize',22)
% title('40~\rm{kHz}', 'interpreter','latex','FontSize',22)
% % xticks([10 30 ])
% % xticklabels({'30','40','50','60','70'})
% ylim([6,95])
% xlim([6,95])
% ax=gca
% ax.TickLabelInterpreter = 'latex';
% ax.FontSize=28;
% 
% exportgraphics(gcf,'40_cdist_noisy.png')

%%
% subplot(1,2,2)
% imagesc(squeeze(a3_rec(:,:)),'AlphaData',~isnan(a_rec))
% hold on
% plot(71,63,'r.','MarkerSize',15);
% hold on
% plot(19,32,'r.','MarkerSize',15);
% colorbar
% %caxis([0 2])
% axis square
% %set(gca, 'YDir','normal')
% %ylabel('y (m)','interpreter','latex','FontSize',22)
% xlabel('x (m)','interpreter','latex','FontSize',22)
% title('True $\alpha$', 'interpreter','latex','FontSize',22)
% % xticks([10 30 ])
% % xticklabels({'30','40','50','60','70'})
% ylim([6,95])
% xlim([6,95])
% ax=gca
% ax.TickLabelInterpreter = 'latex';
% ax.FontSize=28;
%% Mean, MAD
% chat = squeeze(sqrt(c3_rec(:,:)));
% me(cfreq-2) = nanmean(chat(:))
% ma(cfreq-2) = mad(chat(:))
% mp(cfreq-2) = ma(cfreq-2)/me(cfreq-2)
%% RMSE, change c, c_rec to alpha, a_rec to calculate the RMSE for alpha
%  rsum=0;
c_rec=sqrt(c3_rec);
%  for i=1:100
%      for j=1:100
%          if(~isnan(c_rec(i,j)))
%             rsum=rsum+(c(i+10,j+10)-c_rec(i,j))^2;
%          end
%      end
%  end
% 
% rmse = sqrt(rsum/(length(Ix)*length(Iy)))
