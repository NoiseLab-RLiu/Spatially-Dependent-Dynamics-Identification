dict = Phi_tensor;%Integral_res;


Phi_tensor = real(dict);
[M,D,N] = size(dict);

a_raw=zeros(D,N);
a_res=zeros(D,N);

%% define s
% For 1d
%vec_constraint = [1,1,1,-1,-1,1];
% For 2d
vec_constraint = [1,1,1,1,-1,-1,-1,-1,1];
%vec_constraint = [1,1,1,1,1,1,1,1,1];
tic
for i=1:N
    %if(sum(sum(Integral_res(:,:,i)))~=0)
        % select dictionary for location i
        Phi_s = squeeze(dict(:,:,i));%squeeze(Phi_tensor(:,:,i));
        n = size(Phi_s,2);
        % normalize dictionary
        Phi_sn=zeros(size(Phi_s,1),size(Phi_s,2));
        for d=1:size(Phi_s,2)
            Phi_sn(:,d) = Phi_s(:,d)./norm(Phi_s(:,d)); % introduce nan column if the column is all 0.
        end
        % append constraint s to the normalized dictionary
        Phi_ext = [Phi_sn;vec_constraint];
        tar = zeros(size(Phi_s,1)+1,1);
        tar(end) = 1;
       % set lambda
        lam=0.2*max(abs(Phi_ext'*tar))/length(tar);
        % lasso
        [a_ext_tmp, s] = lasso(Phi_ext,tar,'Lambda',lam,'Intercept',false,'RelTol',1e-8,'MaxIter',10^6);
        a_raw(:,i) = a_ext_tmp;
    %end
end
toc


%%

% figure
% stem(Cond)

figure
imagesc((log10(abs((a_raw)))))
%imagesc((((abs(a_raw)))))
% cbh = colorbar('XTickLabel',{'$\leq -5$','-4','-3','-2','-1','0'}, ...
%                'XTick', -5:0,'TickLabelInterpreter','latex')
ylabel('Index of the entry $i$ in ${\bar{\mathbf{a}}}_n$','interpreter','latex')
xlabel('Location index $n$','interpreter','latex')
%name = strcat('$log_{10}|{\bar{\mathbf{a}}}_n(i)|$ for $\space$',strcat(num2str(10*cfreq),' kHz'));
%title('$log_{10}|{\bar{\mathbf{a}}}_n(i)|$','interpreter','latex')%
%title(name,'interpreter','latex')%
yticklabels({'1','2','3','4','5','6','7','8','9'})
%xticks([1,200,400,600,784])
xticks([1,3000,6000, 8100])
%xticks([1,20,40,60,80,97])
ax=gca
ax.FontSize=20;
colorbar
%caxis([0,0.4])
caxis([-5 0])
colormap(hot)
set(gca,'TickLabelInterpreter','latex')

%% Thresholding
epsilon = 0.01;
indicator = zeros(D,N);
for i=1:N
   %  if(sum(sum(Integral_res(:,:,i)))~=0)
        [maxv,maxi] = max(abs(a_raw(:,i)));
        ind=find(abs(a_raw(:,i))>=epsilon*maxv);%maxv);//working Th rel0.05
        indicator(ind,i) = 1;
  %   end
end
