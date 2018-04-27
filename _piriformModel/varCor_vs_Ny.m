
tmpTime = clock;
rng(tmpTime(end)); %seed random number generator

No = 100; %200;
sA = .3;%0.05;
sx = 0.1;

NxList = [10,100,1000];
syList = 0.05; %0.05;%[0.05,0.1,0.2]; %
Nylist = round(logspace( 2,5,10)); %round(logspace( 3,7,20));%
Nylist = Nylist(1:end-1);

v = zeros(length(Nylist),length(NxList),length(syList));
c = zeros(length(Nylist),length(NxList),length(syList));
c2 = zeros(length(Nylist),length(NxList),length(syList));

for mj = 1:length(Nylist)
    display(mj/length(Nylist))
    
    for sj = 1:length(NxList)
        for zj = 1:length(syList)
            
            
            Ny = Nylist(mj);
            Nx = NxList(sj);
            sY = syList(zj);
            
            
            A = sparse(randn(Ny,Nx).*(rand(Ny,Nx)<sA));
            A2 = sparse(randn(Ny,Nx).*(rand(Ny,Nx)<sA));
            x = randn(Nx,No).*(rand(Nx,No)<sx);
            
            y = A*x;
            tmp = sort(reshape(y,Ny*No,1),'descend');
            theta = tmp(round(Ny*No*sY));
            y = sparse(y.*(y>theta));%y.*
            
            y2 = A2*x;
            tmp = sort(reshape(y2,Ny*No,1),'descend');
            theta = tmp(round(Ny*No*sY));
            y2 = sparse(y2.*(y2>theta));%y.*
            
            %             for j=1:Ny
            %
            %                 [~,loc] = sort(y(j,:),'descend');
            %                 y(j,loc(round(.1*No)+1:end))=0;
            %                 y(j,loc(1:round(.1*No)))=1;
            %
            %                 [~,loc] = sort(y2(j,:),'descend');
            %                 y2(j,loc(round(.1*No)+1:end))=0;
            %                 y2(j,loc(1:round(.1*No)))=1;
            %
            %             end
            
            for j=1:No
                x(:,j) = x(:,j)-mean(x(:,j));
                x(:,j) = x(:,j)/norm(x(:,j));
                y(:,j) = y(:,j)-mean(y(:,j));
                y(:,j) = y(:,j)/norm(y(:,j));
                y2(:,j) = y2(:,j)-mean(y2(:,j));
                y2(:,j) = y2(:,j)/norm(y2(:,j));
            end
            
            ctr = 0;
            inCorr = zeros( (No^2-No)/2 , 1);
            outCorr = zeros( (No^2-No)/2 , 1);
            outCorr2 = zeros( (No^2-No)/2 , 1);
            
            for i=1:No
                for j=1:i-1
                    ctr = ctr + 1;
                    inCorr(ctr) = x(:,i)'*x(:,j);
                    outCorr(ctr) = y(:,i)'*y(:,j);
                    outCorr2(ctr) = y2(:,i)'*y2(:,j);
                end
            end
            
            isOk = ~isnan(inCorr) & ~isnan(outCorr) & ~isnan(outCorr2);
            inCorr = inCorr(isOk);
            outCorr = outCorr(isOk);
            outCorr2 = outCorr2(isOk);
            
            v(mj,sj,zj) = outCorr'*outCorr; %var(outCorr);
            tmp = corrcoef( inCorr, outCorr );
            c(mj,sj,zj) = tmp(2); %inCorr'*outCorr/(norm(outCorr)*norm(inCorr)); %
            tmp = corrcoef( outCorr, outCorr2 );
            c2(mj,sj,zj) = tmp(2); %outCorr2'*outCorr/(norm(outCorr)*norm(outCorr2)); %
            
            %save('simpleCorrVarVsNy.mat','Nylist','sxList','syList','c','v','c2');

        end
    end
end

%save('simpleCorrVarVsNy.mat','Nylist','sxList','syList','c','v','c2');

figure; set(gcf,'color','w')
semilogx(Nylist,squeeze(v(:,end,:)),'linewidth',1);
box off
xlabel('N_y','fontsize',16)
ylabel('variance in output correlation','fontsize',16)
legend({'S_y=0.05','S_y=0.1','S_y=0.15','S_y=0.2','S_y=0.25','S_y=0.3'})

figure; set(gcf,'color','w')
semilogx(Nylist,squeeze(c(:,1,:)),'linewidth',1);
box off
xlabel('N_y','fontsize',16)
ylabel('correlation between input/output correlation','fontsize',16)
legend({'S_y=0.05','S_y=0.1','S_y=0.15','S_y=0.2','S_y=0.25','S_y=0.3'})

figure; set(gcf,'color','w')
%semilogx(Nylist,squeeze(c2(:,1,:)),'linewidth',1);
semilogx(Nylist,c2,'linewidth',1);
box off
xlabel('N_y','fontsize',16)
ylabel('correlation between two mice','fontsize',16)
%legend({'S_y=0.05','S_y=0.1','S_y=0.15','S_y=0.2','S_y=0.25','S_y=0.3'})
