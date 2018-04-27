
Nx = 1000;
No = 50;
sx = 1;%0.01;%0.1;
sA = 1;%0.3;
sY = 0.05;

Nylist = round(logspace( 3,6,10));
Nylist = Nylist(1:end-1);
v = zeros(size(Nylist));
c = zeros(size(Nylist));

for mj = 1:length(Nylist)
    display(mj/length(Nylist))
    
    Ny = Nylist(mj);
    
    A = randn(Ny,Nx).*(rand(Ny,Nx)<sA);
    x = randn(Nx,No).*(rand(Nx,No)<sx);
    
    y = A*x;
    %tmp = sort(reshape(y,Ny*No,1),'descend');
    %theta = tmp(round(Ny*No*sY));
    %y = y.*(y>theta);
    %     y = 1*(y>1.25);
    
    for j=1:Ny
        [~,loc] = sort(y(j,:),'descend');
        y(j,loc(round(sY*No)+1:end))=0;
        y(j,loc(1:round(sY*No)))=1;
    end
    
    for j=1:No
        x(:,j) = x(:,j)-mean(x(:,j));
        x(:,j) = x(:,j)/norm(x(:,j));
        y(:,j) = y(:,j)-mean(y(:,j));
        y(:,j) = y(:,j)/norm(y(:,j));
    end
    
    ctr = 0;
    inCorr = zeros( (No^2-No)/2 , 1);
    outCorr = zeros( (No^2-No)/2 , 1);
    
    for i=1:No
        for j=1:i-1
            ctr = ctr + 1;
            inCorr(ctr) = x(:,i)'*x(:,j);
            outCorr(ctr) = y(:,i)'*y(:,j);
        end
    end
    
    v(mj) = var(outCorr);
    tmp = corrcoef( inCorr, outCorr );
    c(mj) = tmp(2);
end

figure; set(gcf,'color','w')
semilogx(Nylist,v,'.','markersize',20);
box off
xlabel('N_y','fontsize',16)
ylabel('variance in output correlation','fontsize',16)

figure; set(gcf,'color','w')
semilogx(Nylist,c,'.','markersize',20);
box off
xlabel('N_y','fontsize',16)
ylabel('correlation between input/output correlation','fontsize',16)

% figure; set(gcf,'color','w')
% semilogx(Mlist,c2,'.','markersize',20);
% box off
% xlabel('N_y','fontsize',16)
% ylabel('correlation between output correlation of 2 mice','fontsize',16)

