




%%
P = GraphModel.Hi(CurLabel(GraphModel.Hi)==3 & CurLabel(GraphModel.Hj)==4);
Q = GraphModel.Hj(CurLabel(GraphModel.Hi)==3 & CurLabel(GraphModel.Hj)==4);
dif = (E1test(P) - E1test(Q))/(0.1047*1.6);
edges = [-1.5 -1.5:0.01:0.4 0.4];
xlim([-1.0 0.3]);
yyaxis left
histogram(dif,edges);
%histogram(dif);

x = edges';
z = exp(-x.^2 /2);
yyaxis right
xlabel('I_p - I_q / \sigma')
plot(x,z);
%%
 x = GraphModel.Z01(:);
 y =GraphModel.H01(:);
 x = x(GraphModel.K01<0);
 y = y(GraphModel.K01<0);
 
 %edges = [-3 -3:0.1:3 3];
 scatter(x,y,'.');
 %xlim([-50,50])

 %histogram(x,edges);

%%
fplot(@(x) exp(-x^2 /2),[0,5],'r')
hold on
fplot(@(x) 1.0,[-5,0],'r')
hold off
grid on
ylim([0 1.0])
yl=ylabel('$f\{\frac{I_p - I_q}{\sigma (L_p,L_q)}\}$');
set(yl,'Interpreter','Latex');

xl=xlabel('$\frac{I_p - I_q}{\sigma (L_p,L_q)}$');
set(xl,'Interpreter','Latex');
%ylabel('$\frac{1}{1}$')
%xlabel('I_p-I_q / \sigma(L_p,L_q)')
%%
outdash = Output- Output2;
imagesc(outdash(:,:,1)');
save_raw(Output,'C:\Users\yourb\Desktop\out.raw','*single');
%%
%テスト正解
proany = zeros(4,1);
any = zeros(4,2);
sigany = zeros(2,2,4);
any = zeros(4,2);
for aa = 1:4
tmp = resultte==aa;
any(aa,1) = mean(E1test(tmp));
any(aa,2) = mean(E4test(tmp));
sigany(:,:,aa) = cov(E1test(tmp),E4test(tmp));

tmp = tmp(:);
proany(aa) = sum(tmp);

end
proany = proany./sum(proany); 
    disp(any);
    disp(sigany);
    disp(proany);

    
%%
t=2;
ttmp = resultte==5;
ttm1 = E1test(ttmp);
ttm2 = E4test(ttmp);
scatter(ttm1,ttm2,'.');
hold on
tmp = resultte==t;
tm1 = E1test(tmp);
tm2 = E4test(tmp);
scatter(tm1,tm2,'.');
%xlim([-0.2,2.0])
%ylim([-0.2,2.0])
hold on



scatter(any(t,1),any(t,2),'k','d');
scatter(Mu(t,1),Mu(t,2),'g','x');
scatter(S.mu(t,1),S.mu(t,2),'r','s');
xlabel('E1')
ylabel('E4')
%xlim([-0.5,3])
%ylim([-0.5,2.5])
C = {'r','g','b','y'};
for i = 1:1
    for j =1:3
        protmp =plot_gaussian_ellipsoid(Mu(i,:),Sig(:,:,i),j);
        protmp.Color = C{1};
        protmp =plot_gaussian_ellipsoid(any(i,:),sigany(:,:,i),j);
        protmp.Color = C{2};
    end
end
hold off
%save_raw(rere,'C:\Users\yourb\Desktop\aa.raw','*single');
%%
imagesc(E1(:,:,6)');
colormap gray
caxis([0,0.6])
axis equal tight off

rectangle('Position',[260,160,20,20],'FaceColor','none','EdgeColor','r','LineWidth',5)
rectangle('Position',[202,234,20,20],'FaceColor','none','EdgeColor','r','LineWidth',5)
%rectangle('Position',[134,276,10,10],'FaceColor','none','EdgeColor','g','LineWidth',1)
%rectangle('Position',[95,220,20,20],'FaceColor','none','EdgeColor','g','LineWidth',1)
%%
roi1E1 = E1(260:279,160:179,6:10); %water
roi2E1 = E1(202:221,234:253,6:10); %fat
roi3E1 = E1(134:143,276:285,35:39); %bone
roi4E1 = E1(95:114,220:239,35:39); %background

roi1E4 = E3(260:279,160:179,6:10); %water
roi2E4 = E3(202:221,234:253,6:10); %fat
roi3E4 = E3(134:143,276:285,35:39); %bone
roi4E4 = E3(95:114,220:239,35:39); %background

roi1E1 = roi1E1(:); roi2E1 = roi2E1(:); roi3E1 = roi3E1(:); roi4E1 = roi4E1(:);
roi1E4 = roi1E4(:); roi2E4 = roi2E4(:); roi3E4 = roi3E4(:); roi4E4 = roi4E4(:);
%%
scatter(roi1E1,roi1E4,'.');
hold on
scatter(roi2E1,roi2E4,'.');
hold off
axis equal tight
xlabel('Energy1')
ylabel('Energy3')
xlim([0.1,0.5])
ylim([0.1,0.5])
legend('Muscle','Fat')
%%
hold on
scatter(roi2E1,roi2E4,'.');
scatter(roi3E1,roi3E4,'.');
scatter(roi4E1,roi4E4,'.');
scatter(Mu(:,1),Mu(:,2),'r','x');
scatter(S.mu(:,1),S.mu(:,2),'r','s');
scatter(any(:,1),any(:,2),'k','d');
axis equal tight
xlabel('E1')
ylabel('E4')
%xlim([-0.2,2.0])
%ylim([-0.2,2.0])
%legend('water','fat','bone','bg')
hold off

C = {'r','g','b','y'};
for i = 1:4
    for j =1:3
        protmp =plot_gaussian_ellipsoid(Mu(i,:),Sig(:,:,i),j);
        protmp.Color = C{i};
    end
end

%%
tm1 = E1(:,:,35); tm2 = E4(:,:,35);
tm1 = tm1(:); tm2 = tm2(:);
scatter(tm1,tm2,'.');
C = {'r','g','b','y'};
for i = 1:4
    for j =1:2
        protmp =plot_gaussian_ellipsoid(Mu(i,:),Sig(:,:,i),j);
        protmp.Color = C{i};
    end
end
axis equal tight
xlabel('E1')
ylabel('E4')
xlim([-0.5,2.5])
ylim([-0.5,2.5])

%%

imagesc(E1(:,:,37)');
colormap gray
caxis([0,0.6])
axis equal tight
rectangle('Position',[170,150,130,130],'FaceColor','none','EdgeColor','g','LineWidth',1)
figure;

testroi1 = E1(170:299,150:279,35:39);
testroi1 = testroi1(:);
testroi2 = E4(170:299,150:279,35:39);
testroi2 = testroi2(:);
scatter(testroi1,testroi2,'.');

C = {'r','g','b','y'};
for i = 3:4
    for j =1:3
        protmp =plot_gaussian_ellipsoid(Mu(i,:),Sig(:,:,i),j);
        protmp.Color = C{i};
    end
end
axis equal tight
xlabel('E1')
ylabel('E4')
xlim([0.1,0.45])
ylim([0.1,0.45])




%%
%plot EM
for c= 1:15
Options = statset('MaxIter',c);
GMModel = fitgmdist([E1test,E4test],6,'start',S,'Options',Options);
Mu = GMModel.mu;
Sig = GMModel.Sigma;
proportion = GMModel.ComponentProportion;
C = {'r','g','b','y','b','r'};
%scatter(E1x(1:10:end),E2x(1:10:end),10,'.')
clf;
for i = 1:6
    for j =1:2
%h = fcontour(@(x,y)mvnpdf([x,y],Mu(i,:),Sig(:,:,i)),[-0.1 0.7 -0.1 0.7]);
protmp =plot_gaussian_ellipsoid(Mu(i,:),Sig(:,:,i),j);
protmp.Color = C{i};
    end
end
drawnow;
pause(0.5);
end



%%
y =275;
rangex = [140,190];
rangey = [y-25,y+25];
test1 = E1ori(:,y,34);
test2 = E1ori(:,y,34);

%{
subplot(3,4,[1,6])
imagesc(E1ori(:,:,36)');
%imagesc(Output(:,:,2)');
%colormap(map);
caxis([0,0.6])
colormap gray
axis equal tight
rectangle('Position',[1,y-0.5,436,1],'FaceColor','none','EdgeColor','r',...
    'LineWidth',1)
xlim(rangex);
ylim(rangey);

%}
subplot(3,4,[3,8])
imagesc(Output(:,:,2)');
%colormap gray
%caxis([0,0.6])
axis equal tight
colormap(map);
 rectangle('Position',[1,y-0.5,436,1],'FaceColor','none','EdgeColor','r',...
    'LineWidth',1)
xlim(rangex);
ylim(rangey);

%{
M1 = movmean(test1,3);
subplot(3,4,[9,10])
h= plot([test1]);
h(1).Color = 'g';
xlim(rangex);
ylim([-0.1,0.7]);
legend({'Gray value'});
%}

M2 = movmean(test2,3);
subplot(3,4,[11,12])
h= plot([test2]);
h(1).Color = 'g';
xlim(rangex);
ylim([-0.1,0.7]);
legend({'Gray value'});


%%
imagesc(Output(:,:,2)');
colormap(map);
axis equal tight off
%caxis([0 0.6]);
%colormap gray
 rectangle('Position',[140,250,50,50],'FaceColor','none','EdgeColor','r',...
    'LineWidth',2)

%%
temp = zeros(sizte);
temp(mmask) = resultte;
imagesc(temp(:,:,2)');
colormap(map);
 rectangle('Position',[140,250,50,50],'FaceColor','none','EdgeColor','r',...
    'LineWidth',2)
axis equal tight off
