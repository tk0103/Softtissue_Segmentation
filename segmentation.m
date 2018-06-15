%Input
E1ori = load_raw('C:\Users\yourb\Desktop\NZdata\softtissue\proposed_bilateral_E1.raw','*single');
E2ori = load_raw('C:\Users\yourb\Desktop\NZdata\softtissue\proposed_bilateral_E2.raw','*single');
E3ori = load_raw('C:\Users\yourb\Desktop\NZdata\softtissue\proposed_bilateral_E3.raw','*single');
E4ori = load_raw('C:\Users\yourb\Desktop\NZdata\softtissue\proposed_bilateral_E4.raw','*single');
%resultori = load_raw('C:\Users\yourb\Desktop\mars.raw','*uint8');
resultori = load_raw('C:\Users\yourb\Desktop\NZdata\softtissue\mars_label.raw','*uint8');

siz = [436,436,63];
sizte = [436,436,5];
E1ori = reshape(E1ori,siz); E2ori = reshape(E2ori,siz); E3ori = reshape(E3ori,siz); E4ori = reshape(E4ori,siz);
resultori = reshape(resultori,siz);

%mask
cx = 215.5; %cx = 218.5;
cy = 220.5; %cy = 218.5;

r = 23/2.0/0.0930111;
[xx,yy,zz] = ndgrid(1:siz(1),1:siz(2),1:siz(3));
mask = (xx-cx).^2 + (yy-cy).^2 <= r^2;

E1 = zeros(siz); E2 = zeros(siz); E3 = zeros(siz); E4 = zeros(siz); result = zeros(siz);
E1(mask) =  E1ori(mask); E2(mask) = E2ori(mask); E3(mask) = E3ori(mask); E4(mask) = E4ori(mask); result(mask) = resultori(mask);

order = [3,4,2,1,5];
result(result==0) = 5;
result = order(result);


%テスト
mmask = mask(:,:,1:5);
E1test = E1(:,:,35:39);
E1test = E1test(mmask);
E4test = E4(:,:,35:39);
E4test = E4test(mmask);
resultte = result(:,:,35:39);
resultte = resultte(mmask);

%学習
E1train = E1(:,:,6:10);
E1train = E1train(mmask);
E4train = E4(:,:,6:10);
E4train = E4train(mmask);
resulttr = result(:,:,6:10);
resulttr = resulttr(mmask);
m = length(E1train);

%Colormap
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];

%%


%result(result==5) = 0;
imagesc(Output(:,:,5)');
%colormap gray
colormap(map);
axis tight equal off
%colormap gray
%caxis([0 0.5])
%save_raw(E1,'C:\Users\yourb\Desktop\Proposed_bilateral_E1mask.raw','*single');
%%
%EM推定
class = 4;
feature = 2;

EMini = EMinitialValue(E1train,E4train,resulttr,class,feature);
GMModel = fitgmdist([E1test,E4test],class,'start',EMini,'Options',statset('Display','iter','MaxIter',100));
Mu = GMModel.mu;
Sig = GMModel.Sigma;
proportion = GMModel.ComponentProportion;

y1 = proportion(1) * mvnpdf([E1test,E4test],Mu(1,:),Sig(:,:,1));
y2 = proportion(2) * mvnpdf([E1test,E4test],Mu(2,:),Sig(:,:,2));
y3 = proportion(3) * mvnpdf([E1test,E4test],Mu(3,:),Sig(:,:,3));
y4 = proportion(4) * mvnpdf([E1test,E4test],Mu(4,:),Sig(:,:,4));
F = [y1,y2,y3,y4];

p_l = proportion;
PP = bsxfun(@times,F,p_l);
p_x = sum(PP,2);
PP = bsxfun(@rdivide,PP,p_x);

%~引数のプレースホルダー
[~,L] = max(PP,[],2);

I = zeros(sizte);
I(mmask) = L;

%#1表示
figure;
label = I(:,:,1);
imagesc(label');
axis equal tight off
colormap(map)

PP1 = zeros(sizte);
PP1(mmask) = PP(:,1);

figure;
subplot(2,2,1)
imagesc(PP1(:,:,1)');
axis equal tight
colormap hot

PP2 = zeros(sizte);
PP2(mmask) = PP(:,2);

subplot(2,2,2)
imagesc(PP2(:,:,1)');
axis equal tight
colormap hot

PP3 = zeros(sizte);
PP3(mmask) = PP(:,3);

subplot(2,2,3)
imagesc(PP3(:,:,1)');
axis equal tight
colormap hot

PP4 = zeros(sizte);
PP4(mmask) = PP(:,4);

subplot(2,2,4)
imagesc(PP4(:,:,1)');
axis equal tight
colormap hot


RP1 = -log(PP1+eps);
RP2 = -log(PP2+eps);
RP3 = -log(PP3+eps);    
RP4 = -log(PP4+eps);
RP = [RP1(mmask),RP2(mmask),RP3(mmask),RP4(mmask)];


%%
imagesc(RP3(:,:,1)');
axis equal tight
%%

label = I(:,:,1);
imagesc(label');
axis equal tight off
colormap default
%%
[lambda,h] =ndgrid(0.2:0.1:1.3,0.2:0.1:1.9);
lambda = lambda(:);
h = h(:);
%%
sumJI = zeros(216,1);
%for n = 101:200
n =174;
GraphModel = CreateFullyConnectedGraphWithMask(mmask);
Sigmat =  abs(bsxfun(@minus,GMModel.mu(:,1),GMModel.mu(:,1)'))*h(n) +eye(4);
%Kmat = [1,1,1,1;0,1,0,0;0,1,1,1;0,1,0,1];
%Kmat = [1,1,1,1;1,1,1,1;1,1,1,1;1,1,1,1];
Kmat = [0,-1,-1,-1;1,0,1,1;1,-1,0,-1;1,-1,1,0];

CurLabel = ones(m,1);
PreLabel = zeros(m,1);
Output = zeros(sizte);
flag = 0;
%%
while(flag ~=1)
    for  aa = 1:4
        PropLabel = zeros(m,1)+aa;

        GraphModel = SetTWeights(GraphModel,RP,CurLabel,PropLabel,lambda(n),m);
        GraphModel = SetNWeights(GraphModel,E1test,CurLabel,PropLabel,Sigmat,Kmat);
        [lowerBound, labels] = qpboMex([GraphModel.Vs,GraphModel.Vt],[GraphModel.Hi,GraphModel.Hj,GraphModel.H00,GraphModel.H01,GraphModel.H10,GraphModel.H11]);
        labels = logical(labels);
        CurLabel(labels) = PropLabel(labels);
        
        Eunary = sum(GraphModel.Vs);
        Epairwise = sum(GraphModel.H00);
        E = Eunary + Epairwise;
        
    end
   % pause;
    %    disp(E); 
        if CurLabel == PreLabel
            flag = 1;
            Output(mmask) = CurLabel;
            imagesc(Output(:,:,1)'); 
            axis equal tight off
        end
        PreLabel = CurLabel;        
end
%%
JI = CalcuJI(Output,result(:,:,35:39),class);
sumJI(n) = sum(JI);
disp(n);
disp(sum(JI));
%end

%%
disp(JI);
disp(sum(JI));
%%
yy = [JIMAP; JI;];
bar(yy);


























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