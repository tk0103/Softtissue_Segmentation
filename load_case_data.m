%load_path
fileID = fopen('InputPath_Soft.txt');
C = textscan(fileID,'%s');
fclose(fileID);
InputPath = C{1,1}; 
InputPath = cell2mat(InputPath);

fileID = fopen('CasePath_Soft.txt');
C = textscan(fileID,'%s');
fclose(fileID);
CasePath = C{1,1};

%load_data
E1ori = load_raw([InputPath CasePath{1,:} '.raw'],'*single');
E2ori = load_raw([InputPath CasePath{2,:} '.raw'],'*single');
E3ori = load_raw([InputPath CasePath{3,:} '.raw'],'*single');
E4ori = load_raw([InputPath CasePath{4,:} '.raw'],'*single');
resultori = load_raw([InputPath CasePath{5,:} '.raw'],'*uint8');

siz = [436,436,63];
sizte = [436,436,5];
E1ori = reshape(E1ori,siz);
E2ori = reshape(E2ori,siz);
E3ori = reshape(E3ori,siz);  
E4ori = reshape(E4ori,siz);
resultori = reshape(resultori,siz);

%%
%confirm
slice = 30;
imagesc(E1ori(:,:,slice)');
axis tight equal off
colormap gray
caxis([0 0.6])