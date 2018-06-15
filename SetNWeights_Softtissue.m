
function G = SetNWeights_Softtissue(G, Img, CurL, PropL, Sigma, K)
%SETEDGE この関数の概要をここに記述
%   詳細説明をここに記述


L00 = double(CurL(G.Hi) ~= CurL(G.Hj));
L01 = double(CurL(G.Hi) ~= PropL(G.Hj));
L10 = double(PropL(G.Hi) ~= CurL(G.Hj));
L11 = double(PropL(G.Hi) ~= PropL(G.Hj));
S00 = sub2ind(size(Sigma), CurL(G.Hi), CurL(G.Hj));
S01 = sub2ind(size(Sigma), CurL(G.Hi), PropL(G.Hj));
S10 = sub2ind(size(Sigma), PropL(G.Hi), CurL(G.Hj));
S11 = sub2ind(size(Sigma), PropL(G.Hi), PropL(G.Hj));

Z00 = (Img(G.Hi)-Img(G.Hj))./Sigma(S00);
Z01 = (Img(G.Hi)-Img(G.Hj))./Sigma(S01);
Z10 = (Img(G.Hi)-Img(G.Hj))./Sigma(S10);
Z11 = (Img(G.Hi)-Img(G.Hj))./Sigma(S11);


%G.H00(:) = L00 .* exp(-X .* K(S00) ./ (2*Sigma(S00).^2)) ./ G.dist;
%G.H01(:) = L01 .* exp(-X .* K(S01) ./ (2*Sigma(S01).^2)) ./ G.dist;
%G.H10(:) = L10 .* exp(-X .* K(S10) ./ (2*Sigma(S10).^2)) ./ G.dist;
%G.H11(:) = L11 .* exp(-X .* K(S11) ./ (2*Sigma(S11).^2)) ./ G.dist;

G.H00(:) = L00 .* func(Z00,K(S00))./G.dist;
G.H01(:) = L01 .* func(Z01,K(S01))./G.dist;
G.H10(:) = L10 .* func(Z10,K(S10))./G.dist;
G.H11(:) = L11 .* func(Z11,K(S11))./G.dist;


%{
G.F00 =  func(Z00,K(S00));
G.F01 = func(Z01,K(S01));
G.F10 = func(Z10,K(S10));
G.F11 = func(Z11,K(S11));

G.Z00 =  Z00;
G.Z01= Z01;
G.Z10 = Z10;
G.Z11 = Z11;

G.K00 = K(S00);
G.K01 = K(S01);
G.K10 = K(S10);
G.K11 = K(S11);
%}
end

function Output = func(Z,K)   
   KZ = K.*Z;
 %  KZ(KZ<0) = 0;
   Output = exp(-(KZ.^2)./2);
  
end
