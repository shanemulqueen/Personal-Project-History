%procrustes mapping of factors
%currencies= csvread('1W120_exOutliers.csv');
a = size(currencies);
window = 120;
length=a(1)-window;

C_angle=zeros(length,length);
Cr_angle=zeros(length,length);
C1_angle=zeros(length,length);
C2_angle=zeros(length,length);
C3_angle=zeros(length,length);
C4_angle=zeros(length,length);
C5_angle=zeros(length,length);

% % %can comment out from here if PCA factors are already stored
% factors=zeros(171,length*5);
% for i = 1:length
%     currencies2 = currencies([i:i+window-1],[1:end]);
%     [coeff,score,latent,tsquared,explained,mu] = pca(currencies2,'NumComponents',5);
%     factors(1:171,(5*(i-1)+1):(5*(i-1)+5))=coeff(1:171,1:5);
%     
% end
% % %to here

for ii = 1:length
B=factors(1:171,(5*(ii-1)+1):(5*(ii-1)+5)); %choose the stationary matrix
%B=rotatefactors(B,'Method','parsimax','Maxit',500); %comment this line out if not doing a rotation

for i = 1:length
    A=factors(1:171,(5*(i-1)+1):(5*(i-1)+5));
    At=transpose(A);
    C=At*B;
    [U,S,V]=svd(C);
    R=U*transpose(V);
    Ar=A*R;
    Cr_angle(ii,i)=subspace(Ar,B);              %calculate hyperplane angle
    C_angle(ii,i)=subspace(A,B);  
    %C_angle(ii,i)=subspace(Ar(:,1:2),B(:,1:2)) %if only looking at the
    %first two components
    dotss=dot(Ar(1:171,1),B(1:171,1));
    C1_angle(ii,i)=abs(acosd(dotss));
    dotss=dot(Ar(1:171,2),B(1:171,2));
    C2_angle(ii,i)=abs(acosd(dotss));
    dotss=dot(Ar(1:171,3),B(1:171,3));
    C3_angle(ii,i)=abs(acosd(dotss));
    dotss=dot(Ar(1:171,4),B(1:171,4));
    C4_angle(ii,i)=abs(acosd(dotss));
    dotss=dot(Ar(1:171,5),B(1:171,5));
    C5_angle(ii,i)=abs(acosd(dotss));

end
end

%plot subspace angles
%mesh(index,index,Cr_angle)