clear

ord=9;
zho=2^ord;
zhom1=zho-1;

n=prnum(ord);
n=2*n-1;
n=n(1:end-1);
% Vorsicht, wenns nicht klappt am 
% minus rumfummeln
n=-n;

Dnn=corr(n,n);

h=randn(1,length(n));
%h=1:length(n);
h=h(:);

y=myconv(n,h);
Dny=corr(n,y);
N=zeros(zhom1,zhom1);
ii=n(:).';
%nein ii=fliplr(ii); 
for(k=1:length(ii)) 
  N(k,:)=ii;
  ii=[ii(end) ii(1:end-1)]; %linksrum/rechtsrum?
end  

myconv(Dnn,h);
% corr(n,y) so nicht !
%corr(y,n) 

S1=[zeros(1,zho-1) ;eye(zho-1)];
S2=S1.';

H=hadamard(zho);
G=H(2:end,2:end);
Ns=((1-N)/2);

iia=2.^(0:ord-1)*Ns(1:ord,:);
[iic,iib]=sort(iia);
presort=iib;
Nst=Ns(:,iib);
P1=inv(inv(Ns)*Nst);

%P2=N*inv(P1)*inv(G);
%P2=N*inv(P1)*(G'-1)/zho;

P2s=N*P1'*(G'-1)/zho;

iii=Nst(:,2.^(0:ord-1))*(2.^(0:ord-1).');
[iik,iim]=sort(iii);
Gs=Nst(iim,:);
P2=inv(Gs*inv(Nst));
postsort=iii;

%P2=P2s;

nn=fliplr(Dnn.');
MDnn=nn;
for(k=1:length(nn)-1)
    nn=[nn(end) nn(1:end-1) ];
    MDnn= [MDnn;nn];
end;

%iid=P2*S2*H*S1*P1*y;
iif=[0 ;y(presort)];
iig=hadamard_butterfly(iif);
iih=iig(2:end);
%iid=P2*iih;
iid=iih(postsort);

iie=(iid-mean(iid))/zho;
iie=[iie(2:end);iie(1)]; %1 mal rotieren
hr=inv(MDnn)*iid;

sum(abs(hr-h))
%sum(abs(iid/zho-h))
plot([h hr iie])
sum(sum(abs(P2'*N*P1'-G)))
sum(sum(abs(P2'*Ns*P1'-Gs)))
sum(sum(abs(Ns-((1-N)/2))))
sum(sum(abs(Gs-((1-G)/2))))
return


function [d]=prnum(m);
%clear
%m=15

if(m==3)
 d=[-1 -1 -1  1 -1 1 1 1];
 d=(d+1)/2;
end;

if(m==4)
d=zeros(1,2^m);d(1:m)= ones(1,m);
for(n=m+1:2^m) d(n)=xor(d(n-m+1),d(n-m));end;
end;

if(m==5)
d=zeros(1,2^m);d(1:m)= ones(1,m);
for(n=m+1:2^m) d(n)=xor(d(n-m+2),d(n-m));end;
end;

if(m==6)
d=zeros(1,2^m);d(1:m)= ones(1,m);
for(n=m+1:2^m) d(n)=xor(d(n-m+1),d(n-m));end;
end;

if(m==7)
d=zeros(1,2^m);d(1:m)= ones(1,m);
for(n=m+1:2^m) d(n)=xor(d(n-m+1),d(n-m));end;
end;

if(m==8)
d=zeros(1,2^m);d(1:m)= ones(1,m);
for(n=m+1:2^m) 
    d(n)=xor(d(n-m+2),d(n-m));
    d(n)=xor(d(n),d(n-m+3));
    d(n)=xor(d(n),d(n-m+4));
end;
end;


if(m==9)
d=zeros(1,512);d(1:9)= ones(1,9);
for(n=10:512) d(n)=xor(d(n-4),d(n-9));end;
end;

%% Alternative

if(m==91)
d=zeros(1,512);d(1:9)= ones(1,9);
for(n=10:512) d(n)=xor(d(n-5),d(n-9));end;
end;

return
