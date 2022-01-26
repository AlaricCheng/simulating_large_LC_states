function rhox=likelyrho( rho0)
% input : rho0 is the direct sum of 4^n pauli tensor,
% output: rhox is a physical density matrix.

[v,dv]=eig(rho0);
a0=real(diag(dv));
a= sort(a0,'descend' );
b=a;
d=length(a);
reg=0;
for i=d:-1:1
    if a(i)+reg/i <0
        b(i)=0;
        reg=reg+a(i);
    else
        for j=i:-1:1
            b(j)=a(j)+reg/i;
        end
        break
    end
end
[~, ind]=sort(a0);
ind=ind(end:-1:1);
[~, ind]=sort(ind);
b0=b(ind);
rhox=zeros(d);
for i=1:d
   rhox=rhox+b0(i)* v(:,i) * v(:,i)';
end
end