function K=interpDispersion(X,Y,ZZ,TT,k,posPart,time)
n=size(posPart,1);
K=zeros(n,1);
for part=1:n
    
    for i=1:size(k,3)

       ki(i) = interp3( X, Y, TT, squeeze(k(:,:,i,:)), posPart(part,1), posPart(part,2),time, 'linear' );
       hi(i) = interp3( X, Y, TT, squeeze(ZZ(:,:,i,:)), posPart(part,1), posPart(part,2),time, 'linear' );
         
    end

        ki(isnan(ki))=0;
        
        if isnan(sum(hi)) || hi(end)>posPart(part,3) %out of domain or sinking or beaching
            K(part)=0;
        else    
            K(part)=interp1(hi,ki,posPart(part,3));
        end

end 
end