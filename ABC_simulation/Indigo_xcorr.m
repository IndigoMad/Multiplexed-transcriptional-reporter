function[result]=Indigo_xcorr(X,Y,location)
%=================Arguments=================
%X,a matric, one column is one array of data.
%Y,a matric, one column is one array of data.
%column number of X and Y must be the same
%location,location of Y aligning to X(1),in column
%===================Output==================
%result,a matric, one column is one array of normalized cross correlation,
%corresponding to loc.
%===========================================
    [LX,NX]=size(X);
    [LY,NY]=size(Y);
    [nl,~]=size(location);
    if NX~=NY
        error('Column number of X,Y must be the same')
    end
    result=zeros(nl,NX);
    for np=1:NX
        for loc=1:nl
            temX=X(max(location(loc)+1,1):min(LY+location(loc),LX),np);
            temY=Y(max(1,1-location(loc)):min(LY,LX-location(loc)),np);
            check=sqrt(sum(temX.^2).*sum(temY.^2));
            if check==0
                result(loc,np)=nan;
            else
                result(loc,np)=sum(temX.*temY)/check;
            end
        end
    end
end
