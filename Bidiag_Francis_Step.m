function Bi_next = Bidiag_Francis_Step( Bi )
[ m , m ] = size(Bi);                                                                %get dimensions of Bi
G = Givens_rotation([(Bi(1,1)*Bi(1,1)) - (Bi(m,m)*Bi(m,m)); (Bi(1,2) * Bi(1,1))]);   %perform first Given's rotation on Bi  
Bi(1:2, 1:2) = Bi(1:2, 1:2) * G;                                                     %only 2x2 region updated 
for i = 1:m-1                                                                        %repeat this multiple times                                 
    G = Givens_rotation([Bi(i,i); Bi(i+1,i)])';                                      %perform next Given's rotation on Bi
    if(i == m-1)                                                                     %If at the end
        Bi(i:i+1, i:i+1) = G * Bi(i:i+1, i:i+1);                                     %only 2x2 region updated
    elseif(i < m-1)                                                                  %otherwise
        Bi(i:i+1, i:i+2) = G * Bi(i:i+1, i:i+2);                                     %2x3 region is updated
    end
    if (i < m-1)                                                                     %handle bulge above superdiagonal
        G  = Givens_rotation([Bi(i,i+1); Bi(i,i+2)]);                                %perform next Given's rotation on Bi
        Bi(i:i+2,i+1:i+2) = Bi(i:i+2,i+1:i+2) * G;                                   %3x2 region is updated
    end
end
Bi_next = Bi;                                                                        %another iteration bites the dust